import os
import re
import json
import argparse
import pandas as pd
import numpy as np
import pickle as pkl
from pytximport import tximport
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text


def main(args):
    ##############
    ## METADATA ##
    ##############
    # Use input sample IDs and clean up metadata
    metadata = pd.read_csv(args.metadata)
    metadata = metadata.loc[:, ~metadata.columns.str.contains("^Unnamed")]
    metadata.index = metadata["ASAP_sample_id"].astype(str) + "_" + metadata["replicate"].astype(str)

    # Remove samples with missing annotations
    samples_to_keep =  ~(metadata["batch"].isna() | metadata["condition_id"].isna())
    metadata = metadata.loc[samples_to_keep]

    # Remove blacklisted samples
    sample_ids_df = pd.read_csv(args.sample_ids, sep="\t", header=None, usecols=[1])
    sample_ids = sample_ids_df.iloc[:, 0]
    metadata = metadata[metadata.index.isin(sample_ids)]

    # Check condition_id contains only valid entries as per current CDE version: ["PD", "Control", "Prodromal", "Other"]
    metadata["condition_id"] = metadata["condition_id"].str.strip()
    valid_conditions = {"PD", "Control", "Prodromal", "Other"}
    actual_conditions = set(metadata["condition_id"].unique())
    invalid_conditions = actual_conditions - valid_conditions
    if invalid_conditions:
        raise ValueError(f"Invalid condition_id values found: {invalid_conditions}")


    ############
    ## COUNTS ##
    ############
    new_column_names = ["transcript_id", "gene_id"]
    gene_map = pd.read_csv(args.gene_map, names=new_column_names, header=0)
    with open(args.gene_ids_and_names, "r") as file:
        gtf_gene_ids_and_names = json.load(file)

    path = os.getcwd()
    samples = metadata.index.tolist()
    files = [os.path.join(path, f"{sample}_salmon_quant/quant.sf") for sample in samples]
    files_dict = dict(zip(samples, files))

    txi_counts = tximport(
        file_paths=files,
        data_type="salmon",
        transcript_gene_map=gene_map,
    )

    # DESeq2 only takes integers and expects a count df with sample (index) by gene (or AnnData object, txi_counts)
    counts_int = pd.DataFrame(
        np.round(txi_counts.X).astype(int),
        index=txi_counts.obs.index,
        columns=txi_counts.var.index,
    )
    counts_int.index = counts_int.index.map({file_name: sample_name for sample_name, file_name in files_dict.items()})

    # Remove genes with less than 10 read counts total
    genes_to_keep = counts_int.columns[counts_int.sum(axis=0) >= 10]
    counts_int = counts_int[genes_to_keep]

    # Note: Single factor analysis vs. multifactor analysis requires more manual coding
    if metadata["batch"].nunique() > 1:
        design_factors = ["batch", "condition_id"]
    else:
        design_factors = ["condition_id"]
    print(f"Using design factors:\n{design_factors}")

    dds = DeseqDataSet(
        counts=counts_int,
        metadata=metadata,
        design_factors=design_factors, # From PyDESeq2: "UserWarning: Same factor names in the design contain underscores ('_'). They will be converted to hyphens ('-')."
    )

    # Fit dispersions and LFCs
    dds.deseq2()
    with open(f"{args.team_id}.{args.salmon_mode}.dds.pkl", "wb") as f:
        pkl.dump(dds, f)

    # Statistical analysis
    log2_fc_threshold = 1
    padj_threshold = 0.05
    # Run pairwise contrasts for each non-Control condition vs. Control
    test_conditions = [c for c in metadata["condition_id"].unique() if c != "Control"]
    if not test_conditions:
        raise ValueError("No test conditions found (only Control in dataset)")

    all_results = {}
    for condition in test_conditions:
        print(f"Running contrast: {condition} vs Control")
        stat_res = DeseqStats(
            dds,
            contrast=["condition_id", condition, "Control"],
        )
        stat_res.summary()
        results_df = stat_res.results_df
        results_df["gene_name"] = results_df.index.map(gtf_gene_ids_and_names)
        results_df["contrast"] = f"{condition}_vs_Control"
        sig_genes = results_df[(results_df["padj"] < padj_threshold) & (results_df["log2FoldChange"].abs() > log2_fc_threshold)]
        sig_genes.to_csv(f"{args.team_id}.{args.salmon_mode}.{condition}_vs_Control.pydeseq2_significant_genes.csv", index_label="ensembl_gene_id")
        all_results[condition] = results_df


    ###################
    ## VISUALIZATION ##
    ###################
    # Volcano plot per contrast
    for condition, results_df in all_results.items():
        results_df["-log10(padj)"] = -np.log10(results_df["padj"])
        results_df["color"] = np.where(
            (results_df["padj"] < padj_threshold) & (results_df["log2FoldChange"] > log2_fc_threshold), "red",
            np.where(
                (results_df["padj"] < padj_threshold) & (results_df["log2FoldChange"] < -log2_fc_threshold), "blue",
                "grey"
            )
        )
        plt.figure(figsize=(10, 6))
        sns.scatterplot(
            x="log2FoldChange",
            y="-log10(padj)",
            data=results_df,
            hue="color",
            palette={"red": "red", "blue": "blue", "grey": "grey"},
            alpha=0.6,
            edgecolor=None,
            legend=False,
        )
        plt.axhline(y=-np.log10(padj_threshold), color="black", linestyle="--", linewidth=1)
        plt.axvline(x=log2_fc_threshold, color="black", linestyle="--", linewidth=1)
        plt.axvline(x=-log2_fc_threshold, color="black", linestyle="--", linewidth=1)

        top_ten_genes = results_df.nsmallest(10, "padj")
        texts = []
        x_list = []
        y_list = []
        for i, row in top_ten_genes.iterrows():
            text = plt.text(
                row["log2FoldChange"],
                row["-log10(padj)"],
                row["gene_name"],
                ha="center",
                fontsize=8,
                color="black",
            )
            texts.append(text)
            x_list.append(row["log2FoldChange"])
            y_list.append(row["-log10(padj)"])

        adjust_text(
            texts,
            x=x_list,
            y=y_list,
        )

        plt.xlabel("Log2 Fold Change")
        plt.ylabel("-Log10 Adjusted P-value")
        plt.title(f"Volcano Plot: {condition} vs Control")
        plt.savefig(f"{args.team_id}.{args.salmon_mode}.{condition}_vs_Control.volcano_plot.png", dpi=300, bbox_inches="tight")
        plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Differential gene expression analysis using Salmon quantification files with PyDESeq2 by comparing non-controls and controls"
    )
    parser.add_argument(
        "-t",
        "--team-id",
        type=str,
        required=True,
        help="Team ID"
    )
    parser.add_argument(
        "-i",
        "--sample-ids",
        type=str,
        required=True,
        help="Sample IDs in a team"
    )
    parser.add_argument(
        "-m",
        "--metadata",
        type=str,
        required=True,
        help="Table containing all sample information including batch, condition, etc. used for pairwise condition"
    )
    parser.add_argument(
        "-g",
        "--gene-map",
        type=str,
        required=True,
        help="Table containing mapped transcript IDs and gene IDs that must be in this order"
    )
    parser.add_argument(
        "-n",
        "--gene-ids-and-names",
        type=str,
        required=True,
        help="JSON containing mapped gene IDs and gene names"
    )
    parser.add_argument(
        "-s",
        "--salmon-mode",
        type=str,
        required=True,
        help="Salmon mode used to quantify transcripts in order to name outputs"
    )

    args = parser.parse_args()

    main(args)
