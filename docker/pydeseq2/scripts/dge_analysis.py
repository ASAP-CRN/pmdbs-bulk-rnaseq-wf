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


def main(args):
    ##############
    ## METADATA ##
    ##############
    # Remove samples with missing annotations
    metadata = pd.read_csv(args.metadata, index_col="sample_id")
    samples_to_keep =  ~(metadata.batch.isna() | metadata.condition_id.isna())
    metadata = metadata.loc[samples_to_keep]

    # Map condition_id to intervention_id (i.e. grouping conditions to "Case", "Control", and "Other")
    condition_dict = pd.read_csv(args.condition_dict)
    metadata_merged = pd.merge(metadata, condition_dict[["condition_id", "intervention_id"]], on="condition_id", how="left")
    metadata_merged.set_index(metadata.index, inplace=True)


    ############
    ## COUNTS ##
    ############
    new_column_names = ["transcript_id", "gene_id"]
    gene_map = pd.read_csv(args.gene_map, names=new_column_names, header=0)
    with open(args.gene_ids_and_names, "r") as file:
        gtf_gene_ids_and_names = json.load(file)

    path = os.getcwd()
    samples = metadata_merged.index.tolist()
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
    if (metadata_merged.groupby("batch").size() > 1).any():
        design_factors = ["batch", "intervention_id"]
    else:
        design_factors = ["intervention_id"]
    print(f"Using design factors:\n{design_factors}")

    dds = DeseqDataSet(
        counts=counts_int,
        metadata=metadata_merged,
        design_factors=design_factors, # From PyDESeq2: "UserWarning: Same factor names in the design contain underscores ('_'). They will be converted to hyphens ('-')."
    )

    # Fit dispersions and LFCs
    dds.deseq2()
    with open(f"{args.cohort_id}.{args.salmon_mode}.dds.pkl", "wb") as f:
        pkl.dump(dds, f)

    # Statistical analysis
    log2_fc_threshold = 1
    padj_threshold = 0.05
    excluded_samples = metadata_merged[metadata_merged["intervention_id"] == "Other"]
    excluded_samples_filtered = excluded_samples[["batch", "condition_id", "intervention_id"]]
    if not excluded_samples_filtered.empty:
        print(f"[WARNING] Samples and levels excluded for contrast:\n{excluded_samples_filtered}")
    # Ensure "Control" is at the end of the list so it's ref_level
    stat_res = DeseqStats(
        dds,
        contrast=["intervention-id", "Case", "Control"], # Comparing intervention_id (adjusted condition) only but model uses information from design factors
    )
    stat_res.summary()
    results_df = stat_res.results_df
    results_df["gene_name"] = results_df.index.map(gtf_gene_ids_and_names)
    sig_genes = results_df[(results_df["padj"] < padj_threshold) & (results_df["log2FoldChange"].abs() > log2_fc_threshold)]
    sig_genes.to_csv(f"{args.cohort_id}.{args.salmon_mode}.pydeseq2_significant_genes.csv")


    ###################
    ## VISUALIZATION ##
    ###################
    # Volcano plot
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
    for i, row in top_ten_genes.iterrows():
        plt.annotate(
            row["gene_name"],
            (row["log2FoldChange"], row["-log10(padj)"]),
            textcoords="offset points",
            xytext=(0,10),
            ha="center",
        )

    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 Adjusted P-value")
    plt.title("Volcano Plot of PyDESeq2 Results")
    plt.savefig(f"{args.cohort_id}.{args.salmon_mode}.volcano_plot.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Differential gene expression analysis using Salmon quantification files with PyDESeq2"
    )
    parser.add_argument(
        "-c",
        "--cohort-id",
        type=str,
        required=True,
        help="Cohort ID"
    )
    parser.add_argument(
        "-d",
        "--condition-dict",
        type=str,
        required=True,
        help="Table containing condition and intervention IDs used to categorize conditions into broader groups for pairwise condition"
    )
    parser.add_argument(
        "-m",
        "--metadata",
        type=str,
        required=True,
        help="Table containing all sample information including batch, condition, etc."
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
