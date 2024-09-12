import os
import argparse
import pandas as pd
import numpy as np
import pydeseq2
import pytximport
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA


def main(args):
    ##############
    ## METADATA ##
    ##############

    # Remove samples with missing annotations
    metadata = pd.read_csv(args.metadata)
    samples_to_keep =  ~(metadata.batch.isna() | metadata.condition.isna())
    metadata = metadata.loc[samples_to_keep]


    ############
    ## COUNTS ##
    ############

    new_column_names = ["transcript_id", "gene_id"]
    gene_map = pd.read_csv(args.gene_map, names=new_column_names, header=0)
    blacklist_genes = pd.read_csv(args.blacklist_genes)

    path = os.getcwd()
    samples = metadata["ASAP_sample_id"].tolist()
    files = [os.path.join(path, "/", sample, "_salmon_quant/quant.sf") for sample in samples]
    #files_dict = dict(zip(samples, files))

    txi_counts = tximport(
        files,
        "salmon",
        gene_map,
    )

    # Note: Single factor analysis vs. multifactor analysis requires more manual coding
    dds = DeseqDataSet(
        counts=txi_counts,
        clinical=metadata,
        design_factors=["batch", "condition"],
    )

    ## Remove genes with less than 10 read counts total
    dds.counts = dds.counts.loc[dds.counts.sum(axis=1) >= 10]

    # Remove blacklist genes
    dds.counts = dds.counts.loc[~dds.counts.index.isin(blacklist_genes)]

    # Fit dispersions and LFCs
    dds.deseq2()

    # Statistical analysis
    log2_fc_threshold = 1
    padj_threshold = 0.05
    unique_conditions = metadata["condition"].unique().tolist()
    unique_conditions.insert(0, "condition")
    stat_res = DeseqStats(
        dds,
        contrast=unique_conditions, # Comparing condition only but model uses information from both the condition and batch variables
    )
    results_df = stat_res.results_df
    sig_genes = results_df[(results_df["padj"] < padj_threshold) & (results_df["log2FoldChange"].abs() > log2_fc_threshold)]
    sig_genes.to_csv(f"{args.cohort_id}.{args.salmon_mode}.pydeseq2_significant_genes.csv")


    ###################
    ## VISUALIZATION ##
    ###################

    # PCA plot
    normalized_counts = dds.normalized_counts
    # PCA expects samples as rows, genes as columns
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(normalized_counts.T)
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
    pca_df["condition"] = metadata["condition"].values
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x="PC1",
        y="PC2",
        hue="condition", # TODO
        data=pca_df,
        s=100,
        palette="Set2",
        alpha=0.7,
    )
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
    plt.title("PCA of Normalized Counts")
    plt.legend(title="Condition", loc="upper right")
    plt.savefig(f"{args.cohort_id}.{args.salmon_mode}.pca_plot.png", dpi=300, bbox_inches="tight")

    plt.close("all")

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
    )
    plt.axhline(y=-np.log10(padj_threshold), color="black", linestyle="--", linewidth=1)
    plt.axvline(x=log2_fc_threshold, color="black", linestyle="--", linewidth=1)
    plt.axvline(x=-log2_fc_threshold, color="black", linestyle="--", linewidth=1)

    top_ten_genes = results_df.nsmallest(10, "padj")
    for i, row in top_ten_genes.iterrows():
        plt.annotate(
            row["gene"],
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
        "-c"
        "--cohort-id",
        type=str,
        required=True,
        help="Cohort ID"
    )
    parser.add_argument(
        "-m"
        "--metadata",
        type=str,
        required=True,
        help="Table containing all sample information including batch, condition, etc."
    )
    parser.add_argument(
        "-g"
        "--gene-map",
        type=str,
        required=True,
        help="Table containing mapped transcript IDs and gene IDs that must be in this order"
    )
    parser.add_argument(
        "-b"
        "--blacklist-genes",
        type=str,
        required=True,
        help="BED file containing the ENCODE Blacklist genes"
    )
    parser.add_argument(
        "-s"
        "--salmon-mode",
        type=str,
        required=True,
        help="Salmon mode used to quantify transcripts in order to name outputs"
    )

    args = parser.parse_args()

    main(args)
