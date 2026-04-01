import argparse
import os
import re
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer


def extract_contrast(filename):
    """Extract contrast name (e.g., 'PD_vs_Control') from a DEG CSV filename.

    Expected filename pattern: {team_id}.{salmon_mode}.{contrast}.pydeseq2_significant_genes.csv
    """
    basename = os.path.basename(filename)
    # Remove the suffix and split to get the contrast portion
    match = re.search(r"\.(\w+_vs_Control)\.pydeseq2_significant_genes\.csv$", basename)
    if match:
        return match.group(1)
    raise ValueError(f"Could not extract contrast name from filename: {basename}")


def main(args):
    ##########################################################
    ## OVERLAPPING DEGS ONLY FOR CROSS TEAM COHORT ANALYSIS ##
    ##########################################################
    if args.n_teams > 1:
        # Group DEG files by contrast, then find overlapping genes per contrast across teams
        contrast_files = {}
        for file in args.degs:
            contrast = extract_contrast(file)
            if contrast not in contrast_files:
                contrast_files[contrast] = []
            contrast_files[contrast].append(file)

        for contrast, files in contrast_files.items():
            dfs = []
            for file in files:
                df = pd.read_csv(file)
                # Infer team_id from filename: {team_id}.{salmon_mode}.{contrast}...
                basename = os.path.basename(file)
                team_id = basename.split(f".{args.salmon_mode}.")[0]
                df["team_id"] = team_id
                dfs.append(df)
            combined_degs = pd.concat(dfs, ignore_index=True)
            combined_degs.set_index(combined_degs.columns[0], inplace=True)
            grouped = combined_degs.groupby("dataset_id").apply(lambda x: set(x.index), include_groups=False)
            gene_dataset_counts = pd.Series(
                {gene: sum(gene in dataset_genes for dataset_genes in grouped) for gene in combined_degs.index.unique()}
            )
            common_degs = sorted(gene_dataset_counts[gene_dataset_counts >= 2].index)
            common_degs_df = combined_degs.loc[combined_degs.index[combined_degs.index.isin(common_degs)]]
            common_degs_df = common_degs_df.reset_index().drop_duplicates()
            common_degs_df.to_csv(
                f"{args.cohort_id}.{args.salmon_mode}.{contrast}.overlapping_significant_genes.csv",
                index=False
            )
            print(f"Found {len(common_degs)} overlapping DEGs for {contrast} across {len(files)} datasets (present in ≥2 datasets)")

            gene_team_overlap_df = pd.DataFrame(
                {dataset: combined_degs[combined_degs["dataset_id"] == dataset].index for dataset in grouped.index},
                index=common_degs
            )
            gene_team_overlap_df = gene_team_overlap_df.apply(lambda col: gene_team_overlap_df.index.isin(col.dropna()), axis=0).astype(bool)
            gene_team_overlap_df["n_datasets"] = gene_team_overlap_df.sum(axis=1)
            gene_team_overlap_df.index.name = "gene_id"
            gene_team_overlap_df.to_csv(f"{args.cohort_id}.{args.salmon_mode}.{contrast}.overlapping_significant_genes_by_dataset.csv")


    ###################
    ## VISUALIZATION ##
    ###################
    # PCA plot
    all_dds_objects = []
    for file in args.dds_object:
        with open(file, "rb") as f:
            dds = pkl.load(f)
            all_dds_objects.append(dds)

    all_metadata = []
    all_normalized_counts = []
    for dds, team_id in zip(all_dds_objects, args.project_ids):
        metadata = dds.obs
        metadata["team_id"] = team_id
        all_metadata.append(metadata)
        normalized_counts = dds.layers["normed_counts"].T
        sample_ids = dds.obs.index
        sample_ids.name = None
        normalized_counts_df = pd.DataFrame(normalized_counts, index=dds.var.index, columns=sample_ids)
        all_normalized_counts.append(normalized_counts_df)
    combined_metadata = pd.concat(all_metadata)
    combined_normalized_counts = pd.concat(all_normalized_counts)

    # Impute missing values with the mean of each feature
    imputer = SimpleImputer(strategy="mean")
    combined_normalized_counts_imputed = imputer.fit_transform(combined_normalized_counts)

    # PCA expects samples as rows, genes as columns
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(combined_normalized_counts_imputed.T)
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
    pca_df["team_id"] = combined_metadata["team_id"].values
    pca_df["condition_id"] = combined_metadata["condition_id"].values
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x="PC1",
        y="PC2",
        hue="condition_id",
        style="team_id",
        data=pca_df,
        s=100,
        palette="Set2",
        alpha=0.7,
    )
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
    plt.title("PCA of Normalized Counts by Condition and Team")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.savefig(f"{args.cohort_id}.{args.salmon_mode}.pca_plot.png", dpi=300, bbox_inches="tight")
    plt.close("all")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Identify overlapping PyDESeq2 differentially expressed genes and generate plots for data visualization"
    )
    parser.add_argument(
        "-c",
        "--cohort-id",
        type=str,
        required=True,
        help="Cohort ID"
    )
    parser.add_argument(
        "-p",
        "--project-ids",
        type=str,
        nargs='+',
        required=True,
        help="Project IDs/team names"
    )
    parser.add_argument(
        "-n",
        "--n-teams",
        type=int,
        required=True,
        help="Number of teams"
    )
    parser.add_argument(
        "-g",
        "--degs",
        type=str,
        nargs='+',
        required=True,
        help="Table containing PyDESeq2 differentially expressed genes for each team"
    )
    parser.add_argument(
        "-d",
        "--dds-object",
        type=str,
        nargs='+',
        required=True,
        help="Pkl file containing the filtered DeSeqDataSet object"
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
