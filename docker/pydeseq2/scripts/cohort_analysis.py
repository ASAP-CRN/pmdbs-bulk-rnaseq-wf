import argparse
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA


def modify_tables(file_list, team_ids):
    dfs = []
    for file, team_id in zip(file_list, team_ids):
        df = pd.read_csv(file)
        df["team_id"] = team_id
        dfs.append(df)
    combined_dfs = pd.concat(dfs, ignore_index=True)
    return combined_dfs


def main(args):
    ##########################################################
    ## OVERLAPPING DEGS ONLY FOR CROSS TEAM COHORT ANALYSIS ##
    ##########################################################
    if args.n_teams > 1:
        combined_degs = modify_tables(args.degs, args.project_ids)
        combined_degs.set_index(combined_degs.columns[0], inplace=True)
        grouped = combined_degs.groupby("team_id").apply(lambda x: set(x.index))
        common_degs = set.intersection(*grouped)
        common_degs_df = combined_degs.loc[combined_degs.index[combined_degs.index.isin(common_degs)]]
        common_degs_df.to_csv(f"{args.cohort_id}.{args.salmon_mode}.overlapping_significant_genes.csv")


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

    # PCA expects samples as rows, genes as columns
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(combined_normalized_counts.T)
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
    pca_df["team_id"] = combined_metadata["team_id"].values
    pca_df["intervention_id"] = combined_metadata["intervention_id"].values
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x="PC1",
        y="PC2",
        hue="intervention_id",
        style="team_id",
        data=pca_df,
        s=100,
        palette="Set2",
        alpha=0.7,
    )
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
    plt.title("PCA of Normalized Counts by Intervention and Team")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.savefig(f"{args.cohort_id}.{args.salmon_mode}.intervention_pca_plot.png", dpi=300, bbox_inches="tight")
    plt.close("all")

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
    plt.savefig(f"{args.cohort_id}.{args.salmon_mode}.condition_pca_plot.png", dpi=300, bbox_inches="tight")
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
