import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO


def DS1_append_deepsea_results(path_to_deepsea_tsv):
    all_results = []
    for arquivo in glob.glob(path_to_deepsea_tsv):
        df = pd.read_csv(arquivo, sep="\t")
        all_results.append(df)
    combined_df = pd.concat(all_results, ignore_index=True)

    combined_df["genome"] = combined_df["Name"].str.split("_").str[0]
    combined_df = combined_df.rename(columns={"Name": "id"})

    combined_df = combined_df.rename(
        columns={
            "Class": "class_deepsea",
            "Prob": "prob_deepsea",
            "genome": "genome_deepsea",
        }
    )

    return combined_df


##################################################################################


def DS2_append_faa(path_to_all_faa):
    faa_list = []
    path = os.path.join(path_to_all_faa, "*.faa")
    
    for file in glob.glob(path):
        for seq_record in SeqIO.parse(file, "fasta"):
            faa_list.append(
                {
                    "id": seq_record.id,
                    "sequence": str(seq_record.seq),
                    "description": seq_record.description,
                    "source_file": os.path.basename(file),
                }
            )

    faa_df = pd.DataFrame(faa_list)
    faa_df = faa_df.rename(
        columns={
            "description": "description_cds",
            "source_file": "source_file_cds",
        }
    )
    return faa_df


##################################################################################


def DS3_merge_deepsea_df(combined_df, faa_df):

    df_faa_with_deepsea = combined_df.merge(faa_df, on="id", how="inner")
    resistance_proteins = df_faa_with_deepsea[
        df_faa_with_deepsea["class_deepsea"] != "NonR"
    ]
    resistance_proteins = resistance_proteins[
        resistance_proteins["prob_deepsea"] > 0.7
    ]
    resistance_proteins["most_frequent_class"] = resistance_proteins.groupby(
        "source_file_cds"
    )["class_deepsea"].transform(lambda x: x.mode()[0])

    return resistance_proteins


##################################################################################


def DS4_resistance_class_distribution(resistance_proteins, output_dir="deepsea_images"):

    count = resistance_proteins["class_deepsea"].value_counts()
    bars = plt.bar(count.index, count.values, color="skyblue")
    for bar in bars:
        tamanho = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            tamanho,
            str(tamanho),
            ha="center",
            va="bottom",
        )
    plt.xlabel("Resistance Classes")
    plt.ylabel("Proteins Count")
    plt.title("Resistance Classes Distribution in Metagenome Proteins")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    output_path = os.path.join(output_dir, "resistance_classes_in_metagenome_proteins.png")

    plt.savefig(output_path, dpi=300)
    plt.close()


##################################################################################


def DS5_resistance_proteins_in_metagenomes(resistance_proteins, output_dir="deepsea_images"):

    table = resistance_proteins.pivot_table(
        index="genome_deepsea",
        columns="class_deepsea",
        values="id",
        aggfunc="count",
        fill_value=0,
    )

    total_per_class = table.sum(axis=0)
    total = total_per_class.sum()
    percentual_classes = (total_per_class / total * 100).round(1)

    table.plot(kind="bar", stacked=True, figsize=(12, 6), colormap="tab20")

    plt.xlabel("Metagenomes")
    plt.ylabel("Resistance Proteins Count")
    plt.legend(
        title="Resistance Classes", bbox_to_anchor=(1.05, 1), loc="upper left"
    )
    text_label = "\n".join(
        [
            f"{classe}: {value:.1f}%"
            for classe, value in percentual_classes.items()
        ]
    )
    plt.gcf().text(
        0.84,
        0.4,
        text_label,
        fontsize=10,
        va="center",
        bbox=dict(facecolor="white", edgecolor="black")
    )
    plt.tight_layout()

    output_path = os.path.join(output_dir, "resistance_proteins_in_metagenomes.png")

    plt.savefig(output_path, dpi=300)
    plt.close()


##################################################################################


def DS6_parse_antismash_gbks(input_directory):
    # input_directory: A folder containing all gbks from antiSMASH
    all_rows = []

    for filename in os.listdir(input_directory):
        if filename.endswith(".gbk"):
            file_path = os.path.join(input_directory, filename)

            try:
                for record in SeqIO.parse(file_path, "genbank"):

                    taxonomy = record.annotations.get("taxonomy", ["NaN"])
                    taxon_str = (
                        "; ".join(taxonomy)
                        if isinstance(taxonomy, list)
                        else str(taxonomy)
                    )

                    region_data = {
                        "product": "NaN",
                        "region_number": "NaN",
                        "taxon_lineage": taxon_str,
                    }

                    for feat in record.features:
                        if feat.type in ["region", "protocluster"]:
                            products = feat.qualifiers.get("product", ["NaN"])
                            region_data["product"] = (
                                ", ".join(products)                              
                                if isinstance(products, list)
                                else products
                            )

                            reg_num = feat.qualifiers.get(
                                "region_number",
                                feat.qualifiers.get(
                                    "protocluster_number", ["NaN"]
                                ),
                            )
                            region_data["region_number"] = (
                                reg_num[0]
                                if isinstance(reg_num, list)
                                else reg_num
                            )

                    for feature in record.features:
                        if feature.type == "CDS":

                            locus_tag = feature.qualifiers.get(
                                "locus_tag", ["NaN"]
                            )[0]
                            gene_kind = feature.qualifiers.get(
                                "gene_kind", ["NaN"]
                            )[0]
                            
                            # ---> 1. EXTRAI A SEQUÊNCIA DE AMINOÁCIDOS AQUI <---
                            translation = feature.qualifiers.get(
                                "translation", ["NaN"]
                            )[0]

                            gene_funcs = feature.qualifiers.get(
                                "gene_functions", ["NaN"]
                            )
                            gene_functions_str = (
                                "; ".join(gene_funcs)
                                if isinstance(gene_funcs, list)
                                else gene_funcs
                            )

                            row = {
                                "record_id": record.id,
                                "feature_type": feature.type,
                                "location": str(feature.location),
                                "product": region_data["product"],
                                "region_number": region_data["region_number"],
                                "locus_tag": locus_tag,
                                "gene_kind": gene_kind,
                                "gene_functions": gene_functions_str,
                                "sequence": translation,
                                "eggnog-coverage": "NaN",
                                "ipr-coverage": "NaN",
                                "taxon-lineage": region_data["taxon_lineage"],
                                "last-update": record.annotations.get(
                                    "date", "NaN"
                                ),
                            }
                            all_rows.append(row)

            except Exception as e:
                print(f"Error in file {filename}: {e}")

    return pd.DataFrame(all_rows)


##################################################################################


def DS7_filtering_antismash_gbks(df_bgc_features):
    df_bgc_features["product"] = df_bgc_features.groupby("record_id")[
        "product"
    ].transform(lambda x: x.ffill().bfill())
    df_bgc_features = df_bgc_features[df_bgc_features["feature_type"] == "CDS"]

    return df_bgc_features


##################################################################################


def DS8_merging_dataframes(df_bgc_features, resistance_proteins):

    df_bgc_resistance_proteins = df_bgc_features.merge(
        resistance_proteins, on="sequence", how="inner"
    )

    return df_bgc_resistance_proteins


##################################################################################


def DS9_download_resistance_proteins_info(df_bgc_resistance_proteins):
    df_bgc_resistance_proteins = df_bgc_resistance_proteins.dropna(
        axis=1, how="all"
    )
    df_bgc_resistance_proteins.to_csv(
        "resistance_proteins_in_bgc.csv", index=False
    )

    return df_bgc_resistance_proteins


##################################################################################
def DS10_deepsea_heatmap(df_bgc_resistance_proteins, output_dir="deepsea_images"):
    df_filtered = df_bgc_resistance_proteins[
        [
            "record_id",
            "genome_deepsea",
            "product",
            "class_deepsea",
            "gene_kind",
        ]
    ]

    table = (
        df_filtered.groupby(["product", "class_deepsea"])
        .size()
        .unstack(fill_value=0)
    )

    plt.figure(figsize=(10, 6))
    sns.heatmap(table, annot=True, cmap="YlGnBu")
    plt.title("Resistance Classes in BGCs heatmap")
    plt.xlabel("Resistance Classes")
    plt.ylabel("BGC Class")
    plt.tight_layout()

    output_path = os.path.join(output_dir, "resistance_proteins_heatmap.png")

    plt.savefig(output_path, dpi=300)
    plt.close()
