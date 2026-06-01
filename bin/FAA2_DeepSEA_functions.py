import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO


def DS1_append_deepsea_results(path_to_deepsea_tsv):
    """
    Reads and concatenates multiple DeepSEA TSV output files into a single DataFrame.

    Extracts the genome name from the 'Name' column and standardizes column names 
    for downstream processing.

    Args:
        path_to_deepsea_tsv (str): Search pattern for DeepSEA TSV files (e.g., "data/*.tsv").

    Returns:
        pd.DataFrame: A combined DataFrame containing all DeepSEA results with 
    """
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
    """
    Parses all FASTA amino acid (.faa) files in a directory and compiles them.
    The .faa files must have the aminoacid sequence for all proteins in the query
    (meta)genomes. 

    Extracts the sequence ID, the raw amino acid sequence, description, and the 
    source file name for every record found.

    Args:
        path_to_all_faa (str): Path to the directory containing the .faa files.

    Returns:
        pd.DataFrame: A DataFrame with sequences and metadata ('id', 'sequence', 
            'description_cds', 'source_file_cds').
    """
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
    """
    Merges DeepSEA predictions with sequence data (from .faa files) and filters 
    for resistance proteins.

    Performs an inner join on sequence IDs, removes non-resistance ('NonR') 
    classifications, and filters for predictions with high probability (> 0.7). 

    Args:
        combined_df (pd.DataFrame): DataFrame containing DeepSEA results (from DS1).
        faa_df (pd.DataFrame): DataFrame containing parsed FASTA sequences (from DS2).

    Returns:
        pd.DataFrame: A filtered DataFrame containing only high-probability 
            antibiotic resistance proteins.
    """

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
    """
    Generates and saves a bar chart showing the distribution of resistance classes.

    Args:
        resistance_proteins (pd.DataFrame): Filtered DataFrame of resistance proteins.
        output_dir (str, optional): Directory to save the plot. Defaults to "deepsea_images".

    Returns:
        None
        
    Saves:
        resistance_classes_in_metagenome_proteins.png in the specified output directory.
    """
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
    """
    Generates and saves a stacked bar chart of resistance proteins per metagenome.

    Displays absolute counts and annotates the overall percentage of each 
    resistance class across the entire dataset.

    Args:
        resistance_proteins (pd.DataFrame): Filtered DataFrame of resistance proteins.
        output_dir (str, optional): Directory to save the plot. Defaults to "deepsea_images".

    Returns:
        None
        
    Saves:
        resistance_proteins_in_metagenomes.png in the specified output directory.
    """

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
    """
    Parses antiSMASH GenBank (.gbk) files to extract BGC and CDS features.

    Iterates over all .gbk files in the directory, capturing global taxonomy, 
    identifying BGC regions/protoclusters, and extracting detailed annotations 
    (sequence, locus tags, product types) for all Coding Sequences (CDS).

    Args:
        input_directory (str): Path to the folder containing antiSMASH .gbk files.

    Returns:
        pd.DataFrame: A DataFrame where each row represents a parsed feature 
            (mostly CDS) mapped to its broader BGC region context.
    """

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
    """
    Filters antiSMASH features to retain only Coding Sequences (CDS).

    Propagates BGC 'product' classifications across missing values within the same 
    record ID using forward and backward filling, ensuring all CDS within a region 
    inherit the region's BGC type.

    Args:
        df_bgc_features (pd.DataFrame): Raw parsed antiSMASH DataFrame (from DS6).

    Returns:
        pd.DataFrame: A filtered DataFrame containing exclusively CDS features.
    """

    df_bgc_features["product"] = df_bgc_features.groupby("record_id")[
        "product"
    ].transform(lambda x: x.ffill().bfill())
    df_bgc_features = df_bgc_features[df_bgc_features["feature_type"] == "CDS"]

    return df_bgc_features


##################################################################################


def DS8_merging_dataframes(df_bgc_features, resistance_proteins):
    """
    Cross-references antiSMASH BGC annotations with DeepSEA resistance predictions.

    Performs an inner join based on the exact amino acid 'sequence' to identify 
    which predicted resistance proteins are physically located within 
    Biosynthetic Gene Clusters (BGCs).

    Args:
        df_bgc_features (pd.DataFrame): Filtered antiSMASH CDS DataFrame (from DS7).
        resistance_proteins (pd.DataFrame): Filtered DeepSEA DataFrame (from DS3).

    Returns:
        pd.DataFrame: A merged DataFrame combining BGC genomic context with 
            antibiotic resistance classifications.
    """

    df_bgc_resistance_proteins = df_bgc_features.merge(
        resistance_proteins, on="sequence", how="inner"
    )

    return df_bgc_resistance_proteins


##################################################################################


def DS9_download_resistance_proteins_info(df_bgc_resistance_proteins):
    """
    Cleans up the final joined DataFrame and exports it to a CSV file.

    Drops any columns that are entirely composed of missing values (NaN) 
    before exporting.

    Args:
        df_bgc_resistance_proteins (pd.DataFrame): Merged BGC-Resistance DataFrame.

    Returns:
        pd.DataFrame: The cleaned DataFrame.
        
    Saves:
        resistance_proteins_in_bgc.csv in the current working directory.
    """
    
    df_bgc_resistance_proteins = df_bgc_resistance_proteins.dropna(
        axis=1, how="all"
    )
    df_bgc_resistance_proteins.to_csv(
        "resistance_proteins_in_bgc.csv", index=False
    )

    return df_bgc_resistance_proteins


##################################################################################
def DS10_deepsea_heatmap(df_bgc_resistance_proteins, output_dir="deepsea_images"):
    """
    Generates and saves a heatmap correlating BGC product types with Resistance Classes.

    Creates a 2D matrix counting the co-occurrences of specific secondary metabolite 
    classes (BGC product) and specific antibiotic resistance categories.

    Args:
        df_bgc_resistance_proteins (pd.DataFrame): Merged BGC-Resistance DataFrame.
        output_dir (str, optional): Directory to save the plot. Defaults to "deepsea_images".

    Returns:
        None
        
    Saves:
        resistance_proteins_heatmap.png in the specified output directory.
    """
    
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
