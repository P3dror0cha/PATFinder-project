import pandas as pd
import re
import os
from Bio import SeqIO

def product_from_gbk_files(gbk_files):
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

    for file_path in gbk_files:
            filename = os.path.basename(file_path) 

            try:
                for record in SeqIO.parse(file_path, "genbank"):
                    for feat in record.features:
                        
                        if feat.type in ["region", "protocluster"]:
                            product_str = ", ".join(feat.qualifiers.get("product", ["NaN"]))
                            
                            all_rows.append({
                                "record_id": os.path.splitext(filename)[0].replace(".region", "_"),
                                "query_bgc_class": product_str
                            })

            except Exception as e:
                print(f"Error in file {filename}: {e}")

    return pd.DataFrame(all_rows)

def uniting_poem_deepsea_results(df_bgc_correlation, poem_result, deepsea_result):
    '''
    Merges POEM and DeepSea results into a single DataFrame.

    This function standardizes the identification column names in both DataFrames
    to 'bgc_id' (renaming 'metagenome_id' from POEM and 'record_id' from DeepSea) 
    and performs an outer merge based on this column.

    Args:
        poem_result (pd.DataFrame): DataFrame containing the POEM results.
        deepsea_result (pd.DataFrame): DataFrame containing the DeepSea results.

    Returns:
        pd.DataFrame: A unified DataFrame containing both results.
    '''
    df_bgc_correlation = df_bgc_correlation.rename(columns={"padronized_name": "bgc_id"})
    poem_result = poem_result.rename(columns={"metagenome_id": "bgc_id"})
    deepsea_result = deepsea_result.rename(columns={"record_id": "bgc_id"})

    df_united = pd.merge(
        poem_result, 
        deepsea_result, 
        on="bgc_id", 
        how="outer"
        )

    df_united = pd.merge(
        df_united,
        df_bgc_correlation,
        on="bgc_id",
        how="outer"
    )

    return df_united

def uniting_kofam_diamond_results(kofam_results, diamond_results):
    '''
    Merges KOfam and Diamond results into a single DataFrame.

    Note: The KOfam table also has the BIG-SCAPE information.
    This function standardizes the gene identification columns to 'gene_id' 
    (renaming 'gene name' from KOfam and 'qseqid' from Diamond) and performs 
    an outer merge based on this column.

    Args:
        kofam_results (pd.DataFrame): DataFrame containing the KOfam results.
        diamond_results (pd.DataFrame): DataFrame containing the Diamond results.

    Returns:
        pd.DataFrame: A unified DataFrame containing the gene annotations.
    '''

    kofam_results = kofam_results.rename(columns={"BGC_ID": "bgc_id"})
    #diamond_results["bgc_id"] = diamond_results["id"].str.extract(r'(sample_\d+_\d+)')

    df_united = pd.merge(
        kofam_results, 
        diamond_results, 
        on="bgc_id", 
        how="outer"
        )
    
    return df_united

def uniting_all_infos(poem_deepsea, kofam_diamond):
    '''
    Consolidates BGC information (POEM/DeepSea) with gene annotations (KOfam/Diamond).

    This function creates a 'bgc_id' column in the gene DataFrame by extracting 
    the prefix from the 'gene_id' column (capturing all text before the first '|' 
    character via Regular Expression) and strips any trailing whitespaces. Finally, 
    it performs an outer merge with the BGC DataFrame.

    Args:
        poem_deepsea (pd.DataFrame): Consolidated DataFrame of POEM and DeepSea results.
        kofam_diamond (pd.DataFrame): Consolidated DataFrame of KOfam and Diamond results.

    Returns:
        pd.DataFrame: The final DataFrame integrating all analyzed information.
    '''

    kofam_diamond['bgc_id'] = kofam_diamond['bgc_id'].str.strip()

    df_united = pd.merge(
        poem_deepsea, 
        kofam_diamond, 
        how="outer", 
        on="bgc_id"
        )

    return df_united

def sorting_concat_columns(df_united):
    '''
    Consolidates BGC information (POEM/DeepSea) with gene annotations (KOfam/Diamond).

    This function creates a 'bgc_id' column in the gene DataFrame by extracting 
    the prefix from the 'gene_id' column (capturing all text before the first '|' 
    character via Regular Expression) and strips any trailing whitespaces. Finally, 
    it performs an outer merge with the BGC DataFrame.

    Args:
        df_united (pd.DataFrame): Output from the function "uniting_all_infos".

    Returns:
        pd.DataFrame: The final DataFrame with final results filtered.
    '''


    df_united = df_united.rename(columns={
        "GBK_b": "MIBIG_bgc", 
        "weights": "BIG_SCAPE_weights",
        "class": "MIBIG_class",
        "compound_name": "MIBIG_compound_name"})

    columns_order = [
        "bgc_id",
        "original_name",
        "query_bgc_class",
        "MIBIG_bgc",
        "MIBIG_class",
        "MIBIG_compound_name",
        "distance",
        "jaccard",
        "adjacency",
        "dss",
        "BIG_SCAPE_weights",
        "deepsea_class",
        "deepsea_prob",
        "deepsea_hits_sequence",
        "deepsea_cds_description",
        "deepsea_hits_count",
        "uniprot_id",
        "uniprot_evalue",
        "uniprot_proteins",
        "uniprot_genes",
        "uniprot_hits_count",
        "KO",
        "filtered_pathways",
        "strand",
        "coordinates",
        "amrfinder_element_symbol",
        "amrfinder_class",
        "amrfinder_subclass",
        "amrfinder_hits_count"
        ]
    
    df_united = df_united[columns_order]

    

    df_united = df_united.drop_duplicates()

    return df_united

def create_bgc_class_correlation(df_bgc_product, ids_correlation_path):
    ids_corr = pd.read_csv(
        ids_correlation_path,
        sep="\t",
        dtype=str
    )

    ids_corr.columns = (
        ids_corr.columns
        .str.replace("\ufeff", "", regex=False)
        .str.strip()
    )

    sample_to_original = dict(
        zip(
            ids_corr["Padronized_ID"].str.strip(),
            ids_corr["Original_ID"].str.strip()
        )
    )

    sorted_samples = sorted(
        sample_to_original.keys(),
        key=len,
        reverse=True
    )

    columns = [
        "original_name",
        "padronized_name",
        "query_bgc_class"
    ]
    rows = []

    for _, row in df_bgc_product.iterrows():
        record_id = str(row["record_id"]).strip()

        matched_sample = next(
            (
                sample_id
                for sample_id in sorted_samples
                if record_id == sample_id
                or record_id.startswith(f"{sample_id}_")
            ),
            None
        )

        suffix = record_id[len(matched_sample):]

        rows.append({
            "original_name": sample_to_original[matched_sample] + suffix,
            "padronized_name": record_id,
            "query_bgc_class": row["query_bgc_class"]
        })

    return pd.DataFrame(rows, columns=columns)


def uniting_amrfinder_results(df_concat, df_amrfinder):
    df_united = pd.merge(
        df_concat,
        df_amrfinder,
        how="outer",
        on="bgc_id"
    )
    return df_united