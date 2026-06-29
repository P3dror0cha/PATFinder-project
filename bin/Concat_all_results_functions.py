import pandas as pd
import re

def uniting_poem_deepsea_results(poem_result, deepsea_result):
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

    poem_result = poem_result.rename(columns={"metagenome_id": "bgc_id"})
    deepsea_result = deepsea_result.rename(columns={"record_id": "bgc_id"})

    df_united = pd.merge(
        poem_result, 
        deepsea_result, 
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
    diamond_results["bgc_id"] = diamond_results["id"].str.extract(r'(sample_\d+_\d+)')

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