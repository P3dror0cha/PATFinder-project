#!/usr/bin/env python3

import pandas as pd

def K5_merging_kofam_df(final_df_filtered, df_list, df_kegg):
    """
    Merges KOfam annotations with KEGG pathway information and summarizes results per BGC.

    Parameters
    ----------
    final_df_filtered : pandas.DataFrame
        Processed KOfam DataFrame (e.g., output from K3_kofam_without_best_hits),
        containing at least 'KO', 'BGC_ID', and 'percentual' columns.

    df_list : pandas.DataFrame
        DataFrame with KEGG pathway descriptions, containing 'pathway' and
        'description' columns.

    df_kegg : pandas.DataFrame
        DataFrame mapping KO identifiers to KEGG pathways, containing 'KO'
        and 'pathway' columns.

    Returns
    -------
    pandas.DataFrame
        Aggregated DataFrame grouped by 'BGC_ID', including:
        - pathway : concatenated unique pathway IDs associated with the BGC
        - description : concatenated unique pathway descriptions
        - KO : concatenated unique KO identifiers
        - percentual : percentage of best hits per BGC (carried over)
    """
    df_merge = final_df_filtered.merge(df_kegg, on="KO", how="left")
    df_merge = df_merge.loc[df_merge["pathway"].str.startswith("map")]
    df_merge = df_merge.merge(df_list, on="pathway", how="left")
    
    df_groupby = df_merge.groupby("BGC_ID").agg({
        "pathway": lambda x: ";".join(x.dropna().unique()),
        "description": lambda x: ";".join(x.dropna().unique()),
        "KO": lambda x: ";".join(x.dropna().unique()),
        "percentual": "first"
    }).reset_index()
    
    return df_groupby