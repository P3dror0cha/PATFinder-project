#!/usr/bin/env python3

import pandas as pd
import requests
import re

def K1_kofam_import(kofam_result_list):
    """
    Reads a KOfam result file and converts it into a structured pandas DataFrame.

    Parameters
    ----------
    kofam_result_list : str
        Path to the KOfam output file to be parsed.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing parsed KOfam results with the following columns:
        - best_hits : int (1 if the hit is marked with '*', otherwise 0)
        - query : str (query sequence ID)
        - KO : str (KEGG Orthology identifier)
        - threshold : float or str (score threshold used)
        - score : float or str (alignment score)
        - evalue : float or str (expectation value)
        - description : str (KO description)
        - function : str (functional annotation)
    """
    data = []

    with open(kofam_result_list) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split(maxsplit=6)
            
            if parts[0] == '*':
                flag = 1
                parts = parts[1:]  
            else:
                flag = 0
            
            while len(parts) < 6:
                parts.append(None)
            
            data.append([flag] + parts)

    df = pd.DataFrame(
        data,
        columns=['best_hits','query', 'KO', 'threshold', 'score', 'evalue', 'description', 'function']
    )
    df.to_csv("/home/pedro/PATFinder-project/results_backup/1_df_do_kofam.csv")
    return df

##############################################################################################

def K2_kofam_best_hits(df):
    """
    Filters and extracts high-confidence best hits from a KOfam DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame generated from KOfam results, expected to contain at least
        the columns: 'best_hits', 'evalue', 'score', 'threshold', and 'function'.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame containing only best hits (best_hits == 1) with
        e-value <= 1e-5. Unnecessary columns ('best_hits', 'threshold',
        'score', 'function') are removed.
    """
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["threshold"] = pd.to_numeric(df["threshold"], errors="coerce")

    df = df[df["evalue"] <= 1.0e-5]

    df_best_hits = df[df["best_hits"] == 1].copy()

    df_best_hits.drop(columns=["best_hits", "threshold", "score", "function"], inplace=True)

    df_best_hits.to_csv("/home/pedro/PATFinder-project/results_backup/2_df_best_hits_do_kofam.csv")
    return df_best_hits

##############################################################################################

def K3_kofam_without_best_hits(df, df_best_hits):
    """
    Complements best-hit annotations by selecting representative non-best hits
    and computing quality metrics per BGC.

    Parameters
    ----------
    df : pandas.DataFrame
        Full KOfam results DataFrame containing all hits. Must include columns
        such as 'query', 'best_hits', 'evalue', 'threshold', 'score',
        'description', and 'function'.

    df_best_hits : pandas.DataFrame
        DataFrame containing only high-confidence best hits (e.g., output from
        K2_kofam_best_hits), including at least the 'query' column.

    Returns
    -------
    pandas.DataFrame
        Final DataFrame combining best hits and selected non-best hits with:
        - quality flag (1 for best hits, 0 for others)
        - BGC_ID extracted from query 
        - gene_id derived from query
        - percentual: percentage of best hits per BGC

        The output excludes columns 'query', 'evalue', and 'description'.
    """
    df_without_best_hits = df[~df["query"].isin(df_best_hits["query"])]
    df_without_best_hits["query"].nunique()
    df_without_best_hits = (
    df_without_best_hits
    .sort_values(["best_hits", "evalue"], ascending=[False, True])
    .drop_duplicates("query")
    )
    
    df_without_best_hits.drop(columns={"best_hits", "threshold", "score", "description"}, inplace=True)
    df_without_best_hits.rename(columns={"function": "description", "('query', None)": "query"}, inplace=True)
    
    df_best_hits["quality"] = 1
    df_without_best_hits["quality"] = 0
    final_df = pd.concat([df_best_hits, df_without_best_hits], ignore_index=True)
    final_df = final_df[~final_df["query"].isna() & (final_df["query"] != "NaN")]

    final_df["BGC_ID"] = final_df["query"].str.extract(r"(sample_\d+_\d+)") \
        .fillna("").astype(str).agg("_".join, axis=1)
    final_df["gene_id"] = final_df["query"].apply(
        lambda x: "_".join(x.split("_")[1:]) if isinstance(x, str) else None
    )

    final_df["percentual"] = final_df.groupby("BGC_ID")["quality"].transform("sum") / final_df.groupby("BGC_ID")["BGC_ID"].transform("size") * 100

    final_df_filtered = final_df.drop(columns=["query", "evalue", "description"])
    final_df_filtered.to_csv("/home/pedro/PATFinder-project/results_backup/3_final_df_filtered.csv", index=False)
    return final_df_filtered

##############################################################################################

def K4_KEGG_API_information(url="https://rest.kegg.jp/link/pathway/ko", url_2="https://rest.kegg.jp/list/pathway"):
    """
    Retrieves KEGG pathway mapping and pathway descriptions via the KEGG REST API.

    Parameters
    ----------
    url : str, optional
        KEGG API endpoint that links KO identifiers to pathways.
        Default is "https://rest.kegg.jp/link/pathway/ko".

    url_2 : str, optional
        KEGG API endpoint that provides pathway IDs and their descriptions.
        Default is "https://rest.kegg.jp/list/pathway".

    Returns
    -------
    tuple of pandas.DataFrame
        - df_kegg : DataFrame mapping KO identifiers to pathway IDs with columns:
            'KO' (KEGG Orthology ID) and 'pathway' (pathway ID).
        - df_list : DataFrame containing pathway metadata with columns:
            'pathway' (pathway ID) and 'description' (pathway name/description).
    """
    r = requests.get(url)

    lines = r.text.strip().split("\n")

    data = [line.split("\t") for line in lines if "\t" in line]

    df_kegg = pd.DataFrame(data, columns=["KO", "pathway"])

    df_kegg["KO"] = df_kegg["KO"].str.replace("ko:", "")
    df_kegg["pathway"] = df_kegg["pathway"].str.replace("path:", "")
    
    r_2 = requests.get(url_2)

    lines = r_2.text.strip().split("\n")

    data = [line.split("\t") for line in lines if "\t" in line]

    df_list = pd.DataFrame(data, columns=["pathway", "description"])
    df_kegg.to_csv("/home/pedro/PATFinder-project/results_backup/df_kegg.csv")
    df_list.to_csv("/home/pedro/PATFinder-project/results_backup/df_list.csv")
    return df_kegg, df_list

##############################################################################################

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
    df_merge = df_merge.loc[df_merge["pathway"].str.startswith("map", na=False)]
    df_merge = df_merge.merge(df_list, on="pathway", how="left")
    df_merge.to_csv("/home/pedro/PATFinder-project/results_backup/3.5_df_merge.csv")

    df_groupby = df_merge.groupby("BGC_ID").agg({
        "pathway": lambda x: ";".join(x.dropna().unique()),
        "description": lambda x: ";".join(x.dropna().unique()),
        "KO": lambda x: ";".join(x.dropna().unique()),
        "percentual": "first"
    }).reset_index()
    df_groupby.to_csv("/home/pedro/PATFinder-project/results_backup/4_df_groupby.csv", index=False)
    return df_groupby

##############################################################################################

def K6_uniting_BGC_information(BIGSCAPE_result, df_groupby):
    """
    Integrates antiSMASH BGC annotations with KEGG/KOfam-derived pathway information.

    Parameters
    ----------
    antiSMASH_result : str
        Path to the antiSMASH output file (CSV format) containing BGC annotations,
        including a 'BGC_ID' column.

    df_groupby : pandas.DataFrame
        DataFrame summarizing KEGG/KOfam information per BGC (e.g., output from
        K5_merging_kofam_df), containing 'BGC_ID' and associated annotation columns.

    Returns
    -------
    pandas.DataFrame
        Merged DataFrame combining antiSMASH annotations with KEGG pathway,
        KO, and percentual information for each BGC.
    """  
    df_bgc = pd.read_csv(BIGSCAPE_result)
    df_bgc.to_csv("/home/pedro/PATFinder-project/results_backup/df_bgc_antes_do_merge.csv", index=False)
    df_groupby.to_csv("/home/pedro/PATFinder-project/results_backup/df_groupby.csv", index=False)
    df_bgc = pd.merge(df_bgc, df_groupby, on="BGC_ID", how="left")

    return df_bgc

##############################################################################################

def K7_kegg_ids_filtering(x):
    """
    Filters KEGG pathway IDs to retain only those associated with antibiotic
    biosynthesis or related pathways.

    Parameters
    ----------
    x : str or NaN
        String containing one or more KEGG pathway IDs (e.g., "map00940;map00945"),
        or NaN.

    Returns
    -------
    str or None
        A semicolon-separated string of filtered KEGG pathway IDs that match
        the predefined set of interest (e.g., antibiotic biosynthesis pathways).
        Returns None if no matching IDs are found or if input is NaN.
    """

    antibiotic_biosynthesis_ids = [    
    "00940", "00945", "00941", "00944", "00942", "00943",
    "00946", "00901", "00403", "00950", "00960", "00996",
    "00232", "00965", "00966", "00402", "00311", "00332",
    "00261", "00331", "00521", "00524", "00525", "00401",
    "00404", "00405", "00333", "00254", "00975", "00998",
    "00999", "00997", "01501", "01502", "01503", "00522", 
    "01051", "01059", "00253", "01053", "01055"]

    kegg_set = set(antibiotic_biosynthesis_ids)
    
    if pd.isna(x):
        return None
    
    ids = re.findall(r"map(\d{5})", x)
    
    filtered_ids = [f"map{id_}" for id_ in ids if id_ in kegg_set]

    return ";".join(filtered_ids) if filtered_ids else None

##############################################################################################

def K8_uniting_all_information(df_bgc, output_path_raw_csv, output_path_filtered_csv):
    """
    Finalizes BGC annotation data by filtering KEGG pathways, renaming columns,
    and exporting both raw and simplified results to CSV files.

    Parameters
    ----------
    df_bgc : pandas.DataFrame
        DataFrame containing merged BGC, KOfam, and KEGG information, including
        a 'pathway' column and 'percentual' metric.

    output_path_raw_csv : str
        File path to save the complete (raw) annotated DataFrame.

    output_path_filtered_csv : str
        File path to save the filtered DataFrame with reduced columns.

    Returns
    -------
    tuple of pandas.DataFrame
        - df_bgc : Full DataFrame with added 'filtered_pathways' column and
          renamed 'percentual' to 'KO_quality_percentual'.
        - df_bgc_filtered : Simplified DataFrame with selected columns removed.
    """
    df_bgc_filtered = df_bgc.drop(columns=["pathway", "description", "KO_quality_percentual", "strand"])
    df_bgc.to_csv(output_path_raw_csv, index=False)
    df_bgc_filtered.to_csv(output_path_filtered_csv, index=False)
    
    return df_bgc, df_bgc_filtered