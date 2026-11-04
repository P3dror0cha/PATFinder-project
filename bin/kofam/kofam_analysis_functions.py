#!/usr/bin/env python3

import pandas as pd
import requests
import re

def K1_kofam_import(kofam_result_list):
    """
    Reads a KOfam result file and converts it into a structured pandas DataFrame.

    Parses the raw text output from KOfam, identifying best hits marked with an 
    asterisk ('*') and standardizing the columns for downstream analysis.

    Args:
        kofam_result_list (str): Path to the KOfam output text file to be parsed.

    Returns:
        pd.DataFrame: A DataFrame containing parsed KOfam results with columns:
            'best_hits', 'query', 'KO', 'threshold', 'score', 'evalue', 
            'description', and 'function'.
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
    #df.to_csv("/home/pedro/PATFinder-project/results_backup/1_df_do_kofam.csv")
    return df

##############################################################################################

def K2_kofam_best_hits(df):
    """
    Filters and extracts high-confidence best hits from a KOfam DataFrame.

    Converts relevant columns to numeric types and filters for records marked 
    as best hits with an e-value <= 1e-5. Drops intermediate columns to streamline 
    the dataset.

    Args:
        df (pd.DataFrame): DataFrame generated from KOfam results (output of K1).

    Returns:
        pd.DataFrame: A filtered DataFrame containing only high-confidence best hits.
    """

    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["threshold"] = pd.to_numeric(df["threshold"], errors="coerce")

    df = df[df["evalue"] <= 1.0e-5]

    df_best_hits = df[df["best_hits"] == 1].copy()

    df_best_hits.drop(columns=["best_hits", "threshold", "score", "function"], inplace=True)

    #df_best_hits.to_csv("/home/pedro/PATFinder-project/results_backup/2_df_best_hits_do_kofam.csv")
    return df_best_hits

##############################################################################################

def K3_kofam_without_best_hits(df, df_best_hits):
    """
    Complements best-hit annotations by selecting non-best hits and computing metrics.

    Identifies query sequences missing from the best hits DataFrame, selects their 
    most significant alignments, and merges them back with the best hits. It also 
    calculates a 'percentual' metric representing the percentage of best hits per BGC.

    Args:
        df (pd.DataFrame): Full KOfam results DataFrame containing all hits.
        df_best_hits (pd.DataFrame): DataFrame containing only high-confidence best hits.

    Returns:
        pd.DataFrame: Final combined DataFrame with a 'quality' flag (1 for best, 
            0 for others), extracted 'BGC_ID', and calculated 'percentual'.
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

    final_df["metagenome_ID"] = final_df["query"].str.extract(r"([a-zA-Z]+_\d+)") \
        .fillna("").astype(str).agg("_".join, axis=1)
    final_df["BGC_ID"] = final_df["query"].str.extract(r"^([a-zA-Z]+_\d+(?:_\d+)*)") \
        .fillna("").astype(str).agg("_".join, axis=1)

    final_df["percentual"] = final_df.groupby("query")["quality"].transform("mean") * 100

    final_df_filtered = final_df.drop(columns=["evalue", "description"])
    #final_df_filtered.to_csv("/home/pedro/PATFinder-project/results_backup/3_final_df_filtered.csv", index=False)
    return final_df_filtered

##############################################################################################

def K4_KEGG_API_information(url="https://rest.kegg.jp/link/pathway/ko", url_2="https://rest.kegg.jp/list/pathway"):
    """
    Retrieves KEGG pathway mapping and descriptions via the KEGG REST API.

    Fetches data from KEGG endpoints to map KEGG Orthology (KO) identifiers 
    to their respective pathway IDs, and retrieves the descriptive names for 
    each pathway.

    Args:
        url (str, optional): KEGG API endpoint linking KOs to pathways. 
            Defaults to "https://rest.kegg.jp/link/pathway/ko".
        url_2 (str, optional): KEGG API endpoint for pathway descriptions. 
            Defaults to "https://rest.kegg.jp/list/pathway".

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames:
            - df_kegg: Maps KO IDs to pathway IDs.
            - df_list: Maps pathway IDs to pathway descriptions.
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
    #df_kegg.to_csv("/home/pedro/PATFinder-project/results_backup/df_kegg.csv")
    #df_list.to_csv("/home/pedro/PATFinder-project/results_backup/df_list.csv")
    return df_kegg, df_list

##############################################################################################

def K5_merging_kofam_df(final_df_filtered, df_list, df_kegg):
    """
    Merges KOfam annotations with KEGG pathway information and summarizes per BGC.

    Joins the filtered KOfam results with KEGG API data based on KO and pathway IDs.
    Aggregates the data by BGC ID, concatenating unique pathways, descriptions, 
    and KOs into single strings per cluster.

    Args:
        final_df_filtered (pd.DataFrame): Processed KOfam DataFrame (output of K3).
        df_list (pd.DataFrame): KEGG pathway descriptions (output of K4).
        df_kegg (pd.DataFrame): KO to pathway mapping (output of K4).

    Returns:
        pd.DataFrame: An aggregated DataFrame grouped by 'BGC_ID'.
    """
    df_merge = final_df_filtered.merge(df_kegg, on="KO", how="left")
    df_merge = df_merge.loc[df_merge["pathway"].str.startswith("map", na=False)]
    df_merge = df_merge.merge(df_list, on="pathway", how="left")
    #df_merge.to_csv("/home/pedro/PATFinder-project/results_backup/3.5_df_merge.csv")

    df_groupby = df_merge.groupby("BGC_ID").agg({
        "pathway": lambda x: ";".join(x.dropna().unique()),
        "description": lambda x: ";".join(x.dropna().unique()),
        "KO": lambda x: ";".join(x.dropna().unique()),
        "percentual": "first"
    }).reset_index()
    #df_groupby.to_csv("/home/pedro/PATFinder-project/results_backup/4_df_groupby.csv", index=False)
    return df_groupby

##############################################################################################

def K6_uniting_BGC_information(BIGSCAPE_result, df_groupby):
    """
    Integrates antiSMASH/BIG-SCAPE BGC annotations with KEGG pathway data.

    Reads a BIG-SCAPE TSV output file, standardizes the BGC ID column, and 
    performs a left join to attach the aggregated KEGG and KOfam metrics.

    Args:
        BIGSCAPE_result (str): Path to the BIG-SCAPE output file (TSV/CSV format).
        df_groupby (pd.DataFrame): DataFrame summarizing KEGG information per BGC (output of K5).

    Returns:
        pd.DataFrame: A merged DataFrame combining structural BGC annotations 
            with functional KEGG annotations.
    """

    df_bgc = pd.read_csv(BIGSCAPE_result, sep='\t')
    df_bgc.columns = df_bgc.columns.str.strip()
    df_bgc.rename(columns={"GBK_a": "BGC_ID"}, inplace=True) 
    #df_bgc.to_csv("/home/pedro/PATFinder-project/results_backup/df_bgc_antes_do_merge.csv", index=False)
    #df_groupby.to_csv("/home/pedro/PATFinder-project/results_backup/df_groupby.csv", index=False)
    df_bgc = pd.merge(df_bgc, df_groupby, on="BGC_ID", how="left")

    return df_bgc

##############################################################################################

def K7_kegg_ids_filtering(x):
    """
    Filters KEGG pathway IDs to retain only those related to antibiotic biosynthesis.

    Scans a string of concatenated KEGG pathway IDs and cross-references them 
    against a predefined, hardcoded list of pathways associated with antibiotic 
    and secondary metabolite biosynthesis.

    Args:
        x (str | float): A string of concatenated KEGG pathway IDs (e.g., "map00940;map00945"), 
            or NaN.

    Returns:
        str | None: A semicolon-separated string of filtered KEGG pathway IDs, 
            or None if the input is NaN or no matches are found.
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
    Finalizes BGC annotation data by exporting raw and simplified DataFrames.

    Creates a filtered version of the final BGC DataFrame by dropping detailed 
    pathway and quality score columns, and exports both the raw and filtered 
    datasets to specified CSV paths.

    Args:
        df_bgc (pd.DataFrame): Merged BGC and KEGG DataFrame (output of K6/K7 processing).
        output_path_raw_csv (str): File path to save the complete (raw) DataFrame.
        output_path_filtered_csv (str): File path to save the simplified DataFrame.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing:
            - df_bgc: The original, complete DataFrame.
            - df_bgc_filtered: A simplified DataFrame with select columns removed.
    """

    df_bgc_filtered = df_bgc.drop(columns=["pathway", "description", "KO_quality_percentual"])
    df_bgc.to_csv(output_path_raw_csv, index=False)
    df_bgc_filtered.to_csv(output_path_filtered_csv, index=False)
    #df_bgc_filtered.to_csv("/home/pedro/PATFinder-project/results_backup/8_df_bgc_filtered.csv", index=False)
    
    return df_bgc, df_bgc_filtered