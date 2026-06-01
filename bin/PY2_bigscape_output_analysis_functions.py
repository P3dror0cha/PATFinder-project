#!/usr/bin/env python3

import pandas as pd
from itertools import product
import requests
import tarfile
import os
import glob
import json


#'./bigscape/*full.network'
def BS1_filtering_bigscape_results(output_dir):
    """
    Filters and preprocesses BiG-SCAPE full_network results.

    Reads a BiG-SCAPE full network file and applies filtering steps to remove 
    self-comparisons (BGC vs BGC, sample vs sample). It normalizes GBK identifiers 
    by removing region suffixes and ensures consistent column ordering so that 
    non-BGC sample entries are always placed in the 'GBK_a' position and 
    reference BGCs in 'GBK_b'.

    Args:
        output_dir (str): Path to the BiG-SCAPE full network file (tab-separated).

    Returns:
        pd.DataFrame: A cleaned and filtered DataFrame ready for downstream analysis.
    """

    df = pd.read_csv(output_dir, sep='\t', header=None)

    df_filtered = df[~(df[1].str.contains("BGC") & df[6].str.contains("BGC"))]
    df_filtered = df_filtered[~(df_filtered[1].str.contains("sample") & df_filtered[6].str.contains("sample"))]
    df_filtered.columns = df_filtered.iloc[0]
    df_filtered = df_filtered[1:]
    df_filtered = df_filtered.drop(columns=['Record_a', 'Record_Type_a', 'Record_Number_a', 'Record_b', 'Record_Type_b', 'Record_Number_b', "alignment_mode", "extend_strategy"])
    df_filtered["GBK_a"] = df_filtered["GBK_a"].str.replace(r"\.region\d+", "", regex=True)
    df_filtered["GBK_b"] = df_filtered["GBK_b"].str.replace(r"\.region\d+", "", regex=True)
    mask = df_filtered["GBK_a"].str.startswith("BGC") & ~df_filtered["GBK_b"].str.startswith("BGC")
    df_filtered.loc[mask, ["GBK_a", "GBK_b", "ORF_coords_a", "ORF_coords_b"]] = df_filtered.loc[mask, ["GBK_b", "GBK_a", "ORF_coords_b", "ORF_coords_a"]].values
    return df_filtered

##############################################################################################

def BS2_reference_bgcs():
    """
    Downloads and parses reference BGCs from the MIBiG database.

    Fetches the MIBiG JSON dataset tarball, extracts its contents to the 
    current working directory, and parses the individual JSON files to extract 
    metadata such as BGC accession ID, biosynthetic class, and compound name.

    Returns:
        pd.DataFrame: A deduplicated DataFrame containing MIBiG reference 
            BGC metadata.
            
    Note:
        The URL for the MIBiG database (version 4.0) and the download/extraction 
        paths are hardcoded to the current working directory (`os.getcwd()`).
    """
    url = "https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz"

    output_dir = os.getcwd()
    tar_path = os.path.join(output_dir, "mibig_json_4.0.tar.gz")

    os.makedirs(output_dir, exist_ok=True)

    with requests.get(url) as r:
        r.raise_for_status()
        with open(tar_path, "wb") as f:
            f.write(r.content)

    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=output_dir)

    extract_path = os.path.join(output_dir, "mibig_json_4.0")
    files = glob.glob(os.path.join(extract_path, "*.json"))

    all_rows = []

    for file in files:
        with open(file) as f:
            data = json.load(f)

        bgc_id = data.get("accession")

        classes = data.get("biosynthesis", {}).get("classes", [])
        compounds = data.get("compounds", [])

        classes = classes or [{}]
        compounds = compounds or [{}]

        for cls, compound in product(classes, compounds):
            all_rows.append({
                "bgc_id": bgc_id,
                "class": cls.get("class"),
                "subclass": cls.get("subclass"),
                "cyclases": cls.get("cyclases"),
                "compound_name": compound.get("name"),
                "compound_structure": compound.get("structure")
            })

    df_ref = pd.DataFrame(all_rows)
    df_ref = df_ref.drop_duplicates(subset="bgc_id", keep="first")

    return df_ref

##############################################################################################

def BS3_uniting_reference_with_sample_BGCs(df_ref, df_filtered):
    """
    Merges sample BGC networks with MIBiG reference annotations.

    Performs a left join between the filtered BiG-SCAPE network and the MIBiG 
    reference dataset based on the reference BGC ID. It then converts the 
    alignment distance to a numeric type and filters the dataset to retain 
    only the best hit (minimum distance) for each sample BGC.

    Args:
        df_ref (pd.DataFrame): DataFrame containing MIBiG reference BGCs (from BS2).
        df_filtered (pd.DataFrame): Filtered BiG-SCAPE network DataFrame (from BS1).

    Returns:
        pd.DataFrame: A merged DataFrame containing only the top MIBiG reference 
            hit for each queried sample BGC.
    """

    df_ref_filtered = df_ref[['bgc_id', 'class', 'compound_name']]
    df_ref_filtered = df_ref_filtered.rename(columns={'bgc_id': 'GBK_b'})
    df_final = pd.merge(df_filtered, df_ref_filtered, on='GBK_b', how='left')
    df_final['distance'] = pd.to_numeric(df_final['distance'], errors='coerce')
    df_best_result = df_final.loc[df_final.groupby('GBK_a')['distance'].idxmin()]

    return df_best_result