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
    Filter and preprocess BiG-SCAPE full_network results.

    This function reads a BiG-SCAPE full network file (typically matching the pattern
    './bigscape/output_files/*full.network') and applies a series of filtering and cleaning steps
    to prepare the data for downstream analysis.

    Steps performed:
    1. Load the full network file as a tab-separated table.
    2. Remove self-comparisons between BGCs (e.g., BGC vs BGC) and between MGYG entries.
    3. Use the first row as the header and drop it from the dataset.
    4. Remove non-essential columns related to record metadata and alignment parameters.
    5. Normalize GBK identifiers by removing region suffixes (e.g., '.region001').
    6. Ensure consistent ordering of comparisons:
    - If one entry is a reference BGC and the other is not, swap columns so that
        non-BGC entries are consistently placed in the same position.
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
    Download, extract, and parse reference biosynthetic gene clusters (BGCs)
    from the MIBiG database, returning a structured pandas DataFrame.
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

    df_ref_filtered = df_ref[['bgc_id', 'class', 'compound_name']]
    df_ref_filtered = df_ref_filtered.rename(columns={'bgc_id': 'GBK_b'})
    df_final = pd.merge(df_filtered, df_ref_filtered, on='GBK_b', how='left')
    df_final['distance'] = pd.to_numeric(df_final['distance'], errors='coerce')
    df_best_result = df_final.loc[df_final.groupby('GBK_a')['distance'].idxmin()]

    return df_best_result