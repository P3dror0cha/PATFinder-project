#!/usr/bin/env python3

import pandas as pd

#'./bigscape/*full.network'
def filtering_bigscape_results(output_dir):
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
    df_filtered = df_filtered[~(df_filtered[1].str.contains("MGYG") & df_filtered[6].str.contains("MGYG"))]
    df_filtered.columns = df_filtered.iloc[0]
    df_filtered = df_filtered[1:]
    df_filtered = df_filtered.drop(columns=['Record_a', 'Record_Type_a', 'Record_Number_a', 'Record_b', 'Record_Type_b', 'Record_Number_b', "alignment_mode", "extend_strategy"])
    df_filtered["GBK_a"] = df_filtered["GBK_a"].str.replace(r"\.region\d+", "", regex=True)
    df_filtered["GBK_b"] = df_filtered["GBK_b"].str.replace(r"\.region\d+", "", regex=True)
    mask = df_filtered["GBK_a"].str.startswith("BGC") & ~df_filtered["GBK_b"].str.startswith("BGC")
    df_filtered.loc[mask, ["GBK_a", "GBK_b", "ORF_coords_a", "ORF_coords_b"]] = df_filtered.loc[mask, ["GBK_b", "GBK_a", "ORF_coords_b", "ORF_coords_a"]].values
    return df_filtered