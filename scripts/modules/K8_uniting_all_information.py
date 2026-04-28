#!/usr/bin/env python3

import pandas as pd
import requests
import re

from bgc_pathways_analysis.python_functions.K7_kegg_ids_filtering import K7_kegg_ids_filtering

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
    df_bgc["filtered_pathways"] = df_bgc["pathway"].apply(K7_kegg_ids_filtering)
    df_bgc.rename(columns={"percentual": "KO_quality_percentual"}, inplace=True)
    df_bgc_filtered = df_bgc.drop(columns=["pathway", "description", "KO_quality_percentual", "strand"])
    df_bgc.to_csv(output_path_raw_csv, index=False)
    df_bgc_filtered.to_csv(output_path_filtered_csv, index=False)
    
    return df_bgc, df_bgc_filtered