#!/usr/bin/env python3

import pandas as pd

def K6_uniting_BGC_information(antiSMASH_result, df_groupby):
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
    df_bgc = pd.read_csv(antiSMASH_result)
    df_bgc = pd.merge(df_bgc, df_groupby, on="BGC_ID", how="left")

    return df_bgc

