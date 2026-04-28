#!/usr/bin/env python3

import pandas as pd

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

    return df_best_hits

