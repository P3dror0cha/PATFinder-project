#!/usr/bin/env python3

import pandas as pd

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
        - BGC_ID extracted from query (e.g., MGYGxxxx_ctgXX)
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

    final_df["BGC_ID"] = final_df["query"].str.extract(r"(MGYG\d+).*?ctg(\d+)") \
        .agg("_".join, axis=1)
    final_df["gene_id"] = final_df["query"].apply(
        lambda x: "_".join(x.split("_")[1:]) if isinstance(x, str) else None
    )

    final_df["percentual"] = final_df.groupby("BGC_ID")["quality"].transform("sum") / final_df.groupby("BGC_ID")["BGC_ID"].transform("size") * 100

    final_df_filtered = final_df.drop(columns=["query", "evalue", "description"])

    return final_df_filtered
