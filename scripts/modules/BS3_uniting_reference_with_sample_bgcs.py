#!/usr/bin/env python3
import pandas as pd

def uniting_reference_with_sample_BGCs(df_ref, df_filtered):

    df_ref_filtered = df_ref[['bgc_id', 'class', 'compound_name']]
    df_ref_filtered = df_ref_filtered.rename(columns={'bgc_id': 'GBK_b'})
    df_final = pd.merge(df_filtered, df_ref_filtered, on='GBK_b', how='left')
    df_final['distance'] = pd.to_numeric(df_final['distance'], errors='coerce')
    df_best_result = df_final.loc[df_final.groupby('GBK_a')['distance'].idxmin()]

    return df_best_result