#!/usr/bin/env python3

import pandas as pd
import requests
import re

from K1_kofam_import import K1_kofam_import
from K2_kofam_best_hits import K2_kofam_best_hits
from K3_kofam_without_best_hits import K3_kofam_without_best_hits
from K4_KEGG_API_information import K4_KEGG_API_information
from K5_merging_kofam_df import K5_merging_kofam_df
from K6_uniting_BGC_information import K6_uniting_BGC_information
from K7_kegg_ids_filtering import K7_kegg_ids_filtering
from K8_uniting_all_information import K8_uniting_all_information

def K9_kofam_pipeline(kofam_result_list):

    # filtered_bigscape_results are obtained after BIG-SCAPE output python functions.
    kofam_results = K1_kofam_import(kofam_result_list)
    kofam_best_hits = K2_kofam_best_hits(kofam_results)
    final_df_filtered = K3_kofam_without_best_hits(kofam_results, kofam_best_hits)
    df_kegg, df_list = K4_KEGG_API_information()
    df_groupby = K5_merging_kofam_df(final_df_filtered, df_list, df_kegg)
    df_bgc = K6_uniting_BGC_information(filtered_bigscape_results, df_groupby)
    list_of_pathways = K7_kegg_ids_filtering()
    df_bgc, df_bgc_filtered = K8_uniting_all_information(df_bgc, output_path_raw_csv, output_path_filtered_csv)