import sys
import argparse
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO

from FAA2_DeepSEA_functions import (
    DS1_append_deepsea_results,
    DS2_append_faa,
    DS3_merge_deepsea_df,
    DS4_resistance_class_distribution,
    DS5_resistance_proteins_in_metagenomes,
    DS6_parse_antismash_gbks,
    DS7_filtering_antismash_gbks,
    DS8_merging_dataframes,
    DS9_download_resistance_proteins_info,
    DS10_deepsea_heatmap
)

def DeepSEA_process(path_to_deepsea_tsv, path_to_all_faa, antismash_output, output_image_dir="deepsea_images"):

    os.makedirs(output_image_dir, exist_ok=True)

    deepsea_raw_df = DS1_append_deepsea_results(path_to_deepsea_tsv)
    all_faa_df = DS2_append_faa(path_to_all_faa)
    resistance_proteins_df = DS3_merge_deepsea_df(deepsea_raw_df, all_faa_df)

    DS4_resistance_class_distribution(resistance_proteins_df, output_image_dir)
    DS5_resistance_proteins_in_metagenomes(resistance_proteins_df, output_image_dir)

    df_bgc_features = DS6_parse_antismash_gbks(antismash_output)
    df_bgc_features = DS7_filtering_antismash_gbks(df_bgc_features)

    df_bgc_resistance_proteins = DS8_merging_dataframes(df_bgc_features, resistance_proteins_df)
    df_bgc_resistance_proteins_filtered = DS9_download_resistance_proteins_info(df_bgc_resistance_proteins)

    DS10_deepsea_heatmap(df_bgc_resistance_proteins, output_image_dir)

    return df_bgc_resistance_proteins_filtered

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Process DeepSEA data and AntiSMASH GBKs")
    parser.add_argument("--path_to_deepsea_tsv", required=True, help="Path to TSV file")
    parser.add_argument("--path_to_all_faa", required=True, help="Path to FAA files")
    parser.add_argument("--antismash_output", required=True, help="Path to AntiSMASH output")
    
    args = parser.parse_args()

    final_deepsea_result = DeepSEA_process(
        args.path_to_deepsea_tsv, 
        args.path_to_all_faa, 
        args.antismash_output
    )

    final_deepsea_result.to_csv("deepsea_final_merged.csv", index=False)