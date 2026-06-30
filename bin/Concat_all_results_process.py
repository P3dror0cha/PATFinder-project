import pandas as pd
import argparse
import re

from Concat_all_results_functions import (
    product_from_gbk_files,
    uniting_poem_deepsea_results,
    uniting_kofam_diamond_results,
    uniting_all_infos
)

def concat_process(gbk_files, poem_results, deepsea_results, kofam_results, diamond_results):

    df_poem = pd.read_csv(poem_results)
    df_deepsea = pd.read_csv(deepsea_results)
    df_kofam = pd.read_csv(kofam_results)
    df_diamond = pd.read_csv(diamond_results)

    df_bgc_product = product_from_gbk_files(gbk_files)

    df_poem_deepsea = uniting_poem_deepsea_results(df_bgc_product, df_poem, df_deepsea)
    df_kofam_diamond = uniting_kofam_diamond_results(df_kofam, df_diamond)

    df_final_concat = uniting_all_infos(df_poem_deepsea, df_kofam_diamond)

    return df_final_concat

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process used to unite the results produced by the pipeline")
    parser.add_argument("-g", "--gbk", nargs="+", required=True, help="List of gbk files")
    parser.add_argument("-p", "--poem", required=True, help="Path to POEM results")
    parser.add_argument("-s", "--deepsea", required=True, help="Path to DeepSea results")
    parser.add_argument("-k", "--kofam", required=True, help="Path to Kofam results")
    parser.add_argument("-d", "--diamond", required=True, help="Path to Diamond and Uniprot results")
    parser.add_argument("-o", "--output", required=True, help="Path to save the final concatenated dataframe")

    args = parser.parse_args()

    df_final = concat_process(
        gbk_files=args.gbk,
        poem_results=args.poem,
        deepsea_results=args.deepsea,
        kofam_results=args.kofam,
        diamond_results=args.diamond
    )
    df_final.to_csv(args.output, index=False)

