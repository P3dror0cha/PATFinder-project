import pandas as pd
import argparse

from AMRFinder_functions import (
    parse_amrfinder_result,
    clustering_results_by_id
)

def AMRFinder_process(amrfinder_output, output_csv):
    df_annotated = parse_amrfinder_result(amrfinder_output)
    df_annotated = clustering_results_by_id(df_annotated)
    df_annotated.to_csv(output_csv, index=False)
    return df_annotated

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process AMRFinderPlus output")
    parser.add_argument("--amrfinder_output", required=True, help="Path to AMRFinderPlus TSV output")
    parser.add_argument("--output_csv", default="amrfinder_annotated_results.csv", help="Name of the output CSV file")

    args = parser.parse_args()
    AMRFinder_process(args.amrfinder_output, args.output_csv)