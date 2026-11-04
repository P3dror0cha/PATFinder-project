import pandas as pd
import argparse

from Diamond_functions import (
    parse_diamond_result,
    clustering_results_by_id
)

def Diamond_process(diamond_output_file, antibiotic_biosynthesis_proteins_tsv_path, antibiotic_resistance_proteins_tsv_path, output_csv="diamond_annotated_results.csv",):

    df_annotated = parse_diamond_result(diamond_output_file, antibiotic_biosynthesis_proteins_tsv_path, antibiotic_resistance_proteins_tsv_path)
    df_annotated = clustering_results_by_id(df_annotated)
    df_annotated.to_csv(output_csv, index=False)

    return df_annotated

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process DIAMOND output and annotate with UniProt information")
    parser.add_argument("--diamond_output", required=True, help="Path to DIAMOND output file (tab-separated)")
    parser.add_argument("--biosynthesis_uniprot_proteins_path", required=True, help="Path to the file containing the uniprot proteins associated with biosynthesis of antibiotics")
    parser.add_argument("--resistance_uniprot_proteins_path", required=True, help="Path to the file containing the uniprot proteins associated with resistance to antibiotics")
    parser.add_argument("--output_csv", default="diamond_annotated_results.csv", help="Name of the output CSV file")

    args = parser.parse_args()

    annotated_df = Diamond_process(args.diamond_output, args.biosynthesis_uniprot_proteins_path, args.resistance_uniprot_proteins_path, args.output_csv,)
    print(f"Annotated DIAMOND results saved to: {args.output_csv}")