import os
import argparse
from Bio import SeqIO

def extract_bgc_to_fasta(gbk_files: list, output_fasta: str) -> int:
    """
    Extracts Biosynthetic Gene Cluster (BGC) regions from GenBank (.gbk) files into a FASTA.
    These files are the normal output of antiSMASH.

    Iterates through a list of .gbk, searching for features specifically 
    annotated as 'region' or 'cluster' (typical BGC annotations). It extracts the 
    nucleotide sequences corresponding to these features and writes them 
    to a single FASTA file.

    Args:
        gbk_files (list[str]): A list of file paths to the .gbk files (Put all in one folder and use the path).
        output_fasta (str): The name or path for the resulting output FASTA file.

    Returns:
        int: The total number of files provided in the input list (0 if empty).
        
    Note:
        - Skips missing files and prints a warning message.
        - Overwrites the `output_fasta` file if it already exists.
    """
    if not gbk_files:
        print("No input files provided.")
        return 0

    with open(output_fasta, "w") as out_file:
        for gbk_path in gbk_files:
         
            if not os.path.exists(gbk_path):
                print(f"File not found -> {gbk_path}")
                continue

            base_name = os.path.basename(gbk_path)

            for record in SeqIO.parse(gbk_path, "genbank"):
                for feature in record.features:
                    if feature.type in ["region", "cluster"]:
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        seq = record.seq[start:end]

                        header = f"{base_name}"

                        out_file.write(f">{header}\n")
                        out_file.write(f"{str(seq)}\n")

    print(f"Processed files: {len(gbk_files)}")
    print(f"Generated file: {output_fasta}")
    
    return len(gbk_files)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extracts BGC regions from GenBank (.gbk) files and converts them to FASTA."
    )
    
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        nargs='+', 
        help="List of input .gbk files."
    )
    
    parser.add_argument(
        "-o", "--output", 
        default="BGCS_fna_sequences.fasta", 
        help="Name of the output FASTA file."
    )
    
    args = parser.parse_args()
    
    extract_bgc_to_fasta(args.input, args.output)