import pandas as pd

def parse_diamond_result(diamond_output_file, antibiotic_biosynthesis_proteins_tsv_path, antibiotic_resistance_proteins_tsv_path):
    """
    Parses a DIAMOND output file and annotates hits with UniProt database info.

    Reads a standard tabular DIAMOND alignment output, filters for the best hit 
    (lowest e-value) per query sequence. It then merges these filtered results 
    with local UniProt TSV databases containing antibiotic biosynthesis and 
    resistance protein data ().

    Args:
        diamond_output_file (str): Path to the tabular DIAMOND output file.

    Returns:
        pd.DataFrame: A merged DataFrame containing the top DIAMOND hits along 
            with their corresponding UniProt annotations.
    """

    df = pd.read_csv(diamond_output_file, sep="\t", header=None)
    
    df.columns = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    df_lowest_evalue = df.loc[df.groupby('qseqid')['evalue'].idxmin()]
    df_filtered = df_lowest_evalue.drop(columns=["length", "qstart", "qend", "sstart", "send"])
    
    df_proteins_biosynthesis_in_uniprot = pd.read_csv("/home/pedro/pedro_nfs/antismash/BIG-SCAPE/uniprotkb_antibiotic_biosynthesis_2026_05_27.tsv", sep="\t")
    df_proteins_resistance_in_uniprot = pd.read_csv("/home/pedro/pedro_nfs/antismash/BIG-SCAPE/uniprotkb_antibiotic_resistance_2026_05_27.tsv", sep="\t")
    
    df_concat = pd.concat([df_proteins_biosynthesis_in_uniprot, df_proteins_resistance_in_uniprot], ignore_index=True)
    df_united = pd.merge(df_filtered, df_concat, left_on="sseqid", right_on="Entry", how="left")
    
    return df_united