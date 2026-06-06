import pandas as pd

def parse_diamond_result(diamond_output_file, antibiotic_biosynthesis_proteins_tsv_path="./uniprot_proteins/uniprotkb_antibiotic_biosynthesis_2026_05_27", antibiotic_resistance_proteins_tsv_path="./uniprot_proteins/uniprotkb_antibiotic_resistance_2026_05_27.tsv"):
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
    
    df_proteins_biosynthesis_in_uniprot = pd.read_csv(antibiotic_biosynthesis_proteins_tsv_path, sep="\t")
    df_proteins_resistance_in_uniprot = pd.read_csv(antibiotic_resistance_proteins_tsv_path, sep="\t")
    
    df_concat = pd.concat([df_proteins_biosynthesis_in_uniprot, df_proteins_resistance_in_uniprot], ignore_index=True)
    df_concat.rename(columns={"Entry": "sseqid"}, inplace=True)
    df_merged = pd.merge(df_lowest_evalue, df_concat, on='sseqid', how='left')

    df_final_alignment = df_merged[["qseqid", "sseqid", "evalue", "Protein names", "Gene Names"]]
    df_final_alignment = df_final_alignment.rename(columns={
        "qseqid": "id", 
        "sseqid": "uniprot_entry", 
        "Protein names": "uniprot_proteins", 
        "Gene Names": "uniprot_genes"
    })
    
    return df_final_alignment

##################################################################################

