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
    
    df_proteins_biosynthesis_in_uniprot = pd.read_csv(antibiotic_biosynthesis_proteins_tsv_path, sep="\t")
    df_proteins_resistance_in_uniprot = pd.read_csv(antibiotic_resistance_proteins_tsv_path, sep="\t")
    
    df_concat = pd.concat([df_proteins_biosynthesis_in_uniprot, df_proteins_resistance_in_uniprot], ignore_index=True)
    df_concat.rename(columns={"Entry Name": "uniprot_id"}, inplace=True)
    df_concat.to_csv("df_concat_of_protein_resistance_and_biosynthesis.csv")

    df_lowest_evalue["uniprot_id"] = df_lowest_evalue["sseqid"].str.extract(r'([^|]+)$')
    df_lowest_evalue.to_csv("df_lowest_evalue.csv")

    df_merged = pd.merge(df_lowest_evalue, df_concat, on='uniprot_id', how='left')
    df_merged["protein_id"] = df_merged["qseqid"].str.extract(r'([^|]+)$')
    df_merged["bgc_id"] = df_merged["qseqid"].str.extract(r'^([^|]+)')
    df_merged.to_csv("df_merged.csv")

    df_final_alignment = df_merged[["bgc_id", "protein_id", "uniprot_id", "evalue", "Protein names", "Gene Names"]]

    df_final_alignment = df_final_alignment.rename(columns={
        "evalue": "uniprot_evalue", 
        "Protein names": "uniprot_proteins", 
        "Gene Names": "uniprot_genes"
    })
    
    return df_final_alignment

##################################################################################

def clustering_results_by_id(df_final_alignment):
    """
    Groups and aggregates UniProt alignment results by Biosynthetic Gene Cluster (BGC) ID.

    Performs a named aggregation on the alignment DataFrame, grouping by 'bgc_id'.
    It concatenates specific UniProt metadata (IDs, e-values, protein names, and genes) 
    into comma-separated strings for each cluster and calculates the total hit count.

    Args:
        df_final_alignment (pd.DataFrame): DataFrame containing raw sequence alignment 
            results with UniProt annotations.

    Returns:
        pd.DataFrame: An aggregated DataFrame with one row per 'bgc_id', listing 
            grouped UniProt information and the total count of hits per cluster.
    """
    df_final_alignment = df_final_alignment.groupby('bgc_id').agg(
    uniprot_id=('uniprot_id', lambda x: ', '.join(x.astype(str))),
    uniprot_evalue=('uniprot_evalue', lambda x: ', '.join(x.astype(str))),
    uniprot_proteins=('uniprot_proteins', lambda x: ', '.join(x.astype(str))),
    uniprot_genes=('uniprot_genes', lambda x: ', '.join(x.astype(str))),
    uniprot_hits_count=('uniprot_id', 'size')  
    ).reset_index()

    return df_final_alignment