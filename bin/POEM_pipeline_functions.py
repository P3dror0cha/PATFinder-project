import pandas as pd
import re

def POEM1_parse_operon_string(result_from_POEM):
    """
    Parses a complex operon string to extract metagenome ID, strand, and coordinates.

    Splits a string of linked genes by directional arrows ('-->' or '<--') and uses 
    regular expressions to extract metadata for each gene. It then reconstructs 
    the coordinates into a unified, direction-aware string.

    Args:
        result_from_POEM (str): Raw operon string containing gene information. 
            Expected to contain substrings matching the pattern: 
            '|{strand}|{start}|{end}$${metagenome_id}'.

    Returns:
        pd.Series: A Series containing three parsed elements, designed to be used 
            with DataFrame.apply():
            - metagenome (str): The unique metagenome identifier (e.g., 'sample_123_456').
            - strand (str): The DNA strand orientation ('+' or '-').
            - coord_string (str): A concatenated string of start/end coordinates 
              linked by arrows indicating transcription direction (e.g., '100/500 --> 600/900').
    """

    genes = re.split(r'-->|<--', str(result_from_POEM))

    metagenome = None
    strand = None
    coords = []

    for g in genes:
        m = re.search(r'\|([+-])\|(\d+)\|(\d+)\$\$(sample_\d+_\d+)\.region(\d+)', g)
        #m = re.search(r'\|([+-])\|(\d+)\|(\d+)\$\$(sample_\d+_\d+)', g)
        if m:
            strand = m.group(1)      
            start = m.group(2)        
            end = m.group(3)         
            metagenome = m.group(4) 
            region = m.group(5)  

            metagenome = f"{metagenome}_{region}"
            coords.append(f"{start}/{end}")

       
    if strand == "+":
        coord_string = " --> ".join(coords)
    else:
        coord_string = " <-- ".join(coords)

    return pd.Series([metagenome, strand, coord_string])

###################################################################################

def POEM2_process_poem_table(table_path, sep="\t"):
    """
    Reads a tabular file and extracts genomic metadata from the 'gene_id' column.

    Loads the data into a pandas DataFrame, standardizes column names, and applies 
    the `POEM1_parse_operon_string` function to generate new metadata columns 
    based on the raw operon strings.
    
    Args:
        table_path (str): Path to the input CSV/TSV file.
        sep (str, optional): The delimiter string used in the file. Defaults to "\t".
        
    Returns:
        pd.DataFrame: The processed DataFrame containing three columns: 
            'metagenome_id', 'strand', and 'coordinates'.
        
    Raises:
        ValueError: If the required 'gene_id' column is missing from the input file.
    """

    df = pd.read_csv(table_path, sep=sep)
    df.columns = df.columns.str.strip().str.lower()

    if "gene_id" not in df.columns:
        raise ValueError(f"Column 'gene_id' not found. Available columns: {list(df.columns)}")

    df[["metagenome_id", "strand", "coordinates"]] = df["gene_id"].apply(POEM1_parse_operon_string)

    return df

###################################################################################

def POEM3_clustering_results(df):
    """
    Aggregate the coordinate results by metagenome id and strand.
    
    Args:
        poem_dataframe (str): Path to the output of the "POEM2_process_poem_table".
        
    Returns:
        pd.DataFrame: The processed DataFrame containing three columns: 
            'metagenome_id', 'strand', and 'coordinates'.
    """

    df = df.groupby(['metagenome_id', 'strand'], as_index=False)['coordinates'].agg(' / '.join)

    return df
