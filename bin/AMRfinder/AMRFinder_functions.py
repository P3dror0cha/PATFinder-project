import pandas as pd

def parse_amrfinder_result(amrfinder_tsv):
    """
    Reads AMRFinderPlus TSV output and extracts bgc_id from Protein id.

    The Protein id column contains the full FASTA header from EXTRACTING_GBKS_CDS,
    which has the format: bgc_id|feature_id_start-end
    This function extracts bgc_id (everything before the '|') and keeps the
    medium-detail columns: Element symbol, Class, Subclass.

    Args:
        amrfinder_tsv (str): Path to AMRFinderPlus TSV output file.

    Returns:
        pd.DataFrame: DataFrame with bgc_id and AMRFinder annotation columns.
    """
    df = pd.read_csv(amrfinder_tsv, sep="\t")
    df["bgc_id"] = df["Protein id"].str.extract(r'^([^|]+)')
    df_final = df[["bgc_id", "Element symbol", "Class", "Subclass"]]
    return df_final

##############################################################################################

def clustering_results_by_id(df):
    """
    Groups AMRFinder results by bgc_id and aggregates into semicolon-joined strings.

    Args:
        df (pd.DataFrame): Output from parse_amrfinder_result.

    Returns:
        pd.DataFrame: One row per bgc_id with aggregated AMR annotations.
    """
    df_grouped = df.groupby("bgc_id").agg(
        amrfinder_element_symbol=("Element symbol", lambda x: "; ".join(x.astype(str))),
        amrfinder_class=("Class", lambda x: "; ".join(x.astype(str))),
        amrfinder_subclass=("Subclass", lambda x: "; ".join(x.astype(str))),
        amrfinder_hits_count=("Element symbol", "size")
    ).reset_index()
    return df_grouped