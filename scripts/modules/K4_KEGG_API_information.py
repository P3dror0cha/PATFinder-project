#!/usr/bin/env python3

import pandas as pd
import requests

def K4_KEGG_API_information(url="https://rest.kegg.jp/link/pathway/ko", url_2="https://rest.kegg.jp/list/pathway"):
    """
    Retrieves KEGG pathway mapping and pathway descriptions via the KEGG REST API.

    Parameters
    ----------
    url : str, optional
        KEGG API endpoint that links KO identifiers to pathways.
        Default is "https://rest.kegg.jp/link/pathway/ko".

    url_2 : str, optional
        KEGG API endpoint that provides pathway IDs and their descriptions.
        Default is "https://rest.kegg.jp/list/pathway".

    Returns
    -------
    tuple of pandas.DataFrame
        - df_kegg : DataFrame mapping KO identifiers to pathway IDs with columns:
            'KO' (KEGG Orthology ID) and 'pathway' (pathway ID).
        - df_list : DataFrame containing pathway metadata with columns:
            'pathway' (pathway ID) and 'description' (pathway name/description).
    """
    r = requests.get(url)

    lines = r.text.strip().split("\n")

    data = [line.split("\t") for line in lines if "\t" in line]

    df_kegg = pd.DataFrame(data, columns=["KO", "pathway"])

    df_kegg["KO"] = df_kegg["KO"].str.replace("ko:", "")
    df_kegg["pathway"] = df_kegg["pathway"].str.replace("path:", "")
    
    r_2 = requests.get(url_2)

    lines = r_2.text.strip().split("\n")

    data = [line.split("\t") for line in lines if "\t" in line]

    df_list = pd.DataFrame(data, columns=["pathway", "description"])

    return df_kegg, df_list

