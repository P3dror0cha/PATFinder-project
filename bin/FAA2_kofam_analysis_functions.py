#!/usr/bin/env python3

import pandas as pd
import requests
import re

def K1_kofam_import(kofam_result_list):
    """
    Reads a KOfam result file and converts it into a structured pandas DataFrame.

    Parameters
    ----------
    kofam_result_list : str
        Path to the KOfam output file to be parsed.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing parsed KOfam results with the following columns:
        - best_hits : int (1 if the hit is marked with '*', otherwise 0)
        - query : str (query sequence ID)
        - KO : str (KEGG Orthology identifier)
        - threshold : float or str (score threshold used)
        - score : float or str (alignment score)
        - evalue : float or str (expectation value)
        - description : str (KO description)
        - function : str (functional annotation)
    """
    data = []

    with open(kofam_result_list) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split(maxsplit=6)
            
            if parts[0] == '*':
                flag = 1
                parts = parts[1:]  
            else:
                flag = 0
            
            while len(parts) < 6:
                parts.append(None)
            
            data.append([flag] + parts)

    df = pd.DataFrame(
        data,
        columns=['best_hits','query', 'KO', 'threshold', 'score', 'evalue', 'description', 'function']
    )

    return df

