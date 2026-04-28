#!/usr/bin/env python3

import pandas as pd
import re

def K7_kegg_ids_filtering(x):
    """
    Filters KEGG pathway IDs to retain only those associated with antibiotic
    biosynthesis or related pathways.

    Parameters
    ----------
    x : str or NaN
        String containing one or more KEGG pathway IDs (e.g., "map00940;map00945"),
        or NaN.

    Returns
    -------
    str or None
        A semicolon-separated string of filtered KEGG pathway IDs that match
        the predefined set of interest (e.g., antibiotic biosynthesis pathways).
        Returns None if no matching IDs are found or if input is NaN.
    """

    antibiotic_biosynthesis_ids = [    
    "00940", "00945", "00941", "00944", "00942", "00943",
    "00946", "00901", "00403", "00950", "00960", "00996",
    "00232", "00965", "00966", "00402", "00311", "00332",
    "00261", "00331", "00521", "00524", "00525", "00401",
    "00404", "00405", "00333", "00254", "00975", "00998",
    "00999", "00997", "01501", "01502", "01503", "00522", 
    "01051", "01059", "00253", "01053", "01055"]

    kegg_set = set(antibiotic_biosynthesis_ids)
    
    if pd.isna(x):
        return None
    
    ids = re.findall(r"map(\d{5})", x)
    
    filtered_ids = [f"map{id_}" for id_ in ids if id_ in kegg_set]

    return ";".join(filtered_ids) if filtered_ids else None