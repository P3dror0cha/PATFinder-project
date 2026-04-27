#!/usr/bin/env python3

import requests
import tarfile
import os
import glob
import json
import pandas as pd
from itertools import product

def reference_bgcs():
    """
    Download, extract, and parse reference biosynthetic gene clusters (BGCs)
    from the MIBiG database, returning a structured pandas DataFrame.
    """
    url = "https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz"

    output_dir = os.getcwd()
    tar_path = os.path.join(output_dir, "mibig_json_4.0.tar.gz")

    os.makedirs(output_dir, exist_ok=True)

    with requests.get(url) as r:
        r.raise_for_status()
        with open(tar_path, "wb") as f:
            f.write(r.content)

    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=output_dir)

    extract_path = os.path.join(output_dir, "mibig_json_4.0")
    files = glob.glob(os.path.join(extract_path, "*.json"))

    all_rows = []

    for file in files:
        with open(file) as f:
            data = json.load(f)

        bgc_id = data.get("accession")

        classes = data.get("biosynthesis", {}).get("classes", [])
        compounds = data.get("compounds", [])

        classes = classes or [{}]
        compounds = compounds or [{}]

        for cls, compound in product(classes, compounds):
            all_rows.append({
                "bgc_id": bgc_id,
                "class": cls.get("class"),
                "subclass": cls.get("subclass"),
                "cyclases": cls.get("cyclases"),
                "compound_name": compound.get("name"),
                "compound_structure": compound.get("structure")
            })

    df_ref = pd.DataFrame(all_rows)
    df_ref = df_ref.drop_duplicates(subset="bgc_id", keep="first")

    return df_ref