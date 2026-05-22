# Parsing CDS sequences from .gbk files in /home/pedro/antismash/resultados/all_BGCs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import glob
import os


def extract_cds_from_gbk(gbk_files):

    results = []

    for gbk in glob.glob(gbk_files):
        for record in SeqIO.parse(gbk, "genbank"):

            for feature in record.features:

                if feature.type not in ["CDS", "aSDomain"]:
                    continue

                qualifiers = feature.qualifiers

                feature_id = (
                    qualifiers.get("locus_tag", [None])[0]
                    or qualifiers.get("protein_id", [None])[0]
                    or qualifiers.get("gene", [None])[0]
                    or qualifiers.get("ID", [None])[0]
                    or "unknown"
                )

                if feature.type == "CDS":
                    if "translation" in qualifiers:
                        sequence = qualifiers["translation"][0]
                    else:
                        sequence = feature.extract(record.seq).translate(to_stop=True)

                elif feature.type == "aSDomain":
                    sequence = str(feature.extract(record.seq))

                results.append({
                    "gbk_file": os.path.basename(gbk),
                    "record_id": record.id,
                    "feature_type": feature.type,
                    "feature_id": feature_id,
                    "start": int(feature.location.start),
                    "end": int(feature.location.end),
                    "strand": feature.location.strand,
                    "sequence": str(sequence)
                })
    df = pd.DataFrame(results)
    df.to_csv("all_gbks_cds_sequences.csv", index=False)
    
    return df

def filtering_cds(df):
    df["start_end"] = df["start"].astype(str) + "-" + df["end"].astype(str)
    df["unique_id"] = df["feature_id"] + "_" + df["start_end"]
    df = df.drop_duplicates(subset="unique_id")
    
    df = df.reindex(columns=[
    "unique_id",
    "feature_type",
    "start_end",
    "sequence"
    ])
    
    records = [
    SeqRecord(
        Seq(seq),
        id=cds_id,
        description=""
    )
    for cds_id, seq in zip(df["unique_id"], df["sequence"])
    ]

    SeqIO.write(records, "antismash_cds.faa", "fasta")
    
    return df

