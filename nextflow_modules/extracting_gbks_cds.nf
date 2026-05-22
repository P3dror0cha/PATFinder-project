#!/usr/bin/env nextflow

process EXTRACTING_GBKS_CDS {

    publishDir "results/gbks_cds_sequences", mode: 'copy'

    input:
    path gbk_files

    output:
    path "BGCS_cds_sequences.fasta", emit: bgc_sequences

    script:
    """
    python3 ${projectDir}/bin/extracting_cds_sequences_from_gbks.py \
        --input ${gbk_files} \
        --output BGCS_cds_sequences.fasta
    """
}