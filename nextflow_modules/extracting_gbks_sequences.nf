#!/usr/bin/env nextflow

process EXTRACTING_GBKS_SEQUENCES {

    publishDir "results/gbks_fna_sequence", mode: 'copy'

    input:
    path gbk_files

    output:
    path "BGCS_fna_sequences.fasta", emit: bgc_sequences

    script:
    """
    python3 ${projectDir}/bin/extracting_fna_sequence_from_gbks.py \
        --input ${gbk_files} \
        --output BGCS_fna_sequences.fasta
    """
}


