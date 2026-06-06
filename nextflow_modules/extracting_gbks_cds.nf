#!/usr/bin/env nextflow

process EXTRACTING_GBKS_CDS {

    publishDir "results/gbks_cds_sequences", mode: 'copy'

    input:
    path gbk_files

    output:
    path "BGCS_cds_sequences.faa", emit: bgc_sequences

    script:
    """
    python3 ${projectDir}/bin/extracting_gbks_cds.py \
        --input "*.gbk" \
        --output BGCS_cds_sequences.faa
    """
}