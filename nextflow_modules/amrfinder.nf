#!/usr/bin/env nextflow

process AMRFINDER {
    label 'amrfinder'
    publishDir "results/amrfinder", mode: 'copy'

    input:
    path faa_files_from_gbks
    path db_ready

    output:
    path "amrfinder_proteins.faa", emit: amrfinder_proteins
    path "amrfinder_results.tsv", emit: amrfinder_result

    script: 
    """
    amrfinder \\
        -p ${faa_files_from_gbks} \\
        --protein_output amrfinder_proteins.faa \\
        --output amrfinder_results.tsv
    """
}