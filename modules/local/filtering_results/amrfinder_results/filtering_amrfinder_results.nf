#!/usr/bin/env nextflow

process FILTERING_AMRFINDER_RESULTS {
    label 'filtering_amrfinder_results'
    publishDir "results/amrfinder", mode: 'copy'

    input:
    path amrfinder_tsv_output

    output:
    path "amrfinder_annotated_results.csv", emit: amrfinder_filtered_result

    script:
    """
    python3 ${projectDir}/bin/AMRFinder_process.py \
        --amrfinder_output ${amrfinder_tsv_output} \
        --output_csv amrfinder_annotated_results.csv
    """
}