#!/usr/bin/env nextflow

process FILTERING_BIGSCAPE_RESULTS {

    label 'filtering_bigscape_results'
    publishDir "results", mode: 'copy'

    input:
    path bigscape_fullnetwork

    output:
    path "filtered_bigscape_results.tsv", emit: filtered_bigscape_results

script:
"""
python3 ${projectDir}/bin/BS4_processing_bigscape_results.py \
    --fullnetwork_file_path ${bigscape_fullnetwork} \
    --output_file filtered_bigscape_results.tsv
"""
}
