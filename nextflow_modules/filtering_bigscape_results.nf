#!/usr/bin/env nextflow

process FILTERING_BIGSCAPE_RESULTS {

    label 'bigscape_results'
    publishDir "results/bigscape_results", mode: 'copy'

    input:
    path bigscape_fullnetwork

    output:
    path "filtered_bigscape_results.tsv", emit: filtered_bigscape_results

script:
    """
    python3 ${projectDir}/bin/PY2_bigscape_output_analysis_process.py \
        --fullnetwork_file_path ${bigscape_fullnetwork} \
        --output_file filtered_bigscape_results.tsv
    """
    }
