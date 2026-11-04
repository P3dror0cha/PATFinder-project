#!/usr/bin/env nextflow

process FILTERING_BIGSCAPE_RESULTS {

    label 'bigscape_results'
    publishDir "results/bigscape_results", mode: 'copy'

    input:
    path bigscape_fullnetwork
    path ids_correlation

    output:
    path "filtered_bigscape_results.tsv", emit: filtered_bigscape_results

script:
    """
    python3 ${projectDir}/bin/bigscape_output_analysis_process.py \
        --fullnetwork_file_path ${bigscape_fullnetwork} \
        --genome_mapping ${ids_correlation} \
        --output_file filtered_bigscape_results.tsv
    """
    }
