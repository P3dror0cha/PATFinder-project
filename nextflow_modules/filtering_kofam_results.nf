#!/usr/bin/env nextflow

process FILTERING_KOFAM_RESULTS {
    label 'filtering_kofam_results'
    
    publishDir "results/kofam", mode: 'copy'

    input:
    path kofam_output
    path filtered_bigscape_results

    output:
    path "kofam_bigscape_raw_results.csv", emit: kofam_bigscape_raw_results
    path "kofam_bigscape_filtered_results.csv", emit: kofam_bigscape_filtered_results

    script:
    """
    python3 ${projectDir}/bin/kofam_analysis_process.py \\
        --kofam_results ${kofam_output} \\
        --bgc_dir ${filtered_bigscape_results} \\
        --out_raw kofam_bigscape_raw_results.csv \\
        --out_filtered kofam_bigscape_filtered_results.csv
    """
}