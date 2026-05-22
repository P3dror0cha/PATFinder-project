#!/usr/bin/env nextflow

process FILTERING_POEM_RESULTS {
    label 'filtering_poem_results'
    publishDir "results/poem_filtered", mode: 'copy'

    input:
    path operon_file

    output:
    path "*.csv", emit: poem_filtered_csv

    script:
    """
    python3 ${projectDir}/bin/POEM_pipeline_process.py \
        --output_POEM ${operon_file}
    """
}