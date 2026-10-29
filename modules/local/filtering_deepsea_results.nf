#!/usr/bin/env nextflow

process FILTERING_DEEPSEA_RESULTS {
    label 'filtering_deepsea_results'
    publishDir "results/deepsea", mode: 'copy'

    input:
    path path_to_deepsea_tsv
    path "faa_temp_folder/*"
    path "gbks_temp_folder/*"

    output:
    path "deepsea_final_merged.csv", emit: deepsea_csv
    path "deepsea_images/", emit: deepsea_images
    

script: 
    """
    python3 ${projectDir}/bin/DeepSEA_process.py \
        --path_to_deepsea_tsv ${path_to_deepsea_tsv} \
        --path_to_all_faa faa_temp_folder \
        --antismash_output gbks_temp_folder
    """
}
