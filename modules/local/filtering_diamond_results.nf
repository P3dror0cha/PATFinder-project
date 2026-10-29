#!/usr/bin/env nextflow

process FILTERING_DIAMOND_RESULTS {
    label 'filtering_diamond_results'
    publishDir "results/diamond", mode: 'copy'

    input:
    path cds_from_gbks_vs_uniprot

    output:
    path "diamond_annotated_results.csv", emit: diamond_filtered_result

    script: 
    """
    python3 ${projectDir}/bin/Diamond_process.py \
        --diamond_output ${cds_from_gbks_vs_uniprot} \
        --biosynthesis_uniprot_proteins_path ${projectDir}/uniprot_proteins/uniprotkb_antibiotic_biosynthesis_2026_05_27.tsv \
        --resistance_uniprot_proteins_path ${projectDir}/uniprot_proteins/uniprotkb_antibiotic_resistance_2026_05_27.tsv \
        --output_csv diamond_annotated_results.csv
    """
}
