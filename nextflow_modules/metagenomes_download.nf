#!/usr/bin/env nextflow

process MGNIFY_DOWNLOAD {
    publishDir "results", mode: 'copy'

    input:
    val max_pages
    val output_prefix

    output:
    path "metagenomes_MGnify/*.faa", emit: faa
    path "metagenomes_MGnify/*.fna", emit: fna

    script:
    """
    python3 ${projectDir}/bin/mgnify_pipeline.py \
        --max_pages ${max_pages} \
        --output_prefix ${output_prefix} 
    """
}
