#!/usr/bin/env nextflow

process MGNIFY_DOWNLOAD {
    publishDir "results", mode: 'copy'

    input:
    val max_pages
    val output_prefix
    path envs_done

    output:
    path "metagenomes_MGnify/*.faa", emit: faa
    path "metagenomes_MGnify/*.fna", emit: fna
    path "${output_prefix}_genome_metadata.csv", emit: metadata

    script:
    """
    python3 ${projectDir}/bin/MGnify_download_process.py \
        --max_pages ${max_pages} \
        --output_prefix ${output_prefix} 
    """
}
