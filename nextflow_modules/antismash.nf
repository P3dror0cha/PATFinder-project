#!/usr/bin/env nextflow

process ANTISMASH {
    label 'antismash'
    publishDir "results/antismash", mode: 'copy'

    input:
    tuple val(meta), path(fasta_file)

    output:
    tuple val(meta), path("${meta.id}/*region*.gbk"), emit: gbk

    script: 
    """
    antismash "$fasta_file" \
        --genefinding-tool prodigal \
        --cb-general \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --tigrfam \
        --output-dir "${meta.id}" 
    """
}
