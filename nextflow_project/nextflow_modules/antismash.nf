#!/usr/bin/env nextflow

process ANTISMASH {
    label 'antismash'
    publishDir "results", mode: 'copy'

    input:
    path fasta_file

    output:
    path "antismash/${fasta_file.baseName}/*.gbk", emit: gbk

    script: 
    """

    mkdir -p antismash
    
    antismash "$fasta_file" \
        --genefinding-tool prodigal \
        --cb-general \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --tigrfam \
        --output-dir "./antismash/${fasta_file.baseName}"
    """
}
