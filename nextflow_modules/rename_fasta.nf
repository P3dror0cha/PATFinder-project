#!/usr/bin/env nextflow

process RENAME_FASTA {
    tag "${meta.id}"

    input:
    tuple val(meta), path(original_fasta)
    path envs_done

    output:
    tuple val(meta), path("${meta.id}.fasta")

    script:
    """ 
    awk '/^>/{ 
        original=\$0; 
        sub(/^>/, "", original); 
        print ">${meta.id}_" ++i " ID_Original:" original; 
        next 
    } {print}' ${original_fasta} > ${meta.id}.fasta
    """
}