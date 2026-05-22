#!/usr/bin/env nextflow

process POEM {
    label 'poem'
    publishDir "results/poem", mode: 'copy'

    input:
    path cds_files

    output:
    path "*"                                   
    path "input.fsa.operon", emit: operon_file 

    script:
    """
    cat ${cds_files.join(' ')} > merged_metagenome_cds.faa

    bash ${projectDir}/POEM_py3k/bin/run_poem.sh \
        -f merged_metagenome_cds.faa \
        -a n \
        -p pka \
        -l y
    """
}