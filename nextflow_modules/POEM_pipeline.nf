#!/usr/bin/env nextflow

process POEM {
    label 'poem'
    
    publishDir "results/poem", mode: 'copy'

    input:
    path fna_files

    output:
    path "${fna_files}_output", emit: poem_results_dir 
    path "${fna_files}_output/input.fsa.operon", emit: operon_file, optional: true 

    script:
    """
    bash ${projectDir}/POEM_py3k/bin/run_poem.sh \\
        -f ${fna_files} \\
        -a n \\
        -p pka \\
        -l y
    """
}