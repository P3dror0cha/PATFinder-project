#!/usr/bin/env nextflow

process KOFAM {
    label 'kofam'
    publishDir "results/kofam", mode: 'copy'

    input:
    path faa_files 

    output:
    path "kofam.out", emit: kofam_output

    script:
    """
    cat ${faa_files} > all_faa_concat.faa

    exec_annotation -o kofam.out \\
        -p ${projectDir}/databases/profiles \\
        -k ${projectDir}/databases/ko_list \\
        all_faa_concat.faa
    """
}