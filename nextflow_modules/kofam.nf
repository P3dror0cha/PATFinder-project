#!/usr/bin/env nextflow

process KOFAM {
    label 'kofam'
    publishDir "results/kofam", mode: 'copy'

    input:
    path faa_files_from_gbks 

    output:
    path "kofam.out", emit: kofam_output

    script:
    """
    cat ${faa_files_from_gbks} > all_faa_concat.faa

    exec_annotation -o kofam.out \\
        -p ${projectDir}/databases/profiles \\
        -k ${projectDir}/databases/ko_list \\
        all_faa_concat.faa
    """
}