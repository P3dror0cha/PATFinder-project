#!/usr/bin/env nextflow

process POEM {
    label 'poem'
    publishDir "results/poem", mode: 'copy'

    input:
    path cds_files

    output:
    path "poem_results"

    script:
    """
    # juntar os faa
    cat ${cds_files.join(' ')} > merged_metagenome_cds.faa

    bash ./bin/run_poem.sh \
        -f merged_metagenome_cds.faa \
        -a n \
        -p pka \
        -l y
    """
}