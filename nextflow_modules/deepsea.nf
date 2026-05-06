#!/usr/bin/env nextflow

process DEEPSEA {
    label 'deepsea'
    publishDir "results", mode: 'copy'

    input:
    path cds_files  

    output:
    path "deepsea.table"
    path "merged_deepsea_results.tsv"

    script:
    """
    cat ${cds_files.join(' ')} > merged_metagenome_cds.faa

    python DeepSEA.py run \
        --input merged_metagenome_cds.faa \
        --outname deepsea.table

    tsv_files=\$(ls *.tsv | grep -v merged_deepsea_results.tsv || true)

    if [ -n "\$tsv_files" ]; then
        cat \$tsv_files > merged_deepsea_results.tsv
    else
        touch merged_deepsea_results.tsv
    fi
    """
}
