#!/usr/bin/env nextflow

process RENAME_GBKS {
publishDir "results/renamed_gbk", mode: 'copy'

    input:
    tuple path(original_name), val(new_name)

    output:
    path("${new_name}"), emit: renamed_gbk

    script:
    """
    mv ${original_name} ${new_name}
    """
}