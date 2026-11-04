#!/usr/bin/env nextflow

process RENAME_FAA {
    tag "${meta.id}"

    input:
    tuple val(meta), path(faa_file)

    output:
    tuple val(meta), path("${meta.id}.faa"), emit: renamed_faa

    script:
    """
    sed -E "s/^>([^_[:space:]]+)_([0-9]+)/>${meta.id}_\\2|\\1_\\2/" ${faa_file} > ${meta.id}.faa
    """
}
