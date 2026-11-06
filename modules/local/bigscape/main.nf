#!/usr/bin/env nextflow

process BIGSCAPE {
    label 'process_medium'
    container ''
    label 'bigscape'
    publishDir "results", mode: 'copy'

    input:
    path bgc_dir
    path pfam_db

    output:
    path "bigscape_output"
    path "bigscape_output/output_files/*full.network", emit: bigscape_fullnetwork

script:
def args = task.ext.args ?: ''

"""
mkdir -p bgcs_from_antismash

for f in ${bgc_dir}; do
    cp "\$f" bgcs_from_antismash/
done

bigscape cluster \
    -i bgcs_from_antismash \
    -o bigscape_output \
    -p ${pfam_db}/Pfam-A.hmm \
    ${args}
"""
}