#!/usr/bin/env nextflow

process DIAMOND {
    label 'process_medium'
    label 'diamond'
    publishDir "results/diamond", mode: 'copy'

    input:
    path faa_files_from_gbks

    output:
    path "cds_from_gbks_vs_uniprot.tsv", emit: diamond_result

    script: 
    def args = task.ext.args ?: ''

    """
    diamond blastp \\
        -q ${faa_files_from_gbks} \\
        --db ${projectDir}/diamond_database/diamond_biosynthesis_database.dmnd \\
        -o cds_from_gbks_vs_uniprot.tsv \\
        -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore \\
        ${args}
    """
}