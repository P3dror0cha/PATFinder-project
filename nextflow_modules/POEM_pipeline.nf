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
    ## Check if nextflow is managed by micromamba or conda and set up the environment accordingly
    if [ "\$CONDA_MANAGEMENT" = "micromamba" ]; then
        conda() {
            micromamba "\$@"
        }
        export -f conda

    elif [ "\$CONDA_MANAGEMENT" = "conda" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
    fi

    bash ${projectDir}/POEM_py3k/bin/run_poem.sh \\
        -f ${fna_files} \\
        -a n \\
        -p pka \\
        -l y
    """
}