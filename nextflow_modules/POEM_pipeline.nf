#!/usr/bin/env nextflow

process POEM {
    label 'poem'
    
    publishDir "results/poem", mode: 'copy'

    input:
    path fna_files
    path cog_db 

    output:
    path "${fna_files}_output", emit: poem_results_dir 
    path "${fna_files}_output/input.fsa.operon", emit: operon_file, optional: true 

    script:
    """
    if [ "\${CONDA_MANAGEMENT:-}" = "micromamba" ]; then
        conda() {
            if [ "\$1" = "activate" ] || [ "\$1" = "deactivate" ] || [ "\$1" = "hook" ]; then
                :
            else
                micromamba "\$@"
            fi
        }
        export -f conda
    elif [ "\${CONDA_MANAGEMENT:-}" = "conda" ]; then
        source \$(conda info --base)/etc/profile.d/conda.sh
    fi

    mkdir -p tmp_poem_env/database
    ln -s \$(pwd)/${cog_db} tmp_poem_env/database/cog
    
    ln -s ${projectDir}/POEM_py3k/bin tmp_poem_env/bin
    ln -s ${projectDir}/POEM_py3k/lib tmp_poem_env/lib
    ln -s ${projectDir}/POEM_py3k/config tmp_poem_env/config

    bash ${projectDir}/POEM_py3k/bin/run_poem.sh \\
        -f ${fna_files} \\
        -a n \\
        -p pka \\
        -l y
    """
}