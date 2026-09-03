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
    if [ -d "/home/mambauser/app" ]; then
        POEM_SOURCE="/home/mambauser/app"
    else
        POEM_SOURCE="${projectDir}/POEM_py3k"
    fi

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

    if ! command -v conda &> /dev/null; then
        conda() { return 0; }
        export -f conda
    fi

    mkdir -p tmp_poem_env/database
    ln -s \$(readlink -f ${cog_db}) tmp_poem_env/database/cog
    
    ln -s \$POEM_SOURCE/bin tmp_poem_env/bin
    ln -s \$POEM_SOURCE/lib tmp_poem_env/lib
    ln -s \$POEM_SOURCE/config tmp_poem_env/config

    export PYTHONPATH=\$(pwd)/tmp_poem_env/lib:\${PYTHONPATH:-}

    bash tmp_poem_env/bin/run_poem.sh \\
        -f ${fna_files} \\
        -a n \\
        -p pro \\
        -l y
    """
}