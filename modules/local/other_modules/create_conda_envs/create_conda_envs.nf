#!/usr/bin/env nextflow

process CREATE_CONDA_ENVS {
    label 'create_conda_envs'

    output:
    path "envs_created.done", emit: done

    script:
    """
    if [ "\$CONDA_MANAGEMENT" = "micromamba" ]; then
        CMD="micromamba"
    elif [ "\$CONDA_MANAGEMENT" = "conda" ]; then
        CMD="conda"
    else
        echo "CONDA_MANAGEMENT not set, skipping environment creation"
        touch envs_created.done
        exit 0
    fi

    for env_yml in ${projectDir}/envs/antismash_bigscape_updated.yml \\
                   ${projectDir}/envs/deepsea_project_updated.yml \\
                   ${projectDir}/envs/diamond_kofam_updated.yml \\
                   ${projectDir}/envs/poem_py3_updated.yml \\
                   ${projectDir}/envs/amrfinder.yml; do

        env_name=\$(grep '^name:' "\$env_yml" | awk '{print \$2}')
        if \$CMD env list | grep -q "\$env_name"; then
            echo "Environment '\$env_name' already exists, skipping."
        else
            echo "Creating environment '\$env_name'..."
            \$CMD env create -f "\$env_yml"
        fi
    done

    touch envs_created.done
    """
}