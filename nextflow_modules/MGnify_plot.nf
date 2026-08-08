process MGNIFY_PLOTS {
    tag "MGnify metadata plots"
    publishDir "results/MGnify_plots", mode: "copy"

    input:
    path metadata_csv
    val output_prefix

    output:
    path "MGnify_plots/*.png", emit: plots

    script:
    """
    python ${projectDir}/bin/MGnify_plot_process.py \
        ${metadata_csv} \
        --output_dir MGnify_plots \
        --output_prefix aquatic
    """
}