process BGC_PLOTS {

    tag "${[gbk_files].flatten().size()} antiSMASH GBKs"

    publishDir "results/bgc_characterization",
        mode: 'copy',
        overwrite: true

    input:
    path gbk_files

    output:
    path "bgc_summary.csv", emit: summary
    path "bgc_cds.csv", emit: cds
    path "bgc_domains.csv", emit: domains
    path "bgc_class_frequency.csv", emit: class_frequency
    path "bgc_functional_composition.csv", emit: functional_composition
    path "gene_map_selection.csv", emit: map_selection
    path "*.png", emit: png_figures
    path "*.svg", emit: svg_figures

    script:
    def gbk_arguments = [gbk_files].flatten().collect { "'${it}'" }.join(' ')

    """
    python3 ${projectDir}/bin/BGC_plots_process.py \
        --gbk ${gbk_arguments} \
        --output-dir . \
        --gene-maps auto \
        --max-maps 6 \
        --max-categories 12 \
        --gene-labels core \
        --formats png svg \
        --dpi 300
    """
}
