"""Command-line process for antiSMASH BGC tables and figures."""

from __future__ import annotations

import argparse
from pathlib import Path

from BGC_plots_functions import (
    build_class_frequency,
    build_functional_composition,
    extract_antismash_tables,
    generate_gene_maps,
    generate_quantitative_plots,
    select_gene_maps,
)


def bgc_plots_process(
    gbk_files: list[str],
    output_dir: str,
    gene_maps: str = "auto",
    map_bgc_ids: list[str] | None = None,
    max_maps: int = 6,
    max_categories: int = 12,
    gene_labels: str = "core",
    formats: tuple[str, ...] = ("png", "svg"),
    dpi: int = 300,
) -> dict[str, object]:
    """Create antiSMASH-derived tables, quantitative plots and optional maps."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    summary, cds, domains = extract_antismash_tables(gbk_files)
    class_frequency = build_class_frequency(summary)
    functional_composition = build_functional_composition(summary, cds)
    map_selection = select_gene_maps(
        summary,
        mode=gene_maps,
        requested_ids=map_bgc_ids,
        max_maps=max_maps,
    )

    table_paths = {
        "summary": output_path / "bgc_summary.csv",
        "cds": output_path / "bgc_cds.csv",
        "domains": output_path / "bgc_domains.csv",
        "class_frequency": output_path / "bgc_class_frequency.csv",
        "functional_composition": output_path / "bgc_functional_composition.csv",
        "map_selection": output_path / "gene_map_selection.csv",
    }
    summary.to_csv(table_paths["summary"], index=False)
    cds.to_csv(table_paths["cds"], index=False)
    domains.to_csv(table_paths["domains"], index=False)
    class_frequency.to_csv(table_paths["class_frequency"], index=False)
    functional_composition.to_csv(table_paths["functional_composition"], index=False)
    map_selection.to_csv(table_paths["map_selection"], index=False)

    figure_paths = generate_quantitative_plots(
        summary,
        cds,
        output_path,
        formats=formats,
        dpi=dpi,
        max_categories=max_categories,
    )
    figure_paths.extend(
        generate_gene_maps(
            summary,
            cds,
            map_selection,
            output_path,
            formats=formats,
            dpi=dpi,
            label_mode=gene_labels,
        )
    )

    return {
        "n_bgcs": len(summary),
        "n_cds": len(cds),
        "n_domains": len(domains),
        "n_gene_maps": len(map_selection),
        "tables": table_paths,
        "figures": figure_paths,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract antiSMASH BGC characteristics from GBKs and generate tables, "
            "quantitative figures and optional gene maps."
        )
    )
    parser.add_argument(
        "-g", "--gbk", nargs="+", required=True, help="List of antiSMASH GBK files"
    )
    parser.add_argument(
        "-o", "--output-dir", default=".", help="Directory for CSV and figure outputs"
    )
    parser.add_argument(
        "--gene-maps",
        choices=("none", "auto", "selected", "all"),
        default="auto",
        help=(
            "Map selection: auto chooses the BGC closest to the median length of "
            "each class; selected uses --map-bgc-ids"
        ),
    )
    parser.add_argument(
        "--map-bgc-ids",
        nargs="*",
        default=[],
        help="BGC IDs used when --gene-maps selected",
    )
    parser.add_argument(
        "--max-maps", type=int, default=6, help="Maximum number of gene maps"
    )
    parser.add_argument(
        "--max-categories",
        type=int,
        default=12,
        help="Maximum number of complete BGC classes displayed in crowded plots",
    )
    parser.add_argument(
        "--gene-labels",
        choices=("none", "core", "all"),
        default="core",
        help="Which CDS labels to display in gene maps",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "svg", "pdf"),
        default=("png", "svg"),
        help="Figure formats",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Raster figure resolution")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.max_maps < 1:
        raise ValueError("--max-maps must be at least 1.")
    if args.max_categories < 1:
        raise ValueError("--max-categories must be at least 1.")
    if args.dpi < 72:
        raise ValueError("--dpi must be at least 72.")

    results = bgc_plots_process(
        gbk_files=args.gbk,
        output_dir=args.output_dir,
        gene_maps=args.gene_maps,
        map_bgc_ids=args.map_bgc_ids,
        max_maps=args.max_maps,
        max_categories=args.max_categories,
        gene_labels=args.gene_labels,
        formats=tuple(args.formats),
        dpi=args.dpi,
    )
    print(
        "BGC characterization completed: "
        f"{results['n_bgcs']} BGCs, {results['n_cds']} CDS, "
        f"{results['n_domains']} domain annotations and "
        f"{results['n_gene_maps']} gene maps."
    )


if __name__ == "__main__":
    main()