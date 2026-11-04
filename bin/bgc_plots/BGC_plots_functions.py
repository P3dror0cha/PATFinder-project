"""Extraction and plotting utilities for antiSMASH GenBank files.

The tables created here deliberately use the antiSMASH ``region`` feature as
the BGC unit.  ``protocluster`` features are not counted as independent BGCs.
"""

from __future__ import annotations

import os
import re
import tempfile
import warnings
from pathlib import Path
from typing import Iterable, Sequence

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "patfinder-matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrow, Patch


FUNCTION_CATEGORIES = [
    "biosynthetic",
    "biosynthetic-additional",
    "regulatory",
    "transport",
    "resistance",
    "other",
]

FUNCTION_COLORS = {
    "biosynthetic": "#D55E00",
    "biosynthetic-additional": "#E69F00",
    "regulatory": "#0072B2",
    "transport": "#009E73",
    "resistance": "#CC79A7",
    "other": "#B3B3B3",
}

SUMMARY_COLUMNS = [
    "bgc_id",
    "source_file",
    "record_id",
    "organism",
    "region_number",
    "bgc_class",
    "bgc_classes",
    "n_classes",
    "is_hybrid",
    "contig_edge",
    "region_start",
    "region_end",
    "bgc_length_bp",
    "bgc_length_kb",
    "n_cds",
    "cds_per_10kb",
    "n_biosynthetic",
    "n_biosynthetic_additional",
    "n_regulatory",
    "n_transport",
    "n_resistance",
    "n_other",
]

CDS_COLUMNS = [
    "bgc_id",
    "gene_id",
    "locus_tag",
    "gene",
    "product",
    "gene_kind",
    "start",
    "end",
    "length_bp",
    "strand",
    "translation_length_aa",
]

DOMAIN_COLUMNS = ["bgc_id", "gene_id", "domain_type", "domain_annotation"]


def _qualifier_values(feature, key: str) -> list[str]:
    """Return a qualifier as a clean list, regardless of its original type."""

    value = feature.qualifiers.get(key, [])
    if value is None:
        return []
    if isinstance(value, str):
        value = [value]
    return [str(item).strip() for item in value if str(item).strip()]


def _first_qualifier(feature, key: str, default: str = "") -> str:
    values = _qualifier_values(feature, key)
    return values[0] if values else default


def _as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _normalise_gene_kind(feature) -> str:
    """Normalise antiSMASH functional labels without inventing annotations."""

    candidates = _qualifier_values(feature, "gene_kind")
    candidates.extend(_qualifier_values(feature, "gene_functions"))
    joined = " ".join(candidates).lower()

    if "biosynthetic-additional" in joined:
        return "biosynthetic-additional"
    if "biosynthetic" in joined:
        return "biosynthetic"
    if "regulat" in joined:
        return "regulatory"
    if "transport" in joined:
        return "transport"
    if "resistance" in joined:
        return "resistance"
    return "other"


def _normalised_bgc_id(file_path: str | Path, region_number: str) -> str:
    """Convert ``sample.region001.gbk`` to ``sample_001``."""

    stem = Path(file_path).stem
    match = re.fullmatch(r"(.+)\.region(\d+)", stem, flags=re.IGNORECASE)
    if match:
        return f"{match.group(1)}_{match.group(2)}"

    clean_region = str(region_number).strip()
    if clean_region.isdigit():
        clean_region = clean_region.zfill(3)
    return f"{stem}_{clean_region}"


def _location_bounds(feature) -> tuple[int, int, int]:
    start = int(feature.location.start)
    end = int(feature.location.end)
    strand = feature.location.strand
    return start, end, int(strand) if strand in (-1, 1) else 0


def extract_antismash_tables(
    gbk_files: Iterable[str | Path],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Extract one BGC row per region, one CDS row per gene and domain rows.

    Parameters
    ----------
    gbk_files:
        antiSMASH GenBank region files.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        ``bgc_summary``, ``bgc_cds`` and ``bgc_domains``.
    """

    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise RuntimeError(
            "Biopython is required to read GBK files. Install it with "
            "'conda install biopython' in the environment used by the process."
        ) from exc

    paths = [Path(path) for path in gbk_files]
    if not paths:
        raise ValueError("No GBK files were supplied.")

    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("GBK files not found: " + ", ".join(missing))

    summary_rows: list[dict] = []
    cds_rows: list[dict] = []
    domain_rows: list[dict] = []

    for file_path in sorted(paths, key=lambda path: str(path)):
        records_found = 0
        for record in SeqIO.parse(str(file_path), "genbank"):
            records_found += 1
            region_features = [
                feature for feature in record.features if feature.type == "region"
            ]
            if not region_features:
                warnings.warn(f"No region feature found in {file_path.name}; skipped.")
                continue

            organism = str(record.annotations.get("organism", ""))
            for region_index, region in enumerate(region_features, start=1):
                region_number = _first_qualifier(
                    region, "region_number", str(region_index).zfill(3)
                )
                bgc_id = _normalised_bgc_id(file_path, region_number)
                region_start, region_end, _ = _location_bounds(region)
                bgc_length_bp = region_end - region_start

                classes = _qualifier_values(region, "product") or ["unknown"]
                classes = list(dict.fromkeys(classes))
                bgc_class = " + ".join(classes)
                contig_edge = _as_bool(_first_qualifier(region, "contig_edge", "false"))

                region_cds = []
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    cds_start, cds_end, strand = _location_bounds(feature)
                    if cds_start >= region_end or cds_end <= region_start:
                        continue

                    locus_tag = _first_qualifier(feature, "locus_tag")
                    gene_name = _first_qualifier(feature, "gene")
                    gene_id = locus_tag or gene_name or f"{bgc_id}_cds_{len(region_cds) + 1:04d}"
                    product = _first_qualifier(feature, "product", "unannotated protein")
                    gene_kind = _normalise_gene_kind(feature)
                    translation = _first_qualifier(feature, "translation")

                    cds_row = {
                        "bgc_id": bgc_id,
                        "gene_id": gene_id,
                        "locus_tag": locus_tag,
                        "gene": gene_name,
                        "product": product,
                        "gene_kind": gene_kind,
                        "start": cds_start,
                        "end": cds_end,
                        "length_bp": cds_end - cds_start,
                        "strand": strand,
                        "translation_length_aa": len(translation) if translation else np.nan,
                    }
                    region_cds.append(cds_row)
                    cds_rows.append(cds_row)

                    for qualifier in ("sec_met_domain", "PFAM_domain", "NRPS_PKS"):
                        for annotation in _qualifier_values(feature, qualifier):
                            domain_rows.append(
                                {
                                    "bgc_id": bgc_id,
                                    "gene_id": gene_id,
                                    "domain_type": qualifier,
                                    "domain_annotation": annotation,
                                }
                            )

                counts = pd.Series(
                    [row["gene_kind"] for row in region_cds], dtype="object"
                ).value_counts()
                n_cds = len(region_cds)
                summary_rows.append(
                    {
                        "bgc_id": bgc_id,
                        "source_file": file_path.name,
                        "record_id": str(record.id),
                        "organism": organism,
                        "region_number": region_number,
                        "bgc_class": bgc_class,
                        "bgc_classes": ";".join(classes),
                        "n_classes": len(classes),
                        "is_hybrid": len(classes) > 1,
                        "contig_edge": contig_edge,
                        "region_start": region_start,
                        "region_end": region_end,
                        "bgc_length_bp": bgc_length_bp,
                        "bgc_length_kb": bgc_length_bp / 1000.0,
                        "n_cds": n_cds,
                        "cds_per_10kb": (
                            n_cds / (bgc_length_bp / 10_000.0)
                            if bgc_length_bp > 0
                            else np.nan
                        ),
                        "n_biosynthetic": int(counts.get("biosynthetic", 0)),
                        "n_biosynthetic_additional": int(
                            counts.get("biosynthetic-additional", 0)
                        ),
                        "n_regulatory": int(counts.get("regulatory", 0)),
                        "n_transport": int(counts.get("transport", 0)),
                        "n_resistance": int(counts.get("resistance", 0)),
                        "n_other": int(counts.get("other", 0)),
                    }
                )

        if records_found == 0:
            warnings.warn(f"No GenBank record could be parsed from {file_path.name}.")

    summary = pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS)
    cds = pd.DataFrame(cds_rows, columns=CDS_COLUMNS)
    domains = pd.DataFrame(domain_rows, columns=DOMAIN_COLUMNS)

    if summary.empty:
        raise ValueError("No antiSMASH region features were found in the supplied GBKs.")

    duplicated = summary.loc[summary["bgc_id"].duplicated(keep=False), "bgc_id"].unique()
    if len(duplicated):
        raise ValueError(
            "Duplicate bgc_id values were generated: " + ", ".join(map(str, duplicated))
        )

    return summary, cds, domains


def build_class_frequency(summary: pd.DataFrame) -> pd.DataFrame:
    classes = summary[["bgc_id", "bgc_classes"]].copy()
    classes["bgc_class"] = classes["bgc_classes"].str.split(";")
    classes = classes.explode("bgc_class")
    classes["bgc_class"] = classes["bgc_class"].str.strip()
    frequency = (
        classes.groupby("bgc_class", dropna=False)["bgc_id"]
        .nunique()
        .rename("n_bgcs")
        .reset_index()
        .sort_values(["n_bgcs", "bgc_class"], ascending=[False, True])
        .reset_index(drop=True)
    )
    frequency["percentage_bgcs"] = 100 * frequency["n_bgcs"] / summary["bgc_id"].nunique()
    return frequency


def build_functional_composition(
    summary: pd.DataFrame, cds: pd.DataFrame
) -> pd.DataFrame:
    columns = ["bgc_class", "gene_kind", "n_cds", "percentage_within_class"]
    if cds.empty:
        return pd.DataFrame(columns=columns)

    joined = cds.merge(summary[["bgc_id", "bgc_class"]], on="bgc_id", how="left")
    composition = (
        joined.groupby(["bgc_class", "gene_kind"], dropna=False)
        .size()
        .rename("n_cds")
        .reset_index()
    )
    totals = composition.groupby("bgc_class")["n_cds"].transform("sum")
    composition["percentage_within_class"] = 100 * composition["n_cds"] / totals
    return composition[columns].sort_values(["bgc_class", "gene_kind"])


def _save_figure(
    figure: plt.Figure,
    output_stem: Path,
    formats: Sequence[str],
    dpi: int,
) -> list[Path]:
    paths = []
    for extension in formats:
        path = output_stem.with_suffix(f".{extension.lower()}")
        figure.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        paths.append(path)
    plt.close(figure)
    return paths


def _top_categories(summary: pd.DataFrame, max_categories: int) -> list[str]:
    counts = summary["bgc_class"].value_counts()
    return counts.head(max_categories).index.tolist()


def _draw_class_distribution(axis, frequency: pd.DataFrame) -> None:
    plot_data = frequency.sort_values("n_bgcs", ascending=True)
    axis.barh(plot_data["bgc_class"], plot_data["n_bgcs"], color="#4472C4")
    axis.set_xlabel("Número de BGCs")
    axis.set_ylabel("Classe antiSMASH")
    axis.set_title("A. Distribuição das classes")
    axis.grid(axis="x", alpha=0.2)


def _draw_size_by_class(axis, summary: pd.DataFrame, max_categories: int) -> None:
    categories = _top_categories(summary, max_categories)
    values = [
        summary.loc[summary["bgc_class"] == category, "bgc_length_kb"].dropna().to_numpy()
        for category in categories
    ]
    # ``labels`` keeps compatibility with the older matplotlib releases that
    # are commonly bundled in bioinformatics Conda environments.
    axis.boxplot(values, labels=categories, vert=False, showfliers=True)
    axis.set_xlabel("Tamanho do BGC (kb)")
    axis.set_ylabel("Classe completa da região")
    axis.set_title("B. Tamanho por classe")
    axis.grid(axis="x", alpha=0.2)


def _draw_length_vs_cds(axis, summary: pd.DataFrame, max_categories: int) -> None:
    categories = _top_categories(summary, max_categories)
    plot_data = summary.copy()
    plot_data["plot_class"] = np.where(
        plot_data["bgc_class"].isin(categories), plot_data["bgc_class"], "Outras"
    )
    plot_categories = categories + (["Outras"] if "Outras" in plot_data["plot_class"].values else [])
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(plot_categories), 1)))
    for category, color in zip(plot_categories, colors):
        group = plot_data[plot_data["plot_class"] == category]
        axis.scatter(
            group["bgc_length_kb"],
            group["n_cds"],
            s=38,
            alpha=0.75,
            edgecolor="white",
            linewidth=0.4,
            color=color,
            label=category,
        )
    axis.set_xlabel("Tamanho do BGC (kb)")
    axis.set_ylabel("Número de CDS")
    axis.set_title("C. Tamanho × número de CDS")
    axis.grid(alpha=0.2)
    axis.legend(
        fontsize=7,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
    )


def _draw_functional_composition(
    axis, summary: pd.DataFrame, cds: pd.DataFrame, max_categories: int
) -> None:
    categories = _top_categories(summary, max_categories)
    if cds.empty:
        axis.text(0.5, 0.5, "Nenhum CDS encontrado", ha="center", va="center")
        axis.set_axis_off()
        return

    joined = cds.merge(summary[["bgc_id", "bgc_class"]], on="bgc_id", how="left")
    joined = joined[joined["bgc_class"].isin(categories)]
    counts = pd.crosstab(joined["bgc_class"], joined["gene_kind"])
    counts = counts.reindex(index=categories, fill_value=0)
    proportions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0) * 100

    left = np.zeros(len(proportions))
    for gene_kind in FUNCTION_CATEGORIES:
        values = (
            proportions[gene_kind].fillna(0).to_numpy()
            if gene_kind in proportions
            else np.zeros(len(proportions))
        )
        axis.barh(
            proportions.index,
            values,
            left=left,
            color=FUNCTION_COLORS[gene_kind],
            label=gene_kind,
        )
        left += values
    axis.set_xlim(0, 100)
    axis.set_xlabel("CDS da classe (%)")
    axis.set_ylabel("Classe completa da região")
    axis.set_title("D. Composição funcional")
    axis.legend(
        fontsize=7,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
    )
    axis.grid(axis="x", alpha=0.2)


def generate_quantitative_plots(
    summary: pd.DataFrame,
    cds: pd.DataFrame,
    output_dir: str | Path,
    formats: Sequence[str] = ("png", "svg"),
    dpi: int = 300,
    max_categories: int = 12,
) -> list[Path]:
    """Create four individual figures and a combined 2×2 panel."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    frequency = build_class_frequency(summary)
    created: list[Path] = []

    specifications = [
        ("bgc_class_distribution", _draw_class_distribution, (frequency,)),
        ("bgc_size_by_class", _draw_size_by_class, (summary, max_categories)),
        ("bgc_length_vs_cds", _draw_length_vs_cds, (summary, max_categories)),
        (
            "bgc_functional_composition",
            _draw_functional_composition,
            (summary, cds, max_categories),
        ),
    ]
    for name, draw_function, arguments in specifications:
        figure, axis = plt.subplots(figsize=(10, 6))
        draw_function(axis, *arguments)
        figure.tight_layout()
        created.extend(_save_figure(figure, output_dir / name, formats, dpi))

    panel, axes = plt.subplots(2, 2, figsize=(18, 13))
    _draw_class_distribution(axes[0, 0], frequency)
    _draw_size_by_class(axes[0, 1], summary, max_categories)
    _draw_length_vs_cds(axes[1, 0], summary, max_categories)
    _draw_functional_composition(axes[1, 1], summary, cds, max_categories)
    panel.tight_layout()
    created.extend(_save_figure(panel, output_dir / "bgc_characterization_panel", formats, dpi))
    return created


def select_gene_maps(
    summary: pd.DataFrame,
    mode: str = "auto",
    requested_ids: Sequence[str] | None = None,
    max_maps: int = 6,
) -> pd.DataFrame:
    """Select maps explicitly or by proximity to each class median length.

    Automatic selections are *not* GCF medoids.  A BiG-SCAPE distance matrix is
    required for a genuine medoid selection.
    """

    columns = ["bgc_id", "bgc_class", "selection_criterion", "class_median_kb"]
    requested_ids = list(requested_ids or [])
    if mode == "none":
        return pd.DataFrame(columns=columns)

    if mode == "selected":
        if not requested_ids:
            raise ValueError("--gene-maps selected requires --map-bgc-ids.")
        unknown = sorted(set(requested_ids) - set(summary["bgc_id"]))
        if unknown:
            raise ValueError("Unknown BGC IDs requested for maps: " + ", ".join(unknown))
        selected = summary.set_index("bgc_id").loc[requested_ids].reset_index()
        selected["selection_criterion"] = "explicitly_selected"
        selected["class_median_kb"] = np.nan
    elif mode == "all":
        selected = summary.sort_values("bgc_id").head(max_maps).copy()
        selected["selection_criterion"] = "all_limited_by_max_maps"
        selected["class_median_kb"] = np.nan
    elif mode == "auto":
        class_order = summary["bgc_class"].value_counts().index.tolist()
        rows = []
        for bgc_class in class_order:
            group = summary[summary["bgc_class"] == bgc_class].copy()
            median = float(group["bgc_length_kb"].median())
            group["distance_to_median"] = (group["bgc_length_kb"] - median).abs()
            representative = group.sort_values(
                ["distance_to_median", "bgc_id"]
            ).iloc[0].copy()
            representative["selection_criterion"] = "closest_to_class_median_length"
            representative["class_median_kb"] = median
            rows.append(representative)
        selected = pd.DataFrame(rows).head(max_maps)
    else:
        raise ValueError(f"Unsupported gene-map mode: {mode}")

    if len(selected) > max_maps:
        raise ValueError(
            f"{len(selected)} BGCs were selected, but --max-maps is {max_maps}."
        )
    return selected[columns].reset_index(drop=True)


def _safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def plot_gene_map(
    summary_row: pd.Series,
    cds: pd.DataFrame,
    output_dir: str | Path,
    formats: Sequence[str] = ("png", "svg"),
    dpi: int = 300,
    label_mode: str = "core",
) -> list[Path]:
    """Draw a uniform CDS-arrow map for one selected BGC."""

    bgc_id = str(summary_row["bgc_id"])
    genes = cds[cds["bgc_id"] == bgc_id].sort_values(["start", "end"])
    if genes.empty:
        warnings.warn(f"No CDS available for gene map {bgc_id}; skipped.")
        return []

    region_start = int(summary_row["region_start"])
    region_length = int(summary_row["bgc_length_bp"])
    figure_width = max(12.0, min(22.0, region_length / 5000.0))
    figure, axis = plt.subplots(figsize=(figure_width, 3.2))
    axis.hlines(0, 0, region_length / 1000.0, color="#666666", linewidth=0.8)

    present_categories = []
    for _, gene in genes.iterrows():
        start_kb = (int(gene["start"]) - region_start) / 1000.0
        end_kb = (int(gene["end"]) - region_start) / 1000.0
        length_kb = max(end_kb - start_kb, 0.02)
        strand = int(gene["strand"])
        kind = str(gene["gene_kind"])
        if kind not in FUNCTION_COLORS:
            kind = "other"
        present_categories.append(kind)

        if strand == -1:
            x, dx = end_kb, -length_kb
        else:
            x, dx = start_kb, length_kb
        head_length = min(max(length_kb * 0.28, 0.08), length_kb * 0.65)
        arrow = FancyArrow(
            x,
            0,
            dx,
            0,
            width=0.34,
            head_width=0.62,
            head_length=head_length,
            length_includes_head=True,
            facecolor=FUNCTION_COLORS[kind],
            edgecolor="#333333",
            linewidth=0.45,
        )
        axis.add_patch(arrow)

        should_label = label_mode == "all" or (
            label_mode == "core" and kind in {"biosynthetic", "biosynthetic-additional"}
        )
        if should_label:
            axis.text(
                (start_kb + end_kb) / 2,
                0.48,
                str(gene["locus_tag"] or gene["gene_id"]),
                rotation=45,
                ha="left",
                va="bottom",
                fontsize=6,
            )

    axis.set_xlim(-0.5, region_length / 1000.0 + 0.5)
    axis.set_ylim(-0.85, 1.25)
    axis.set_yticks([])
    axis.set_xlabel("Posição relativa no BGC (kb)")
    axis.set_title(f"{bgc_id} — {summary_row['bgc_class']}")
    for spine_name in ("left", "right", "top"):
        axis.spines[spine_name].set_visible(False)
    axis.grid(axis="x", alpha=0.15)

    legend_categories = [
        category for category in FUNCTION_CATEGORIES if category in present_categories
    ]
    handles = [
        Patch(facecolor=FUNCTION_COLORS[category], edgecolor="#333333", label=category)
        for category in legend_categories
    ]
    axis.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7)
    figure.tight_layout()
    return _save_figure(
        figure,
        Path(output_dir) / f"gene_map_{_safe_filename(bgc_id)}",
        formats,
        dpi,
    )


def generate_gene_maps(
    summary: pd.DataFrame,
    cds: pd.DataFrame,
    selection: pd.DataFrame,
    output_dir: str | Path,
    formats: Sequence[str] = ("png", "svg"),
    dpi: int = 300,
    label_mode: str = "core",
) -> list[Path]:
    created: list[Path] = []
    summary_by_id = summary.set_index("bgc_id")
    for bgc_id in selection["bgc_id"].tolist():
        summary_row = summary_by_id.loc[bgc_id].copy()
        summary_row["bgc_id"] = bgc_id
        created.extend(
            plot_gene_map(
                summary_row,
                cds,
                output_dir,
                formats=formats,
                dpi=dpi,
                label_mode=label_mode,
            )
        )
    return created