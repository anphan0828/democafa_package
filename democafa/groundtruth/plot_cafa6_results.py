#!/usr/bin/env python3
"""Plot CAFA6 evaluator weighted F-measure results by subset and GO aspect."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


EVALUATION_FILE = "evaluation_best_f_micro_w.tsv"
SUBSET_DIRS = {
    "NK": "test_cafaeval_NK",
    "LK": "test_cafaeval_LK",
    "PK": "test_cafaeval_PK",
}
SUBSET_LABELS = {
    "NK": "No knowledge",
    "LK": "Limited knowledge",
    "PK": "Partial knowledge",
}
ASPECTS = {
    "biological_process": "Biological process",
    "cellular_component": "Cellular component",
    "molecular_function": "Molecular function",
}


def strip_prediction_suffix(filename: str) -> str:
    """Return a compact prediction label from a filename/path."""
    path = Path(str(filename)).name
    while Path(path).suffix:
        path = Path(path).stem
    return path


def extract_submission_id(filename: str) -> str:
    """Return the submission ID portion of an evaluator filename."""
    return strip_prediction_suffix(filename)


def read_evaluator_scores(result_file: str | Path, subset: str) -> pd.DataFrame:
    """Read one CAFA-evaluator ``evaluation_best_f_micro_w.tsv`` file."""
    result_file = Path(result_file)
    if not result_file.exists():
        raise FileNotFoundError(f"Result file not found: {result_file}")

    df = pd.read_csv(result_file, sep="\t", dtype={"filename": str, "ns": str})
    required = {"filename", "ns", "f_micro_w"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{result_file} is missing columns: {sorted(missing)}")

    df["f_micro_w"] = pd.to_numeric(df["f_micro_w"], errors="coerce").fillna(0.0)
    scores = (
        df.groupby(["filename", "ns"], as_index=False)
        .agg(f_micro_w=("f_micro_w", "mean"))
        .assign(
            subset=subset,
            subset_label=SUBSET_LABELS.get(subset, subset),
            submission_id=lambda frame: frame["filename"].map(extract_submission_id),
            prediction_label=lambda frame: frame["filename"].map(strip_prediction_suffix),
        )
    )
    return scores


def read_score_metadata(metadata_file: str | Path) -> dict[str, str]:
    """Read score metadata CSV and return ``SubmissionId`` to ``TeamName``."""
    metadata = pd.read_csv(metadata_file, dtype={"SubmissionId": str, "TeamName": str})
    required = {"SubmissionId", "TeamName"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"{metadata_file} is missing columns: {sorted(missing)}")

    metadata = metadata.dropna(subset=["SubmissionId"]).copy()
    metadata["SubmissionId"] = metadata["SubmissionId"].astype(str)
    metadata["TeamName"] = metadata["TeamName"].fillna("").astype(str)
    return metadata.drop_duplicates("SubmissionId", keep="last").set_index("SubmissionId")["TeamName"].to_dict()


def apply_score_metadata(scores: pd.DataFrame, metadata_file: str | Path | None = None) -> pd.DataFrame:
    """Replace filename labels with team names when metadata contains a match."""
    if metadata_file is None:
        return scores

    id_to_team = read_score_metadata(metadata_file)
    scores = scores.copy()
    scores["prediction_label"] = scores.apply(
        lambda row: id_to_team.get(str(row["submission_id"]), "") or row["prediction_label"],
        axis=1,
    )
    return scores


def collect_evaluator_scores(data_dir: str | Path, subset_dirs: dict[str, str] | None = None) -> pd.DataFrame:
    """Collect weighted F-measure scores from NK/LK/PK evaluator folders."""
    data_dir = Path(data_dir)
    subset_dirs = subset_dirs or SUBSET_DIRS
    frames = []
    for subset, dirname in subset_dirs.items():
        result_file = data_dir / dirname / EVALUATION_FILE
        frames.append(read_evaluator_scores(result_file, subset=subset))
    return pd.concat(frames, ignore_index=True)


def write_score_table(scores: pd.DataFrame, output_file: str | Path) -> None:
    """Write the long-form score table used for plotting."""
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(output_file, sep="\t", index=False)


def plot_cafa6_results(
    scores: pd.DataFrame,
    output_file: str | Path,
    figsize: tuple[float, float] | None = None,
) -> tuple[plt.Figure, list[list[plt.Axes]]]:
    """Plot a 3x3 bar grid: rows are NK/LK/PK and columns are GO aspects."""
    if figsize is None:
        n_predictions = max(scores["filename"].nunique(), 1)
        figsize = (max(10.0, 0.75 * n_predictions * 3), 13.0)

    subsets = list(SUBSET_DIRS.keys())
    aspects = list(ASPECTS.keys())
    cmap = plt.get_cmap("Dark2")
    color_lookup = {
        filename: cmap(index % cmap.N)
        for index, filename in enumerate(sorted(scores["filename"].unique()))
    }
    y_max = max(1.0, scores["f_micro_w"].max() * 1.4)

    fig, axes = plt.subplots(len(subsets), len(aspects), figsize=figsize, sharey=True, squeeze=False)
    for row_index, subset in enumerate(subsets):
        for col_index, aspect in enumerate(aspects):
            ax = axes[row_index][col_index]
            panel = scores[(scores["subset"] == subset) & (scores["ns"] == aspect)].copy()
            panel = panel.sort_values(["f_micro_w", "prediction_label"], ascending=[False, True]).reset_index(drop=True)
            team_label_position = "above" if subset == "PK" else "inside"
            _plot_result_panel(ax, panel, color_lookup, y_max=y_max, team_label_position=team_label_position)

            ax.set_title(f"{ASPECTS[aspect]}: {SUBSET_LABELS[subset]}", fontsize=12)
            if col_index == 0:
                ax.set_ylabel("f_micro_w")
            else:
                ax.set_ylabel("")
            if row_index == len(subsets) - 1:
                ax.set_xlabel("Teams")
            else:
                ax.set_xlabel("")

    fig.tight_layout()
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    return fig, axes.tolist()


def _plot_result_panel(
    ax: plt.Axes,
    panel: pd.DataFrame,
    color_lookup: dict[str, tuple[float, float, float, float]],
    y_max: float,
    team_label_position: str = "inside",
) -> None:
    if panel.empty:
        ax.text(0.5, 0.5, "No results", ha="center", va="center", transform=ax.transAxes)
        ax.set_ylim(0, y_max)
        ax.set_xticks([])
        return

    x_positions = range(len(panel))
    bar_colors = [color_lookup[filename] for filename in panel["filename"]]
    bars = ax.bar(x_positions, panel["f_micro_w"], color=bar_colors, edgecolor="white", linewidth=0.8)
    ax.set_ylim(0, y_max)
    ax.set_xticks(list(x_positions))
    ax.set_xticklabels([])
    ax.grid(axis="y", color="0.88", linewidth=0.8)
    ax.set_axisbelow(True)

    for bar, row in zip(bars, panel.itertuples(index=False)):
        height = float(row.f_micro_w)
        x = bar.get_x() + bar.get_width() / 2
        if team_label_position == "above":
            ax.text(
                x,
                height + y_max * 0.08,
                row.prediction_label,
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=12,
                color="black",
                clip_on=False,
            )
            ax.text(
                x,
                height + y_max * 0.02,
                f"{height:.3f}",
                ha="center",
                va="bottom",
                fontsize=11,
                color="black",
            )
        else:
            ax.text(
                x,
                height + 0.015,
                f"{height:.3f}",
                ha="center",
                va="bottom",
                fontsize=11,
            )
            ax.text(
                x,
                max(height * 0.5, 0.03),
                row.prediction_label,
                ha="center",
                va="center",
                rotation=90,
                fontsize=12,
                color="black",
                clip_on=True,
            )


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot CAFA6 f_micro_w evaluator results by subset and GO aspect.")
    parser.add_argument(
        "--data-dir",
        "--DATA_DIR",
        dest="data_dir",
        required=True,
        help="Directory containing test_cafaeval_NK, test_cafaeval_LK, and test_cafaeval_PK folders.",
    )
    parser.add_argument("--output-plot", required=True, help="Output plot path, e.g. cafa6_f_micro_w_grid.png.")
    parser.add_argument("--output-table", default=None, help="Optional long-form TSV of plotted scores.")
    parser.add_argument(
        "--score-metadata",
        default=None,
        help="Optional CAFA score metadata CSV with SubmissionId and TeamName columns.",
    )
    parser.add_argument("--width", type=float, default=None, help="Optional figure width in inches.")
    parser.add_argument("--height", type=float, default=None, help="Optional figure height in inches.")
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    scores = collect_evaluator_scores(args.data_dir)
    scores = apply_score_metadata(scores, args.score_metadata)
    if args.output_table:
        write_score_table(scores, args.output_table)

    figsize = (args.width, args.height) if args.width and args.height else None
    plot_cafa6_results(scores, output_file=args.output_plot, figsize=figsize)


if __name__ == "__main__":
    main()
