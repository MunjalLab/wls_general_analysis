#!/usr/bin/env python3
"""Before/after CV tests for the anterior and posterior ablation figures.

Anterior
    The plotted points are the 10 embryos on plotted_pairs in
    anterior_ablation_before_after.xlsx. Before is the 01 spot and after is
    the 05 spot. The figure's printed coefficients of variation, 0.070 and
    0.092, were calculated on those 10 plus the 3 embryos on
    cv_annotation_only. Both sets are tested.

Posterior
    posterior_ablation_before_after.xlsx has the 11 Prism pairs. The figure
    says n = 10 and its printed coefficients of variation, 0.092 and 0.112,
    match the rows marked in_figure_n10 = yes. That set is the primary test.
    All 11 pairs are tested as well, because the Prism paired t-test on the
    figure used 11.

The percent change is (CV_after - CV_before) / CV_before, with the sample
standard deviation. Embryos are resampled with each embryo's before and after
angle kept together. The angle change is tested with an exact Wilcoxon
signed-rank test and a paired t-test.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from haase_before_after_test import group_summary, percent_cv_change
from project_paths import measurement_dir, results_dir

ANTERIOR_FILE = "anterior_ablation_before_after.xlsx"
POSTERIOR_FILE = "posterior_ablation_before_after.xlsx"


def _as_pairs(frame: pd.DataFrame, group: str) -> pd.DataFrame:
    paired = pd.DataFrame(
        {
            "group": group,
            "embryo": frame["embryo"].astype(str) if "embryo" in frame.columns else frame["prism_row"].astype(str),
            "before": frame["before_angle_degrees"].to_numpy(dtype=float),
            "after": frame["after_angle_degrees"].to_numpy(dtype=float),
        }
    )
    if not np.all(np.isfinite(paired[["before", "after"]].to_numpy())):
        raise ValueError(f"{group} has a missing angle.")
    return paired


def load_sets(data_dir: Path) -> list[tuple[str, pd.DataFrame]]:
    anterior = pd.read_excel(data_dir / ANTERIOR_FILE, sheet_name="plotted_pairs")
    anterior_extra = pd.read_excel(data_dir / ANTERIOR_FILE, sheet_name="cv_annotation_only")
    posterior = pd.read_excel(data_dir / POSTERIOR_FILE, sheet_name="pairs")
    posterior_figure = posterior.loc[posterior["in_figure_n10"].astype(str).str.lower() == "yes"].copy()
    return [
        ("Anterior plotted n=10", _as_pairs(anterior, "Anterior plotted n=10")),
        (
            "Anterior CV labels n=13",
            _as_pairs(pd.concat([anterior, anterior_extra], ignore_index=True), "Anterior CV labels n=13"),
        ),
        ("Posterior figure n=10", _as_pairs(posterior_figure, "Posterior figure n=10")),
        ("Posterior Prism n=11", _as_pairs(posterior, "Posterior Prism n=11")),
    ]


def render(summaries: list[dict[str, float | str]], seed: int) -> str:
    lines = [
        "Ablation posterior and anterior angle, before versus after",
        "Percent CV change resamples embryos and keeps each before/after pair together.",
        f"Bootstrap replicates: {summaries[0]['n_bootstrap']}. Seed: {seed}.",
        "",
    ]
    for summary in summaries:
        lines.extend(
            [
                str(summary["group"]),
                (
                    f"  n = {summary['n_embryos']}   "
                    f"mean angle {summary['mean_before']:.1f} -> {summary['mean_after']:.1f}   "
                    f"mean change {summary['mean_difference_after_minus_before']:.2f} degrees"
                ),
                (
                    f"  Wilcoxon signed-rank p = {summary['wilcoxon_p']:.4g}   "
                    f"paired t p = {summary['paired_t_p']:.4g}"
                ),
                (
                    f"  CV {summary['cv_before']:.4f} -> {summary['cv_after']:.4f}   "
                    f"percent change {100 * summary['percent_cv_change']:.2f}%"
                ),
                (
                    f"  percent CV change 95% bootstrap interval "
                    f"[{100 * summary['percent_cv_change_ci95_low']:.2f}%, "
                    f"{100 * summary['percent_cv_change_ci95_high']:.2f}%]"
                ),
                f"  two-sided p for no CV change = {summary['percent_cv_change_p_display']}",
                "",
            ]
        )
    lines.extend(
        [
            "Anterior plotted is the 10 connected points. Anterior CV labels adds the 3 embryos",
            "used for the printed 0.070 and 0.092. Posterior figure drops Prism row 10 so that",
            "n = 10 and the sample CVs round to the printed 0.092 and 0.112.",
            "",
        ]
    )
    return "\n".join(lines)


def _self_check(data_dir: Path) -> None:
    sets = dict(load_sets(data_dir))
    anterior = sets["Anterior plotted n=10"]
    if len(anterior) != 10:
        raise AssertionError("The anterior plotted sheet does not have 10 embryos.")
    posterior = sets["Posterior figure n=10"]
    if len(posterior) != 10:
        raise AssertionError("The posterior figure set does not have 10 embryos.")
    before = anterior["before"].to_numpy()
    after = anterior["after"].to_numpy()
    if abs(percent_cv_change(before, after) - 0.22487) > 0.001:
        raise AssertionError("Anterior plotted percent CV change is not the expected 22.5%.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Before/after CV tests for the ablation figures.")
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--n-boot", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_boot < 1000:
        raise SystemExit("Use at least 1000 bootstrap replicates.")
    data_dir = args.data_dir or measurement_dir()
    _self_check(data_dir)
    rng = np.random.default_rng(args.seed)
    summaries = [group_summary(paired, group, rng, args.n_boot) for group, paired in load_sets(data_dir)]
    for summary in summaries:
        summary["seed"] = args.seed
    result_path = results_dir() / "ablation_cv_change_results.csv"
    pd.DataFrame(summaries).to_csv(result_path, index=False)
    print(render(summaries, args.seed), end="")
    print(f"Wrote {result_path}")


if __name__ == "__main__":
    main()
