#!/usr/bin/env python3
"""Before/after tests for the DMSO and forskolin posterior-angle experiment.

PosteriorAngle_1hr-forskolin.xlsx is a wide table. Each row is one treatment.
Columns 1-10 are the angles before treatment and columns 11-20 are the angles
1 h later, in the same embryo order. DMSO and forskolin each have 10 pairs.

The questions match the PBS/HAase analysis.

Angle change within a group
    Did the posterior angle move from before to 1 h later? The primary test
    is the exact Wilcoxon signed-rank test. A paired t-test is reported
    beside it.

Angle change between groups
    The number printed above the figure is the equal-variance two-sample t
    test of the paired changes: forskolin's mean change minus DMSO's mean
    change. An exact Mann-Whitney test of those same changes is reported
    beside it.

CV change
    The percent change is (CV_after - CV_before) / CV_before, using the
    sample standard deviation. The figure rounds those CVs to three decimals
    and then computes the percent from the rounded labels (DMSO -33.33%,
    forskolin -43.48%). This script tests the unrounded percent. Embryos are
    resampled with each embryo's before and after angle kept together.

    The comparison between treatments is forskolin's percent change minus
    DMSO's percent change. Each group is resampled on its own.
"""

from __future__ import annotations

import argparse
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from project_paths import measurement_dir, results_dir
from haase_before_after_test import (
    _resample_percent,
    _tail_p,
    coefficient_of_variation,
    group_summary,
    paired_t,
    percent_cv_change,
    student_t_two_sided_p,
    wilcoxon_signed_rank_p,
)

GROUP_ORDER = ("DMSO", "Forskolin")
DATA_FILE = "PosteriorAngle_1hr-forskolin.xlsx"


def load_pairs(path: Path) -> pd.DataFrame:
    sheet = pd.read_excel(path, header=None)
    rows = []
    for row_index in range(1, len(sheet)):
        group = str(sheet.iloc[row_index, 0]).strip()
        before = sheet.iloc[row_index, 1:11].to_numpy(dtype=float)
        after = sheet.iloc[row_index, 11:21].to_numpy(dtype=float)
        if len(before) != 10 or len(after) != 10:
            raise ValueError(f"{group} does not have 10 before and 10 after angles.")
        if not np.all(np.isfinite(before)) or not np.all(np.isfinite(after)):
            raise ValueError(f"{group} has a missing angle.")
        for embryo, (angle_before, angle_after) in enumerate(zip(before, after), start=1):
            rows.append(
                {
                    "group": group,
                    "embryo": embryo,
                    "before": float(angle_before),
                    "after": float(angle_after),
                    "difference": float(angle_after - angle_before),
                }
            )
    paired = pd.DataFrame(rows)
    unknown = sorted(set(paired["group"]) - set(GROUP_ORDER))
    if unknown or set(paired["group"]) != set(GROUP_ORDER):
        raise ValueError(f"Expected rows {GROUP_ORDER}, found {sorted(paired['group'].unique())}.")
    return paired.sort_values(["group", "embryo"]).reset_index(drop=True)


def t_critical(degrees: int, tail: float = 0.025) -> float:
    lo, hi = 0.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if 0.5 * student_t_two_sided_p(mid, degrees) > tail:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def angle_contrast(paired: pd.DataFrame) -> dict[str, float]:
    """Equal-variance t test and exact Mann-Whitney test of the paired changes."""
    changes = {
        group: paired.loc[paired["group"] == group, "difference"].to_numpy(dtype=float)
        for group in GROUP_ORDER
    }
    dmso = changes["DMSO"]
    forskolin = changes["Forskolin"]
    n_dmso = len(dmso)
    n_forskolin = len(forskolin)
    mean_dmso = float(dmso.mean())
    mean_forskolin = float(forskolin.mean())
    estimate = mean_forskolin - mean_dmso
    degrees = n_dmso + n_forskolin - 2
    pooled_variance = (
        (n_dmso - 1) * float(dmso.var(ddof=1)) + (n_forskolin - 1) * float(forskolin.var(ddof=1))
    ) / degrees
    standard_error = math.sqrt(pooled_variance * (1.0 / n_dmso + 1.0 / n_forskolin))
    t_statistic = estimate / standard_error
    p_value = student_t_two_sided_p(t_statistic, degrees)
    critical = t_critical(degrees)
    rank_sum, mann_whitney_p = mann_whitney_exact(forskolin, dmso)
    return {
        "contrast": "Forskolin mean angle change minus DMSO mean angle change",
        "n_dmso": n_dmso,
        "n_forskolin": n_forskolin,
        "mean_change_dmso": mean_dmso,
        "mean_change_forskolin": mean_forskolin,
        "estimate_degrees": estimate,
        "equal_variance_t": t_statistic,
        "degrees_of_freedom": degrees,
        "t_p": p_value,
        "t_ci95_low": estimate - critical * standard_error,
        "t_ci95_high": estimate + critical * standard_error,
        "forskolin_rank_sum": rank_sum,
        "mann_whitney_p": mann_whitney_p,
    }


def mann_whitney_exact(first: np.ndarray, second: np.ndarray) -> tuple[float, float]:
    """Two-sided exact Mann-Whitney p-value for the rank sum of `first`."""
    combined = np.concatenate([np.asarray(first, dtype=float), np.asarray(second, dtype=float)])
    ranks = pd.Series(combined).rank(method="average").to_numpy()
    n_first = len(first)
    observed = float(ranks[:n_first].sum())
    center = float(ranks.sum()) * n_first / len(ranks)
    deviation = abs(observed - center)
    as_extreme = 0
    n_assignments = 0
    for chosen in combinations(range(len(ranks)), n_first):
        n_assignments += 1
        if abs(float(ranks[list(chosen)].sum()) - center) >= deviation - 1e-9:
            as_extreme += 1
    return observed, as_extreme / n_assignments


def cv_contrast(paired: pd.DataFrame, rng: np.random.Generator, n_boot: int) -> dict[str, float | str]:
    groups = {
        group: paired.loc[paired["group"] == group, ["before", "after"]].to_numpy(dtype=float)
        for group in GROUP_ORDER
    }
    observed = percent_cv_change(groups["Forskolin"][:, 0], groups["Forskolin"][:, 1]) - percent_cv_change(
        groups["DMSO"][:, 0], groups["DMSO"][:, 1]
    )
    boot = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        percents = []
        for group in GROUP_ORDER:
            values = groups[group]
            index = rng.integers(0, len(values), size=len(values))
            percents.append(percent_cv_change(values[index, 0], values[index, 1]))
        boot[i] = percents[1] - percents[0]
    low, high = np.quantile(boot, [0.025, 0.975])
    p_value, p_display, minority = _tail_p(observed, boot)
    return {
        "contrast": "Forskolin percent CV change minus DMSO percent CV change",
        "estimate": observed,
        "ci95_low": float(low),
        "ci95_high": float(high),
        "bootstrap_p": p_value,
        "bootstrap_p_display": p_display,
        "n_bootstrap_contrast_ge_0": int(np.sum(boot >= 0.0)),
        "minority_side_count": minority,
        "n_bootstrap": n_boot,
    }


def render(
    summaries: list[dict[str, float | str]],
    angles: dict[str, float],
    contrast: dict[str, float | str],
    seed: int,
) -> str:
    lines = [
        "DMSO and forskolin posterior angle, before versus 1 h after",
        "Columns are paired in embryo order. Each group has 10 embryos.",
        f"Bootstrap replicates: {contrast['n_bootstrap']}. Seed: {seed}.",
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
                    f"  95% t interval for the mean change "
                    f"[{summary['difference_t_ci95_low']:.2f}, {summary['difference_t_ci95_high']:.2f}]"
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
                    f"  sample CV rounded to 3 decimals, as labeled on the figure: "
                    f"{summary['cv_before']:.3f} -> {summary['cv_after']:.3f}"
                ),
                (
                    f"  percent CV change 95% bootstrap interval "
                    f"[{100 * summary['percent_cv_change_ci95_low']:.2f}%, "
                    f"{100 * summary['percent_cv_change_ci95_high']:.2f}%]"
                ),
                (
                    f"  bootstrap replicates with percent change >= 0: "
                    f"{summary['n_bootstrap_percent_ge_0']} of {summary['n_bootstrap']}"
                ),
                f"  two-sided p for no CV change = {summary['percent_cv_change_p_display']}",
                "",
            ]
        )
    lines.extend(
        [
            "Forskolin minus DMSO, mean angle change",
            (
                f"  estimate {angles['estimate_degrees']:.2f} degrees   "
                f"95% t interval [{angles['t_ci95_low']:.2f}, {angles['t_ci95_high']:.2f}]"
            ),
            (
                f"  equal-variance t = {angles['equal_variance_t']:.3f}, "
                f"df = {angles['degrees_of_freedom']}, two-sided p = {angles['t_p']:.4g}"
            ),
            f"  exact Mann-Whitney p = {angles['mann_whitney_p']:.4g}",
            "",
            "Forskolin minus DMSO, percent CV change",
            (
                f"  estimate {100 * contrast['estimate']:.2f} percentage points   "
                f"95% interval [{100 * contrast['ci95_low']:.2f}, {100 * contrast['ci95_high']:.2f}]"
            ),
            (
                f"  bootstrap replicates >= 0: {contrast['n_bootstrap_contrast_ge_0']} "
                f"of {contrast['n_bootstrap']}"
            ),
            f"  two-sided p = {contrast['bootstrap_p_display']}",
            "",
        ]
    )
    return "\n".join(lines)


def _self_check() -> None:
    before = np.array([1.0, 2.0, 3.0, 4.0])
    if abs(percent_cv_change(before, before * 2.0)) > 1e-12:
        raise AssertionError("Scaling every angle by the same factor changed the CV.")
    _, p_value = wilcoxon_signed_rank_p(np.zeros(3), np.array([1.0, 2.0, 3.0]))
    if abs(p_value - 0.25) > 1e-12:
        raise AssertionError(f"Wilcoxon exact check failed: p={p_value}")
    _, mann_p = mann_whitney_exact(np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0]))
    if abs(mann_p - 0.1) > 1e-12:
        raise AssertionError(f"Mann-Whitney exact check failed: p={mann_p}")
    reference = pd.read_excel(measurement_dir() / DATA_FILE, header=None)
    dmso_before = reference.iloc[1, 1:11].to_numpy(dtype=float)
    dmso_after = reference.iloc[1, 11:21].to_numpy(dtype=float)
    forskolin_before = reference.iloc[2, 1:11].to_numpy(dtype=float)
    forskolin_after = reference.iloc[2, 11:21].to_numpy(dtype=float)
    if abs(coefficient_of_variation(dmso_before) - 0.075) > 0.0005:
        raise AssertionError("DMSO before CV does not match the figure label 0.075.")
    if abs(coefficient_of_variation(forskolin_before) - 0.115) > 0.0005:
        raise AssertionError("Forskolin before CV does not match the figure label 0.115.")
    dmso_change = dmso_after - dmso_before
    forskolin_change = forskolin_after - forskolin_before
    estimate = float(forskolin_change.mean() - dmso_change.mean())
    degrees = 18
    pooled = (
        9 * float(dmso_change.var(ddof=1)) + 9 * float(forskolin_change.var(ddof=1))
    ) / degrees
    t_statistic = estimate / math.sqrt(pooled * (0.1 + 0.1))
    if abs(student_t_two_sided_p(t_statistic, degrees) - 0.01416) > 5e-5:
        raise AssertionError("Between-group angle t-test does not match the figure p-value 0.0142.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Paired before/after tests for DMSO and forskolin angles.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory containing PosteriorAngle_1hr-forskolin.xlsx.",
    )
    parser.add_argument("--n-boot", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_boot < 1000:
        raise SystemExit("Use at least 1000 bootstrap replicates.")
    _self_check()
    data_dir = args.data_dir or measurement_dir()
    path = data_dir / DATA_FILE
    if not path.is_file():
        raise SystemExit(f"Missing data file: {path}")
    paired = load_pairs(path)
    rng = np.random.default_rng(args.seed)
    summaries = [group_summary(paired, group, rng, args.n_boot) for group in GROUP_ORDER]
    angles = angle_contrast(paired)
    contrast = cv_contrast(paired, rng, args.n_boot)
    for summary in summaries:
        summary["seed"] = args.seed
    angles["seed"] = args.seed
    contrast["seed"] = args.seed
    paired_path = results_dir() / "forskolin_paired_angles.csv"
    result_path = results_dir() / "forskolin_cv_change_results.csv"
    angle_path = results_dir() / "forskolin_angle_contrast.csv"
    contrast_path = results_dir() / "forskolin_cv_change_contrast.csv"
    paired.to_csv(paired_path, index=False)
    pd.DataFrame(summaries).to_csv(result_path, index=False)
    pd.DataFrame([angles]).to_csv(angle_path, index=False)
    pd.DataFrame([contrast]).to_csv(contrast_path, index=False)
    print(render(summaries, angles, contrast, args.seed), end="")
    print(f"Wrote {paired_path}")
    print(f"Wrote {result_path}")
    print(f"Wrote {angle_path}")
    print(f"Wrote {contrast_path}")


if __name__ == "__main__":
    main()
