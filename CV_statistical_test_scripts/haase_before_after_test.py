#!/usr/bin/env python3
"""Before/after tests for the PBS and HAase posterior-angle experiment.

Each embryo is imaged twice at 60 hpf. The file name encodes the pair:
01-spot is before treatment and 03-spot is 0.5 h later. Embryos that lack
one of those two images are left out. The two imaging dates are pooled,
which gives the figure's sample sizes: 10 PBS embryos and 9 HAase embryos.
One HAase embryo on 2026-03-05 (e02) has no after image.

Two questions are kept separate.

Angle change within a group
    Did the posterior angle itself move from before to after? The embryos
    are paired, and the samples are small, so the primary test is the exact
    Wilcoxon signed-rank test. A paired t-test is reported beside it.

CV change
    The figure's percent change is

        percent = (CV_after - CV_before) / CV_before

    CV is the sample standard deviation divided by the mean. Dividing the
    standard deviation by n instead of n - 1 rescales both CVs by the same
    factor, so the percent change is identical. The printed figure rounds
    the divisor-n CVs to three decimals and then computes the percent from
    those rounded values (PBS -3.39%, HAase -75.93%). This script tests the
    unrounded percent.

    Embryos are resampled with replacement, and each embryo keeps its before
    and after angle together. That produces a bootstrap distribution of the
    percent change. The 95% interval is the 2.5 and 97.5 percentiles. For a
    negative observed change, the two-sided p-value is twice the fraction of
    bootstrap percents that are greater than or equal to zero.

Difference between groups
    The comparison statistic is

        percent_HAase - percent_PBS

    Each group is resampled on its own. The same percentile interval and
    two-sided tail probability are used for the null hypothesis that the
    two percent changes are equal.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from project_paths import measurement_dir, results_dir

FILE_PATTERN = re.compile(
    r"(?P<date>\d{8})_mNG_bung_(?P<group>HAase|PBS)_60hpf_(?P<embryo>e\d+)-(?P<spot>01|03)-Spot"
)
GROUP_ORDER = ("PBS", "HAase")
DATA_FILES = (
    "posteriorAngles_20260305_HAase_60hpf_01.csv",
    "posteriorAngles_20260306_HAase_60hpf_01.csv",
)


def coefficient_of_variation(values: np.ndarray, ddof: int = 1) -> float:
    values = np.asarray(values, dtype=float)
    mean = float(values.mean())
    if mean == 0.0:
        raise ValueError("Coefficient of variation is undefined when the mean is 0.")
    return float(values.std(ddof=ddof) / mean)


def percent_cv_change(before: np.ndarray, after: np.ndarray) -> float:
    baseline = coefficient_of_variation(before)
    if baseline == 0.0:
        raise ValueError("Percent CV change is undefined when the before CV is 0.")
    return (coefficient_of_variation(after) - baseline) / baseline


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function."""
    max_iter = 200
    epsilon = 3e-14
    tiny = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < epsilon:
            break
    return h


def regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_beta = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log(1.0 - x)
    )
    beta = math.exp(log_beta)
    if x < (a + 1.0) / (a + b + 2.0):
        return beta * _betacf(a, b, x) / a
    return 1.0 - beta * _betacf(b, a, 1.0 - x) / b


def student_t_two_sided_p(t_statistic: float, degrees: int) -> float:
    x = degrees / (degrees + t_statistic**2)
    return regularized_incomplete_beta(degrees / 2.0, 0.5, x)


def paired_t(before: np.ndarray, after: np.ndarray) -> tuple[float, float, float, float, float]:
    """Mean after-minus-before, its t statistic, two-sided p, and 95% CI."""
    difference = np.asarray(after, dtype=float) - np.asarray(before, dtype=float)
    n = len(difference)
    mean = float(difference.mean())
    sd = float(difference.std(ddof=1))
    se = sd / math.sqrt(n)
    t_statistic = mean / se
    degrees = n - 1
    p_value = student_t_two_sided_p(t_statistic, degrees)
    # A short numerical inverse of the t cdf supplies the 95% critical value.
    target = 0.025
    lo, hi = 0.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        tail = 0.5 * student_t_two_sided_p(mid, degrees)
        if tail > target:
            lo = mid
        else:
            hi = mid
    critical = 0.5 * (lo + hi)
    return mean, t_statistic, p_value, mean - critical * se, mean + critical * se


def wilcoxon_signed_rank_p(before: np.ndarray, after: np.ndarray) -> tuple[float, float]:
    """Exact two-sided Wilcoxon signed-rank p-value.

    Returns the smaller signed-rank sum and the p-value. Differences of zero
    are dropped. Average ranks are used if absolute differences tie.
    """
    difference = np.asarray(after, dtype=float) - np.asarray(before, dtype=float)
    difference = difference[difference != 0.0]
    n = len(difference)
    if n == 0:
        raise ValueError("Every paired difference is zero.")
    if n > 20:
        raise ValueError("The exact Wilcoxon test in this script supports at most 20 pairs.")
    ranks = pd.Series(np.abs(difference)).rank(method="average").to_numpy()
    positive = float(np.sum(ranks[difference > 0.0]))
    total = float(ranks.sum())
    smaller = min(positive, total - positive)
    as_extreme = 0
    for mask in range(1 << n):
        signed = 0.0
        for bit in range(n):
            if mask & (1 << bit):
                signed += ranks[bit]
        if min(signed, total - signed) <= smaller + 1e-9:
            as_extreme += 1
    return smaller, as_extreme / (1 << n)


def load_pairs(paths: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = [pd.read_csv(path, usecols=["Angle", "File"]) for path in paths]
    measurements = pd.concat(frames, ignore_index=True)
    parsed = measurements["File"].str.extract(FILE_PATTERN)
    if parsed.isna().any().any():
        bad = measurements.loc[parsed.isna().any(axis=1), "File"].tolist()
        raise ValueError(f"Could not parse embryo, group, or spot from: {bad}")
    measurements = pd.concat([measurements, parsed], axis=1)
    wide = measurements.pivot_table(
        index=["date", "group", "embryo"],
        columns="spot",
        values="Angle",
        aggfunc="mean",
    )
    missing = wide[wide["01"].isna() | wide["03"].isna()].reset_index()
    paired = (
        wide.dropna()
        .rename(columns={"01": "before", "03": "after"})
        .reset_index()
        .sort_values(["group", "date", "embryo"])
        .reset_index(drop=True)
    )
    paired["difference"] = paired["after"] - paired["before"]
    return paired, missing


def _resample_percent(before: np.ndarray, after: np.ndarray, rng: np.random.Generator, n_boot: int) -> np.ndarray:
    n = len(before)
    draws = rng.integers(0, n, size=(n_boot, n))
    percents = np.empty(n_boot, dtype=float)
    for i, index in enumerate(draws):
        percents[i] = percent_cv_change(before[index], after[index])
    if not np.all(np.isfinite(percents)):
        raise RuntimeError("A bootstrap replicate had an undefined percent CV change.")
    return percents


def _tail_p(observed: float, replicates: np.ndarray) -> tuple[float, str, int]:
    """Two-sided percentile p-value for H0: statistic = 0."""
    n_boot = len(replicates)
    n_positive = int(np.sum(replicates >= 0.0))
    n_negative = int(np.sum(replicates <= 0.0))
    minority = n_positive if observed < 0.0 else n_negative
    if observed == 0.0:
        minority = min(n_positive, n_negative)
    if minority == 0:
        return 2.0 / n_boot, f"< {2.0 / n_boot:.1e}", minority
    p_value = min(1.0, 2.0 * minority / n_boot)
    return p_value, f"{p_value:.4g}", minority


def group_summary(paired: pd.DataFrame, group: str, rng: np.random.Generator, n_boot: int) -> dict[str, float | str]:
    subset = paired.loc[paired["group"] == group]
    before = subset["before"].to_numpy(dtype=float)
    after = subset["after"].to_numpy(dtype=float)
    mean_difference, t_statistic, t_p, ci_low, ci_high = paired_t(before, after)
    wilcoxon_statistic, wilcoxon_p = wilcoxon_signed_rank_p(before, after)
    cv_before = coefficient_of_variation(before)
    cv_after = coefficient_of_variation(after)
    cv_before_figure = coefficient_of_variation(before, ddof=0)
    cv_after_figure = coefficient_of_variation(after, ddof=0)
    percent = (cv_after - cv_before) / cv_before
    boot = _resample_percent(before, after, rng, n_boot)
    low, high = np.quantile(boot, [0.025, 0.975])
    p_value, p_display, minority = _tail_p(percent, boot)
    return {
        "group": group,
        "n_embryos": len(subset),
        "mean_before": float(before.mean()),
        "mean_after": float(after.mean()),
        "mean_difference_after_minus_before": mean_difference,
        "difference_t_ci95_low": ci_low,
        "difference_t_ci95_high": ci_high,
        "paired_t": t_statistic,
        "paired_t_p": t_p,
        "wilcoxon_signed_rank": wilcoxon_statistic,
        "wilcoxon_p": wilcoxon_p,
        "cv_before": cv_before,
        "cv_after": cv_after,
        "cv_before_divisor_n": cv_before_figure,
        "cv_after_divisor_n": cv_after_figure,
        "percent_cv_change": percent,
        "percent_cv_change_ci95_low": float(low),
        "percent_cv_change_ci95_high": float(high),
        "percent_cv_change_p": p_value,
        "percent_cv_change_p_display": p_display,
        "n_bootstrap_percent_ge_0": int(np.sum(boot >= 0.0)),
        "minority_side_count": minority,
        "n_bootstrap": n_boot,
    }


def contrast_summary(
    paired: pd.DataFrame,
    rng: np.random.Generator,
    n_boot: int,
) -> dict[str, float | str]:
    groups = {
        group: paired.loc[paired["group"] == group, ["before", "after"]].to_numpy(dtype=float)
        for group in GROUP_ORDER
    }
    observed = percent_cv_change(groups["HAase"][:, 0], groups["HAase"][:, 1]) - percent_cv_change(
        groups["PBS"][:, 0], groups["PBS"][:, 1]
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
        "contrast": "HAase percent CV change minus PBS percent CV change",
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
    paired: pd.DataFrame,
    missing: pd.DataFrame,
    summaries: list[dict[str, float | str]],
    contrast: dict[str, float | str],
    seed: int,
) -> str:
    lines = [
        "PBS and HAase posterior angle, before versus 0.5 h after",
        "01-spot is before. 03-spot is after. Dates 2026-03-05 and 2026-03-06 are pooled.",
        f"Bootstrap replicates: {contrast['n_bootstrap']}. Seed: {seed}.",
        "",
    ]
    if len(missing):
        lines.append("Embryos omitted because one time point is missing:")
        for row in missing.itertuples(index=False):
            lines.append(f"  {row.date} {row.group} {row.embryo}")
        lines.append("")
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
                    f"  divisor-n CV, as labeled on the figure: "
                    f"{summary['cv_before_divisor_n']:.3f} -> {summary['cv_after_divisor_n']:.3f}"
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
            "HAase minus PBS, percent CV change",
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
            "Percent CV change by imaging date, descriptive only",
        ]
    )
    for (date, group), subset in paired.groupby(["date", "group"], sort=True):
        percent = percent_cv_change(subset["before"].to_numpy(), subset["after"].to_numpy())
        lines.append(f"  {date} {group:5}  n={len(subset):2d}  percent change {100 * percent:.2f}%")
    lines.append("")
    return "\n".join(lines)


def _self_check() -> None:
    if abs(student_t_two_sided_p(1.775938714757421, 9) - 0.1094759) > 1e-6:
        raise AssertionError("Paired t p-value does not match the reference PBS result.")
    if abs(student_t_two_sided_p(3.2452572331642084, 8) - 0.01178447) > 1e-6:
        raise AssertionError("Paired t p-value does not match the reference HAase result.")
    _, p_value = wilcoxon_signed_rank_p(np.array([0.0, 0.0, 0.0]), np.array([1.0, 2.0, 3.0]))
    if abs(p_value - 0.25) > 1e-12:
        raise AssertionError(f"Wilcoxon exact check failed: p={p_value}")
    before = np.array([1.0, 2.0, 3.0, 4.0])
    after = before * 2.0
    if abs(percent_cv_change(before, after)) > 1e-12:
        raise AssertionError("Scaling every angle by the same factor changed the CV.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Paired before/after tests for PBS and HAase angles.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory containing the two posteriorAngles HAase CSVs.",
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
    paths = [data_dir / name for name in DATA_FILES]
    missing_files = [path for path in paths if not path.is_file()]
    if missing_files:
        raise SystemExit(f"Missing data files: {missing_files}")
    paired, omitted = load_pairs(paths)
    rng = np.random.default_rng(args.seed)
    summaries = [group_summary(paired, group, rng, args.n_boot) for group in GROUP_ORDER]
    contrast = contrast_summary(paired, rng, args.n_boot)
    for summary in summaries:
        summary["seed"] = args.seed
    contrast["seed"] = args.seed
    paired_path = results_dir() / "haase_paired_angles.csv"
    result_path = results_dir() / "haase_cv_change_results.csv"
    paired.to_csv(paired_path, index=False)
    pd.DataFrame(summaries).to_csv(result_path, index=False)
    contrast_path = results_dir() / "haase_cv_change_contrast.csv"
    pd.DataFrame([contrast]).to_csv(contrast_path, index=False)
    print(render(paired, omitted, summaries, contrast, args.seed), end="")
    print(f"Wrote {paired_path}")
    print(f"Wrote {result_path}")
    print(f"Wrote {contrast_path}")


if __name__ == "__main__":
    main()
