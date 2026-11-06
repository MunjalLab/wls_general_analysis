#!/usr/bin/env python3
"""Bootstrap test of the pillar-angle CV against lumen volume.

The volume bins are the ones in the Fig. 3C' plotting script
2026-06-09_plot_allpts_lumenVol_pillarAngle_alltime.R:

    binsize = 90
    breaks = seq(from = 315, to = 1350, by = binsize)

R's cut() keeps volumes from 315 through 1305 pL. The first interval is
closed on the left and the rest are right-closed, so the bins are

    [315, 405], (405, 495], ..., (1215, 1305]

and the x-axis labels are the bin centers 360, 450, ..., 1260. Embryos
outside that range are omitted, as in the R script. The same breaks are
used for anterior, posterior, and ventral angles. Posterior angles are the
non-missing P-Angle rows in WT_AllAspectRatiosandAngles.csv, which match
Fig. 3C'. Anterior and ventral angles are the WT_Anterior and WT_Ventral
tables in this folder; those tables match the A-Angle and V-Angle columns
of the same source file.

The test then follows the binned-CV regression:

1. In each occupied bin, CV = sample standard deviation (divisor n - 1) / mean.
   Bins with fewer than two embryos, or with zero standard deviation, are
   dropped, matching the R summary filter.
2. Fit the weighted line CV_i = beta0 + beta1 * V_i. V_i is the bin center,
   the value printed on the figure axis. Each bin is weighted by its number
   of embryos.
3. Resample embryos with replacement, recompute the bin CVs, and refit.
   Repeat 10,000 times. The replicates are the bootstrap distribution of beta1.

The hypothesis is that the slope is zero: angle CV does not change with
lumen volume. The 95% interval is the 2.5 and 97.5 percentiles of that
distribution. The two-sided p-value is twice the fraction of bootstrap
slopes that fall on the other side of zero from the observed slope. With
B replicates and K of them across zero,

    p = min(1, 2 * K / B).

If K is 0, the p-value is reported as less than 2/B.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from project_paths import measurement_dir, results_dir

BIN_WIDTH_PL = 90.0
BIN_BREAKS_PL = np.arange(315.0, 1350.0, BIN_WIDTH_PL)
MIN_BIN_COUNT = 2
PAPER_TIMEPOINTS = ("54hpf", "60hpf", "66hpf", "72hpf")
PAPER_DATA = Path(
    "/Users/jiacheng/Library/CloudStorage/OneDrive-DukeUniversity/"
    "Wls paper/Fig3/WT_AllAspectRatiosandAngles.csv"
)
ASPECTS = (
    ("Anterior", "WT_AnteriorAspectRatiosandAngles.csv", "Angle"),
    ("Posterior", None, "P-Angle"),
    ("Ventral", "WT_VentralAspectRatiosandAngles.csv", "Angle"),
)


def coefficient_of_variation(values: np.ndarray) -> float:
    mean = float(values.mean())
    if mean == 0.0:
        raise ValueError("Coefficient of variation is undefined when the mean is 0.")
    return float(values.std(ddof=1) / mean)


def weighted_line(volume: np.ndarray, response: np.ndarray, weight: np.ndarray) -> tuple[float, float, float]:
    """Weighted least squares for response = beta0 + beta1 * volume.

    Returns the intercept, the slope, and the weighted R-squared.
    """
    volume = np.asarray(volume, dtype=float)
    response = np.asarray(response, dtype=float)
    weight = np.asarray(weight, dtype=float)
    if len(volume) < 2:
        raise ValueError("The regression needs at least two bins.")
    total = float(weight.sum())
    if total <= 0.0:
        raise ValueError("Weights must sum to a positive number.")
    volume_mean = float(np.sum(weight * volume) / total)
    response_mean = float(np.sum(weight * response) / total)
    volume_centered = volume - volume_mean
    response_centered = response - response_mean
    denominator = float(np.sum(weight * volume_centered**2))
    if denominator == 0.0:
        raise ValueError("Volume does not vary across bins, so the slope is undefined.")
    slope = float(np.sum(weight * volume_centered * response_centered) / denominator)
    intercept = response_mean - slope * volume_mean
    fitted = intercept + slope * volume
    residual = float(np.sum(weight * (response - fitted) ** 2))
    total_ss = float(np.sum(weight * response_centered**2))
    r_squared = 1.0 - residual / total_ss if total_ss > 0.0 else np.nan
    return intercept, slope, r_squared


def bin_index(volume: np.ndarray) -> np.ndarray:
    """Bin index from R's cut(breaks, right = TRUE, include.lowest = TRUE)."""
    volume = np.asarray(volume, dtype=float)
    breaks = BIN_BREAKS_PL
    index = np.searchsorted(breaks, volume, side="left") - 1
    index[volume == breaks[0]] = 0
    outside = (volume < breaks[0]) | (volume > breaks[-1]) | ~np.isfinite(volume)
    index[outside] = -1
    return index.astype(int)


def bin_summary(angle: np.ndarray, volume: np.ndarray) -> pd.DataFrame:
    """One row per paper volume bin with at least two embryos and sd > 0."""
    assigned = bin_index(volume)
    rows = []
    for bin_id in range(len(BIN_BREAKS_PL) - 1):
        selected = assigned == bin_id
        count = int(selected.sum())
        if count < MIN_BIN_COUNT:
            continue
        angles = angle[selected]
        sd = float(angles.std(ddof=1))
        if sd == 0.0:
            continue
        low = float(BIN_BREAKS_PL[bin_id])
        high = float(BIN_BREAKS_PL[bin_id + 1])
        rows.append(
            {
                "bin_low_pl": low,
                "bin_high_pl": high,
                "bin_center_pl": (low + high) / 2.0,
                "interval": f"[{low:.0f}, {high:.0f}]" if bin_id == 0 else f"({low:.0f}, {high:.0f}]",
                "mean_volume_pl": float(volume[selected].mean()),
                "n_embryos": count,
                "mean_angle": float(angles.mean()),
                "sd_angle": sd,
                "cv": coefficient_of_variation(angles),
            }
        )
    if len(rows) < 2:
        raise ValueError("Fewer than two bins have enough embryos to fit a line.")
    return pd.DataFrame(rows)


def slope_from_bins(bins: pd.DataFrame) -> tuple[float, float, float]:
    return weighted_line(
        bins["bin_center_pl"].to_numpy(),
        bins["cv"].to_numpy(),
        bins["n_embryos"].to_numpy(dtype=float),
    )


def bootstrap_slopes(
    angle: np.ndarray,
    volume: np.ndarray,
    rng: np.random.Generator,
    n_boot: int,
) -> np.ndarray:
    """Refit beta1 after resampling embryos with replacement."""
    count = len(angle)
    slopes = np.empty(n_boot, dtype=float)
    for replicate in range(n_boot):
        drawn = rng.integers(0, count, size=count)
        slopes[replicate] = slope_from_bins(bin_summary(angle[drawn], volume[drawn]))[1]
    if not np.all(np.isfinite(slopes)):
        raise RuntimeError("A bootstrap replicate produced an undefined slope.")
    return slopes


def two_sided_p(n_positive: int, n_negative: int, n_bootstrap: int) -> tuple[float, str]:
    """Percentile p-value for H0: slope = 0.

    ``n_positive`` counts bootstrap slopes >= 0 and ``n_negative`` counts
    slopes <= 0. The p-value is twice the smaller of those fractions.
    When every replicate has the same sign, the displayed value is < 2/B.
    """
    minority = min(n_positive, n_negative)
    if minority == 0:
        return 2.0 / n_bootstrap, f"< {2.0 / n_bootstrap:.1e}"
    p_value = min(1.0, 2.0 * minority / n_bootstrap)
    return p_value, f"{p_value:.4g}"


def load_angles(path: Path, angle_column: str) -> pd.DataFrame:
    """Load one aspect's angles and the lumen volume of each embryo.

    The embryo id may be stored as Embryo or File. A wide table is limited
    to 54, 60, 66, and 72 hpf, matching the plotting script. Rows with a
    missing angle or volume are dropped.
    """
    table = pd.read_csv(path)
    if "Volume (pl)" not in table.columns:
        raise ValueError(f"{path.name} needs a Volume (pl) column.")
    if angle_column not in table.columns:
        raise ValueError(f"{path.name} needs a {angle_column} column.")
    id_column = "Embryo" if "Embryo" in table.columns else "File" if "File" in table.columns else None
    if id_column is None:
        raise ValueError(f"{path.name} needs an Embryo or File column.")
    selected = table.loc[:, [id_column, "Volume (pl)", angle_column]].rename(
        columns={id_column: "Embryo", angle_column: "Angle"}
    )
    selected = selected.loc[selected["Angle"].notna() & selected["Volume (pl)"].notna()].copy()
    if angle_column != "Angle":
        timepoint = selected["Embryo"].astype(str).str.extract(r"(\d+hpf)", expand=False)
        selected = selected.loc[timepoint.isin(PAPER_TIMEPOINTS)]
    selected = selected.reset_index(drop=True)
    if selected["Embryo"].duplicated().any():
        raise ValueError(f"{path.name} has duplicate embryo names, so the bootstrap unit is ambiguous.")
    if len(selected) < 2:
        raise ValueError(f"{path.name} has fewer than two angles in {angle_column}.")
    return selected


def angles_inside_bins(table: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Keep embryos inside the paper volume range. Return the count omitted."""
    inside = bin_index(table["Volume (pl)"].to_numpy(dtype=float)) >= 0
    return table.loc[inside].reset_index(drop=True), int((~inside).sum())


def run(
    aspect: str,
    path: Path,
    angle_column: str,
    n_boot: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, float | str], np.ndarray]:
    loaded = load_angles(path, angle_column)
    table, n_outside = angles_inside_bins(loaded)
    angle = table["Angle"].to_numpy(dtype=float)
    volume = table["Volume (pl)"].to_numpy(dtype=float)
    bins = bin_summary(angle, volume)
    bins.insert(0, "aspect", aspect)
    intercept, slope, r_squared = slope_from_bins(bins)
    slopes = bootstrap_slopes(angle, volume, np.random.default_rng(seed), n_boot)
    low, high = np.quantile(slopes, [0.025, 0.975])
    n_positive = int(np.sum(slopes >= 0.0))
    n_negative = int(np.sum(slopes <= 0.0))
    p_value, p_display = two_sided_p(n_positive, n_negative, n_boot)
    summary = {
        "aspect": aspect,
        "source": path.name,
        "n_embryos": len(table),
        "n_outside_bins": n_outside,
        "n_bins": len(bins),
        "bin_width_pl": BIN_WIDTH_PL,
        "intercept": intercept,
        "slope_per_pl": slope,
        "slope_per_90pl": slope * BIN_WIDTH_PL,
        "slope_per_100pl": slope * 100.0,
        "weighted_r_squared": r_squared,
        "ci95_low_per_pl": float(low),
        "ci95_high_per_pl": float(high),
        "ci95_low_per_90pl": float(low * BIN_WIDTH_PL),
        "ci95_high_per_90pl": float(high * BIN_WIDTH_PL),
        "bootstrap_p": p_value,
        "bootstrap_p_display": p_display,
        "n_bootstrap_slopes_ge_0": n_positive,
        "n_bootstrap_slopes_le_0": n_negative,
        "n_bootstrap": n_boot,
        "seed": seed,
    }
    return bins, summary, slopes


def render(bins: pd.DataFrame, summary: dict[str, float | str]) -> str:
    lines = [
        f"{summary['aspect']} pillar-angle CV versus lumen volume",
        "Bins match the paper: 90 pL, breaks at 315, 405, ..., 1305 pL.",
        "CV = sample standard deviation (divisor n - 1) / mean.",
        "Weighted line: CV = intercept + slope * bin center.",
        "Weight of each bin = number of embryos in that bin.",
        (
            f"Source: {summary['source']}. Embryos inside the bins: {summary['n_embryos']}. "
            f"Omitted outside 315-1305 pL: {summary['n_outside_bins']}."
        ),
        (
            f"Bins used: {summary['n_bins']}. "
            f"Bootstrap replicates: {summary['n_bootstrap']}. Seed: {summary['seed']}."
        ),
        "",
        "Bins",
    ]
    for row in bins.itertuples(index=False):
        lines.append(
            f"  {row.interval:>16}  center={row.bin_center_pl:6.0f}  "
            f"n={row.n_embryos:3d}  mean angle={row.mean_angle:6.1f}  CV={row.cv:.4f}"
        )
    lines.extend(
        [
            "",
            f"Intercept = {summary['intercept']:.4f}",
            (
                f"Slope = {summary['slope_per_pl']:.6g} per pL "
                f"({summary['slope_per_90pl']:.4f} per 90 pL)"
            ),
            (
                f"95% bootstrap interval = [{summary['ci95_low_per_pl']:.6g}, "
                f"{summary['ci95_high_per_pl']:.6g}] per pL"
            ),
            (
                f"                       = [{summary['ci95_low_per_90pl']:.4f}, "
                f"{summary['ci95_high_per_90pl']:.4f}] per 90 pL"
            ),
            f"Weighted R-squared = {summary['weighted_r_squared']:.3f}",
            (
                f"Bootstrap slopes >= 0: {summary['n_bootstrap_slopes_ge_0']} "
                f"of {summary['n_bootstrap']}"
            ),
            f"Two-sided p = {summary['bootstrap_p_display']}",
            "",
        ]
    )
    return "\n".join(lines)


def _self_check() -> None:
    intercept, slope, r_squared = weighted_line(
        np.array([0.0, 1.0, 2.0]),
        np.array([1.0, 3.0, 5.0]),
        np.array([1.0, 1.0, 1.0]),
    )
    if abs(intercept - 1.0) > 1e-12 or abs(slope - 2.0) > 1e-12 or abs(r_squared - 1.0) > 1e-12:
        raise AssertionError("Weighted line failed on a perfect line.")
    # Weighting the high endpoint more heavily must steepen this slope.
    _, equal_weight_slope, _ = weighted_line(
        np.array([0.0, 1.0, 2.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([1.0, 1.0, 1.0]),
    )
    _, heavy_slope, _ = weighted_line(
        np.array([0.0, 1.0, 2.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([1.0, 1.0, 100.0]),
    )
    if abs(equal_weight_slope - 0.5) > 1e-12:
        raise AssertionError("Unweighted slope check failed.")
    if not heavy_slope > equal_weight_slope:
        raise AssertionError("Weights did not pull the slope.")
    angles = np.array([1.0, 2.0, 3.0, 4.0])
    if abs(coefficient_of_variation(angles) - angles.std(ddof=1) / angles.mean()) > 1e-12:
        raise AssertionError("CV definition check failed.")
    probe = np.array([315.0, 405.0, 405.1, 1305.0, 1305.1, 300.0])
    found = bin_index(probe)
    if found.tolist() != [0, 0, 1, len(BIN_BREAKS_PL) - 2, -1, -1]:
        raise AssertionError(f"Bin edges do not match R cut(): {found.tolist()}")
    centers = (BIN_BREAKS_PL[:-1] + BIN_BREAKS_PL[1:]) / 2.0
    if not np.array_equal(centers, np.arange(360.0, 1260.0 + 1.0, 90.0)):
        raise AssertionError(f"Bin centers are not the figure labels: {centers.tolist()}")


def aspect_sources(directory: Path) -> list[tuple[str, Path, str]]:
    sources = []
    for aspect, filename, angle_column in ASPECTS:
        path = PAPER_DATA if filename is None else directory / filename
        sources.append((aspect, path, angle_column))
    return sources


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Bootstrap the slope of pillar-angle CV against lumen volume."
    )
    parser.add_argument(
        "--aspect",
        choices=[aspect for aspect, _, _ in ASPECTS],
        default=None,
        help="Run one aspect. The default runs anterior, posterior, and ventral.",
    )
    parser.add_argument("--n-boot", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for the result CSVs. Defaults to statistical_test_results.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_boot < 1000:
        raise SystemExit("Use at least 1000 bootstrap replicates.")
    _self_check()
    output_dir = args.output_dir or results_dir()
    sources = aspect_sources(measurement_dir())
    if args.aspect is not None:
        sources = [item for item in sources if item[0] == args.aspect]
    bin_frames = []
    summaries = []
    for aspect, path, angle_column in sources:
        if not path.is_file():
            raise SystemExit(f"Missing {aspect} data: {path}")
        bins, summary, _slopes = run(aspect, path, angle_column, args.n_boot, args.seed)
        bin_frames.append(bins)
        summaries.append(summary)
        print(render(bins, summary), end="")
        print()
    bin_path = output_dir / "cv_by_volume_bin.csv"
    slope_path = output_dir / "cv_volume_slope_bootstrap.csv"
    pd.concat(bin_frames, ignore_index=True).to_csv(bin_path, index=False)
    pd.DataFrame(summaries).to_csv(slope_path, index=False)
    print(f"Wrote {bin_path}")
    print(f"Wrote {slope_path}")


if __name__ == "__main__":
    main()
