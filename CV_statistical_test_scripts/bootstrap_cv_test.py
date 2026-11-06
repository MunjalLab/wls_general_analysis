#!/usr/bin/env python3
"""Bootstrap test of pillar-angle coefficients of variation.

The annotations on wls_pillarAngle.pdf are the coefficient of variation of
pillar angle for each region (anterior, posterior, ventral) and genotype
(tyrosinase, wls):

    CV = s / mean

where s is the sample standard deviation with divisor n - 1. This script
reads the three angle tables, reproduces those values, and tests

    H0: CV_wls = CV_tyrosinase

in each region. The test statistic is the difference

    CV_wls - CV_tyrosinase

so a positive value means the wls angles are more variable relative to
their mean.

Bootstrap test
--------------
The p-value is a two-sided nonparametric bootstrap test. The pooled CV

    tau = sum((n_g - 1) * CV_g) / sum(n_g - 1)

is the common value under H0. Each genotype is then rescaled,

    x' = mean(x) + (x - mean(x)) * (tau * mean(x) / s)

so the rescaled sample has the original mean and the original shape of
deviations, but CV exactly equal to tau. Bootstrap samples are drawn with
replacement from those rescaled observations (Davison and Hinkley, Bootstrap
Methods and their Application, 1997, chapter 4).

The statistic is the studentized difference

    Z = (CV_wls - CV_tyrosinase) / sqrt(SE_wls^2 + SE_tyrosinase^2)

where each SE is the jackknife standard error of that sample's CV. The
two-sided p-value is

    (1 + #{|Z*| >= |Z observed|}) / (B + 1).

The difference is studentized because, with 10 and 11 embryos, the sampling
distribution of a CV is skewed and its width depends on the CV itself. An
unstudentized bootstrap of the raw difference does not hold its nominal
false-positive rate in that setting; the jackknife standard error corrects it.

Confidence intervals
--------------------
The 95% interval for the CV difference is a bootstrap-t interval from the
same studentized statistic, resampling the original observations without
imposing H0. Individual CVs also get percentile and BCa 95% intervals
(Efron, Better bootstrap confidence intervals, JASA 1987).

Companion test
--------------
The Feltz and Miller (1996) statistic is also reported. It is asymptotically
chi-squared with 1 degree of freedom when the angles are normal, which is a
stronger assumption than the bootstrap uses. With two groups it is

    tau = sum(m_g * CV_g) / sum(m_g),   m_g = n_g - 1
    D = (sum(m_g * CV_g^2) - tau^2 * sum(m_g)) / (tau^2 * (0.5 + tau^2))

Holm adjustment is applied across the three regions for each test. The
regions are measured on the same embryos, so the three tests are not
independent and the Holm p-values are a conservative family-wise correction.
Both imaging dates in the tables are pooled, matching the figure. The script
also prints each date separately. Those splits have only five or six embryos,
so they are descriptive and are not a second test.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

from project_paths import measurement_dir, results_dir

REGIONS = ("Anterior", "Posterior", "Ventral")
FILES = {
    "Anterior": "anteriorAngleswls_01.csv",
    "Posterior": "posteriorAngleswls_01.csv",
    "Ventral": "ventralAngleswls_01.csv",
}
GENOTYPES = ("tyrosinase", "wls")


def cv(values: np.ndarray) -> float:
    """Sample coefficient of variation, divisor n - 1 in the standard deviation."""
    values = np.asarray(values, dtype=float)
    mean = float(values.mean())
    if mean == 0.0:
        raise ValueError("Coefficient of variation is undefined when the mean is 0.")
    return float(values.std(ddof=1) / mean)


def pooled_cv(first: np.ndarray, second: np.ndarray) -> float:
    """Degrees-of-freedom weighted mean of the two sample CVs."""
    weights = np.array([len(first) - 1, len(second) - 1], dtype=float)
    estimates = np.array([cv(first), cv(second)], dtype=float)
    return float(np.sum(weights * estimates) / np.sum(weights))


def rescale_to_cv(values: np.ndarray, target: float) -> np.ndarray:
    """Keep the mean and the shape of deviations; set the CV to ``target``."""
    values = np.asarray(values, dtype=float)
    mean = float(values.mean())
    scale = float(values.std(ddof=1))
    if mean == 0.0 or scale == 0.0:
        raise ValueError("Cannot rescale a sample with mean 0 or standard deviation 0.")
    return mean + (values - mean) * (target * mean / scale)


def feltz_miller(cvs: list[float], sizes: list[int]) -> tuple[float, float]:
    """Feltz-Miller asymptotic test of equal CVs.

    Returns the chi-squared statistic and its upper-tail p-value on
    k - 1 degrees of freedom. For the two-group comparisons in this
    folder, k - 1 = 1. The three-group check below exercises df = 2,
    where the survival function is exp(-x / 2).
    """
    cvs_arr = np.asarray(cvs, dtype=float)
    degrees = np.asarray(sizes, dtype=float) - 1.0
    if np.any(degrees <= 0):
        raise ValueError("Each group needs at least two observations.")
    total = float(degrees.sum())
    tau = float(np.sum(degrees * cvs_arr) / total)
    statistic = float(
        (np.sum(degrees * cvs_arr**2) - (tau**2) * total) / (tau**2 * (0.5 + tau**2))
    )
    df = len(cvs) - 1
    if df == 1:
        p_value = math.erfc(math.sqrt(max(statistic, 0.0) / 2.0))
    elif df == 2:
        p_value = math.exp(-max(statistic, 0.0) / 2.0)
    else:
        raise ValueError("This script only evaluates the chi-squared tail for 1 or 2 df.")
    return statistic, p_value


def holm(p_values: list[float]) -> list[float]:
    """Monotone Holm adjustment. Order is the same as the input."""
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    running = 0.0
    families = len(values)
    for rank, index in enumerate(order):
        step = min(1.0, (families - rank) * float(values[index]))
        running = max(running, step)
        adjusted[index] = running
    return adjusted.tolist()


def norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def norm_ppf(p: float) -> float:
    """Inverse standard normal CDF by bisection. Used only for BCa endpoints."""
    if not 0.0 < p < 1.0:
        raise ValueError(f"Probability must be strictly between 0 and 1, got {p}.")
    lo, hi = -20.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if norm_cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _acceleration(jackknife: np.ndarray) -> float:
    centered = float(jackknife.mean()) - jackknife
    denominator = float(np.sum(centered**2) ** 1.5)
    if denominator == 0.0:
        return 0.0
    return float(np.sum(centered**3) / (6.0 * denominator))


def bca_interval(
    estimate: float,
    replicates: np.ndarray,
    jackknife: np.ndarray,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """BCa interval from bootstrap replicates and leave-one-out statistics."""
    replicates = np.asarray(replicates, dtype=float)
    count = replicates.size
    proportion = float(np.mean(replicates < estimate))
    proportion = min(max(proportion, 1.0 / (count + 1.0)), count / (count + 1.0))
    bias = norm_ppf(proportion)
    accel = _acceleration(np.asarray(jackknife, dtype=float))

    def adjusted(normal_quantile: float) -> float:
        denominator = 1.0 - accel * (bias + normal_quantile)
        if denominator <= 1e-12:
            return 0.0 if normal_quantile < 0.0 else 1.0
        probability = norm_cdf(bias + (bias + normal_quantile) / denominator)
        return min(max(probability, 0.0), 1.0)

    lower_p = adjusted(norm_ppf(alpha / 2.0))
    upper_p = adjusted(norm_ppf(1.0 - alpha / 2.0))
    if upper_p < lower_p:
        lower_p, upper_p = upper_p, lower_p
    lower, upper = np.quantile(replicates, [lower_p, upper_p])
    return float(lower), float(upper)


def _resample(values: np.ndarray, rng: np.random.Generator, n_boot: int) -> np.ndarray:
    return values[rng.integers(0, len(values), size=(n_boot, len(values)))]


def _cvs_from_draws(draws: np.ndarray) -> np.ndarray:
    means = draws.mean(axis=1)
    if np.any(means == 0.0):
        raise ValueError("A bootstrap sample had mean 0, so its CV is undefined.")
    return draws.std(axis=1, ddof=1) / means


def _bootstrap_cvs(values: np.ndarray, rng: np.random.Generator, n_boot: int) -> np.ndarray:
    return _cvs_from_draws(_resample(values, rng, n_boot))


def jackknife_se_cv(draws: np.ndarray) -> np.ndarray:
    """Jackknife standard error of the sample CV.

    ``draws`` has shape (n,) or (n_replicates, n). A one-dimensional input
    returns a scalar; a matrix returns one standard error per row.
    """
    draws = np.asarray(draws, dtype=float)
    single = draws.ndim == 1
    if single:
        draws = draws[None, :]
    if draws.shape[1] < 3:
        raise ValueError("The jackknife standard error of the CV needs at least 3 observations.")
    leave = np.empty_like(draws)
    for index in range(draws.shape[1]):
        rest = np.delete(draws, index, axis=1)
        leave[:, index] = rest.std(axis=1, ddof=1) / rest.mean(axis=1)
    centered = leave - leave.mean(axis=1, keepdims=True)
    standard_error = np.sqrt((draws.shape[1] - 1) / draws.shape[1] * np.sum(centered**2, axis=1))
    return float(standard_error[0]) if single else standard_error


def studentized_difference(tyrosinase_draws: np.ndarray, wls_draws: np.ndarray) -> np.ndarray:
    """(CV_wls - CV_tyrosinase) / jackknife SE, for one sample or a batch of resamples."""
    difference = _cvs_from_draws(wls_draws) - _cvs_from_draws(tyrosinase_draws)
    standard_error = np.sqrt(
        jackknife_se_cv(tyrosinase_draws) ** 2 + jackknife_se_cv(wls_draws) ** 2
    )
    return difference / standard_error


def _jackknife_cv(values: np.ndarray) -> np.ndarray:
    return np.array([cv(np.delete(values, i)) for i in range(len(values))], dtype=float)


def _jackknife_difference(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """Leave-one-out values of CV(second) - CV(first). ``first`` is tyrosinase."""
    cv_first = cv(first)
    cv_second = cv(second)
    from_first = np.array(
        [cv_second - cv(np.delete(first, i)) for i in range(len(first))],
        dtype=float,
    )
    from_second = np.array(
        [cv(np.delete(second, j)) - cv_first for j in range(len(second))],
        dtype=float,
    )
    return np.concatenate([from_first, from_second])


def compare_cvs(
    tyrosinase: np.ndarray,
    wls: np.ndarray,
    rng: np.random.Generator,
    n_boot: int,
) -> dict[str, float]:
    """Bootstrap comparison of CV_wls and CV_tyrosinase for one region."""
    if n_boot < 1000:
        raise ValueError("Use at least 1000 bootstrap replicates.")
    cv_tyrosinase = cv(tyrosinase)
    cv_wls = cv(wls)
    difference = cv_wls - cv_tyrosinase
    ratio = cv_wls / cv_tyrosinase
    common = pooled_cv(tyrosinase, wls)

    tyrosinase_draws = _resample(tyrosinase, rng, n_boot)
    wls_draws = _resample(wls, rng, n_boot)
    tyrosinase_boot = _cvs_from_draws(tyrosinase_draws)
    wls_boot = _cvs_from_draws(wls_draws)
    difference_boot = wls_boot - tyrosinase_boot
    ratio_boot = wls_boot / tyrosinase_boot
    observed_statistic = float(
        studentized_difference(tyrosinase[None, :], wls[None, :])[0]
    )
    bootstrap_t = (difference_boot - difference) / np.sqrt(
        jackknife_se_cv(tyrosinase_draws) ** 2 + jackknife_se_cv(wls_draws) ** 2
    )

    null_tyrosinase = rescale_to_cv(tyrosinase, common)
    null_wls = rescale_to_cv(wls, common)
    null_statistic = studentized_difference(
        _resample(null_tyrosinase, rng, n_boot),
        _resample(null_wls, rng, n_boot),
    )
    if (
        np.isfinite(null_statistic).mean() < 0.99
        or np.isfinite(bootstrap_t).mean() < 0.99
    ):
        raise RuntimeError("More than 1% of bootstrap replicates had a zero jackknife SE.")
    null_statistic = null_statistic[np.isfinite(null_statistic)]
    bootstrap_t = bootstrap_t[np.isfinite(bootstrap_t)]
    exceed = int(np.sum(np.abs(null_statistic) >= abs(observed_statistic)))
    bootstrap_p = (exceed + 1) / (null_statistic.size + 1)
    # Bootstrap-t interval: theta is inside
    # diff - se * t_(1-a/2)  to  diff - se * t_(a/2).
    t_low, t_high = np.quantile(bootstrap_t, [0.025, 0.975])
    observed_se = float(np.sqrt(jackknife_se_cv(tyrosinase) ** 2 + jackknife_se_cv(wls) ** 2))
    bootstrap_t_low = float(difference - observed_se * t_high)
    bootstrap_t_high = float(difference - observed_se * t_low)

    statistic, asymptotic_p = feltz_miller(
        [cv_tyrosinase, cv_wls],
        [len(tyrosinase), len(wls)],
    )
    diff_bca = bca_interval(difference, difference_boot, _jackknife_difference(tyrosinase, wls))
    tyrosinase_bca = bca_interval(cv_tyrosinase, tyrosinase_boot, _jackknife_cv(tyrosinase))
    wls_bca = bca_interval(cv_wls, wls_boot, _jackknife_cv(wls))
    diff_percentile = np.quantile(difference_boot, [0.025, 0.975])
    ratio_percentile = np.quantile(ratio_boot, [0.025, 0.975])

    return {
        "n_tyrosinase": len(tyrosinase),
        "n_wls": len(wls),
        "mean_tyrosinase": float(tyrosinase.mean()),
        "mean_wls": float(wls.mean()),
        "sd_tyrosinase": float(tyrosinase.std(ddof=1)),
        "sd_wls": float(wls.std(ddof=1)),
        "cv_tyrosinase": cv_tyrosinase,
        "cv_wls": cv_wls,
        "cv_tyrosinase_bca95_low": tyrosinase_bca[0],
        "cv_tyrosinase_bca95_high": tyrosinase_bca[1],
        "cv_wls_bca95_low": wls_bca[0],
        "cv_wls_bca95_high": wls_bca[1],
        "pooled_cv_under_null": common,
        "cv_difference_wls_minus_tyrosinase": difference,
        "studentized_difference": observed_statistic,
        "difference_bootstrap_t95_low": bootstrap_t_low,
        "difference_bootstrap_t95_high": bootstrap_t_high,
        "difference_percentile95_low": float(diff_percentile[0]),
        "difference_percentile95_high": float(diff_percentile[1]),
        "difference_bca95_low": diff_bca[0],
        "difference_bca95_high": diff_bca[1],
        "cv_ratio_wls_over_tyrosinase": ratio,
        "ratio_percentile95_low": float(ratio_percentile[0]),
        "ratio_percentile95_high": float(ratio_percentile[1]),
        "bootstrap_p": float(bootstrap_p),
        "feltz_miller_statistic": statistic,
        "feltz_miller_p": float(asymptotic_p),
        "null_statistic_median": float(np.median(null_statistic)),
    }


def genotype_of(filename: str) -> str:
    lowered = filename.lower()
    if "tyrosinase" in lowered:
        return "tyrosinase"
    if "wls" in lowered:
        return "wls"
    raise ValueError(f"No genotype (wls or tyrosinase) found in filename: {filename}")


def load_region(path: Path) -> dict[str, np.ndarray]:
    table = pd.read_csv(path)
    expected = {"Timepoint", "Angle", "File"}
    missing = expected.difference(table.columns)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")
    grouped = {name: [] for name in GENOTYPES}
    for angle, filename in zip(table["Angle"], table["File"]):
        grouped[genotype_of(str(filename))].append(float(angle))
    arrays = {name: np.asarray(values, dtype=float) for name, values in grouped.items()}
    if any(len(values) < 3 for values in arrays.values()):
        raise ValueError(f"{path.name} needs at least three embryos in each genotype.")
    return arrays


def load_all(data_dir: Path) -> dict[str, dict[str, np.ndarray]]:
    loaded = {}
    for region in REGIONS:
        path = data_dir / FILES[region]
        if not path.is_file():
            raise FileNotFoundError(path)
        loaded[region] = load_region(path)
    return loaded


def run(data_dir: Path, n_boot: int, seed: int) -> pd.DataFrame:
    measurements = load_all(data_dir)
    rng = np.random.default_rng(seed)
    rows = []
    for region in REGIONS:
        result = compare_cvs(
            measurements[region]["tyrosinase"],
            measurements[region]["wls"],
            rng,
            n_boot,
        )
        result["region"] = region
        result["n_bootstrap"] = n_boot
        result["seed"] = seed
        rows.append(result)
    frame = pd.DataFrame(rows)
    frame["bootstrap_p_holm"] = holm(frame["bootstrap_p"].tolist())
    frame["feltz_miller_p_holm"] = holm(frame["feltz_miller_p"].tolist())
    leading = ["region", "n_tyrosinase", "n_wls", "cv_tyrosinase", "cv_wls"]
    ordered = leading + [column for column in frame.columns if column not in leading]
    return frame[ordered]


def imaging_date(filename: str) -> str:
    return filename.split("_", 1)[0]


def summarize_by_date(data_dir: Path) -> pd.DataFrame:
    """CV within each region, genotype, and imaging date. Descriptive only."""
    rows = []
    for region in REGIONS:
        table = pd.read_csv(data_dir / FILES[region])
        table["genotype"] = table["File"].map(lambda name: genotype_of(str(name)))
        table["imaging_date"] = table["File"].map(lambda name: imaging_date(str(name)))
        for (genotype, date), subset in table.groupby(["genotype", "imaging_date"], sort=False):
            values = subset["Angle"].to_numpy(dtype=float)
            rows.append(
                {
                    "region": region,
                    "genotype": genotype,
                    "imaging_date": date,
                    "n": int(len(values)),
                    "mean": float(values.mean()),
                    "sd": float(values.std(ddof=1)),
                    "cv": cv(values),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["region", "imaging_date", "genotype"], kind="mergesort"
    ).reset_index(drop=True)


def _format_p(p_value: float) -> str:
    if p_value < 1e-4:
        return f"{p_value:.2e}"
    return f"{p_value:.4f}"


def render(frame: pd.DataFrame) -> str:
    lines = [
        "Pillar-angle coefficient of variation: wls versus tyrosinase",
        "CV = sample standard deviation (divisor n - 1) / mean.",
        (
            f"Bootstrap replicates: {int(frame['n_bootstrap'].iloc[0])}. "
            f"Seed: {int(frame['seed'].iloc[0])}."
        ),
        "Difference = CV_wls - CV_tyrosinase. Positive means wls is more variable.",
        "p-values are two-sided. Holm adjustment is across the three regions,",
        "which share the same embryos, so that adjustment is conservative.",
        "",
    ]
    for row in frame.itertuples(index=False):
        lines.extend(
            [
                row.region,
                (
                    f"  tyrosinase  n={row.n_tyrosinase:2d}  "
                    f"mean={row.mean_tyrosinase:.2f}  sd={row.sd_tyrosinase:.2f}  "
                    f"CV={row.cv_tyrosinase:.4f}  "
                    f"95% BCa CI [{row.cv_tyrosinase_bca95_low:.4f}, "
                    f"{row.cv_tyrosinase_bca95_high:.4f}]"
                ),
                (
                    f"  wls         n={row.n_wls:2d}  "
                    f"mean={row.mean_wls:.2f}  sd={row.sd_wls:.2f}  "
                    f"CV={row.cv_wls:.4f}  "
                    f"95% BCa CI [{row.cv_wls_bca95_low:.4f}, {row.cv_wls_bca95_high:.4f}]"
                ),
                (
                    f"  difference = {row.cv_difference_wls_minus_tyrosinase:.4f}  "
                    f"95% bootstrap-t CI [{row.difference_bootstrap_t95_low:.4f}, "
                    f"{row.difference_bootstrap_t95_high:.4f}]"
                ),
                (
                    f"  ratio (wls / tyrosinase) = {row.cv_ratio_wls_over_tyrosinase:.3f}  "
                    f"95% percentile CI [{row.ratio_percentile95_low:.3f}, "
                    f"{row.ratio_percentile95_high:.3f}]"
                ),
                (
                    f"  studentized bootstrap p = {_format_p(row.bootstrap_p)}   "
                    f"Holm p = {_format_p(row.bootstrap_p_holm)}"
                ),
                (
                    f"  Feltz-Miller D = {row.feltz_miller_statistic:.3f}   "
                    f"p = {_format_p(row.feltz_miller_p)}   "
                    f"Holm p = {_format_p(row.feltz_miller_p_holm)}"
                ),
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def render_dates(by_date: pd.DataFrame) -> str:
    lines = [
        "",
        "Coefficient of variation by imaging date",
        "Descriptive only. Each date has 5 or 6 embryos, so these are not separate tests.",
        "",
    ]
    for region in REGIONS:
        lines.append(region)
        subset = by_date[by_date["region"] == region]
        for row in subset.itertuples(index=False):
            lines.append(
                f"  {row.genotype:12} {row.imaging_date:10}  n={row.n:2d}  "
                f"mean={row.mean:.2f}  sd={row.sd:.2f}  CV={row.cv:.4f}"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _self_check() -> None:
    """Lock the CV definition, the null rescaling, and the Feltz-Miller formula."""
    values = np.array([1.0, 2.0, 3.0, 4.0])
    expected = values.std(ddof=1) / values.mean()
    if abs(cv(values) - expected) > 1e-12:
        raise AssertionError("CV does not match sample SD / mean.")
    rescaled = rescale_to_cv(values, 0.2)
    if abs(cv(rescaled) - 0.2) > 1e-12 or abs(rescaled.mean() - values.mean()) > 1e-12:
        raise AssertionError("Rescaling did not preserve the mean and set the CV.")
    # cvequality vignette, Feltz and Miller's three-assay example:
    # statistic 5.53045, chi-squared df = 2, p = 0.06296.
    statistic, p_value = feltz_miller([0.090, 0.462, 0.340], [5, 5, 5])
    if abs(statistic - 5.5304524) > 1e-4 or abs(p_value - 0.0629619) > 1e-4:
        raise AssertionError(f"Feltz-Miller check failed: D={statistic}, p={p_value}")
    adjusted = holm([0.01, 0.04, 0.03])
    if not np.allclose(adjusted, [0.03, 0.06, 0.06]):
        raise AssertionError(f"Holm check failed: {adjusted}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Bootstrap test of pillar-angle CVs, wls versus tyrosinase."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory containing the three angle CSV files.",
    )
    parser.add_argument("--n-boot", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="CSV path for the results. Defaults to cv_bootstrap_results.csv in statistical_test_results.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _self_check()
    data_dir = args.data_dir or measurement_dir()
    frame = run(data_dir, args.n_boot, args.seed)
    by_date = summarize_by_date(data_dir)
    output = args.output or (results_dir() / "cv_bootstrap_results.csv")
    date_output = output.with_name("cv_by_imaging_date.csv")
    frame.to_csv(output, index=False)
    by_date.to_csv(date_output, index=False)
    print(render(frame), end="")
    print(render_dates(by_date), end="")
    print(f"Wrote {output}")
    print(f"Wrote {date_output}")


if __name__ == "__main__":
    main()
