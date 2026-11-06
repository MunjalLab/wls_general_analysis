"""Locations of the measurement tables, test results, and the presentation."""

from pathlib import Path


def project_root() -> Path:
    here = Path(__file__).resolve().parent
    if here.name == "statistical_test_scripts":
        return here.parent
    return here


def measurement_dir() -> Path:
    candidate = project_root() / "measurement_values"
    return candidate if candidate.is_dir() else project_root()


def results_dir() -> Path:
    candidate = project_root() / "statistical_test_results"
    return candidate if candidate.is_dir() else project_root()


def presentation_dir() -> Path:
    candidate = project_root() / "presentation"
    return candidate if candidate.is_dir() else project_root()
