"""conf/hc_genes.yaml: the clade table (classifier_clade per clade). Pure Python, no Nextflow."""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml

TABLE = Path(__file__).resolve().parent.parent / "conf" / "hc_genes.yaml"


@pytest.fixture(scope="module")
def clades() -> dict:
    return yaml.safe_load(TABLE.read_text())["clades"]


@pytest.mark.parametrize("clade, classifier_clade", [
    ("Mesangiospermae", "embryophyta"),
    ("Fungi", "fungi"),
    ("Vertebrata", "vertebrates"),
    ("Mammalia", "vertebrates"),
    ("Insecta", "insects"),
    ("Lepidoptera", "insects"),
])
def test_classifier_clade(clades: dict, clade: str, classifier_clade: str) -> None:
    assert clades[clade]["classifier_clade"] == classifier_clade


@pytest.mark.parametrize("clade", ["Chlorophyta", "Bacillariophyta"])
def test_clades_unseen_by_the_classifier_are_absent(clades: dict, clade: str) -> None:
    assert clade.lower() not in {k.lower() for k in clades}


def test_no_drusilla_model_key(clades: dict) -> None:
    """Drusilla runs one model for all clades (params.drusilla.model); the table does not choose it."""
    assert not [k for k, v in clades.items() if "drusilla_model" in (v or {})]
    assert "drusilla_model" not in yaml.safe_load(TABLE.read_text())
