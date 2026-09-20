"""Unit tests for the launcher. No Nextflow needed."""
from __future__ import annotations

import os
import stat
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from paludamentum import launcher

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"


def make_args(tmp_path: Path, params: dict, **overrides) -> SimpleNamespace:
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump(params))
    values = dict(
        params_yaml=str(params_file), nf_config="local", nextflow_args=[], work_dir=None,
        profile=None, nextflow_bin="nextflow", resume=False, check_tools=False,
        skip_singularity_check=True, dry_run=False,
    )
    values.update(overrides)
    return SimpleNamespace(**values)


def fake_executable(directory: Path, name: str, body: str = "exit 0") -> Path:
    path = directory / name
    path.write_text(f"#!/bin/sh\n{body}\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


@pytest.fixture
def fake_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """PATH with fake nextflow, java (version 17) and python3."""
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    fake_executable(bindir, "nextflow")
    fake_executable(bindir, "java", "echo 'openjdk version \"17.0.1\" 2021-10-19' >&2")
    fake_executable(bindir, "python3")
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}/usr/bin{os.pathsep}/bin")
    return bindir


# ---------------------------------------------------------------- layout

def test_pipeline_paths_point_into_the_checkout():
    root, main_nf, base_config = launcher.pipeline_paths()
    assert root == ROOT
    assert main_nf == ROOT / "main.nf" and main_nf.is_file()
    assert base_config == ROOT / "conf" / "base.config" and base_config.is_file()


def test_pipeline_paths_root_override(tmp_path: Path):
    root, main_nf, base_config = launcher.pipeline_paths(tmp_path)
    assert (root, main_nf, base_config) == (tmp_path, tmp_path / "main.nf", tmp_path / "conf" / "base.config")


def test_bin_scripts_are_executable():
    scripts = [p for p in (ROOT / "bin").iterdir() if p.is_file()]
    assert scripts
    assert not [p.name for p in scripts if not os.access(p, os.X_OK)]


def test_base_config_references_existing_scoring_matrix():
    text = (ROOT / "conf" / "base.config").read_text()
    assert '"${projectDir}/conf/blosum62.csv"' in text
    assert (ROOT / "conf" / "blosum62.csv").is_file()


# ---------------------------------------------------------------- resolve_nf_config

@pytest.mark.parametrize("value", [
    "slurm_generic", "slurm_generic.config", "conf/slurm_generic.config",
    str(ROOT / "conf" / "slurm_generic.config"),
])
def test_resolve_nf_config_accepts_names_and_paths(value: str, tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # 'conf/...' does not exist relative to here
    assert launcher.resolve_nf_config(value) == ROOT / "conf" / "slurm_generic.config"


def test_resolve_nf_config_prefers_existing_user_file(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    mine = tmp_path / "local.config"
    mine.write_text("// mine\n")
    assert launcher.resolve_nf_config("local.config") == mine


def test_resolve_nf_config_missing(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit, match="Nextflow config not found"):
        launcher.resolve_nf_config("no_such_cluster")


# ---------------------------------------------------------------- input validation

def test_validate_input_data_ok(tmp_path: Path):
    params = {
        "genome": str(DATA / "tiny.fa"),
        "proteins": [str(DATA / "tiny_proteins.faa")],
        "rnaseq_paired": str(DATA / "reads_{1,2}.fastq"),
        "isoseq": str(DATA / "iso*.fastq"),
    }
    assert launcher.validate_input_data(params, tmp_path / "params.yaml") == []


def test_validate_input_data_reports_problems(tmp_path: Path):
    errors = launcher.validate_input_data({}, tmp_path / "params.yaml")
    assert any("'genome'" in e for e in errors)
    errors = launcher.validate_input_data(
        {"genome": str(DATA / "tiny.fa"), "rnaseq_single": str(DATA / "nothing*.fq")},
        tmp_path / "params.yaml",
    )
    assert any("No files matched" in e for e in errors)


def test_relative_paths_are_resolved_against_the_launch_directory(tmp_path: Path, monkeypatch):
    """Nextflow resolves relative paths against the launch directory; so does the validation."""
    launch = tmp_path / "launch"
    launch.mkdir()
    (launch / "genome.fa").write_text(">s\nACGT\n")
    (launch / "r_1.fq").write_text("")
    (launch / "r_2.fq").write_text("")
    elsewhere = tmp_path / "out"
    elsewhere.mkdir()
    monkeypatch.chdir(launch)
    params = {"genome": "genome.fa", "rnaseq_paired": "r_{1,2}.fq"}
    assert launcher.validate_input_data(params, elsewhere / "params.yaml") == []
    monkeypatch.chdir(elsewhere)
    assert launcher.validate_input_data(params, launch / "params.yaml")


def test_expand_braces():
    assert launcher.expand_braces("a_{1,2}.fq") == ["a_1.fq", "a_2.fq"]


# ---------------------------------------------------------------- launching

def test_command_assembly(tmp_path: Path, fake_path: Path, monkeypatch):
    calls = []
    real_run = launcher.subprocess.run

    def spy(cmd, *a, **kw):
        if Path(str(cmd[0])).name == "nextflow":
            calls.append(list(cmd))
            return SimpleNamespace(returncode=0)
        return real_run(cmd, *a, **kw)  # the java -version probe

    monkeypatch.setattr(launcher.subprocess, "run", spy)
    args = make_args(
        tmp_path, {"genome": str(DATA / "tiny.fa")},
        resume=True, profile="test", work_dir=str(tmp_path / "w"),
        nextflow_args=["--", "-with-dag", "dag.html"],
    )
    launcher.run_nextflow_pipeline(args, genefinder="tiberius")

    assert len(calls) == 1
    cmd = calls[0]
    assert cmd[:3] == ["nextflow", "run", str(ROOT / "main.nf")]
    assert cmd[cmd.index("-params-file") + 1] == args.params_yaml
    configs = [cmd[i + 1] for i, part in enumerate(cmd) if part == "-c"]
    assert configs == [str(ROOT / "conf" / "base.config"), str(ROOT / "conf" / "local.config")]
    assert cmd[-7:] == ["-profile", "test", "-resume", "-work-dir", str(tmp_path / "w"), "-with-dag", "dag.html"]


def test_base_config_is_not_passed_twice(tmp_path: Path, fake_path: Path, monkeypatch):
    calls = []
    real_run = launcher.subprocess.run

    def spy(cmd, *a, **kw):
        if Path(str(cmd[0])).name == "nextflow":
            calls.append(list(cmd))
            return SimpleNamespace(returncode=0)
        return real_run(cmd, *a, **kw)

    monkeypatch.setattr(launcher.subprocess, "run", spy)
    args = make_args(tmp_path, {"genome": str(DATA / "tiny.fa")}, nf_config="base.config")
    launcher.run_nextflow_pipeline(args)
    assert calls[0].count("-c") == 1


def test_dry_run_does_not_start_nextflow(tmp_path: Path, fake_path: Path, monkeypatch):
    real_run = launcher.subprocess.run

    def guard(cmd, *a, **kw):
        assert Path(str(cmd[0])).name != "nextflow", "nextflow must not be started in a dry run"
        return real_run(cmd, *a, **kw)

    monkeypatch.setattr(launcher.subprocess, "run", guard)
    args = make_args(tmp_path, {"genome": str(DATA / "tiny.fa")}, dry_run=True)
    launcher.run_nextflow_pipeline(args)


def test_nextflow_failure_is_propagated(tmp_path: Path, fake_path: Path):
    fake_executable(fake_path, "nextflow", "exit 3")
    args = make_args(tmp_path, {"genome": str(DATA / "tiny.fa")})
    with pytest.raises(SystemExit) as excinfo:
        launcher.run_nextflow_pipeline(args)
    assert excinfo.value.code == 3


def test_missing_genome_fails_validation(tmp_path: Path, fake_path: Path):
    args = make_args(tmp_path, {"genome": str(tmp_path / "missing.fa")})
    with pytest.raises(SystemExit, match="Input validation failed"):
        launcher.run_nextflow_pipeline(args)


def test_unknown_genefinder(tmp_path: Path):
    args = make_args(tmp_path, {"genome": str(DATA / "tiny.fa")})
    with pytest.raises(SystemExit, match="Unknown gene finder"):
        launcher.run_nextflow_pipeline(args, genefinder="augustus")


@pytest.mark.parametrize("finder,cli", sorted(launcher.GENEFINDER_CLI.items()))
def test_check_tools_requires_the_genefinder_cli(finder: str, cli: str, tmp_path: Path, fake_path: Path):
    params = {finder: {"run": True}}
    errors, _ = launcher.validate_executables(
        params, nextflow_bin="nextflow", check_tool_binaries=True,
        skip_singularity_check=True, repo_root=ROOT,
    )
    assert any(f"({cli})" in e for e in errors)
    fake_executable(fake_path, cli)
    errors, _ = launcher.validate_executables(
        params, nextflow_bin="nextflow", check_tool_binaries=True,
        skip_singularity_check=True, repo_root=ROOT,
    )
    assert not any(f"({cli})" in e for e in errors)


def test_old_java_is_rejected(tmp_path: Path, fake_path: Path):
    fake_executable(fake_path, "java", "echo 'java version \"1.8.0_292\"' >&2")
    errors, _ = launcher.validate_executables(
        {}, nextflow_bin="nextflow", check_tool_binaries=False,
        skip_singularity_check=True, repo_root=ROOT,
    )
    assert any("11+" in e for e in errors)


# ---------------------------------------------------------------- params helpers

@pytest.mark.parametrize("finder", sorted(launcher.GENEFINDER_CLI))
def test_default_params_switch_on_one_gene_finder(finder: str):
    params = launcher.default_params(finder)
    assert params["genefinder"] == finder
    assert params["outdir"] == f"{finder}_results"
    assert {name: params[name]["run"] for name in launcher.GENEFINDER_CLI} == {
        name: name == finder for name in launcher.GENEFINDER_CLI
    }
    assert Path(params["scoring_matrix"]) == ROOT / "conf" / "blosum62.csv"


def test_default_params_unknown_gene_finder():
    with pytest.raises(SystemExit, match="Unknown gene finder"):
        launcher.default_params("augustus")


def test_merge_params_is_recursive_and_ignores_none():
    base = {"threads": 48, "vipsania": {"run": True, "model": None}, "proteins": None}
    launcher.merge_params(base, {"threads": None, "vipsania": {"model": "Fungi"}, "proteins": ["p.faa"]})
    assert base == {"threads": 48, "vipsania": {"run": True, "model": "Fungi"}, "proteins": ["p.faa"]}


def test_write_params_yaml(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    params = launcher.default_params("vipsania")
    path = launcher.write_params_yaml(params)
    assert path == tmp_path / "vipsania_results" / "params.yaml"
    written = yaml.safe_load(path.read_text())
    assert written["outdir"] == str(tmp_path / "vipsania_results")
    assert written["vipsania"] == {"run": True}
