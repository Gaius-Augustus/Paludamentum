"""Unit tests for the launcher. No Nextflow needed."""
from __future__ import annotations

import os
import stat
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from paludamentum import cli, launcher

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
TIBERIUS = ROOT / "tiberius"
needs_tiberius = pytest.mark.skipif(
    not (TIBERIUS / "tiberius.py").is_file(), reason="Tiberius submodule not checked out"
)


def make_args(tmp_path: Path, params: dict, **overrides) -> SimpleNamespace:
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump(params))
    values = dict(
        params_yaml=str(params_file), nf_config="local", nextflow_args=[], work_dir=None,
        profile=None, nextflow_bin="nextflow", resume=False, check_tools=False,
        skip_singularity_check=True, dry_run=False, genefinder=None,
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
    """Scripts (with a shebang) are executable; imported modules are not."""
    files = [p for p in (ROOT / "bin").iterdir() if p.is_file()]
    assert files
    for path in files:
        has_shebang = path.read_bytes().startswith(b"#!")
        assert os.access(path, os.X_OK) == has_shebang, path.name


def test_pipeline_root_from_environment(tmp_path: Path, monkeypatch):
    """pip install . copies only the launcher; PALUDAMENTUM_ROOT names the checkout then."""
    monkeypatch.setenv(launcher.ROOT_ENV, str(ROOT))
    assert launcher.pipeline_paths()[0] == ROOT
    monkeypatch.setenv(launcher.ROOT_ENV, str(tmp_path))
    with pytest.raises(SystemExit, match="not a Paludamentum checkout"):
        launcher.pipeline_paths()


def test_installed_launcher_without_checkout_explains(tmp_path: Path, monkeypatch):
    fake_site = tmp_path / "site-packages" / "paludamentum"
    fake_site.mkdir(parents=True)
    monkeypatch.delenv(launcher.ROOT_ENV, raising=False)
    monkeypatch.setattr(launcher, "__file__", str(fake_site / "launcher.py"))
    with pytest.raises(SystemExit, match="pip install -e"):
        launcher.pipeline_paths()


def test_base_config_references_existing_scoring_matrix():
    text = (ROOT / "conf" / "base.config").read_text()
    assert '"${projectDir}/conf/blosum62.csv"' in text
    assert (ROOT / "conf" / "blosum62.csv").is_file()


# ---------------------------------------------------------------- submodules

def test_submodules_are_declared():
    text = (ROOT / ".gitmodules").read_text()
    for name, path in launcher.SUBMODULES.items():
        assert f"path = {path}\n" in text, name
        assert launcher.submodule_root(name) == ROOT / path


def test_container_tags_of_base_config():
    tags = launcher.container_tags(ROOT / "conf" / "base.config")
    assert set(tags) >= {"tiberius", "vipsania", "drusilla"}
    assert tags["drusilla"] != "latest"


@pytest.mark.skipif(
    not all(launcher.submodule_version(n) for n in launcher.SUBMODULES),
    reason="submodules not checked out",
)
def test_submodule_versions_match_the_image_tags():
    """The submodule pins the version whose image conf/base.config runs."""
    assert launcher.version_mismatches() == []


def test_version_mismatch_is_reported(tmp_path: Path):
    (tmp_path / "conf").mkdir()
    (tmp_path / "conf" / "base.config").write_text('container = "docker://x/tiberius:9.9.9"\n')
    (tmp_path / "tiberius").mkdir()
    (tmp_path / "tiberius" / "pyproject.toml").write_text('[project]\nversion = "2.0.7"\n')
    problems = launcher.version_mismatches(tmp_path)
    assert len(problems) == 1 and "2.0.7" in problems[0] and "9.9.9" in problems[0]


@needs_tiberius
def test_nextflow_env_appends_a_callable_tiberius_py(tmp_path: Path, monkeypatch):
    """tiberius.py is not executable in the Tiberius repository: a wrapper makes it callable."""
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    env = launcher.nextflow_env()
    bin_dir = Path(env["PATH"].split(os.pathsep)[-1])
    assert env["PATH"].startswith(os.environ["PATH"])
    script = bin_dir / "tiberius.py"
    assert os.access(script, os.X_OK)
    if bin_dir != TIBERIUS:
        assert str(TIBERIUS / "tiberius.py") in script.read_text()


def test_nextflow_env_without_checkout(tmp_path: Path):
    assert launcher.nextflow_env(tmp_path)["PATH"] == os.environ["PATH"]


# ---------------------------------------------------------------- resolve_model_cfg

@needs_tiberius
@pytest.mark.parametrize("value", ["diatoms", "diatoms.yaml", "diatoms.yml"])
def test_resolve_model_cfg_by_name(value: str, tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert launcher.resolve_model_cfg(value) == TIBERIUS / "model_cfg" / "diatoms.yaml"


def test_resolve_model_cfg_path(tmp_path: Path):
    mine = tmp_path / "my_model.yaml"
    mine.write_text("target_species: x\n")
    assert launcher.resolve_model_cfg(str(mine), tmp_path) == mine


@needs_tiberius
def test_resolve_model_cfg_unknown_name_lists_the_names():
    with pytest.raises(SystemExit, match="Available names: .*diatoms"):
        launcher.resolve_model_cfg("no_such_clade")


def test_resolve_model_cfg_without_submodule(tmp_path: Path):
    with pytest.raises(SystemExit, match="git submodule update --init tiberius"):
        launcher.resolve_model_cfg("diatoms", tmp_path)


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


def test_resolve_nf_config_mistyped_path_does_not_fall_back(tmp_path: Path, monkeypatch):
    """~/x/slurm_generic.config must not silently become conf/slurm_generic.config."""
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit, match="Nextflow config not found"):
        launcher.resolve_nf_config(str(tmp_path / "x" / "slurm_generic.config"))


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


def test_validate_gene_finder_files(tmp_path: Path):
    genome = str(DATA / "tiny.fa")
    errors = launcher.validate_input_data(
        {"genome": genome, "tiberius": {"run": True, "result": str(tmp_path / "no.gtf")},
         "vipsania": {"model_dir": str(tmp_path / "nodir")}, "rnaseq_bam": [str(tmp_path / "no.bam")]},
        tmp_path / "params.yaml")
    assert any("Tiberius result missing" in e for e in errors)
    assert any("Vipsania model_dir missing" in e for e in errors)
    assert any("RNA-Seq BAM missing" in e for e in errors)
    (tmp_path / "models").mkdir()
    ok = launcher.validate_input_data(
        {"genome": genome, "tiberius": {"run": True, "result": genome},
         "vipsania": {"model_dir": str(tmp_path / "models")}, "rnaseq_bam": genome},
        tmp_path / "params.yaml")
    assert ok == []


def test_paths_with_whitespace_are_rejected(tmp_path: Path):
    genome = tmp_path / "my genome.fa"
    genome.write_text(">s\nACGT\n")
    errors = launcher.validate_input_data({"genome": str(genome)}, tmp_path / "params.yaml")
    assert any("whitespace" in e for e in errors)


def test_non_path_value_is_an_error_not_a_traceback(tmp_path: Path):
    with pytest.raises(SystemExit, match="'proteins'"):
        launcher.validate_input_data({"genome": str(DATA / "tiny.fa"), "proteins": 42}, tmp_path / "params.yaml")


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
    # tmp_path as root: no submodule checkout that could satisfy the check
    errors, _ = launcher.validate_executables(
        params, nextflow_bin="nextflow", check_tool_binaries=True,
        skip_singularity_check=True, repo_root=tmp_path,
    )
    assert any(f"({cli})" in e for e in errors)
    fake_executable(fake_path, cli)
    errors, _ = launcher.validate_executables(
        params, nextflow_bin="nextflow", check_tool_binaries=True,
        skip_singularity_check=True, repo_root=tmp_path,
    )
    assert not any(f"({cli})" in e for e in errors)


@needs_tiberius
def test_check_tools_accepts_tiberius_from_the_checkout(fake_path: Path):
    """Without an installed tiberius.py, the submodule's copy is used (appended to PATH)."""
    errors, _ = launcher.validate_executables(
        {"tiberius": {"run": True}}, nextflow_bin="nextflow", check_tool_binaries=True,
        skip_singularity_check=True, repo_root=ROOT,
    )
    assert not any("(tiberius.py)" in e for e in errors)


def test_old_java_is_rejected(tmp_path: Path, fake_path: Path):
    fake_executable(fake_path, "java", "echo 'java version \"1.8.0_292\"' >&2")
    errors, _ = launcher.validate_executables(
        {}, nextflow_bin="nextflow", check_tool_binaries=False,
        skip_singularity_check=True, repo_root=ROOT,
    )
    assert any("17+" in e for e in errors)


def test_apptainer_satisfies_the_container_check(tmp_path: Path, fake_path: Path):
    errors, _ = launcher.validate_executables(
        {}, nextflow_bin="nextflow", check_tool_binaries=False,
        skip_singularity_check=False, repo_root=ROOT,
    )
    assert any("singularity or apptainer" in e for e in errors)
    fake_executable(fake_path, "apptainer")
    errors, _ = launcher.validate_executables(
        {}, nextflow_bin="nextflow", check_tool_binaries=False,
        skip_singularity_check=False, repo_root=ROOT,
    )
    assert errors == []


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


# ---------------------------------------------------------------- build_params

def cli_args(*argv: str):
    return cli.build_parser().parse_args(list(argv))


@needs_tiberius
def test_build_params_tiberius_from_the_command_line(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    finder, path = launcher.build_params(cli_args(
        "--genome", str(DATA / "tiny.fa"), "--model_cfg", "diatoms",
        "--proteins", str(DATA / "tiny_proteins.faa"), "--rnaseq_paired", "r_{1,2}.fq",
        "--max_parallel", "1",
    ))
    assert finder == "tiberius"
    assert path == tmp_path / "tiberius_results" / "params.yaml"
    params = yaml.safe_load(path.read_text())
    assert params["genefinder"] == "tiberius"
    assert params["tiberius"] == {
        "run": True, "model_cfg": str(TIBERIUS / "model_cfg" / "diatoms.yaml"), "max_parallel": 1,
    }
    assert params["vipsania"] == {"run": False}
    assert params["proteins"] == [str(DATA / "tiny_proteins.faa")]
    # one value is a glob, not a list; relative paths are written absolute
    assert params["rnaseq_paired"] == str(tmp_path / "r_{1,2}.fq")
    assert params["outdir"] == str(tmp_path / "tiberius_results")


def test_build_params_vipsania_from_the_command_line(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    finder, path = launcher.build_params(cli_args(
        "--genome", str(DATA / "tiny.fa"), "--model", "Fungi", "--finetune", "--finetune_epochs", "1",
    ))
    assert finder == "vipsania"   # inferred from --model
    params = yaml.safe_load(path.read_text())
    assert params["vipsania"] == {"run": True, "model": "Fungi", "finetune": True, "finetune_epochs": 1}
    assert params["tiberius"] == {"run": False}
    assert params["outdir"] == str(tmp_path / "vipsania_results")


def test_build_params_merges_file_and_command_line(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump({
        "genome": str(DATA / "tiny.fa"), "threads": 8, "outdir": "out",
        "vipsania": {"run": True, "model": "Fungi", "batch_size": 4},
    }))
    finder, path = launcher.build_params(cli_args("-p", str(params_file), "--threads", "16", "--model", "Insecta"))
    assert finder == "vipsania"   # inferred from the block with run: true
    assert path == tmp_path / "out" / "params.yaml"
    params = yaml.safe_load(path.read_text())
    assert params["threads"] == 16
    assert params["vipsania"] == {"run": True, "model": "Insecta", "batch_size": 4}


def test_build_params_explicit_genefinder_and_result(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    finder, path = launcher.build_params(cli_args(
        "--genefinder", "tiberius", "--genome", str(DATA / "tiny.fa"), "--result", "old.gtf",
    ))
    assert finder == "tiberius"
    params = yaml.safe_load(path.read_text())
    assert params["tiberius"] == {"run": True, "result": str(tmp_path / "old.gtf")}   # no model_cfg needed


def test_genefinder_option_overrides_the_params_file(tmp_path: Path, monkeypatch):
    """--genefinder vipsania with a Tiberius params file runs Vipsania, and only Vipsania."""
    monkeypatch.chdir(tmp_path)
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump({
        "genome": str(DATA / "tiny.fa"), "genefinder": "tiberius",
        "tiberius": {"run": True, "model_cfg": str(DATA / "tiny.fa")},
    }))
    finder, path = launcher.build_params(cli_args("-p", str(params_file), "--genefinder", "vipsania", "--model", "Insecta"))
    assert finder == "vipsania"
    params = yaml.safe_load(path.read_text())
    assert params["genefinder"] == "vipsania"
    assert params["vipsania"] == {"run": True, "model": "Insecta"}
    assert params["tiberius"]["run"] is False


def test_genefinder_option_switches_a_disabled_block_on(tmp_path: Path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump({
        "genome": str(DATA / "tiny.fa"), "vipsania": {"run": False, "model": "Fungi"},
    }))
    _, path = launcher.build_params(cli_args("-p", str(params_file), "--genefinder", "vipsania"))
    assert yaml.safe_load(path.read_text())["vipsania"]["run"] is True
    assert "switches vipsania.run on" in capsys.readouterr().out


def test_two_enabled_blocks_without_a_choice_are_an_error(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump({
        "genome": str(DATA / "tiny.fa"),
        "tiberius": {"run": True, "model_cfg": str(DATA / "tiny.fa")}, "vipsania": {"run": True, "model": "Fungi"},
    }))
    with pytest.raises(SystemExit, match="More than one gene finder"):
        launcher.build_params(cli_args("-p", str(params_file)))


def test_tilde_and_variables_are_expanded_in_the_written_params(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("MYDATA", str(DATA))
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump({
        "genome": "$MYDATA/tiny.fa", "proteins": ["~/p.faa"], "tiberius": {"run": True, "result": "~/old.gtf"},
        "rnaseq_paired": [["~/a_1.fq", "~/a_2.fq"]],
    }))
    _, path = launcher.build_params(cli_args("-p", str(params_file)))
    params = yaml.safe_load(path.read_text())
    assert params["genome"] == str(DATA / "tiny.fa")
    assert params["proteins"] == [str(tmp_path / "p.faa")]
    assert params["tiberius"]["result"] == str(tmp_path / "old.gtf")
    assert params["rnaseq_paired"] == [[str(tmp_path / "a_1.fq"), str(tmp_path / "a_2.fq")]]


def test_zero_is_a_command_line_value(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _, path = launcher.build_params(cli_args(
        "--genefinder", "vipsania", "--genome", str(DATA / "tiny.fa"), "--model", "Fungi", "--threads", "0"))
    params = yaml.safe_load(path.read_text())
    assert params["threads"] == 0
    assert launcher.validate_input_data(params, path)   # and the validation rejects it


def test_three_paired_files_on_the_command_line_are_rejected(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit, match="rnaseq_paired"):
        launcher.build_params(cli_args("--genome", str(DATA / "tiny.fa"), "--model", "Fungi",
                                       "--rnaseq_paired", "a_1.fq", "a_2.fq", "b_1.fq"))


def test_mode_has_choices():
    with pytest.raises(SystemExit):
        cli_args("--genome", "g.fa", "--mode", "rnaseqq")


def test_build_params_evidence_only(tmp_path: Path, monkeypatch):
    """A params file can switch the gene finder off; then no model is required."""
    monkeypatch.chdir(tmp_path)
    params_file = tmp_path / "params.yaml"
    params_file.write_text(yaml.safe_dump({"genome": str(DATA / "tiny.fa"), "tiberius": {"run": False}}))
    finder, path = launcher.build_params(cli_args("-p", str(params_file)))
    assert finder == "tiberius"
    assert yaml.safe_load(path.read_text())["tiberius"] == {"run": False}


@pytest.mark.parametrize("argv,message", [
    (["--model_cfg", "diatoms"], "genome is required"),
    (["--genome", "g.fa"], "Tiberius needs a model configuration"),
    (["--genome", "g.fa", "--genefinder", "vipsania"], "Vipsania needs a model"),
    (["-p", "missing.yaml"], "Params YAML not found"),
])
def test_build_params_reports_missing_values(argv, message, tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit, match=message):
        launcher.build_params(cli_args(*argv))


@needs_tiberius
def test_cli_dry_run_writes_params_and_validates(tmp_path: Path, fake_path: Path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    cli.main([
        "--nf_config", "local", "--genome", str(DATA / "tiny.fa"), "--model_cfg", "diatoms",
        "--dry_run", "--skip_singularity_check",
    ])
    out = capsys.readouterr().out
    assert "Gene finder: tiberius" in out and "Dry run requested" in out
    assert (tmp_path / "tiberius_results" / "params.yaml").is_file()


def test_cli_forwards_nextflow_args(tmp_path: Path, fake_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calls = []
    real_run = launcher.subprocess.run

    def spy(cmd, *a, **kw):
        if Path(str(cmd[0])).name == "nextflow":
            calls.append((list(cmd), kw.get("env")))
            return SimpleNamespace(returncode=0)
        return real_run(cmd, *a, **kw)

    monkeypatch.setattr(launcher.subprocess, "run", spy)
    cli.main(["--genefinder", "vipsania", "--genome", str(DATA / "tiny.fa"), "--model", "Fungi",
              "--skip_singularity_check", "--", "-stub-run"])
    cmd, env = calls[0]
    assert cmd[-1] == "-stub-run"
    assert cmd[cmd.index("-params-file") + 1] == str(tmp_path / "vipsania_results" / "params.yaml")
    assert env is not None and "PATH" in env


# ---------------------------------------------------------------- version

def test_version_is_the_same_everywhere():
    """paludamentum/__init__.py is the single source; the other places must follow it."""
    import re
    version = launcher.__name__ and __import__("paludamentum").__version__
    manifest = re.search(r"version\s*=\s*'([^']+)'", (ROOT / "nextflow.config").read_text()).group(1)
    assert manifest == version
    assert f"**Status (v{version}).**" in (ROOT / "README.md").read_text()
    assert f"\nversion: {version}\n" in (ROOT / "CITATION.cff").read_text()
