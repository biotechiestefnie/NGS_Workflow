"""
unit tests for stage_registry.resolve_stage

covers every branch resolve_stage can return: success, unknown stage,
unimplemented stage, missing param(s) & invalid yaml. run with:

pytest tests/test_stage_registry.py
"""

import textwrap

import pytest

from stage_registry import get_nested, resolve_stage


FULL_CONFIG = """
sample:
  accession: "test_accession"

reference:
  name: "test_reference"
  url: "https://example.com/test_reference.fa.gz"

trimming:
  adapter_type: "test_adapter"
  adapter_url: "https://example.com/test_adapter.fa"

decontamination:
  spike_ins: "test_contaminant"
  
variant_calling:
  ploidy: 1
  min_base_quality: 20
  min_map_quality: 20
  min_variant_quality: 30
  min_depth: 10

threads: 1
"""

NO_REFERENCE_CONFIG = """
sample:
  accession: "test_accession"

trimming:
  adapter_type: "test_adapter"
  adapter_url: "https://example.com/test_adapter.fa"

decontamination:
  spike_ins: "test_contaminant"

threads: 1
"""

# verify dotted yaml path resolves correct nested value
def test_get_nested_returns_value():
    config = {
        "sample": {
            "accession": "test_accession"
        }
    }

    assert get_nested(config, "sample.accession") == "test_accession"

# verify missing dotted path returns none for downstream validation
def test_get_nested_returns_none_for_missing_key():
    config = {
        "sample": {
            "accession": "test_accession"
        }
    }

    assert get_nested(config, "sample.reference_name") is None
    
    
# create temporary config file from yaml string
@pytest.fixture
def write_config(tmp_path):
    def _write(contents):
        config_path = tmp_path / "config.yaml"
        # remove python source indentation before writing yaml
        config_path.write_text(textwrap.dedent(contents))
        return str(config_path)
    # provide config path for tests
    return _write

# test for when reference genome not required for stage run
def test_qc_stage_does_not_require_reference_params(write_config):
    config_path = write_config(NO_REFERENCE_CONFIG)

    lines = resolve_stage(config_path, "qc")
    
    # verify qc stage skips reference resource requirements
    assert "need_reference=0" in lines
    assert "snakemake_target=qc" in lines
    assert not any(line.startswith("MISSING:") for line in lines)

# confirm alignment stage fails without reference genome
def test_align_stage_requires_reference_params(write_config):
    config_path = write_config(NO_REFERENCE_CONFIG)

    lines = resolve_stage(config_path, "align")

    assert len(lines) == 1
    assert lines[0].startswith("MISSING:")
    assert "reference.name" in lines[0]
    assert "reference.url" in lines[0]

# confirm that running entire pipeline requires all params
def test_all_stage_succeeds_with_full_config(write_config):
    config_path = write_config(FULL_CONFIG)

    lines = resolve_stage(config_path, "all")

    assert "need_reference=1" in lines
    assert "need_bwa_index=1" in lines
    assert "snakemake_target=all" in lines

# ensure execution command contains correct staging option
def test_unknown_stage_name(write_config):
    config_path = write_config(FULL_CONFIG)

    lines = resolve_stage(config_path, "bogus")

    assert lines[0].startswith("UNKNOWN_STAGE:bogus:")

# test for stages not yet created
def test_planned_but_unimplemented_stage(write_config):
    config_path = write_config(FULL_CONFIG)

    lines = resolve_stage(config_path, "amr")

    assert lines == ["UNIMPLEMENTED:amr"]

# verify malformed yaml param entry returns invalid status
def test_invalid_yaml_reported_cleanly(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text("sample:\n  accession: [unclosed")

    lines = resolve_stage(str(config_path), "qc")

    assert lines[0].startswith("INVALID:")

# verify shell-special characters receive safe bash quoting
def test_values_are_shell_quoted(write_config):
    config_path = write_config(
        FULL_CONFIG.replace(
            'adapter_url: "https://example.com/test_adapter.fa"',
            'adapter_url: "https://example.com/test adapter.fa"',
        )
    )

    lines = resolve_stage(config_path, "qc")

    adapter_line = next(
        line
        for line in lines
        if line.startswith("adapter_url=")
    )

    assert (
        adapter_line
        == "adapter_url='https://example.com/test adapter.fa'"
    )
    
# confirm variant-calling stage resolves to underscore snakemake target
def test_variant_calling_stage_succeeds_with_full_config(write_config):
    config_path = write_config(FULL_CONFIG)

    lines = resolve_stage(config_path, "variant-calling")

    assert "snakemake_target=variant_calling" in lines
    assert "ploidy=1" in lines
    assert "need_reference=1" in lines


# confirm variant-calling stage fails when variant params are absent
def test_variant_calling_stage_requires_variant_params(write_config):
    config_path = write_config(
        FULL_CONFIG.replace("ploidy: 1", "")
    )

    lines = resolve_stage(config_path, "variant-calling")

    assert len(lines) == 1
    assert lines[0].startswith("MISSING:")
    assert "variant_calling.ploidy" in lines[0]
    