
"""
unit tests for stage_registry

covers helper functions, stage validation, parameter validation,
resource resolution, shell quoting, invalid yaml handling,
unknown stages, unimplemented stages, success paths.
run with:

pytest scripts/test_stage_registry.py
"""

# import packages and modules
import shlex
import sys

import yaml

# dotted yaml path mapped to resulting bash variable name in dict, shared across stages
PARAM_VAR_MAP = {
    "sample.accession": "sample",
    "sample.reference_name": "reference_name",
    "sample.reference_url": "reference_url",
    "trimming.adapter_type": "adapter_type",
    "trimming.adapter_url": "adapter_url",
    "decontamination.spike_ins": "spike_ins",
    "threads": "threads",
}

# manual override for raw read urls/checksums from ena api
OPTIONAL_PARAM_VAR_MAP = {
    "sample.fastq_url_r1": "fastq_url_r1",
    "sample.fastq_url_r2": "fastq_url_r2",
    "sample.fastq_md5_r1": "fastq_md5_r1",
    "sample.fastq_md5_r2": "fastq_md5_r2",
}

QC_PARAMS = [
    "sample.accession",
    "trimming.adapter_type",
    "trimming.adapter_url",
    "decontamination.spike_ins",
    "threads",
]
ALIGN_PARAMS = QC_PARAMS + ["sample.reference_name", "sample.reference_url"]

# stage name -> required params, required resource categories & snakemake target
STAGE_REGISTRY = {
    "qc": {
        "target": "qc",
        "params": QC_PARAMS,
        "resources": {"reads": True, "adapters": True, "reference": False, "bwa_index": False},
    },
    "align": {
        "target": "align",
        "params": ALIGN_PARAMS,
        "resources": {"reads": True, "adapters": True, "reference": True, "bwa_index": True},
    },
    "all": {
        "target": "all",
        "params": ALIGN_PARAMS,
        "resources": {"reads": True, "adapters": True, "reference": True, "bwa_index": True},
    },
    # target none marks stage planned but not built into snakefile
    "variant-calling": {"target": None, "params": [], "resources": {}},
    "assembly": {"target": None, "params": [], "resources": {}},
    "annotation": {"target": None, "params": [], "resources": {}},
    "coverage": {"target": None, "params": [], "resources": {}},
    "amr": {"target": None, "params": [], "resources": {}},
    "taxonomy": {"target": None, "params": [], "resources": {}},
}


def get_nested(config, dotted_key):
    """
    Resolves dotted yaml path to nested value. Shared by required and optional
    parameter validation
    :param config: (dict) parsed yaml configuration structure
    :param dotted_key: (str) dot-delimited yaml path to resolve
    :return: returns resolved value when full path exists. returns none when any path
             component missing
    """
    
    # start traversal from root config dictionary
    value = config
    # walk nested yaml structure using each component of dotted path
    for part in dotted_key.split("."):
        if isinstance(value, dict) and part in value:
            value = value[part]
        else:
            return None
    return value

def resolve_stage(config_path, stage):
    """
    Validates stage name & config params required for chosen stage to execute
    :param config_path: (str) path to yaml config file
    :param stage: (str) name of stage to execute (ie qc, align, all)
    :returns: (list) status/result lines
                     success produces bash-ready key=value lines (config vars,
                     need_<resource> flags, snakemake_target). Failure produces
                     single prefixed status line: unknown_stage, unimplemented,
                     missing or invalid, each carrying enough detail for caller
                     to report
    """
    
    # notice for execution stage not recognized
    if stage not in STAGE_REGISTRY:
        valid_stages = ",".join(sorted(STAGE_REGISTRY))
        return [f"UNKNOWN_STAGE:{stage}:{valid_stages}"]
    
    # load stage definition containing target, required params, resource flags
    stage_spec = STAGE_REGISTRY[stage]
    # notice for specified stage not built
    if stage_spec["target"] is None:
        return [f"UNIMPLEMENTED:{stage}"]
    
    # read and parse yaml config for stage validation
    try:
        with open(config_path) as handle:
            config = yaml.safe_load(handle) or {}
    # notice for unreadable file or invalid yaml syntax
    except Exception as exc:
        return [f"INVALID:{exc}"]
    
    # collect missing required params and resolved bash variable values
    missing = []
    resolved = {}
    
    # helper fxn to validate only params required for selected stage
    for dotted_key in stage_spec["params"]:
        value = get_nested(config, dotted_key)
        
        if value in (None, ""):
            missing.append(dotted_key)
        # store resolved value using bash variable name mapping
        else:
            resolved[PARAM_VAR_MAP[dotted_key]] = value
    
    if missing:
        return ["MISSING:" + ",".join(missing)]
    
    # collect optional read download overrides
    optional_resolved = {}
    # resolve optional parameters without triggering validation failures
    for dotted_key, var_name in OPTIONAL_PARAM_VAR_MAP.items():
        value = get_nested(config, dotted_key)
        # missing optional fields remain empty rather than failing validation
        optional_resolved[var_name] = (
            value if value not in (None, "") else ""
        )

    # assemble bash-ready output
    lines = [f"{var_name}={shlex.quote(str(value))}" for var_name, value in resolved.items()]
    lines += [f"{var_name}={shlex.quote(str(value))}" for var_name, value in optional_resolved.items()]
    lines += [f"need_{name}={1 if needed else 0}" for name, needed in stage_spec["resources"].items()]
    lines.append(f"snakemake_target={shlex.quote(stage_spec['target'])}")
    return lines

# program entry point for command line execution
def main():
    # receive config path and stage from shell script
    config_path, stage = sys.argv[1], sys.argv[2]
    # emit validation results for bash consumption
    for line in resolve_stage(config_path, stage):
        print(line)

# execute only when invoked directly, not when imported for testing
if __name__ == "__main__":
    main()
