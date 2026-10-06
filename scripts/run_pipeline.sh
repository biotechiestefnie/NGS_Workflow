#!/bin/bash

# ngs workflow launcher and resource manager
#
# responsibilities: validate execution environment, validate config through stage_registry,
# retrieve sequencing reads and reference resources, verify download integrity through md5
# checksums, provision stage-specific dependencies, and execute snakemake workflows through
# stage targets
#
# supported stages: -qc, -align, -all
#
# future expansion stages: variant calling, assembly, annotation, coverage, amr, and taxonomy

# stop execution on failure or unspecified variables
# flags: -e exits if command nonzero exit, -u unset vars as errors, -q quit if command fails
set -euo pipefail

# resolve working directory for execution call
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
cd "${PROJECT_ROOT}"

# configuration file location
CONFIG_PATH="config/config.yaml"
# stage registry file location
STAGE_REGISTRY_SCRIPT="scripts/stage_registry.py"

# default stage runs full pipeline when no --stage flag input
stage="all"

# parse cli execution args: supports --stage <name>, --stage=<name>, -h/--help
while [ $# -gt 0 ]; do
    case "$1" in
        --stage)
            stage="$2"
            shift 2
            ;;
        --stage=*)
            stage="${1#--stage=}"
            shift
            ;;
        -h|--help)
            echo "Usage: ./run_pipeline.sh [--stage <stage_name>]"
            echo "Implemented Stages: qc, align, all"
            echo "Planned Stages: variant-calling, assembly, annotation, coverage, amr, taxonomy"
            exit 0
            ;;
        *)
            echo "Error: Unrecognized Argument '$1'."
            echo "Usage: ./run_pipeline.sh [--stage <stage_name>]"
            exit 1
            ;;
    esac
done

# lightweight check for environment activation/dependency availability
if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "${CONDA_DEFAULT_ENV}" != "bioinfo" ]; then
    echo "Error: Please Activate the Required Conda Environment 'bioinfo' Before Executing the Program."
    echo "Terminal Command: conda activate bioinfo"
    exit 1
fi

# keeps checksum verification portable across operating systems
compute_md5() {
    local file="$1"

    if command -v md5sum >/dev/null 2>&1; then  # linux check
        md5sum "${file}" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then  # macos check
        md5 -q "${file}"
    else
        echo "Error: No Md5 Utility Found on This System (Expected 'md5sum' or 'md5')." >&2
        return 1
    fi
}

# retry wrapper shared by network fetches: 3 attempts, centralize retry so downloads share
# identical robustness rather than duplicating ad hoc error
# handling per call site
download_with_retry() {
    local url="$1"
    local dest="$2"
    local attempt=1
    local max_attempts=3

    while [ "${attempt}" -le "${max_attempts}" ]; do
        if wget -q -O "${dest}" "${url}"; then
            return 0
        fi
        echo "[warn] Download Attempt ${attempt} of ${max_attempts} Failed for '${url}'. Retrying..." >&2
        attempt=$((attempt + 1))
        sleep 5
    done

    echo "Error: Unable to Download '${url}' After ${max_attempts} Attempts." >&2
    return 1
}

# retrieve read download locations and checksum metadata from ena. retries api failures
# instead of terminating workflow d/t temp service issues
fetch_ena_metadata() {
    local accession="$1"
    local attempt=1
    local max_attempts=3
    local metadata
    local url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"

    while [ "${attempt}" -le "${max_attempts}" ]; do
        if metadata=$(curl -fsSL "${url}" 2>/dev/null | tail -n 1); then
            echo "${metadata}"
            return 0
        fi
        echo "[warn] ENA Metadata Query Attempt ${attempt} of ${max_attempts} Failed for Accession '${accession}'. Retrying..." >&2
        attempt=$((attempt + 1))
        sleep 5
    done

    echo "Error: Unable to Reach ENA Portal API for Accession '${accession}' After ${max_attempts} Attempts." >&2
    return 1
}

# validate download integrity against md5 checksum. prevents corrupted or incomplete
# files from entering workflow
verify_checksum() {
    local file="$1"
    local expected_md5="$2"
    local actual_md5

    actual_md5=$(compute_md5 "${file}") || return 1
    if [ "${actual_md5}" != "${expected_md5}" ]; then
        echo "Error: Checksum Mismatch for '${file}'."
        echo "Expected: ${expected_md5}"
        echo "Got: ${actual_md5}"
        echo "File May Be Corrupted or Incomplete. Delete File and Re-Run Program."
        return 1
    fi
}

# resolve stage requirements, resource flags and snakemake target through centralized
# registry validation layer
config_output=$(python3 "${STAGE_REGISTRY_SCRIPT}" "${CONFIG_PATH}" "${stage}")

# handle validation failures returned by stage registry
case "${config_output}" in
    UNKNOWN_STAGE:*)
        rest="${config_output#UNKNOWN_STAGE:}"
        bad_stage="${rest%%:*}"
        valid_stages="${rest#*:}"
        echo "Error: Unrecognized Stage '${bad_stage}'."
        echo "Valid Stages: ${valid_stages//,/, }"
        exit 1
        ;;
    UNIMPLEMENTED:*)
        bad_stage="${config_output#UNIMPLEMENTED:}"
        echo "Error: Stage '${bad_stage}' Is Not Yet Implemented in Snakefile."
        echo "Currently Available Stages: qc, align, all."
        exit 1
        ;;
    MISSING:*)
        missing_params="${config_output#MISSING:}"
        echo "Error: Config File Is Missing Required Parameter(s) for Stage '${stage}': ${missing_params//,/, }"
        echo "Return to Config File and Implement Your Desired Specifications for This Run."
        exit 1
        ;;
    INVALID:*)
        echo "Error: Config File Could Not Be Parsed - ${config_output#INVALID:}"
        echo "Check Yaml Syntax in Config File and Re-Run Program."
        exit 1
        ;;
esac

# declare variables populated by registry output. prevents undefined variable warnings
declare sample reference_name reference_url adapter_type adapter_url \
    threads need_reads need_adapters need_reference need_bwa_index snakemake_target \
    fastq_url_r1 fastq_url_r2 fastq_md5_r1 fastq_md5_r2

# load params from config & state registry files
eval "${config_output}"

# start runtime counter
start_time=$SECONDS

# create project subdirectories if not present
mkdir -p "raw"
mkdir -p "resources/reference"
mkdir -p "resources/adapters"
mkdir -p "logs"
mkdir -p "benchmarks"

echo "[info] Selected Stage: ${stage}"

# provision paired-end reads required by selected stage. reuse local files when present,
# otherwise resolve sources from config or ena
if [ "${need_reads}" = "1" ]; then
    echo "[info] Target Dataset: ${sample}"

    if [ -f "raw/${sample}_1.fastq.gz" ] && [ -f "raw/${sample}_2.fastq.gz" ]; then
        echo "[info] Skipping: Compressed FASTQ Files Already Present"

    elif [ -f "raw/${sample}_1.fastq" ] && [ -f "raw/${sample}_2.fastq" ]; then
        echo "[info] Compressing Existing FASTQ Files..."
        gzip "raw/${sample}_1.fastq"
        gzip "raw/${sample}_2.fastq"

    elif [ -n "${fastq_url_r1}" ] && [ -n "${fastq_url_r2}" ]; then
        echo "[info] Using Manually Specified FASTQ URLs From Config File..."

        echo "[info] Downloading Forward Reads..."
        download_with_retry "${fastq_url_r1}" "raw/${sample}_1.fastq.gz" || exit 1
        if [ -n "${fastq_md5_r1}" ]; then
            verify_checksum "raw/${sample}_1.fastq.gz" "${fastq_md5_r1}" || exit 1
        else
            echo "[warn] No Checksum Provided for Forward Reads (sample.fastq_md5_r1). Skipping Verification."
        fi

        echo "[info] Downloading Reverse Reads..."
        download_with_retry "${fastq_url_r2}" "raw/${sample}_2.fastq.gz" || exit 1
        if [ -n "${fastq_md5_r2}" ]; then
            verify_checksum "raw/${sample}_2.fastq.gz" "${fastq_md5_r2}" || exit 1
        else
            echo "[warn] No Checksum Provided for Reverse Reads (sample.fastq_md5_r2). Skipping Verification."
        fi

    else
        echo "[info] Querying ENA Portal API for Read Metadata..."
        read_metadata=$(fetch_ena_metadata "${sample}") || exit 1

        if [ -z "${read_metadata}" ]; then
            echo "Error: ENA Portal API Returned No Metadata for Accession '${sample}'."
            echo "Confirm Accession Is Correct in Config File."
            exit 1
        fi

        fastq_urls=$(echo "${read_metadata}" | cut -f1)
        fastq_md5s=$(echo "${read_metadata}" | cut -f2)
        url_r1="https://$(echo "${fastq_urls}" | cut -d';' -f1)"
        url_r2="https://$(echo "${fastq_urls}" | cut -d';' -f2)"
        md5_r1=$(echo "${fastq_md5s}" | cut -d';' -f1)
        md5_r2=$(echo "${fastq_md5s}" | cut -d';' -f2)

        echo "[info] Downloading Forward Reads..."
        download_with_retry "${url_r1}" "raw/${sample}_1.fastq.gz" || exit 1
        verify_checksum "raw/${sample}_1.fastq.gz" "${md5_r1}" || exit 1

        echo "[info] Downloading Reverse Reads..."
        download_with_retry "${url_r2}" "raw/${sample}_2.fastq.gz" || exit 1
        verify_checksum "raw/${sample}_2.fastq.gz" "${md5_r2}" || exit 1
    fi
fi

# provision reference genome required by selected stage, skip if already present
if [ "${need_reference}" = "1" ]; then
    if [ -f "resources/reference/${reference_name}.fa" ]; then
        echo "[info] Skipping: Reference Genome Already Present"

    else
        if [ ! -f "resources/reference/${reference_name}.fa.gz" ]; then
            echo "[info] Downloading Reference Genome..."
            download_with_retry "${reference_url}" "resources/reference/${reference_name}.fa.gz" || exit 1
        fi

        echo "[info] Decompressing Reference Genome..."
        gunzip -k "resources/reference/${reference_name}.fa.gz"
    fi
fi

# download sequence adapters according to ngs run for trimming, only if qc executed & file
# not already present
if [ "${need_adapters}" = "1" ]; then
    if [ ! -f "resources/adapters/${adapter_type}.fa" ]; then
        echo "[info] Downloading ${adapter_type} Adapter Fasta..."
        download_with_retry "${adapter_url}" "resources/adapters/${adapter_type}.fa" || exit 1

    else
        echo "[info] Adapter Fasta File Already Present"
    fi
fi

# build bwa index for reference genome when absent, only when selected stage align or all
if [ "${need_bwa_index}" = "1" ]; then
    if [ ! -f "resources/reference/${reference_name}.fa.bwt" ]; then
        echo "[info] Building BWA Index..."
        bwa index "resources/reference/${reference_name}.fa" > "logs/bwa_index_${sample}.log" 2>&1 || {
            # error advising inspect log for building fail details
            echo "Error: BWA Index Build Failed. See logs/bwa_index_${sample}.log for Details."
            exit 1
        }

    else
        echo "[info] Skipping: BWA Index Already Present"
    fi
fi

# dry run to confirm workflow integrity for selected stage target
echo "[info] Running Validation Dry Run..."
snakemake -n --rerun-incomplete "${snakemake_target}"

# execute specified workflow through selected stage target
echo "[info] Executing Live Multi-Core Workflow..."
snakemake \
    --cores "${threads}" \
    --printshellcmds \
    --rerun-incomplete \
    "${snakemake_target}"

# calculate total workflow runtime sans file downloads, ref indexing
duration=$((SECONDS-start_time))
minutes=$((duration/60))
seconds=$((duration%60))

# completion status, total runtime, dashboard location to terminal
echo "Pipeline Completed Successfully for Stage '${stage}'!"
echo "Total Workflow Runtime: ${minutes}m ${seconds}s"

if [ "${stage}" = "all" ]; then
    echo "MultiQC Dashboard Located in 'deliverables/' Folder"
fi
