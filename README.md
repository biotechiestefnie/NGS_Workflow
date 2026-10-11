# NGS Workflow: Automated, Reproducible Paired-End QC, Alignment & Variant Calling

**Author:** Stefanie Moreno
**Last updated:** October 11, 2026

## Overview

NGS Workflow is a stage-driven next-generation sequencing (NGS) pipeline for paired-end
Illumina runs. It automates read retrieval, preprocessing, quality control, contaminant
removal, alignment, duplicate removal, variant calling, and reporting. Every stage can be
run on its own or as a full workflow with a single command.

* **Bash:** orchestration, argument parsing, resource download, checksum verification, stage execution
* **Python:** validation, configuration validation and stage resolution through a central stage registry
* **Snakemake:** workflow engine, dependency tracking, rule execution, dry-run validation
* **Conda / micromamba:** environment management, pinned tool versions
* **pytest:** testing, unit tests covering every branch of stage resolution

The pipeline:

* Downloads paired-end reads from user-specified URLs and verifies MD5 checksums when provided
* Runs baseline and post-processing FastQC for before/after comparison
* Trims adapters and low-quality bases (Trimmomatic) and removes spike-in contaminants (BBDuk)
* Aligns cleaned reads to a reference genome (BWA-MEM) and removes PCR duplicates (samtools)
* Calls haploid variants (bcftools) and applies hard filters on quality and depth
* Produces one MultiQC dashboard per stage plus a unified dashboard for the full workflow
* Skips any download or index that already exists, and reports total runtime on completion
* Writes a log file and a benchmark file for every rule

**Supported platforms:** macOS and Linux. Windows users can run it through WSL. Docker
support is planned.

## Execution Parameters

All behavior is driven by `config/config.yaml`. The values below are the example run
included in this repository.

### Sequencing run

* **Accession:** `SRR25083113`
* **Organism:** *Escherichia coli* (ENA study: *E. coli* O157:H7 str. EDL933 genome sequencing and assembly)
* **Platform:** Illumina HiSeq 2000
* **Library layout:** Paired-end
* **Source:** ENA FASTQ URLs set in `sample.fastq_url_r1` and `sample.fastq_url_r2`

### Reference genome

* **Strain:** *E. coli* O157:H7 str. EDL933
* **RefSeq assembly:** `GCF_000006665.1` (ASM666v1)
* **Sequences:** chromosome `NC_002655.2` and plasmid pO157 `NC_007414.1`
* **`reference.name`:** `ecoli_o157h7_edl933` (label used for files in `resources/reference/`)
* **`reference.url`:** NCBI RefSeq genomic FASTA, gzip compressed

The reference is matched to the sample's strain. Aligning this run to the *E. coli* K-12
reference instead gave a 75% mapped-read rate and tens of thousands of variant calls,
while the EDL933 reference maps 97% of reads. Changing the strain means editing only the
two `reference` values.

### Threads

* **`threads`:** `10`

The value is passed to Snakemake as `--cores` and to every multithreaded tool (BWA,
samtools, bcftools). Set it to the number of cores on your machine.

### Adapters

* **`trimming.adapter_type`:** `TruSeq3-PE`
* **`trimming.adapter_url`:** Trimmomatic adapter FASTA, downloaded to `resources/adapters/`

### Spike-ins

* **`decontamination.spike_ins`:** `phix` (BBDuk built-in reference keyword)

Other BBDuk keywords such as `lambda`, `adapters` and `artifacts` can be used instead.

### Variant calling

* **`variant_calling.ploidy`** = `1`: haploid calling, bacteria have one copy of the chromosome
* **`variant_calling.min_base_quality`** = `20`: minimum base quality counted in the pileup
* **`variant_calling.min_map_quality`** = `20`: minimum mapping quality for a read to be used
* **`variant_calling.min_variant_quality`** = `30`: hard filter on the QUAL field
* **`variant_calling.min_depth`** = `10`: hard filter on read depth at a variant site

### Fixed tool settings

These are set in the Snakefile and are not yet configurable.

* **Trimmomatic:** `ILLUMINACLIP:<adapters>:2:30:10`, `LEADING:3`, `TRAILING:3`, `SLIDINGWINDOW:4:15`, `MINLEN:36`
* **BBDuk:** `k=31`, `hdist=1`
* **BWA-MEM:** `-M`, read group with `PL:illumina`
* **Duplicate removal:** `samtools fixmate -m`, then `samtools markdup -r`

### Example configuration

```yaml
sample:
  accession: "SRR25083113"
  # manual read URLs, md5 checksums are optional and verified when provided
  fastq_url_r1: "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR250/013/SRR25083113/SRR25083113_1.fastq.gz"
  fastq_url_r2: "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR250/013/SRR25083113/SRR25083113_2.fastq.gz"
  fastq_md5_r1: ""
  fastq_md5_r2: ""

reference:
  name: "ecoli_o157h7_edl933"
  url: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/665/GCF_000006665.1_ASM666v1/GCF_000006665.1_ASM666v1_genomic.fna.gz"

trimming:
  adapter_type: "TruSeq3-PE"
  adapter_url: "https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa"

decontamination:
  spike_ins: "phix"

variant_calling:
  ploidy: 1
  min_base_quality: 20
  min_map_quality: 20
  min_variant_quality: 30
  min_depth: 10

threads: 10
```

## Stages

The launcher exposes stages and not individual Snakemake rules.

* **`qc`**: baseline FastQC, trimming, decontamination, post-QC FastQC (reference required: no)
* **`align`**: `qc` + BWA alignment, duplicate removal, BAM index, flagstat, alignment statistics (reference required: yes)
* **`variant-calling`**: `align` + variant calling, filtering, VCF indexing, variant statistics (reference required: yes)
* **`all`**: every implemented stage plus the unified MultiQC dashboard (reference required: yes)

Planned stages: `assembly`, `annotation`, `coverage`, `amr`, `taxonomy`.

```text
raw FASTQ -> baseline_fastqc -> trim -> decon -> processed_fastqc
                                          |
                                          v
                  bwa_align -> dedup -> index_bam -> flagstat / stats
                                          |
                                          v
        faidx_reference -> variant_call -> index_vcf -> filter_variants
                                                              |
                                                              v
                                                variant_stats -> dashboards
```

## Running the Pipeline

```bash
# quality control only
bash scripts/run_pipeline.sh --stage qc

# quality control + alignment
bash scripts/run_pipeline.sh --stage align

# quality control + alignment + variant calling
bash scripts/run_pipeline.sh --stage variant-calling

# full workflow (default when no stage is given)
bash scripts/run_pipeline.sh
bash scripts/run_pipeline.sh --stage all

# list stages
bash scripts/run_pipeline.sh --help
```

Each run validates the configuration, downloads only what is missing, performs a Snakemake
dry run to verify file setup and dependencies present, then executes the workflow.

## Environment Setup

```bash
conda env create -f envs/env.yaml
conda activate bioinfo
```

The launcher checks that the `bioinfo` environment is active before running.

## Outputs

* **`results/fastqc/`:** baseline and post-processing FastQC reports
* **`results/alignment/`:** deduplicated BAM, BAM index, duplicate metrics, flagstat
* **`results/variants/`:** raw and filtered VCFs with indexes, variant statistics
* **`results/deliverables/`:** stage dashboards and statistics
* **`logs/`:** one log per rule
* **`benchmarks/`:** one runtime record per rule (see note below)

Deliverables by stage:

```text
qc -> qc_dashboard.html
align -> alignment_dashboard.html, <sample>_stats.txt
variant-calling -> variant_dashboard.html
all -> multiqc_dashboard.html, <sample>_stats.txt
```

On macOS, Snakemake cannot collect CPU and memory statistics, so some benchmark files
contain runtimes only.

## Directory Structure

```text
NGS_Workflow/
├── config/
│   └── config.yaml
├── envs/
│   └── env.yaml
├── scripts/
│   ├── run_pipeline.sh
│   └── stage_registry.py
├── tests/
│   └── test_stage_registry.py
├── raw/
├── resources/
│   ├── reference/
│   └── adapters/
├── results/
│   ├── fastqc/
│   ├── trimmed/
│   ├── processed/
│   ├── alignment/
│   ├── variants/
│   └── deliverables/
├── logs/
├── benchmarks/
├── pyproject.toml
└── Snakefile
```

## Design

**Separation of concerns.** The Bash launcher handles environment checks, downloads and
checksums. The Python stage registry validates the stage and configuration and returns
shell-safe variables. Snakemake handles the bioinformatics workflow. Each part can be tested 
and changed on its own.

**Validation before execution.** Unknown stages, unimplemented stages, missing parameters
and unreadable YAML are caught by the registry before Snakemake starts, with a clear error
message and no stack trace.

**Reproducibility.** Inputs, reference, adapters, contaminants and thresholds all live in
one configuration file. The reference genome is a tracked Snakemake input, so changing it
re-runs alignment and every downstream step instead of reusing stale results.

**Extensibility.** A new stage needs one entry in the stage registry and its Snakemake
rules.

## Testing

```bash
python -m pytest -v
```

The tests cover stage validation, parameter validation, unknown and unimplemented stages,
invalid YAML, resource resolution, the variant-calling stage, and shell-safe quoting.

## Troubleshooting

* **Environment error:** `conda activate bioinfo`
* **`MISSING:` error:** add the named parameter to `config/config.yaml`
* **`UNKNOWN_STAGE` error:** run the script with `--help` for valid stages
* **Download keeps failing:** check the URL and your connection, downloads retry three times
* **Checksum mismatch:** delete `raw/*.fastq.gz` and rerun
* **Workflow problems:** run `snakemake -n <target>` to inspect the dependency graph

## Roadmap

* Docker image for fully portable execution
* Assembly (SPAdes), annotation (Prokka or Bakta), coverage analysis (mosdepth)
* AMR detection (AMRFinderPlus), taxonomic classification (Kraken2)
* Minimum allele-fraction filter for variant calls
* Per-stage tool requirement checks in place of the environment-name check
