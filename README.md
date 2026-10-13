# NGS Workflow: Automated, Reproducible Paired-End Quality Control, Alignment, & Variant Calling

author: Stefanie Moreno
last updated: October 11, 2026

## Overview:

NGS Workflow is a stage-driven next-generation sequencing (NGS) pipeline for paired-end
Illumina runs. It automates read retrieval, preprocessing, quality control, contaminant
removal, alignment, duplicate removal, variant calling, and reporting. Every stage can be
run on its own or as a full workflow with a single command.


## Programming Components:

* Bash: Orchestration, argument parsing, resource download, checksum verification, stage execution

* Python: Validation, configuration validation and stage resolution through a central stage registry

* Snakemake: Workflow engine, dependency tracking, rule execution, dry-run validation

* Conda / Micromamba: Environment management, pinned tool versions

* Pytest: Testing, unit tests covering every branch of stage resolution


## The Pipeline:

* Auto-queries ENA for raw datafiles and md5 checksums based on user-apecified accession in the configuration file 

  -----OR-----

  *Note:* If after 3 automated attempts, raw data files are not successfully retrieved from ENA portal, the user can elect to uncomment
  the manual override in the configuration file and enter the specified URLs, which skips the ENA portal api lookup when both urls
  present. This design choice was explicitly implemented to safeguard against possible down ENA api or a broken accession metadata
  record. MD5 checksums are also verified when the URLs are directly provided in manual override.  
  
* Runs baseline and post-processing FastQC for before/after read quality comparison
  
* Trims adapters and low-quality bases with Trimmomatic

* Removes spike-in contaminants with BBMac.BBDuk
  
* Aligns processed reads to the corresponding indexed reference genome with BWA-MEM

* Removes PCR duplicates with samtools

* Conducts overview and comprehensive alignment statistics with flagstat and samtools stats
  
* Calls haploid variants with bcftools, applying parameterized hard filters on quality and depth
  
* Produces one MultiQC dashboard per stage or a unified dashboard for the full workflow
  
* Skips any download or index that already exists, and reports total runtime on completion
  
* Writes log and benchmark files for every rule
  

## Supported Platforms: 

MacOS and Linux. Windows users can run it through WSL. Docker support is planned as an upcoming implementation


## Execution Parameters:

All behavior is driven by `config/config.yaml`. The parameterized arguments are listed below:

* sample.accession -> Accession identifying the paired-end NGS run to be analyzed

**MANUAL OVERRIDE PARAMETERS- COMMENT OUT IF USING AUTO-ENA API RETRIEVAL BY ACCESSION**

* sample.fastq_url_r1 -> Manual override link for forward reads file

* sample.fastq_url_r2 -> Manual override link for reverse reads file

* sample.fastq_md5_r1 -> Manual override link for forward reads file md5 checksum

* sample.fastq_md5_r2 -> Manual override link for reverse reads file md5 checksum
  *Note:* Link for md5 checksums can be left as "" to omit during manual override

* reference.name -> Naming convention for specified reference genome (somewhat arbitrary)

* reference.url -> Link to reference genome file

* trimming.adapter_type -> Targeted adapter name according to specifications of NGS run

* trimming.adapter_url -> Link for targeted adapter file

* decontamination.spike_ins -> Name of targeted spike-in contaminant (BBDuk built-in reference keyword)
  *Note:* Other BBDuk keywords such as `lambda`, `adapters` and `artifacts` can be used depending on decontamination needs

* variant_calling.ploidy -> The number of complete sets of chromosomes in the cell that was sequenced

* variant_calling.min_base_quality ->






## Parameter Specifications and Associated Details for Example Execution

**Sequencing Run**

* Accession: `SRR25083113`
  
* Organism: *Escherichia coli* (ENA study: *E. coli* O157:H7 str. EDL933 genome sequencing and assembly)

* Platform: Illumina HiSeq 2000
  
* Library Layout: Paired-end
  
* Source: ENA FASTQ URLs set in `sample.fastq_url_r1` and `sample.fastq_url_r2` and commented out, most recent execution & push reflect
  ENA auto-query based on accession
  
**Reference Genome**

* Strain: *E. coli* `O157:H7 str. EDL933`

* RefSeq Assembly: `GCF_000006665.1` (ASM666v1)
  
* Sequences: Chromosome `NC_002655.2` and Plasmid pO157 `NC_007414.1`
  
* Reference URL: NCBI RefSeq genomic FASTA, gzip compressed

*Note:* The reference is matched to the sample's strain. Aligning this run to the *E. coli* K-12 reference instead gave a 75% mapped-read rate and tens of thousands of variant calls, while the EDL933 reference maps 97% of reads. Changing the strain means editing only the
two reference arguments.


**Threads**

* Threads: 10

The value is passed to Snakemake as `--cores` and to every multithreaded tool (BWA, samtools, bcftools), and is set according to the capacity specifications of my operating system. Customize to the number of cores on your machine.


**Adapters**

* Trimming Adapter Type: TruSeq3-PE
* Trimming Adapter Url: Trimmomatic adapter FASTA, downloaded to `resources/adapters/`


**Spike-ins**

* Decontamination Spike-Ins: `phix` (BBDuk built-in reference keyword)
  *Note:* Other BBDuk keywords such as `lambda`, `adapters` and `artifacts` can be used depending on decontamination needs 


**Variant Calling**

* Variant Calling Ploidy: Haploid calling, bacteria have one copy of the chromosome
  
* Variant Calling Minimum Base Quality: Minimum base quality counted in the pileup
  
* Variant Calling Minimum Map Quality: Minimum mapping quality for a read to be used
  
* Variant Calling Minimum Variant Quality: Hard filter on the QUAL field
  
* Variant Calling Minimum Depth: Hard filter on read depth at a variant site


**Fixed Tool Settings**

These are set in the Snakefile and are not yet configurable.

* Trimmomatic: `ILLUMINACLIP:<adapters>:2:30:10`, `LEADING:3`, `TRAILING:3`, `SLIDINGWINDOW:4:15`, `MINLEN:36`
  
* BBDuk: `k=31`, `hdist=1`
  
* BWA-MEM: `-M`, read group with `PL:illumina`
  
* Duplicate Removal: `samtools fixmate -m`, then `samtools markdup -r`
  

### Configuration Used for Example Execution

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


## Stages of Execution

The launcher is organized into stages and not individual Snakemake rules.

* Quality Control: Baseline FastQC, trimming, decontamination, post-QC FastQC (reference required: no)
  
* Alignment to Reference Genome: QC + BWA alignment, duplicate removal, BAM index, flagstat, alignment stats (reference required: yes)
  
* Variant Calling: align + variant calling, filtering, VCF indexing, variant statistics (reference required: yes)
  
* All: Every implemented stage plus the unified MultiQC dashboard (reference required: yes)

* Future Stage Implementations: Assembly, Annotation, Coverage Analysis, Circos Visualization, Anti-Microbial Resistance (AMR) Profiling, Taxonomy, Gene Expression Profiling, Comparative Genoics.


**Stage Diagram**

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

Each execution run validates the configuration, downloads only what is missing, performs a Snakemake
dry run to verify file setup and dependencies present, then initiates the workflow.

## Environment Setup

```bash
conda env create -f envs/env.yaml
conda activate bioinfo
```

The launcher checks that the `bioinfo` environment is active before running.

## Outputs

* `results/fastqc/`: Baseline and post-processing FastQC reports
  
* `results/alignment/`: Deduplicated BAM, BAM index, duplicate metrics, flagstat
  
* `results/variants/`: Raw and filtered VCFs with indexes, variant statistics
  
* `results/deliverables/`: Stage dashboards and statistics

* `logs/`: One log per rule
  
* `benchmarks/`: One runtime record per rule. *Note:* On macOS, Snakemake cannot collect CPU and memory statistics, so some benchmark
   files contain runtimes only from example execution.


**Deliverables By Stage:**

* Quality Control -> qc_dashboard.html

* Alignment -> alignment_dashboard.html, <sample>_stats.txt

* Variant Calling -> variant_dashboard.html

* All -> multiqc_dashboard.html, <sample>_stats.txt


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

* Separation of Concerns: The Bash launcher handles environment checks, downloads and checksums. The Python stage registry validates the
  stage and configuration and returns shell-safe variables. Snakemake handles the bioinformatics workflow. Each part can be tested and
  changed on its own.

* Validation Before Execution: Unknown stages, unimplemented stages, missing parameters and unreadable YAML are caught by the registry before Snakemake starts, with a clear error message and no stack trace.

* Reproducibility: Inputs, reference, adapters, contaminants and thresholds all live in one configuration file. The reference genome is a tracked Snakemake input, so changing it re-runs alignment and every downstream step instead of reusing stale results.

* Extensibility: A new stage needs one entry in the stage registry and its Snakemake rule(s), along with configuration for any parameters


## Testing

Execution Call:

```bash
python -m pytest -v
```
The tests cover stage validation, parameter validation, unknown and unimplemented stages, invalid YAML, resource resolution, the variant-calling stage, and shell-safe quoting.


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
