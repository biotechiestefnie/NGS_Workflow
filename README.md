# NGS Workflow: Automated, Reproducible Paired-End Whole Genome Sequencing (WGS) Quality Control, Alignment, & Variant Calling

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
  

## Supported Operating Systems: 

MacOS and Linux. Windows users can run it through WSL. Docker support is planned as an upcoming addition


## Execution Parameters:

All parameterization behavior is driven by `config/config.yaml`. The arguments that the user can adjust are listed below:

* `{sample.accession}` -> Accession identifying the NGS paired-end run to be analyzed

**MANUAL OVERRIDE PARAMETERS- COMMENT OUT IF USING AUTO-ENA API RETRIEVAL BY ACCESSION**

* `{sample.fastq_url_r1}` -> Manual override link for forward reads file

* `{sample.fastq_url_r2}` -> Manual override link for reverse reads file

* `{sample.fastq_md5_r1} -> Manual override link for forward reads file md5 checksum

* `{sample.fastq_md5_r2}` -> Manual override link for reverse reads file md5 checksum
  *Note:* Link for md5 checksums can be left as "" to omit during manual override

* `{reference.name}` -> Naming convention for specified reference genome (somewhat arbitrary)

* `{reference.url} -> Link to reference genome file

* `{trimming.adapter_type}` -> Targeted adapter name according to specifications of NGS run
  

* ~{trimming.adapter_url}` -> Link for targeted adapter file

* `{decontamination.spike_ins} -> Name of targeted spike-in contaminant (BBDuk built-in reference keyword)
  *Note:* Other BBDuk keywords such as `lambda`, `adapters` and `artifacts` can be used depending on decontamination needs

* `{variant_calling.ploidy}` -> The number of complete sets of chromosomes in the cell that was sequenced

* `{variant_calling.min_base_quality}` -> Minimum base quality counted in the pileup

* `{variant_calling.min_map_quality} -> Minimum mapping quality for a read to be used
  
* `{variant_calling.min_variant_quality}`-> Hard filter on the QUAL field
  
* `{variant_calling.min_depth}`-> Hard filter on read depth at a variant site

* `{threads}`: The value is passed to Snakemake as `--cores` and to every multithreaded tool (BWA, samtools, bcftools), and is set according to the capabilities of my operating system (Macbook M2Pro). The multiple threads of a given process may be executed in parallel (via multithreading capabilities), sharing resources such as memory. In particular, the threads of a process share its executable code and the values of its dynamically allocated variables and non-thread-local global variables at any given time. This results in significant enhancements in execution speed for the multithreaded tools. Customize to the number of available cores you intend to allocate on your machine.


## Parameters Used as Inputs with Associated Details for Example Execution

**Sequencing Run**

* Accession: SRR25083113
  
* Organism: *Escherichia coli* (ENA study: *E. coli* O157:H7 str. EDL933 genome sequencing and assembly)

* Platform: Illumina HiSeq 2000
  
* Library Layout: Paired-end
  
* Source: Automated ENA API extraction
  *Note:* Automating extraction of the raw files by accession simplifies the data collection process for the user and minimizes the
  parameters required for input, as well as the potential for entry errors or inability to locate desired datafiles.
  
**Reference Genome**

* Strain: *Escherichia coli* O157:H7 str. EDL933

* RefSeq Assembly: GCF_000006665.1 (ASM666v1)
  
* Sequences: Chromosome NC_002655.2 and Plasmid pO157 NC_007414.1
  
* Reference URL: g-zipped fasta file from NCBI RefSeq genomic FASTA link:
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/665/GCF_000006665.1_ASM666v1/GCF_000006665.1_ASM666v1_genomic.fna.gz"
  *Note:* The reference is matched to the sample's strain. Aligning this run to the *E. coli* K-12 reference instead gave a 75% mapped
   read rate and tens of thousands of variant calls, while the EDL933 reference maps 97% of reads. Changing the strain means editing only
   the two reference arguments.

**Trimming**

* Trimming Adapter Type: TruSeq3-PE
  *Note:* TruSeq adapters are synthetic DNA sequences about 80 bases long that are attached to DNA or RNA fragments during Illumina NGS library preparation. They contain the mandatory P5 and P7 terminal sequences that anchor library molecules to the Illumina sequencing flow cell. TruSeq adapters also contain barcodes that incorporate 6 to 10bp index sequences (i7 and i5), allowing multiple samples to be pooled and sequenced together. Furthermore, they provide the required binding sites for sequencing and amplification primers, such as PCR Primer 1.0 and 2.0. Because E. coli has a relatively small genome (~4.6 Mb) and genomic runs usually yield high coverage with varying fragment sizes, adapter trimming is essential to prevent alignment failures. Using TruSeq3-PE or TruSeq3-PE-2 via Trimmomatic is the standard approach for trimming Illumina TruSeq adapters from paired-end data.
   
* Trimming Adapter Url: Trimmomatic adapter FASTA, downloaded to the resources folder from
  https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa

**Decontamination**

* Spike-Ins: phix

**Variant Calling**

* variant_calling.ploidy: 1
  *Note:* *E. coli* is a haploid bacterium with a single circular chromosome, and therefore, each position should carry one allele.
  Setting ploidy to 1 stops the caller from modeling heterozygous genotypes (like 0/1), which would be biologically meaningless here and
  would inflate false-positive calls from sequencing errors or mismapped reads.

* variant_calling.min_base_quality: 20 (Phred)
  *Note:* 20 Q20 means a 1% chance the base call is wrong. Filtering below this removes the noisiest bases (Illumina quality drops toward read ends) so they do not contribute to allele counts and trigger spurious SNPs.

* variant_calling.min_map_quality: 20
  *Note:* MAPQ 30 corresponds to roughly a 0.1% chance the read is misaligned. This discards ambiguous alignments, typically reads from
  repeats (IS elements, rRNA operons, REP sequences) in the genome, which are a major source of false variants in bacterial data. 

* variant_calling.min_variant_quality: 30
  *NOTE:* This is the Phred-scaled confidence that a variant exists at the site at all, so Q30 keeps calls with about 99.9% confidence.
  This is a common, conventional threshold that filters weak, marginal calls after the caller has pooled all the evidence.

* min_depth: 10
  *Note:* With a haploid genome, enough independent reads are needed to distinguish a real variant from random error and to trust that
  the consensus allele is correct. About 10x is a widely used floor for bacterial isolates, and it is often applied alongside an allele
  fraction cutoff (e.g., ≥ 80-90% of reads supporting the alt).

**Threads**

* Threads: 10 (Wit a high-end processor (like an 12-core CPU), assigning 10 threads means the OS can assign one thread to each individual
  core. All 10 threads will execute instructions at the exact same physical instant. This is called parallelism, and it dramatically
  speeds up the execution of the multi-threaded tools.)
  

**Fixed Tool Settings**

These are set in the Snakefile and are not yet configurable.

* Trimmomatic: ILLUMINACLIP:<adapters>:2:30:10, LEADING:3, TRAILING:3, SLIDINGWINDOW:4:15, MINLEN:36
  
* BBDuk: k=31, hdist=1
  
* BWA-MEM: -M, read group with PL:illumina
  
* Duplicate Removal: samtools fixmate -m, then samtools markdup -r
  

### Argument Imput Locations Within the Configuration File

```yaml

sample:
  accession: 

  # Manual Override: read URLs, md5 checksums are optional and verified when provided
  # fastq_url_r1: 
  # fastq_url_r2: 
  # fastq_md5_r1: ""
  # fastq_md5_r2: ""

reference:
  name: 
  url: 

trimming:
  adapter_type: 
  adapter_url:

decontamination:
  spike_ins: 

variant_calling:
  ploidy: 
  min_base_quality: 
  min_map_quality: 
  min_variant_quality: 
  min_depth: 

threads: 

```


## Stages of Execution

The launcher is organized into stages and not individual Snakemake rules.

* Quality Control: Baseline FastQC, trimming, decontamination, post-processing FastQC (reference required: no)
  
* Alignment to Reference Genome: Quality control + BWA alignment, duplicate removal, BAM indexing of aligned reads, flagstat overview
  alignment statistics, comprehensive alignment statistics (reference required: yes)
  
* Variant Calling: align plus variant calling, filtering, VCF indexing, variant statistics (reference required: yes)
  
* All: Every implemented stage plus the consolidated MultiQC dashboard (reference required: yes)

* Future Stage Implementations: Assembly, Annotation, Coverage Analysis, Circos Visualization, Anti-Microbial Resistance (AMR) Profiling, Taxonomy, Gene Expression Profiling, Comparative Genoics.


**Flow of Data Through Pipeline**

```text

raw FASTQ.gz -> baseline_fastqc -> trim -> decon -> processed_fastqc -> quality control multiqc_dashboard
                                             |
                                             v
                                          bwa align -> dedup -> index_bam -> flagstat / stats -> alignment multiqc_dashboard
                                                                     |
                                                                     v
                                                                faidx_reference -> variant_call -> index_vcf -> filter_variants
                                                                                                                      |
                                                                                                                      v
                                                                                                               variant_stats -> Full WF
                                                                                                                                 --or--
                                                                                                                                Variant
                                                                                                                         multiqc_dashboard
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

Fully automated from the Bash file, each execution run validates the configuration, downloads only parameter files that are missing, performs a Snakemake dry run to verify file setup and dependencies present, then initiates the workflow.

## Environment Setup

```bash
conda env create -f envs/env.yaml
conda activate bioinfo
```

The launcher checks that the bioinfo environment is active before running.

## Outputs

* results/fastqc/: Baseline and post-processing FastQC reports
  
* results/alignment/: Deduplicated BAM, BAM index, duplicate metrics, flagstat
  
* results/variants/: Raw and filtered VCFs with indexes, variant statistics
  
* results/deliverables/: Stage dashboards and statistics

* logs/: One log per rule
  
* benchmarks/: One runtime record per rule.
  *Note:* On macOS, Snakemake cannot collect CPU and memory statistics, so some benchmark files contain runtimes only from example
   execution.


**Deliverables By Stage:**

* Quality Control -> `qc_dashboard.html`

* Alignment -> `alignment_dashboard.html`, `<sample>_stats.txt`

* Variant Calling -> `variant_dashboard.html`

* All -> `multiqc_dashboard.html`, `<sample>_stats.txt`


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

* Separation of Concerns: The Bash launcher handles environment checks, downloads and checksums. It also initiates the dry-run and upon successful completion, executes the specified workflow. The Python stage registry validates the stage and configuration and returns shell-safe variables. Snakemake handles the bioinformatics workflow. Each part can be tested and changed on its own.

* Validation Before Execution: Unknown stages, unimplemented stages, missing parameters and unreadable YAML are caught by the registry before Snakemake starts, with a clear error message and no stack trace.

* Reproducibility: Run accession, reference, adapters, contaminants, variant thresholds, and OS specifications all live in one configuration file for easy replacement of arguments to analyze a different run. The reference genome is a tracked Snakemake input, so changing it re-runs alignment and every downstream step instead of reusing stale results.

* Extensibility: All future stages need only one entry in the stage registry and its Snakemake rule(s), along with configuration for any
  parameters


## Testing

Execution Call:

```bash
python -m pytest -v
```
The tests cover stage validation, parameter validation, unknown and unimplemented stages, invalid YAML, resource resolution, the variant-calling stage, and shell-safe quoting.


## Troubleshooting

* Environment Error: `conda activate bioinfo`
 
* `MISSING:` Error: Add the named parameter to config/config.yaml
  
* `UNKNOWN_STAGE` error: Run the script with `--help` for valid stages
* 
* Automated Download by Accession Keeps Failing:
  
  1. Check the accession and your connection, downloads auto-retry three times
   
  2. Check terminal to review details
     
  3. If still failing, uncomment the Manual Override urls in the configuration file under sample and input the direct links to your files
     of interest. Example execution links are included to demonstrate input requirements. M5 checksums can be input along with the direct
     links, or they can be left as is with "" to skip checks
     
* Checksum Mismatch: Delete raw fastq files and attempt rerun- may contain residual files
  
* Workflow Problems: Run `snakemake -n <target>` to inspect the dependency graph. If unable to identify problem, please contact author on
  GitHub asap.


## Future Additions

* Docker image for fully portable execution
  
* Assembly (SPAdes)

* Annotation (Prokka or Bakta)

* Coverage Analysis (mosdepth)

* Circos Visualization (pyCirclize)
  
* Anti-Microbial Resistance Profiling (AMRFinderPlus)

* Taxonomic Classification (Kraken2)

* Per-stage tool requirement checks in place of the environment-name check
