# project args
configfile: "config/config.yaml"

# project configuration values shared across workflow rules
# sequencing run param
sample = config["sample"]["accession"]
# reference genome param
reference_name = config["reference"]["name"]
# threads param for multi-thread tools
num_threads = config["threads"]
# trimming param
adapter_type = config["trimming"]["adapter_type"]
# decontamination param
spike_ins = config["decontamination"]["spike_ins"]
# variant calling params
ploidy = config["variant_calling"]["ploidy"]
min_base_quality = config["variant_calling"]["min_base_quality"]
min_map_quality = config["variant_calling"]["min_map_quality"]
min_variant_quality = config["variant_calling"]["min_variant_quality"]
min_depth = config["variant_calling"]["min_depth"]
min_allele_fraction = config["variant_calling"]["min_allele_fraction"]

# shared reference path, reused by every rule
reference = f"resources/reference/{reference_name}.fa"

# stage target set defining files indicating each stage success
# aggregate target for quality-control workflow only
rule qc:
    input:
        f"results/deliverables/qc_dashboard.html"

# aggregate target for alignment workflow only
rule align:
    input:
        f"results/deliverables/alignment_dashboard.html",
        f"results/deliverables/{sample}_stats.txt"

# aggregate target for variant calling stage only
rule variant_calling:
    input:
        f"results/variants/{sample}_filtered.vcf.gz.tbi",
        f"results/variants/{sample}_variant_stats.txt",
        "results/deliverables/variant_dashboard.html"

# aggregate target representing entire workflow
rule all:
    input:
        f"results/deliverables/{sample}_stats.txt",
        f"results/deliverables/multiqc_dashboard.html",
        f"results/variants/{sample}_filtered.vcf.gz",
        f"results/variants/{sample}_variant_stats.txt"

# rule 1: baseline fastqc raw seq paired end run files
rule baseline_fastqc:
    input:
        r1=f"raw/{sample}_1.fastq.gz",
        r2=f"raw/{sample}_2.fastq.gz"
    output:
        html1=f"results/fastqc/baseline/{sample}_1_fastqc.html",
        zip1=f"results/fastqc/baseline/{sample}_1_fastqc.zip",
        html2=f"results/fastqc/baseline/{sample}_2_fastqc.html",
        zip2=f"results/fastqc/baseline/{sample}_2_fastqc.zip"
    benchmark:
        f"benchmarks/baseline_fastqc_{sample}.txt"
    log:
        f"logs/baseline_fastqc_{sample}.log"
    threads: 2
    shell:
        """
        fastqc \
            -t {threads} \
            {input.r1} \
            {input.r2} \
            --outdir results/fastqc/baseline \
            > {log} 2>&1
        """

# rule 2: trim adapters & low quality bases
rule trim:
    input:
        r1=f"raw/{sample}_1.fastq.gz",
        r2=f"raw/{sample}_2.fastq.gz"
    output:
        r1=f"results/trimmed/{sample}_1_trimmed.fastq.gz",
        r2=f"results/trimmed/{sample}_2_trimmed.fastq.gz",
        u1=f"results/trimmed/{sample}_1_unpaired.fastq.gz",
        u2=f"results/trimmed/{sample}_2_unpaired.fastq.gz"
    params:
        adapters=f"resources/adapters/{adapter_type}.fa"
    benchmark:
        f"benchmarks/trim_{sample}.txt"
    log:
        f"logs/trim_{sample}.log"
    shell:
        """
        trimmomatic PE \
            {input.r1} \
            {input.r2} \
            {output.r1} \
            {output.u1} \
            {output.r2} \
            {output.u2} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            > {log} 2>&1
        """

# rule 3: remove spike_in contaminants from reads
rule decon:
    input:
        r1=f"results/trimmed/{sample}_1_trimmed.fastq.gz",
        r2=f"results/trimmed/{sample}_2_trimmed.fastq.gz"
    output:
        r1=f"results/processed/{sample}_1_decon.fastq.gz",
        r2=f"results/processed/{sample}_2_decon.fastq.gz"
    params:
        contaminant=spike_ins
    benchmark:
        f"benchmarks/decon_{sample}.txt"
    log:
        f"logs/decon_{sample}.log"
    shell:
        """
        bbduk.sh \
            in1={input.r1} \
            in2={input.r2} \
            out1={output.r1} \
            out2={output.r2} \
            ref={params.contaminant} \
            k=31 \
            hdist=1 \
            > {log} 2>&1
        """

# rule 4: generate processed fastqc for comparison
rule processed_fastqc:
    input:
        r1=f"results/processed/{sample}_1_decon.fastq.gz",
        r2=f"results/processed/{sample}_2_decon.fastq.gz"
    output:
        html1=f"results/fastqc/postqc/{sample}_1_decon_fastqc.html",
        zip1=f"results/fastqc/postqc/{sample}_1_decon_fastqc.zip",
        html2=f"results/fastqc/postqc/{sample}_2_decon_fastqc.html",
        zip2=f"results/fastqc/postqc/{sample}_2_decon_fastqc.zip"
    benchmark:
        f"benchmarks/processed_fastqc_{sample}.txt"
    log:
        f"logs/processed_fastqc_{sample}.log"
    threads: 2
    shell:
        """
        fastqc \
            -t {threads} \
            {input.r1} \
            {input.r2} \
            --outdir results/fastqc/postqc \
            > {log} 2>&1
        """

# rule 5-qc: consolidate preprocessing and quality-control metrics into stage dashboard
rule qc_dashboard:
    input:
        f"results/fastqc/baseline/{sample}_1_fastqc.zip",
        f"results/fastqc/baseline/{sample}_2_fastqc.zip",
        f"results/fastqc/postqc/{sample}_1_decon_fastqc.zip",
        f"results/fastqc/postqc/{sample}_2_decon_fastqc.zip",
        f"logs/trim_{sample}.log",
        f"logs/decon_{sample}.log"
    output:
        html="results/deliverables/qc_dashboard.html"
    log:
        "logs/qc_multiqc.log"
    shell:
        """
        multiqc results/fastqc logs \
            -o results/qc_multiqc \
            -n qc_dashboard.html \
            --force \
            > {log} 2>&1

        mv results/qc_multiqc/qc_dashboard.html \
           results/deliverables/qc_dashboard.html
        """

# rule 6: align processed reads to indexed reference genome with bwa
rule bwa_align:
    input:
        r1=f"results/processed/{sample}_1_decon.fastq.gz",
        r2=f"results/processed/{sample}_2_decon.fastq.gz",
        ref=reference
    output:
        bam=temp(f"results/alignment/{sample}.bam")
    params:
        rg=f"@RG\\tID:{sample}\\tPL:illumina\\tPU:flowcell1\\tSM:{sample}"
    benchmark:
        f"benchmarks/alignment_{sample}.txt"
    log:
        f"logs/alignment_{sample}.log"
    threads:
        num_threads
    resources:
        mem_mb=4000  # realistic ceiling for bacterial genome
    shell:
        """
        (
        bwa mem \
            -M \
            -R "{params.rg}" \
            -t {threads} \
            {input.ref} \
            {input.r1} \
            {input.r2} | \
        samtools view \
            -@ {threads} \
            -b \
            -o {output.bam} -
        ) > {log} 2>&1
        """

# rule 7: mark & remove pcr duplicates
rule dedup:
    input:
        bam=f"results/alignment/{sample}.bam"
    output:
        # aligned output cleaned to prevent fp's for downstream variant analysis
        bam=f"results/alignment/{sample}_dedup.bam",
        metrics=f"results/alignment/{sample}_markdup_metrics.txt"
    benchmark:
        f"benchmarks/dedup_{sample}.txt"
    log:
        f"logs/dedup_{sample}.log"
    threads:
        num_threads
    shell:
        """
        samtools sort \
            -n \
            -@ {threads} \
            {input.bam} | \
        samtools fixmate \
            -m \
            - \
            - | \
        samtools sort \
            -@ {threads} \
            - | \
        samtools markdup \
            -r \
            -s \
            -f {output.metrics} \
            - \
            {output.bam} \
            2> {log} 
        """

# rule 8: build aligned, deduplicated bam index
rule index_bam:
    input:
        bam=f"results/alignment/{sample}_dedup.bam"
    output:
        bai=f"results/alignment/{sample}_dedup.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

# rule 9: collect overview mapping quality statistics
rule flagstat:
    input:
        bam=f"results/alignment/{sample}_dedup.bam"
    output:
        f"results/alignment/{sample}_flagstat.txt"
    shell:
        """
        samtools flagstat {input.bam} > {output}
        """

# rule 10: collect comprehensive alignment statistics
rule stats:
    input:
        bam=f"results/alignment/{sample}_dedup.bam"
    output:
        f"results/deliverables/{sample}_stats.txt"
    shell:
        """
        samtools stats {input.bam} > {output}
        """

# rule 11-align: aggregate alignment metrics into stage dashboard
rule alignment_dashboard:
    input:
        f"results/alignment/{sample}_markdup_metrics.txt",
        f"results/alignment/{sample}_flagstat.txt",
        f"results/deliverables/{sample}_stats.txt"
    output:
        html="results/deliverables/alignment_dashboard.html"
    log:
        "logs/alignment_multiqc.log"
    shell:
        """
        multiqc results/alignment \
            -o results/alignment_multiqc \
            -n alignment_dashboard.html \
            --force \
            > {log} 2>&1

        mv results/alignment_multiqc/alignment_dashboard.html \
           results/deliverables/alignment_dashboard.html
        """

# rule 11-all: combine all results into unified dashboard
rule multiqc_dashboard:
    input:
        # takes results from all stages to consolidate
        f"results/fastqc/baseline/{sample}_1_fastqc.zip",
        f"results/fastqc/baseline/{sample}_2_fastqc.zip",
        f"results/fastqc/postqc/{sample}_1_decon_fastqc.zip",
        f"results/fastqc/postqc/{sample}_2_decon_fastqc.zip",
        f"results/alignment/{sample}_markdup_metrics.txt",
        f"results/alignment/{sample}_flagstat.txt",
        f"results/deliverables/{sample}_stats.txt",
        f"results/variants/{sample}_variant_stats.txt"
    output:
        # High-value dashboard asset: in deliverables folder
        html=f"results/deliverables/multiqc_dashboard.html"
    log:
        "logs/multiqc.log"
    shell:
        """
        multiqc results \
            -o results/multiqc \
            -n multiqc_dashboard.html \
            --force \
            > {log} 2>&1
        # move dashboard to deliverables folder
        mv results/multiqc/multiqc_dashboard.html \
           results/deliverables/multiqc_dashboard.html
        """

# rule 12: index reference fasta so parallel jobs never race to create it
rule faidx_reference:
    input:
        reference
    output:
        f"{reference}.fai"
    log:
        f"logs/faidx_{reference_name}.log"
    shell:
        "samtools faidx {input} 2> {log}"

# rule 13: call variants with bcftools mpileup + call
rule variant_call:
    input:
        bam=f"results/alignment/{sample}_dedup.bam",
        bai=f"results/alignment/{sample}_dedup.bam.bai",
        ref=reference,
        fai=f"{reference}.fai"
    output:
        vcf=f"results/variants/{sample}_raw.vcf.gz"
    params:
        ploidy=ploidy,
        min_bq=min_base_quality,
        min_mq=min_map_quality
    benchmark:
        f"benchmarks/variant_call_{sample}.txt"
    log:
        f"logs/variant_call_{sample}.log"
    threads:
        num_threads
    shell:
        """
        (
        bcftools mpileup \
            --fasta-ref {input.ref} \
            --min-BQ {params.min_bq} \
            --min-MQ {params.min_mq} \
            --annotate FORMAT/AD,FORMAT/DP \
            --output-type u \
            --threads {threads} \
            {input.bam} | \
        bcftools call \
            --multiallelic-caller \
            --variants-only \
            --ploidy {params.ploidy} \
            --output-type z \
            --threads {threads} \
            --output {output.vcf} 
        ) > {log} 2>&1  
        """


# rule 14: hard-filter variants by quality, depth, & allele fraction
rule filter_variants:
    input:
        vcf=f"results/variants/{sample}_raw.vcf.gz"
    output:
        vcf=f"results/variants/{sample}_filtered.vcf.gz",
        tbi=f"results/variants/{sample}_filtered.vcf.gz.tbi"
    params:
        min_qual=min_variant_quality,
        min_dp=min_depth,
        min_af=min_allele_fraction
    benchmark:
        f"benchmarks/filter_variants_{sample}.txt"
    log:
        f"logs/filter_variants_{sample}.log"
    shell:
        """
        (
        bcftools filter \
            --include 'QUAL>={params.min_qual} && FORMAT/DP>={params.min_dp} && FORMAT/AD[0:1]/FORMAT/DP>={params.min_af}' \
            --output-type z \
            --output {output.vcf} \
            {input.vcf}
        bcftools index --tbi {output.vcf}
        ) > {log} 2>&1
        """

# rule 15: generate variant summary statistics
rule variant_stats:
    input:
        vcf=f"results/variants/{sample}_filtered.vcf.gz"
    output:
        stats=f"results/variants/{sample}_variant_stats.txt"
    log:
        f"logs/variant_stats_{sample}.log"
    shell:
        """
        bcftools stats {input.vcf} \
            > {output.stats} \
            2> {log}
        """

# rule 16: aggregate variant statistics into stage dashboard
rule variant_dashboard:
    input:
        f"results/variants/{sample}_variant_stats.txt"
    output:
        html="results/deliverables/variant_dashboard.html"
    log:
        "logs/variant_multiqc.log"
    shell:
        """
        multiqc results/variants \
            -o results/variant_multiqc \
            -n variant_dashboard.html \
            --force \
            > {log} 2>&1

        mv results/variant_multiqc/variant_dashboard.html \
           results/deliverables/variant_dashboard.html
        """
