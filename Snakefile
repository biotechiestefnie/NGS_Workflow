# project args
configfile: "config/config.yaml"

# project configuration values shared across workflow rules
sample = config["sample"]["accession"]
reference_name = config["sample"]["reference_name"]
num_threads = config["threads"]
adapter_type = config["trimming"]["adapter_type"]
spike_ins = config["decontamination"]["spike_ins"]

# stage target set defining files that indicate each stage success

# aggregate target for quality-control workflow only
rule qc:
    input:
        f"results/deliverables/qc_dashboard.html"

# aggregate target for alignment workflow only
rule align:
    input:
        f"results/deliverables/alignment_dashboard.html",
        f"results/deliverables/{sample}_stats.txt"

# aggregate target representing entire workflow
rule all:
    input:
        f"results/deliverables/{sample}_stats.txt",
        f"results/deliverables/multiqc_dashboard.html"

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

# rule 3: remove spike_in contaminants
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

# rule 4: generate processed fastqc
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

# rule 5-qc: aggregate preprocessing and quality-control metrics into stage dashboard
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

# rule 5: align processed reads to indexed reference genome with bwa
rule bwa_align:
    input:
        r1=f"results/processed/{sample}_1_decon.fastq.gz",
        r2=f"results/processed/{sample}_2_decon.fastq.gz"
    output:
        bam=f"results/alignment/{sample}.bam"
    params:
        reference_path=f"resources/reference/{reference_name}.fa",
        rg=f"@RG\\tID:{sample}\\tPL:illumina\\tPU:flowcell1\\tSM:{sample}"
    benchmark:
        f"benchmarks/alignment_{sample}.txt"
    log:
        f"logs/alignment_{sample}.log"
    threads:
        num_threads
    shell:
        """
        bwa mem \
            -M \
            -R "{params.rg}" \
            -t {threads} \
            {params.reference_path} \
            {input.r1} \
            {input.r2} \
            2> {log} | \
        samtools view \
            -@ {threads} \
            -Sb - \
            > {output.bam}
        """

# rule 6: mark & remove pcr duplicates
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
    # remove first alignment file containing pcr dupes
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
            2> {log} && \
        rm -f {input.bam}
        """

# rule 7: build aligned, deduplicated bam index
rule index_bam:
    input:
        bam=f"results/alignment/{sample}_dedup.bam"
    output:
        bai=f"results/alignment/{sample}_dedup.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

# rule 8: collect overview mapping quality statistics
rule flagstat:
    input:
        bam=f"results/alignment/{sample}_dedup.bam"
    output:
        f"results/alignment/{sample}_flagstat.txt"
    shell:
        """
        samtools flagstat {input.bam} > {output}
        """

# rule 9: collect comprehensive alignment statistics
rule stats:
    input:
        bam=f"results/alignment/{sample}_dedup.bam"
    output:
        f"results/deliverables/{sample}_stats.txt"
    shell:
        """
        samtools stats {input.bam} > {output}
        """

# rule 10-align: aggregate alignment metrics into stage dashboard
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

# rule 10-all: combine all results into unified dashboard
rule multiqc_dashboard:
    input:
        # takes results from qc & alignment stage to consolidate
        f"results/fastqc/baseline/{sample}_1_fastqc.zip",
        f"results/fastqc/baseline/{sample}_2_fastqc.zip",
        f"results/fastqc/postqc/{sample}_1_decon_fastqc.zip",
        f"results/fastqc/postqc/{sample}_2_decon_fastqc.zip",
        f"results/alignment/{sample}_markdup_metrics.txt",
        f"results/alignment/{sample}_flagstat.txt",
        f"results/deliverables/{sample}_stats.txt"
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
