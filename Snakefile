import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.4.4")

##### load config and sample sheets #####
configfile: "src/config.yaml"

validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")
##### target rules #####

rule all:
    input:
        # igvtools
        expand("analysis/igvtools/{samples.sample}.tdf", samples=samples.itertuples()),
        # multiQC
        "qc/multiqc/multiqc_report.html",

rule multiqc:
    input:
        # fastqc
        expand("qc/fastqc/{samples.sample}_fastqc.html", samples=samples.itertuples()),
        expand("qc/fastqc/{samples.sample}_fastqc.zip", samples=samples.itertuples()),
        # trim_galore_se
        expand("trimmed/{samples.sample}.fastq.gz_trimming_report.txt", samples=samples.itertuples()),
        # bwa_mem_qc
        expand("qc/bwa/{samples.sample}.idxstats", samples=samples.itertuples()),
        expand("qc/bwa/{samples.sample}.stats", samples=samples.itertuples()),
        expand("qc/bwa/{samples.sample}.flagstat", samples=samples.itertuples()),
        # deeptools
        "analysis/deeptools/multibamsum.npz",
        "analysis/deeptools/multibamsum.tab",
        "analysis/deeptools/pearsoncor_multibamsum.png",
        "analysis/deeptools/pearsoncor_multibamsum_matrix.txt",
        "analysis/deeptools/multibamsum_cov.png",
        "analysis/deeptools/multibamsum_cov_counts.txt",
        "analysis/deeptools/multibamsum_fingerprint.png",
        "analysis/deeptools/multibamsum_fingerprint_metrics.txt",
        "analysis/deeptools/multibamsum_fingerprint_rawcounts.txt",
        # macs2
        expand("analysis/macs2/{id}_model.r", id=samples.id.unique().tolist()),
        expand("analysis/macs2/{id}_peaks.narrowPeak", id=samples.id.unique().tolist()),
        expand("analysis/macs2/{id}_peaks.xls", id=samples.id.unique().tolist()),
        expand("analysis/macs2/{id}_summits.bed", id=samples.id.unique().tolist()),
    params:
        "qc/fastqc/",
        "trimmed/",
        "qc/bwa/",
        "analysis/deeptools/",
        "analysis/macs2/",
    output:
        # "qc/multiqc/multiqc_data_1/multiqc_macs.txt",
        # "qc/multiqc/multiqc_data/multiqc_samtools_stats.txt",
        # "qc/multiqc/multiqc_data/multiqc_samtools_flagstat.txt",
        # "qc/multiqc/multiqc_data/multiqc_samtools_idxstats.txt",
        # "qc/multiqc/multiqc_data/multiqc_cutadapt.txt",
        # "qc/multiqc/multiqc_data/multiqc_fastqc.txt",
        # "qc/multiqc/multiqc_data/multiqc_data.json",
        # "qc/multiqc/multiqc_data/multiqc_general_stats.txt",
        # "qc/multiqc/multiqc_data/multiqc_sources.txt",
        # "qc/multiqc/multiqc_data/multiqc.log",
        directory("qc/multiqc/multiqc_report_data/"),
        "qc/multiqc/multiqc_report.html",
    conda:
        "envs/multiqc.yaml"
    log:
        "logs/multiqc.log"
    shell:
        """
        multiqc -f {params} \
        -o qc/multiqc/ \
        -n multiqc_report.html
        """

rule symlink:
    output:
        expand("raw_data/{samples.sample}.fastq.gz", samples=samples.itertuples()),
    log:
        "logs/symlink.log"
    shell:
        "Rscript src/symlink.R"


rule fastqc:
    input:
        "raw_data/{sample}.fastq.gz"
    output:
        html="qc/fastqc/{sample}_fastqc.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    params: ""
    log:
        "logs/fastqc.{sample}.log"
    wrapper:
        "0.35.1/bio/fastqc"

rule trim_galore_se:
    input:
        # get_fastq
        "raw_data/{sample}.fastq.gz"
    output:
        "trimmed/{sample}_trimmed.fq.gz",
        "trimmed/{sample}.fastq.gz_trimming_report.txt",
    params:
        extra = "--retain_unpaired -q 20"
    log:
        "logs/trim_galore.{sample}.log"
    wrapper:
        "0.31.1/bio/trim_galore/se"

rule bwa_mem:
    input:
        "trimmed/{sample}_trimmed.fq.gz"
    output:
        bam     = "analysis/bwa/{sample}.bam",
        bai     = "analysis/bwa/{sample}.bam.bai"
    params:
        index = config["ref"]["index"]
    threads: 32
    log:
        "logs/bwa_mem.{sample}.log"
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} -M {params.index} {input} | samtools view -hbu -q 20 -F 4 - | samtools sort -@ {threads} -O BAM -o {output.bam} -
        samtools index {output.bam}
        """

rule bwa_mem_qc:
    input:
        "analysis/bwa/{sample}.bam"
    output:
        idxstats = "qc/bwa/{sample}.idxstats",
        stats = "qc/bwa/{sample}.stats",
        flagstat = "qc/bwa/{sample}.flagstat"
    log:
        "logs/bwa_mem_qc.{sample}.log"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools idxstats {input} > {output.idxstats}
        samtools stats {input} > {output.stats}
        samtools flagstat {input} > {output.flagstat}
        """
rule bam_to_tdf:
    input: "analysis/bwa/{sample}.bam"
    output: "analysis/igvtools/{sample}.tdf"
    params:
        igv_genome = "mm10"
    log:
        "logs/bam_to_tdf.{sample}.log"
    conda:
        "envs/igvtools.yaml"
    shell:
        """
        igvtools count \
        --minMapQuality 20 \
        -z 5 \
        -w 25 \
        -e 225 \
        {input} \
        {output} \
        {params.igv_genome}
        """

rule deeptools_summary:
    input:
        expand("analysis/bwa/{samples.sample}.bam", samples=samples.itertuples()),
    output:
        sum     = "analysis/deeptools/multibamsum.npz",
        counts  = "analysis/deeptools/multibamsum.tab"
    params:
        labels=samples["sample"].tolist()
    threads: 32
    log:
        "logs/deeptools_summary.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        """
        multiBamSummary bins \
        -p {threads} \
        -b {input} \
        --minMappingQuality 20 \
        -e 225 \
        --ignoreDuplicates \
        --centerReads \
        --labels {params.labels} \
        -out {output.sum} \
        --outRawCounts {output.counts}
        """
rule deeptools_correlation:
    input: "analysis/deeptools/multibamsum.npz"
    output:
        fig     = "analysis/deeptools/pearsoncor_multibamsum.png",
        matrix  = "analysis/deeptools/pearsoncor_multibamsum_matrix.txt"
    conda:
        "envs/deeptools.yaml"
    log:
        "logs/deeptools_correlation.log"
    shell:
        """
        plotCorrelation \
        --corData {input} \
        --plotFile {output.fig} \
        --outFileCorMatrix {output.matrix} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --skipZeros \
        --plotTitle "Pearson Correlations of BWA Alignments" \
        --plotNumbers \
        --colorMap RdYlBu
        """

rule deeptools_coverage:
    input: expand("analysis/bwa/{samples.sample}.bam", samples=samples.itertuples()),
    output:
        fig     = "analysis/deeptools/multibamsum_cov.png",
        counts  = "analysis/deeptools/multibamsum_cov_counts.txt"
    threads: 32
    params:
        labels=samples["sample"].tolist()
    conda:
        "envs/deeptools.yaml"
    log:
        "logs/deeptools_coverage.log"
    shell:
        """
        plotCoverage \
        -p {threads} \
        -b {input} \
        --plotFile {output.fig} \
        --outRawCounts {output.counts} \
        -e 225 \
        --plotTitle "Coverage of BWA Alignments" \
        --labels {params.labels} \
        --minMappingQuality 20 \
        --ignoreDuplicates \
        --skipZeros \
        --centerReads
        """

rule deeptools_fingerprint:
    input: expand("analysis/bwa/{samples.sample}.bam", samples=samples.itertuples()),
    output:
        fig         = "analysis/deeptools/multibamsum_fingerprint.png",
        metrics     = "analysis/deeptools/multibamsum_fingerprint_metrics.txt",
        rawcounts   = "analysis/deeptools/multibamsum_fingerprint_rawcounts.txt"
    threads: 32
    params:
        labels = samples["sample"].tolist()
    conda:
        "envs/deeptools.yaml"
    log:
        "logs/deeptools_fingerprint.log"
    shell:
        """
        plotFingerprint -p {threads} \
        -b {input} \
        --plotFile {output.fig} \
        --outQualityMetrics {output.metrics} \
        --outRawCounts {output.rawcounts} \
        -e 225 \
        --plotTitle "ChIP Enrichment Profile of BWA Alignments" \
        --labels {params.labels} \
        --minMappingQuality 20 \
        --ignoreDuplicates \
        --skipZeros \
        --centerReads
        """

rule macs2_callpeaks:
    input:
        chip    = "analysis/bwa/{id}ChIP.bam",
        input   = "analysis/bwa/{id}Input.bam",
    output:
        "analysis/macs2/{id}_model.r",
        "analysis/macs2/{id}_peaks.narrowPeak",
        "analysis/macs2/{id}_peaks.xls",
        "analysis/macs2/{id}_summits.bed",
    params:
        id          = "{id}",
        organism    = config["macs2"]["gsize"]
    conda:
        "envs/macs2.yaml"
    log:
        "logs/macs2_callpeaks.{id}.log"
    shell:
        """
        macs2 callpeak \
        --name {params.id} \
        --treatment {input.chip} \
        --control {input.input} \
        --gsize {params.organism} \
        --outdir analysis/macs2
        """
