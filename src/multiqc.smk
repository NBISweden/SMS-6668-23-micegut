import os

localrules: multiqc, all

sample_file = config["samples"]

samples = []
with open(sample_file, 'r') as fhin:
    for i, line in enumerate(fhin):
        if i == 0:
            continue
        samples.append(line.rstrip().rsplit()[0])

rule all:
    input:
        "atlas/multiqc/multiqc_report.html",
        expand("atlas/fastqc/{sample}_QC_R{n}_fastqc.zip", 
            sample = samples, n = ["1","2"])

rule fastqc:
    output:
        "atlas/fastqc/{sample}_QC_R{n}_fastqc.zip"
    input:
        "atlas/{sample}/sequence_quality_control/{sample}_QC_R{n}.fastq.gz"
    log:
        "atlas/logs/fastqc/{sample}_R{n}.fastqc.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    resources:
        runtime = 60
    shell:
        """
        mkdir -p $TMPDIR/{wildcards.sample}.fastqc
        fastqc -o $TMPDIR/{wildcards.sample}.fastqc {input} > {log} 2>&1
        mv $TMPDIR/{wildcards.sample}.fastqc/{wildcards.sample}_QC_R{wildcards.n}_fastqc.zip {params.outdir}
        rm -rf $TMPDIR/{wildcards.sample}.fastqc
        """

rule multiqc:
    output:
        "atlas/multiqc/multiqc_report.html"
    input:
        bbduk = expand("atlas/{sample}/logs/QC/quality_filter.log", 
            sample = samples),
        fastqc = expand("atlas/fastqc/{sample}_QC_R{n}_fastq.zip", 
            sample = samples, n = ["1", "2"])
    log:
        "atlas/logs/multiqc/multiqc.log"
    shell:
        """
        multiqc --outdir {params.outdir} {input} > {log} 2>&1
        """
