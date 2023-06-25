import os
import pandas as pd

samples = pd.read_csv(config["sample_list"], sep="\t", index_col=0)
results = config["results_dir"]

def mem_allowed(wildcards, threads):
    mem_per_core = config["mem_per_core"]
    return max(threads * mem_per_core, mem_per_core)

rule all:
    input:
        html=expand(results + "/{sample}/annotation/mmseqs2/{seqTaxDB}.report.html",
                    sample = list(samples.index), seqTaxDB = config["taxonomy"]["seqTaxDB"]),
        tsv=expand(results + "/{sample}/annotation/mmseqs2/{seqTaxDB}.taxonomy.tsv",
                   sample = list(samples.index), seqTaxDB = config["taxonomy"]["seqTaxDB"]),

rule mmseqs2_create_queryDB:
    output:
        queryDB=expand(
            results + "/{{sample}}/annotation/mmseqs2/queryDB{suff}",
            suff=["", "_h", ".index", "_h.index"],
        ),
    input:
        fa=results + "/{sample}/assembly/{sample}_final_contigs.fasta",
    log:
        results + "/{sample}/logs/mmseqs2/createquerydb.log",
    params:
        mem_mib=mem_allowed,
        outdir=lambda wildcards, output: os.path.dirname(output.queryDB[0]),
        tmpdir="$TMPDIR/{sample}.mmseqs2",
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "envs/mmseqs2.yml"
    resources:
        runtime=60 * 2,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"],
    shell:
        """
        mkdir -p {params.tmpdir}
        mmseqs createdb {input.fa} {params.tmpdir}/queryDB --dbtype 2 --shuffle 0 \
            --createdb-mode 1 --write-lookup 0 --id-offset 0 --compressed 0 -v 3 > {log} 2>&1
        mv {params.tmpdir}/queryDB* {params.outdir}
        rm -rf {params.tmpdir}
        """

rule mmseqs2_taxonomy:
    output:
        expand(
            results + "/{{sample}}/annotation/mmseqs2/{{seqTaxDB}}{suff}",
            suff=[".dbtype", ".index", "_aln.dbtype", "_aln.index"],
        ),
    input:
        queryDB=rules.mmseqs2_create_queryDB.output.queryDB[0],
        seqTaxDB=expand(config["database_dir"] + "/mmseqs2/{{seqTaxDB}}{suff}",
                        suff=["","_taxonomy","_mapping","_h.dbtype","_h.index","_h",".lookup",".dbtype",".index"]),
    log:
        results + "/{sample}/logs/mmseqs2/{seqTaxDB}.taxonomy.log",
    params:
        lca_ranks=",".join(config["taxonomy"]["ranks"]),
        out=lambda wildcards, output: os.path.dirname(output[0])+ "/"+ wildcards.seqTaxDB,
        tmpdir="$TMPDIR/mmseqs.taxonomy.{sample}",
        mem_mib=mem_allowed,
    threads: 20
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "envs/mmseqs2.yml"
    resources:
        runtime=1440,
        mem_mib=mem_allowed,
    shell:
        """
        mkdir -p {params.tmpdir}
        mmseqs taxonomy {input.queryDB} {input.seqTaxDB[0]} {params.out} {params.tmpdir} \
            --lca-mode 3 --tax-output-mode 2 --lca-ranks {params.lca_ranks} --tax-lineage 1 \
            --threads {threads} --local-tmp {params.tmpdir} --remove-tmp-files 1 > {log} 2>&1 
        """

rule mmseqs2_createtsv:
    output:
        tsv=results + "/{sample}/annotation/mmseqs2/{seqTaxDB}.taxonomy.tsv",
    input:
        queryDB=rules.mmseqs2_create_queryDB.output.queryDB[0],
        taxRes=rules.mmseqs2_taxonomy.output,
    log:
        results + "/{sample}/logs/mmseqs2/{seqTaxDB}.createtsv.log",
    params:
        taxRes=lambda wildcards, input: os.path.dirname(input.taxRes[0])+ "/" + wildcards.seqTaxDB,
        mem_mib=mem_allowed,
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "envs/mmseqs2.yml"
    threads: 10
    resources:
        runtime=60,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"],
    shell:
        """
        mmseqs createtsv {input.queryDB} {params.taxRes} {output.tsv} \
            --threads {threads} --first-seq-as-repr 0 --target-column 1 \
            --full-header 0 --idx-seq-src 0 --db-output 0 --compressed 0 -v 3
        """

rule mmseqs2_kronareport:
    output:
        html=results + "/{sample}/annotation/mmseqs2/{seqTaxDB}.report.html"
    input:
        taxRes=rules.mmseqs2_taxonomy.output,
        seqTaxDB=expand(config["database_dir"] + "/mmseqs2/{{seqTaxDB}}{suff}",
                        suff=["","_taxonomy","_mapping","_h.dbtype","_h.index","_h",".lookup",".dbtype",".index"]),
    log:
        results + "/{sample}/logs/mmseqs2/{seqTaxDB}.kronareport.log"
    params:
        taxRes=lambda wildcards, input: os.path.dirname(input.taxRes[0])+ "/" + wildcards.seqTaxDB,
        mem_mib=mem_allowed,
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "envs/mmseqs2.yml"
    resources:
        runtime=60,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
    shell:
        """
        mmseqs taxonomyreport {input.seqTaxDB[0]} {params.taxRes} {output.html} --report-mode 1 >{log} 2>&1 
        """
