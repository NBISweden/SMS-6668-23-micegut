import os
import pandas as pd
import glob

results = config["results_dir"]
samples = pd.read_csv(config["sample_list"], sep="\t", index_col=0)
genomes = [os.path.basename(x.replace(".faa", "")) for x in glob.glob(results + "/genomes/annotations/genes/*.faa")]

localrules:
    quantify_taxonomy,
    download_rgi_data

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

def split_ranks(df, ranks):
    uniq_lin = df.lineage.unique()
    d = {}
    for lin in uniq_lin:
        d[lin] = dict(zip(ranks,lin.split(";")))
    return pd.DataFrame(d).T


rule quantify_taxonomy:
    output:
        cov=results + "/{sample}/annotation/mmseqs2/{seqTaxDB}.median_fold.tsv",
        contig_cov= results + "/{sample}/annotation/mmseqs2/{seqTaxDB}.contig.median_fold.tsv"
    input:
        tsv=rules.mmseqs2_createtsv.output.tsv, # e.g. atlas/C1/annotation/mmseqs2/UniRef100.taxonomy.tsv
        cov=results + "/{sample}/assembly/contig_stats/prefilter_coverage_stats.txt" # e.g. atlas/C1/assembly/contig_stats/prefilter_coverage_stats.txt
    params:
        ranks = config["taxonomy"]["ranks"]
    run:
        import pandas as pd
        # Read taxonomic assignments and fill NA values
        taxdf = pd.read_csv(input.tsv, sep="\t", header=None, index_col=0, usecols=[0,8], names=["contig","lineage"])
        taxdf.fillna("uc", inplace=True)
        # Read coverage info for contigs, only storing the median fold values
        covdf = pd.read_csv(input.cov, sep="\t", header=0, index_col=0, usecols=[0,9], names=["contig","Median_fold"])
        # Extract ranks from unique assignments
        lineage_df = split_ranks(taxdf, params.ranks)
        # Merge rank assignments to taxonomy dataframe
        dataframe = pd.merge(taxdf, lineage_df, left_on="lineage", right_index=True)
        # Merge with coverage, taking the outer join
        dataframe = pd.merge(dataframe, covdf, left_index=True, right_index=True, how="outer")
        # Fill NA values
        dataframe.loc[dataframe.lineage != dataframe.lineage, "lineage"] = ";".join(["unknown"]*len(params.ranks))
        dataframe.fillna("unknown", inplace=True)
        # Sum median fold values to lineage
        lineage_cov = dataframe.groupby(["lineage"]+params.ranks).sum(numeric_only=True)
        # Write contig taxonomy + summed median fold values to files
        dataframe.drop("lineage", axis=1).to_csv(output.contig_cov, sep="\t")
        lineage_cov.to_csv(output.cov, sep="\t")

def merge_taxcov(samples, results, seqTaxDB):
    df = pd.DataFrame()
    taxdict = {}
    for sample in samples:
        f = f"{results}/{sample}/annotation/mmseqs2/{seqTaxDB}.median_fold.tsv"
        _df = pd.read_csv(f,sep="\t",index_col=0)
        _df.rename(columns={"Median_fold": sample},inplace=True)
        taxcols = _df.drop(sample,axis=1)
        _taxdict = taxcols.to_dict(orient="index")
        taxdict.update(_taxdict)
        df = pd.merge(df,pd.DataFrame(_df.loc[:, sample]),right_index=True,left_index=True,how="outer")
    return df, taxdict

rule collate_taxcov:
    output:
        tsv=results + "/taxonomy/{seqTaxDB}.median_fold.tsv"
    input:
        expand(results + "/{sample}/annotation/mmseqs2/{{seqTaxDB}}.median_fold.tsv", sample=samples.index)
    run:
        df = pd.DataFrame()
        df, taxdict = merge_taxcov(samples.index, results, wildcards.seqTaxDB)
        dataframe = pd.merge(pd.DataFrame(taxdict).T, df, left_index=True, right_index=True)
        dataframe.fillna(0, inplace=True)
        dataframe.index.name = "lineage"
        dataframe.to_csv(output.tsv, sep="\t")


##### resistance gene identifier #####

rule download_rgi_data:
    output:
        json="resources/card/card.json",
        version="resources/card/card.version"
    log:
        "resources/card/log"
    params:
        tar="resources/card/data.tar.gz",
        dir=lambda w, output: os.path.dirname(output.json)
    shell:
         """
         curl -L -o {params.tar} \
            https://card.mcmaster.ca/latest/data >{log} 2>&1
         tar -C {params.dir} -xf {params.tar} ./card.json
         # Store download date in versionfile
         date > {output.version}
         rm {params.tar}
         """

rule rgi_genecatalog:
    input:
        faa=results+"/Genecatalog/gene_catalog.faa",
        db="resources/card/card.json"
    output:
        json=results+"/Genecatalog/annotations/rgi.out.json",
        txt=results+"/Genecatalog/annotations/rgi.out.txt"
    log:
        results+"/logs/rgi/rgi.log",
    params:
        out=results+"/Genecatalog/annotations/rgi.out",
        settings="-a diamond --local --clean --input_type protein -d wgs",
        tmpdir="$TMPDIR/genecatalog.rgi",
        faa="$TMPDIR/genecatalog.rgi/gene_catalog.faa",
    shadow: "minimal"
    conda:
        "envs/rgi.yml"
    threads: 10
    resources:
        runtime=60*24,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
    shell:
        """
        exec &>{log}
        mkdir -p {params.tmpdir}
        rgi load --card_json {input.db} --local
        sed 's/*//g' {input.faa} > {params.faa}
        rgi main -i {params.faa} -o {params.out} -n {threads} {params.settings}
        rm -r {params.tmpdir}
        """

rule rgi_genomes:
    output:
        json=results+"/genomes/annotations/rgi/genomes.out.json",
        txt=results+"/genomes/annotations/rgi/genomes.out.txt"
    input:
        faa=expand(results+"/genomes/annotations/genes/{genome}.faa", genome = genomes),
        db="resources/card/card.json"
    log:
        results+"/logs/genomes/rgi/genomes.log",
    params:
        out=results+"/genomes/annotations/rgi/genomes.out",
        settings="-a diamond --local --clean --input_type protein -d wgs",
        tmpdir="$TMPDIR/genomes.rgi",
        faa="$TMPDIR/genomes.rgi/genomes.faa",
    shadow: "minimal"
    conda:
        "envs/rgi.yml"
    threads: 10
    resources:
        runtime=60*10,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
    shell:
        """
        exec &>{log}
        mkdir -p {params.tmpdir}
        rgi load --card_json {input.db} --local
        sed 's/*//g' {input.faa} > {params.faa}
        rgi main -i {params.faa} -o {params.out} -n {threads} {params.settings} -d wgs
        rm -r {params.tmpdir}
        """

rule rgi_parse_genecatalog:
    output:
        tsv=results + "/Genecatalog/annotations/rgi.parsed.tsv"
    input:
        txt=rules.rgi_genecatalog.output.txt
    run:
        annot = pd.read_csv(input.txt, sep="\t", index_col=0)
        annot = annot.loc[:,
                ["Model_ID", "AMR Gene Family", "Resistance Mechanism"]]
        annot.loc[:, "Model_ID"] = ["RGI_{}".format(x) for x in annot.Model_ID]
        annot.rename(index=lambda x: x.split(" ")[0], inplace=True)
        annot.to_csv(output.tsv, sep="\t", index=True)
