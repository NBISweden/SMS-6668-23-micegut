import os
import pandas as pd
import glob

results = config["results_dir"]
samples = pd.read_csv(config["sample_list"], sep="\t", index_col=0)
genomes = [os.path.basename(x.replace(".faa", "")) for x in glob.glob(results + "/genomes/annotations/genes/*.faa")]

localrules:
    quantify_taxonomy,
    quantify_taxonomy_GC,
    download_rgi_data,
    rgi_parse_genecatalog,
    rgi_parse_genomes,
    sum_rgi_genecatalog

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

rule mmseqs2_create_queryDB_GC:
    """
    MMseqs2 on Genecatalog
    """
    input:
        fa=results + "/Genecatalog/gene_catalog.faa"
    output:
        queryDB=expand(
            results + "/Genecatalog/annotations/mmseqs2/queryDB{suff}",
            suff=["", "_h", ".index", "_h.index"],
        ),
    log:
        results + "logs/Genecatalog/mmseqs2/createquerydb.log",
    params:
        mem_mib=mem_allowed,
        outdir=lambda wildcards, output: os.path.dirname(output.queryDB[0]),
        tmpdir="$TMPDIR/Genecatalog.mmseqs2",
    envmodules:
        "bioinfo-tools",
        "MMseqs2/14-7e284",
    conda:
        "envs/mmseqs2.yml"
    resources:
        runtime=60 * 10,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"],
    shell:
        """
        mkdir -p {params.tmpdir}
        mmseqs createdb {input.fa} {params.tmpdir}/queryDB --dbtype 1 --shuffle 0 \
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

rule mmseqs2_taxonomy_GC:
    """
    MMseqs2 on Genecatalog
    """
    output:
        expand(
            results + "/Genecatalog/annotation/mmseqs2/{{seqTaxDB}}{suff}",
            suff=[".dbtype", ".index", "_aln.dbtype", "_aln.index"],
        ),
    input:
        queryDB=rules.mmseqs2_create_queryDB_GC.output.queryDB[0],
        seqTaxDB=expand(config["database_dir"] + "/mmseqs2/{{seqTaxDB}}{suff}",
                        suff=["","_taxonomy","_mapping","_h.dbtype","_h.index","_h",".lookup",".dbtype",".index"]),
    log:
        results + "/logs/Genecatalog/mmseqs2/{seqTaxDB}.taxonomy.log",
    params:
        lca_ranks=",".join(config["taxonomy"]["ranks"]),
        out=lambda wildcards, output: os.path.dirname(output[0])+ "/"+ wildcards.seqTaxDB,
        tmpdir="$TMPDIR/mmseqs.taxonomy.GC",
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

rule mmseqs2_createtsv_GC:
    output:
        tsv=results + "/Genecatalog/annotation/mmseqs2/{seqTaxDB}.taxonomy.tsv",
    input:
        queryDB=rules.mmseqs2_create_queryDB_GC.output.queryDB[0],
        taxRes=rules.mmseqs2_taxonomy_GC.output,
    log:
        results + "/logs/Genecatalog/logs/mmseqs2/{seqTaxDB}.createtsv.log",
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
        runtime=60 * 24 * 10,
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

rule mmseqs2_kronareport_GC:
    output:
        html=results + "/Genecatalog/annotation/mmseqs2/{seqTaxDB}.report.html"
    input:
        taxRes=rules.mmseqs2_taxonomy_GC.output,
        seqTaxDB=expand(config["database_dir"] + "/mmseqs2/{{seqTaxDB}}{suff}",
                        suff=["","_taxonomy","_mapping","_h.dbtype","_h.index","_h",".lookup",".dbtype",".index"]),
    log:
        results + "/logs/Genecatalog/mmseqs2/{seqTaxDB}.kronareport.log"
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


def parse_mmseqs2(tsv, cov, ranks, usecols=[0,8]):
    import pandas as pd
    # Read taxonomic assignments and fill NA values
    taxdf = pd.read_csv(tsv, sep="\t", header=None, index_col=0, usecols=usecols, names=["contig","lineage"])
    taxdf.fillna("uc", inplace=True)
    # Read coverage info for contigs, only storing the median fold values
    covdf = pd.read_csv(cov, sep="\t", header=0, index_col=0, usecols=[0,9], names=["contig","Median_fold"])
    # Extract ranks from unique assignments
    lineage_df = split_ranks(taxdf, ranks)
    # Merge rank assignments to taxonomy dataframe
    dataframe = pd.merge(taxdf, lineage_df, left_on="lineage", right_index=True)
    # Merge with coverage, taking the outer join
    dataframe = pd.merge(dataframe, covdf, left_index=True, right_index=True, how="outer")
    # Fill NA values
    dataframe.loc[dataframe.lineage != dataframe.lineage, "lineage"] = ";".join(["unknown"]*len(ranks))
    dataframe.fillna("unknown", inplace=True)
    # Sum median fold values to lineage
    lineage_cov = dataframe.groupby(["lineage"]+ranks).sum(numeric_only=True)
    return lineage_cov, dataframe

def parse_mmseqs2_GC(tsv, cov, ranks, usecols=[0,4]):
    import pandas as pd
    # Read taxonomic assignments and fill NA values
    taxdf = pd.read_csv(tsv, sep="\t", header=None, index_col=0, usecols=usecols, names=["Gene","lineage"])
    taxdf.fillna("uc", inplace=True)
    # Extract ranks from unique assignments
    lineage_df = split_ranks(taxdf, ranks)
    # Merge rank assignments to taxonomy dataframe
    dataframe = pd.merge(taxdf, lineage_df, left_on="lineage", right_index=True)
    # Merge with coverage, taking the outer join
    dataframe = pd.merge(dataframe, cov, left_index=True, right_index=True, how="outer")
    # Fill NA values
    dataframe.loc[dataframe.lineage != dataframe.lineage, "lineage"] = ";".join(["unknown"]*len(ranks))
    dataframe.fillna("unknown", inplace=True)
    # Sum median fold values to lineage
    lineage_cov = dataframe.groupby(["lineage"]+ranks).sum(numeric_only=True)
    return lineage_cov, dataframe

rule quantify_taxonomy:
    output:
        cov=results + "/{sample}/annotation/mmseqs2/{seqTaxDB}.median_fold.tsv",
        contig_cov= results + "/{sample}/annotation/mmseqs2/{seqTaxDB}.contig.median_fold.tsv"
    input:
        tsv=rules.mmseqs2_createtsv.output.tsv, # e.g. atlas/C1/annotation/mmseqs2/UniRef100.taxonomy.tsv
        cov=results + "/{sample}/assembly/contig_stats/postfilter_coverage_stats.txt" # e.g. atlas/C1/assembly/contig_stats/prefilter_coverage_stats.txt
    params:
        ranks = config["taxonomy"]["ranks"]
    run:
        lineage_cov, dataframe = parse_mmseqs2(input.tsv, input.cov, params.ranks)
        # Write contig taxonomy + summed median fold values to files
        dataframe.drop("lineage", axis=1).to_csv(output.contig_cov, sep="\t")
        lineage_cov.to_csv(output.cov, sep="\t")

def read_h5(filename):
    """
    Reads a hdf file with h5py package and returns a data frame
    """
    with h5py.File(filename, 'r') as hdf_file:
        data_matrix = hdf_file['data'][:]
        sample_names = hdf_file['data'].attrs['sample_names'].astype(str)
    return pd.DataFrame(dict(zip(sample_names, data_matrix)))

rule quantify_taxonomy_GC:
    output:
        cov=results + "/Genecatalog/annotation/mmseqs2/{seqTaxDB}.median_fold.tsv",
        gene_cov = results + "/Genecatalog/annotation/mmseqs2/{seqTaxDB}.gene.median_fold.tsv",
    input:
        tsv=rules.mmseqs2_createtsv_GC.output.tsv,
        cov=results + "/Genecatalog/counts/median_coverage.h5",
        gene_coverage_stats = results + "/Genecatalog/counts/gene_coverage_stats.parquet"
    params:
        ranks = config["taxonomy"]["ranks"]
    run:
        import pandas as pd
        import h5py
        cov = read_h5(input.cov)
        gene_coverage_stats = pd.read_parquet(input.gene_coverage_stats)
        cov.index = [x.split(" ")[0] for x in gene_coverage_stats["#Name"]]
        lineage_cov, dataframe = parse_mmseqs2_GC(cov, input.tsv, params.ranks, [0,4])
        # Write contig taxonomy + summed median fold values to files
        dataframe.drop("lineage", axis=1).to_csv(output.gene_cov, sep="\t")
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
        outdir=lambda wildcards, output: os.path.abspath(".") + "/" + os.path.dirname(output.txt),
        settings="-a diamond --local --clean --input_type protein --debug --include_loose --include_nudge",
        tmpdir="$TMPDIR/genecatalog.rgi",
        faa="$TMPDIR/genecatalog.rgi/input.faa",
    conda:
        "envs/rgi.yml"
    threads: 10
    resources:
        runtime=60,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"]
    shell:
        """
        exec &>{log}
        mkdir -p {params.tmpdir}
        cp {input.db} {params.tmpdir}
        sed 's/*//g' {input.faa} > {params.faa}
        cd {params.tmpdir}
        rgi load --card_json card.json --local
        rgi main -i input.faa -o rgi.out -n {threads} {params.settings}
        mv rgi.out* {params.outdir}
        rm -rf {params.tmpdir}
        """

rule rgi_genomes:
    output:
        json=results+"/genomes/annotations/rgi/rgi.out.json",
        txt=results+"/genomes/annotations/rgi/rgi.out.txt"
    input:
        faa=expand(results+"/genomes/annotations/genes/{genome}.faa", genome = genomes),
        db="resources/card/card.json"
    log:
        results+"/logs/genomes/rgi/genomes.log",
    params:
        outdir=lambda wildcards, output: os.path.abspath(".") + "/" + os.path.dirname(output.txt),
        settings="-a diamond --local --clean --input_type protein --debug --include_nudge --include_loose",
        tmpdir="$TMPDIR/genomes.rgi",
        faa="$TMPDIR/genomes.rgi/input.faa",
    conda:
        "envs/rgi.yml"
    threads: 20
    resources:
        runtime=60,
        mem_mib=mem_allowed,
        slurm_account=lambda wildcards: config["slurm_account"],
        constraint="mem256GB"
    shell:
        """
        exec &>{log}
        mkdir -p {params.tmpdir}
        cp {input.db} {params.tmpdir}
        sed 's/*//g' {input.faa} > {params.faa}
        cd {params.tmpdir}
        rgi load --card_json card.json --local
        rgi main -i input.faa -o rgi.out -n {threads} {params.settings}
        mv rgi.out* {params.outdir}
        rm -rf {params.tmpdir}
        """

rule rgi_parse_genecatalog:
    output:
        tsv=results + "/Genecatalog/annotations/rgi.parsed.tsv"
    input:
        txt=rules.rgi_genecatalog.output.txt
    params:
        best_hit_id=40,
        best_hit_bitscore=100,
        percent_length=50
    run:
        annot = pd.read_csv(input.txt, sep="\t", index_col=0)
        annot = annot.loc[(annot.Cut_Off!="Loose")|((annot["Best_Identities"]>=params.best_hit_id)&(annot["Best_Hit_Bitscore"]>=params.best_hit_bitscore)&(annot["Percentage Length of Reference Sequence"]>=params.percent_length))]
        annot = annot.loc[:,
                ["Model_ID", "Best_Hit_ARO", "AMR Gene Family", "Resistance Mechanism", "Cut_Off"]]
        annot.loc[:, "Model_ID"] = ["RGI_{}".format(x) for x in annot.Model_ID]
        annot.rename(index=lambda x: x.split(" ")[0], inplace=True)
        annot.to_csv(output.tsv, sep="\t", index=True)

rule rgi_parse_genomes:
    output:
        tsv=results + "/genomes/annotations/rgi/rgi.parsed.tsv"
    input:
        txt=rules.rgi_genomes.output.txt
    params:
        best_hit_id=40,
        best_hit_bitscore=100,
        percent_length=50
    run:
        annot = pd.read_csv(input.txt, sep="\t", index_col=0)
        annot = annot.loc[(annot.Cut_Off!="Loose")|((annot["Best_Identities"]>=params.best_hit_id)&(annot["Best_Hit_Bitscore"]>=params.best_hit_bitscore)&(annot["Percentage Length of Reference Sequence"]>=params.percent_length))]
        annot = annot.loc[:,
                ["Model_ID", "Best_Hit_ARO", "AMR Gene Family", "Resistance Mechanism", "Cut_Off"]]
        annot.loc[:, "Model_ID"] = ["RGI_{}".format(x) for x in annot.Model_ID]
        annot.rename(index=lambda x: x.split(" ")[0], inplace=True)
        annot.to_csv(output.tsv, sep="\t", index=True)

def read_h5(filename):
    """
    Reads a hdf file with h5py package and returns a data frame
    """
    import h5py
    with h5py.File(filename, 'r') as hdf_file:
        data_matrix = hdf_file['data'][:]
        sample_names = hdf_file['data'].attrs['sample_names'].astype(str)
    return pd.DataFrame(dict(zip(sample_names, data_matrix)))

rule sum_rgi_genecatalog:
    output:
        rgi_model=results + "/Genecatalog/counts/rgi_model.parsed.tsv",
        rgi_aro=results + "/Genecatalog/counts/rgi_aro.parsed.tsv",
        rgi_family=results + "/Genecatalog/counts/rgi_family.parsed.tsv",
        rgi_model_strict = results + "/Genecatalog/counts/rgi_model_strict.parsed.tsv",
        rgi_aro_strict = results + "/Genecatalog/counts/rgi_aro_strict.parsed.tsv",
        rgi_family_strict = results + "/Genecatalog/counts/rgi_family_strict.parsed.tsv"
    input:
        tsv=rules.rgi_parse_genecatalog.output.tsv,
        h5=results + "/Genecatalog/counts/median_coverage.h5",
        parquet=results + "/Genecatalog/counts/gene_coverage_stats.parquet",
    run:
        import pandas as pd
        import h5py
        median_coverage = read_h5(input.h5)
        gene_coverage_stats = pd.read_parquet(input.parquet)
        median_coverage.index = [x.split(" ")[0] for x in gene_coverage_stats["#Name"]]
        annot = pd.read_csv(input.tsv, sep="\t", index_col=0)
        annot_cov = pd.merge(annot, median_coverage, left_index=True, right_index=True)
        # Sum to AMR model
        rgi_model_info = annot.set_index("Model_ID").groupby(level=0).first().loc[:, ["AMR Gene Family", "Resistance Mechanism"]]
        rgi_model_sum = annot_cov.groupby("Model_ID").sum(numeric_only=True)
        rgi_model_sum = pd.merge(rgi_model_info, rgi_model_sum, left_index=True, right_index=True)
        rgi_model_sum.to_csv(output.rgi_model, sep="\t")
        # Sum to ARO term
        rgi_aro_info = annot.set_index("Best_Hit_ARO").groupby(level=0).first().loc[:, ["AMR Gene Family", "Resistance Mechanism"]]
        rgi_aro_sum = annot_cov.groupby("Best_Hit_ARO").sum(numeric_only=True)
        rgi_aro_sum = pd.merge(rgi_aro_info, rgi_aro_sum, left_index=True, right_index=True)
        rgi_aro_sum.to_csv(output.rgi_aro, sep="\t")
        # Sum to AMR family
        rgi_family_info = annot.set_index("AMR Gene Family").groupby(level=0).first().loc[:, ["Resistance Mechanism"]]
        rgi_family_sum = annot_cov.groupby("AMR Gene Family").sum(numeric_only=True)
        rgi_family_sum = pd.merge(rgi_family_info, rgi_family_sum, left_index=True, right_index=True)
        rgi_family_sum.to_csv(output.rgi_family, sep="\t")
        # Sum to AMR model strict
        rgi_model_strictsum = annot_cov.loc[annot_cov.Cut_Off!="Loose"].groupby("Model_ID").sum(numeric_only=True)
        rgi_model_strictsum = pd.merge(rgi_model_info, rgi_model_strictsum, left_index=True, right_index=True)
        rgi_model_strictsum.to_csv(output.rgi_model_strict, sep="\t")
        # Sum to ARO term strict
        rgi_aro_strictsum = annot_cov.loc[annot_cov.Cut_Off!="Loose"].groupby("Best_Hit_ARO").sum(numeric_only=True)
        rgi_aro_strictsum = pd.merge(rgi_aro_info, rgi_aro_strictsum, left_index=True, right_index=True)
        rgi_aro_strictsum.to_csv(output.rgi_aro_strict, sep="\t")
        # Sum to AMR family strict
        rgi_family_strictsum = annot_cov.loc[annot_cov.Cut_Off!="Loose"].groupby("AMR Gene Family").sum(numeric_only=True)
        rgi_family_strictsum = pd.merge(rgi_family_info, rgi_family_strictsum, left_index=True, right_index=True)
        rgi_family_strictsum.to_csv(output.rgi_family_strict, sep="\t")

        