# National Bioinformatics Support - Support #6668 Shotgun metagenomic sequencing - Three-generations microbiome study

## Setup

```bash
mkdir -p envs
mamba env create -f environment.yml -p envs/mice-gut
mamba activate envs/mice-gut
```

Make sample list

```bash
python src/make_sample_list.py > data/sample_list.csv
```

### Setup nf-core/mag pipeline

```bash
module purge
module load uppmax bioinfo-tools Nextflow/22.10.1
NXF_HOME="/proj/snic2020-5-486/nobackup/SMS-23-6668-micegut/.nextflow"
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_SINGULARITY_CACHEDIR="/sw/data/ToolBox/nf-core/"
export NXF_SINGULARITY_LIBRARYDIR="/sw/bioinfo/nf-core-pipelines/latest/rackham/singularity_cache_dir"
nextflow pull nf-core/mag -r 2.3.0
```

### Run pipeline

```bash
nextflow -c conf/custom.config run -params-file conf/mag.config.yml nf-core/mag -r 2.3.0 -resume -profile uppmax --project snic2022-5-350 
```
