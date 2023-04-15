# National Bioinformatics Support - Support #6668 Shotgun metagenomic sequencing - Three-generations microbiome study

## Installation

```bash
mamba env create -f environment.yml
mamba activate envs/mice-gut
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d $CONDA_PREFIX/etc/conda/activate.d
cp src/activate.d/env_vars.sh $CONDA_PREFIX/etc/conda/activate.d/
cp src/deactivate.d/env_vars.sh $CONDA_PREFIX/etc/conda/deactivate.d/
```

## Setup

Make sample list

```bash
python src/make_sample_list.py > data/sample_list.csv
```

### Setup nf-core/mag pipeline

```bash
nextflow pull nf-core/mag -r 2.3.0
```

### Run pipeline

```bash
nextflow -c conf/custom.config run -params-file conf/mag.config.yml nf-core/mag -r 2.3.0 -resume -profile uppmax --project snic2022-5-350 
```
