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
bash src/init_mag.sh
```

### Run pipeline

```bash
nextflow -c conf/custom.config run -params-file conf/mag.config.yml nf-core/mag -r 2.3.0 -resume -profile uppmax --project snic2022-5-350 
```
