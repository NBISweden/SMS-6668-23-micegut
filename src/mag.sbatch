#!/bin/bash -l
#SBATCH -A snic2022-5-350
#SBATCH -p node
#SBATCH -t 240:00:00
#SBATCH -J mag.node


mamba activate mice-gut
basedir="/proj/snic2020-5-486/nobackup/SMS-23-6668-micegut"
params="$basedir/conf/mag.config.yml"
workdir="work"
outdir="mag"
nextflow run -c conf/custom.config -w $TMPDIR/$workdir nf-core/mag -r 2.3.0 -resume -profile uppmax -params-file $params --max_cpus 20