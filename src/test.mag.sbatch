#!/bin/bash -l
#SBATCH -A snic2022-5-350
#SBATCH -p node
#SBATCH -t 24:00:00
#SBATCH -J mag.test.node


mamba activate mice-gut
basedir="/proj/snic2020-5-486/nobackup/SMS-23-6668-micegut"
workdir="work.test"
outdir="mag.test"
mkdir -p $TMPDIR/$workdir
nextflow run nf-core/mag -r 2.3.0 -profile test,singularity -resume -work-dir $TMPDIR/$workdir --max_cpus 20 --outdir $basedir/$outdir
rsync -r $TMPDIR/$workdir $basedir/
