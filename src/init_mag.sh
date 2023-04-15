#!/bin/bash -c
module load bioinfo-tools Nextflow/22.10.1
NXF_HOME="/proj/snic2020-5-486/nobackup/SMS-23-6668-micegut/.nextflow"
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_SINGULARITY_CACHEDIR="/sw/data/ToolBox/nf-core/"
export NXF_SINGULARITY_LIBRARYDIR="/sw/bioinfo/nf-core-pipelines/latest/rackham/singularity_cache_dir"
nextflow pull nf-core/mag -r 2.3.0
