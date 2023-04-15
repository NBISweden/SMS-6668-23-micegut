#!/usr/bin/env python

import pandas as pd
import glob
import sys
from collections import defaultdict

### nf-core/mag example sample list
"""
sample,group,short_reads_1,short_reads_2,long_reads
sample1,0,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz,data/sample1.fastq.gz
sample2,0,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,data/sample2.fastq.gz
sample3,1,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,
"""

failed=["H19","C33"]
basedir="/proj/snic2020-5-486/dbp_gut_microbiome/DataDelivery_2023-01-10_15-17-23_ngisthlm00104/files/P27457"
sample_info=f"{basedir}/00-Reports/O.Karlsson_22_02_sample_info.txt"
sample_groups="/proj/snic2020-5-486/nobackup/SMS-23-6668-micegut/data/sample_groups.csv"
dfgroups = pd.read_csv(sample_groups, index_col=0).to_dict(orient="index")
df = pd.read_csv(sample_info, sep="\t", index_col=1)
sample_dict = {}
sample_num = defaultdict(lambda: 0)
for row in df.iterrows():
	# Fetch assembly group from groupdict
	sample = row[0].split(" ")[0]
	# skip failed samples
	if sample in failed:
		continue
	try:
		group = dfgroups[sample]["group"]
	except KeyError:
		group = "unknown"
	# Replace space with underscore in sample name
	sample = row[0].replace(" ", "_")
	# Add suffix to sample (makes e.g. 'm.c' samples unique)
	suffix = ""
	if sample in sample_num.keys() and sample_num[sample] > 0:
		suffix = f"_{sample_num[sample]+1}"
	sample_num[sample]+=1
	ngi_id = row[1]["NGI ID"]
	R1 = glob.glob(f"{basedir}/{ngi_id}/**/*R1*.fastq.gz", recursive=True)
	R2 = glob.glob(f"{basedir}/{ngi_id}/**/*R2*.fastq.gz", recursive=True)
	if len(R1) > 1 or len(R2) > 1:
		sys.stderr.write(f"ERROR: Multiple fastq files for sample {sample}\n")
		sys.exit(f"{','.join(R1)}\n{','.join(R2)}\n")
	sample_dict[f"{sample}{suffix}"] = {
		"group": group,
		"short_reads_1": R1[0],
		"short_reads_2": R2[0],
		"long_reads": "",
	}

sample_list = pd.DataFrame(sample_dict).T
sample_list.index.name = "sample"
with sys.stdout as fhout:
	sample_list.to_csv(fhout, sep=",")
