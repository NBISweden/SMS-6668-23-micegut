# National Bioinformatics Support

_Support #6668 Shotgun metagenomic sequencing - Three-generations microbiome study_

## Installation

A conda environment file with the software used to run the Atlas workflow is provided as `atlas-env.yml`. To create the environment run:

```bash
conda env create -f atlas-env.yml
```

Then activate the environment with:

```bash
conda activate atlas
```

Patching the atlas workflow installation:

1. The version of the GUNC software used in atlas v2.15.0 was found to be incompatible with version 2 of the python package `pandas`. We therefore manually set the pandas version to `1.5.3` in the `gunc.yaml` file of the Atlas workflow. (Locate this file with `find $CONDA_PREFIX -name gunc.yaml`).
2. In order to get the DRAM software to work with our setup we had to replace the `CM.pm` perl module for the `tRNAscan-SE` program with a custom module obtained from the tRNAscan-SE developers (see GitHub issue [#23](https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/23)). The custom module is found in the GitHub repository linked above under `src/CM.pm`. The non-working module was replaced by activating the specific conda environment containing tRNAscan-SE then running `cp src/CM.pm $CONDA_PREFIX/lib/tRNAscan-SE/tRNAscanSE/CM.pm`.
3. We also found that the DRAM software used by Atlas required a configuration file to be passed to the `DRAM.py annotate` and `DRAM.py distill` with `--config_loc`. This configuration file is also found in the GitHub repository linked above under `resources/DRAM/DRAM.config`. A modified version of the `dram.smk` rules file is availble under `src/dram.smk` and can be used to replace the original file in the Atlas workflow by running `cp src/dram.smk $(find $CONDA_PREFIX -name "dram.smk")`
4. Because atlas version 2.15.0 does not include any filtering of mapped reads we modified the atlas source code to allow passing `minmapq` to the `pileup.sh` script used when calculating coverages. These changes are tracked at [https://github.com/johnne/atlas](https://github.com/johnne/atlas). Following these modifications we set `minimum_map_quality: 10` in the atlas config. In order to use the `minmapq` parameter, replace the following files `genecatalog.smk`, `assemble.smk` and `genomes.smk` from the `src/` directory in this repository with the corresponding files in the atlas workflow installation. The files can be located with `find $CONDA_PREFIX -name "genecatalog.smk"` etc.


## Running Jupyter notebook with R kernel

Several of the R-packages used require a Linux kernel with the AMD 
architecture. This repo contains a Dockerfile at
[envs/Dockerfile](envs/Dockerfile) which can be used to create a Linux based 
Docker image with a conda environment that allows you to run several 
R-packages for statistical analyses in a Jupyter notebook.

To create the image run:

```bash
docker build -f envs/Dockerfile -t stats_linux envs/
```

Then to start a Jupyter server run:

```bash
docker run -p 8888:8888 --platform linux/x86_64 -v $(pwd):/analysis stats_linux
```

This connects port 8888 on the host and inside the container, explicitly 
runs the container on the Linux platform and mounts the current directory on 
the host to the `/analysis` directory inside the container. 

Once the container starts, go to `localhost:8888` you will see something 
similar to this in your terminal:

```bash
[I 06:39:10.737 NotebookApp] Writing notebook server cookie secret to /root/.local/share/jupyter/runtime/notebook_cookie_secret
[I 06:39:11.510 NotebookApp] Serving notebooks from local directory: /analysis
[I 06:39:11.511 NotebookApp] Jupyter Notebook 6.5.4 is running at:
[I 06:39:11.511 NotebookApp] http://f6b8639983e7:8888/?token=19c9352d9291fdfeaa045f6364773711
[I 06:39:11.511 NotebookApp]  or http://127.0.0.1:8888/?token=19c9352d9291fdfeaa045f6364773711
[I 06:39:11.511 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 06:39:11.525 NotebookApp] 
    
    To access the notebook, open this file in a browser:
        file:///root/.local/share/jupyter/runtime/nbserver-8-open.html
    Or copy and paste one of these URLs:
        http://f6b8639983e7:8888/?token=19c9352d9291fdfeaa045f6364773711
     or http://127.0.0.1:8888/?token=19c9352d9291fdfeaa045f6364773711
```

Now in your browser, go to `localhost:8888` and in the password field enter the 
token you saw printed to the terminal when you started the container (in 
this case it would be `19c9352d9291fdfeaa045f6364773711`).