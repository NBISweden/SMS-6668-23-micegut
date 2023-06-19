# National Bioinformatics Support

_Support #6668 Shotgun metagenomic sequencing - Three-generations microbiome study_

- Project report: [doc/project-report.ipynb](doc/project-report.ipynb).

## Running Jupyter notebook with R kernel

Several of the R-packages used require a Linux kernel with the AMD 
architecture. This repo contains a Dockerfile at
[doc/Dockerfile](doc/Dockerfile) which can be used to create a Linux based 
Docker image with a conda environment that allows you to run several 
R-packages for statistical analyses in a Jupyter notebook.

To create the image run:

```bash
docker build -f doc/Dockerfile -t stats_linux doc/
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