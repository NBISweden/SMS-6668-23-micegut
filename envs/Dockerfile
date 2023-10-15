FROM --platform=linux/x86_64 condaforge/mambaforge

ARG ENVNAME="r-statistics"
ARG ENVFILE="r-statistics.yml"

LABEL description = "Image for running data exploration with jupyter on a linux kernel"
MAINTAINER "John Sundh" john.sundh@scilifelab.se

# Use bash as shell
SHELL ["/bin/bash", "--login", "-c"]

# Set workdir
WORKDIR /analysis

# Set timezone
ENV TZ="Europe/Stockholm"
ENV DEBIAN_FRONTEND=noninteractive

# Install package for setting timezone
RUN apt-get update && apt-get install -y tzdata gcc && apt-get clean

# Configure Conda/Mamba
RUN mamba init bash && conda config --set channel_priority strict && \
    conda config --append channels bioconda && \
    conda config --append channels r && \
    conda config --set subdir linux-64

# Install environment
COPY $ENVFILE ./
RUN mamba env create -f $ENVFILE -n $ENVNAME && \
    mamba clean -a

# Set mrsa-workflow environment as active at start-up
RUN echo "source activate $ENVNAME" >> ~/.bashrc

# Add environment to PATH
ENV PATH /opt/conda/envs/$ENVNAME/bin:${PATH}

# Open up port 8888
EXPOSE 8888

# Run notebook
CMD jupyter notebook --allow-root --ip 0.0.0.0 --no-browser
