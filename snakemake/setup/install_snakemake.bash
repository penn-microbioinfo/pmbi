#!/bin/bash

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

source "/home/ubuntu/miniconda3/etc/profile.d/conda.sh"

conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake

