#!/bin/bash


# # conda remove --name HIERVAE --all
# conda create -n HIERVAE python=3.6
# echo "y\n"
# source "~/anaconda3/etc/profile.d/conda.sh"
# conda activate HIERVAE

conda install -c rdkit rdkit=2020
# conda install -c rdkit rdkit=2019.03
conda install -c rdkit rdkit=2018
conda install -c rdkit rdkit=2017

# conda install pytorch==1.10.1 torchvision==0.11.2 torchaudio==0.10.1 cudatoolkit=10.2 -c pytorch

python get_vocab.py --ncpu 16 < ../zdata/chembl/all.txt > vocab.txt


conda activate
conda remove --name HIERVAE --all
conda create -n HIERVAE python=3.6
conda activate HIERVAE