#!/usr/bin/env bash
set -u
set -e
echo "Downloading data archive from Zenodo..."
wget -O data.zip https://zenodo.org/api/records/20208872/files-archive
echo "Uncompressing..."
unzip data.zip
rm data.zip
mkdir -p data
mv yields.tsv *.zip data/  
(cd data && for i in *.zip; do unzip $i; done)
find data -type f -name "*.gz" -exec gunzip {} \;
rm data/*.zip
mkdir -p fasta
./env.sh
echo "Setup Completed"
