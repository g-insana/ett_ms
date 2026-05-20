# ett_ms
[![DOI](https://zenodo.org/badge/1224686730.svg)](https://doi.org/10.5281/zenodo.20313420)

This repository contains the code for the analysis and figures of the manuscript:
# Early terminated transcripts and missing proteins reflect artifacts in bacterial proteomes
[(doi 10.64898/2026.05.19.725897)](https://doi.org/10.64898/2026.05.19.725897)

## Abstract

MMMseqs2 clustering was used to examine the uniformity of proteomes from 20 bacterial species. Clusters with proteins from ≥50% of proteomes typically contain proteins from 95% of the proteomes and capture more than 80% of the proteins in an organism.

Protein clusters are highly uniform in length; across the 20 bacteria, the median cluster has more than 99% of the proteins at the mode length. In contrast to this uniformity, some clusters contain dozens to hundreds of proteins that are considerably shorter (<75%) than the mode-length, and a few clusters include proteins that are >133% the mode length.

Most “outlier” proteins are found in fewer than 10% of clusters, and “high-outlier” clusters are over-represented in a small fraction of proteomes, that often have poor Proteome BUSCO fragment scores.

Short-outlier proteins are artifacts; at least 80% of short-outlier genomes contain mode-length copies of the protein, which were missed because of frame-shifts, termination codons, or initiation codon choice.

As with “short-outlier” proteins, the ∼5% of proteomes missing from the core (50% participation) cluster set encode the missing protein more than 98% of the time. MMseqs2 clustering with 50% participation provides robust sets of core bacterial proteins.

## Contents of the repository
```
analysis_code/      # bash and python scripts to perform the analysis
data/               # location for data files and test data
figures_code/       # datafiles and .R code to recreate the figures in the manuscript
env.sh              # script to setup environment variables to run the analysis
fetch_datasets.sh   # script to download and uncompress the datasets
testrun.sh          # script to run a mini test analysis to ensure requirements are met and code is functional
```

## Requirements

The python scripts under `analysis_code/` do not require any module which is not already included in the standard library,
but you will need to have the binaries from the [FASTA36](https://github.com/wrpearson/fasta36/) suite to run protein vs dna and protein vs protein searches.
The binaries can be downloaded from the [releases tab of the FASTA36 github](https://github.com/wrpearson/fasta36/releases)
or from the [University of Virginia](https://fasta.bioch.virginia.edu/wrpearson/fasta/fasta36).

The R scripts under `figures_code` require the libraries detailed in the [figures code README](figures_code/README.md).

The input data (processed clusters) is available at [Zenodo](https://doi.org/10.5281/zenodo.20208872) or [Figshare](https://doi.org/10.6084/m9.figshare.32301477).
It was produced by [MMseqs2](https://github.com/soedinglab/MMseqs2) via the [ProteomeCluster pipeline](https://github.com/g-insana/ProteomeCluster) [v1.0.0](https://doi.org/10.5281/zenodo.20208647) and
needs to be downloaded and unpacked. A script in the main directory (`fetch_datasets.sh`) can be used for that purpose.

Proteome files in FASTA format also need to be present, if you'd like to re-run the sequence searches. These can be downloaded to `proteomes/` named subfolders under each `OSCODE/` directory, using either the `upid` or the `gca_set_acc` information present in the `OSCODE/OSCODE.proteomes.tsv` files.
For example, to download the sequence file for the proteome with upid `UP000434630` and save it as `proteome_4348828.fa`, the following command could be used:
```
wget -O data/SHIFL/proteomes/proteome_4348828.fa.gz "https://rest.uniprot.org/uniparc/proteome/UP000434630/stream?compressed=true&format=fasta"
gunzip data/SHIFL/proteomes/proteome_4348828.fa.gz
```
A few proteomes have been provided under `data/TEST/proteomes/` to run the test.

## DOCUMENTATION

- [clustering](https://github.com/g-insana/ProteomeCluster/blob/main/README.md)
- [analysis](analysis_code/README.md)
- [figures](figures_code/README.md)


## LINKS

- [bioRxiv preprint](https://doi.org/10.64898/2026.05.19.725897)
- [MMseqs2](https://github.com/soedinglab/MMseqs2) search and clustering suite
- [FASTA36](https://github.com/wrpearson/fasta36) sequence search and comparison software
- [ProteomeCluster](https://github.com/g-insana/ProteomeCluster) clustering pipeline
- input datasets at [Zenodo](https://doi.org/10.5281/zenodo.20208872) or [Figshare](https://doi.org/10.6084/m9.figshare.32301477)

## CITATION

If you find this software useful, please consider citing our [paper](https://doi.org/10.64898/2026.05.19.725897):

``` 
Insana, G., Martin, M.J. & Pearson, W.R.
Early terminated transcripts and missing proteins reflect artifacts in bacterial proteomes
BioRxiv (2026). https://doi.org/10.64898/2026.05.19.725897
```

Bibtex:
```
@article {Insana2026.05.19.725897,
	author = {Insana, Giuseppe and Martin, Maria J. and Pearson, William R.},
	title = {Early terminated transcripts and missing proteins reflect artifacts in bacterial proteomes},
	elocation-id = {2026.05.19.725897},
	year = {2026},
	doi = {10.64898/2026.05.19.725897},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {MMseqs2 clustering was used to examine the uniformity of proteomes from 20 bacterial species. Clusters with proteins from >=\50% of proteomes typically contain proteins from 95\% of the proteomes and capture more than 80\% of the proteins in an organism. Protein clusters are highly uniform in length; across the 20 bacteria, the median cluster has more than 99\% of the proteins at the mode length. In contrast to this uniformity, some clusters contain dozens to hundreds of proteins that are considerably shorter (<75\%) than the mode-length, and a few clusters include proteins that are >133\% the mode length. Most "outlier" proteins are found in fewer than 10\% of clusters, and "high-outlier" clusters are over-represented in a small fraction of proteomes, that often have poor Proteome BUSCO fragment scores. Short-outlier proteins are artifacts; at least 80\% of short-outlier genomes contain mode-length copies of the protein, which were missed because of frame-shifts, termination codons, or initiation codon choice. As with "short-outlier" proteins, the ∼5\% of proteomes missing from the core (50\% participation) cluster set encode the missing protein more than 98\% of the time. MMseqs2 clustering with 50\% participation provides robust sets of core bacterial proteins. Competing Interest StatementThe authors have declared no competing interest. European Molecular Biology Laboratory core funds},
	URL = {https://www.biorxiv.org/content/early/2026/05/19/2026.05.19.725897},
	eprint = {https://www.biorxiv.org/content/early/2026/05/19/2026.05.19.725897.full.pdf},
	journal = {bioRxiv}
}
```


