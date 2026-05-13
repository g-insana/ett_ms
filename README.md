# ett_ms
## Early terminated transcripts and missing proteins reflect artifacts in bacterial proteomes
(code for the analysis and figures of the manuscript)

# Abstract

MMseqs2 clustering was used to examine the uniformity and heterogeneity of proteomes from 20 bacterial species. Using clustering parameters
that required 50% sequence overlap, clusters with proteins from ≥50% of proteomes typically contain proteins from 95% of the proteomes and capture
more than 80% of the proteins in an organism.

Protein clusters are highly uniform in length; across the 20 bacteria, the median cluster has more than 99% of the proteins at the mode length. While protein lengths in clusters are highly uniform, some clusters contain dozens to hundreds of proteins that are considerably shorter (<75%) than the mode length, and a few clusters include proteins that are >133% the mode length.

Most “outlier” proteins are found in fewer than 10% of clusters, and “high-outlier” clusters are over-represented in a small
fraction of proteomes. Short-outlier proteins are artifacts; at least 80% of short-outlier genomes contain mode-length copies of the protein in the
cluster; 40% of short protein artifacts are produced by sequencing errors (frameshifts and termination
codons) while another 40% by initiation codon choice.

High “outlier” clusters are concentrated in a small fraction of proteomes, which often have poor Proteome BUSCO fragment scores. As with “short-outlier” proteins,
the ∼5% of proteomes that are excluded from the core (50% participation) cluster set encode the missing protein more than 98% of the time; these proteins were missed because
of frameshifts in the genome sequence. MMseqs2 clustering with 50% participation provides robust sets of core bacterial proteins.

# Contents of the repository
```
analysis_code/      # bash and python scripts to perform the analysis
figures_code/       # datafiles and .R code to recreate the figures in the manuscript
```

## INSTALLATION
- git clone the repository: 

```git clone https://github.com/g-insana/ett_ms.git``` 

- install requirements (virtual environment is optional but recommended) via pip or conda/mamba:

via pip:
```
cd ett_ms && python3 -m venv venv_ett
source venv_ett/bin/activate
pip3 install -r requirements.txt
```

via conda or mamba:
```
cd ett_ms && mamba create --name ett --file requirements.txt --channel conda-forge
mamba activate ett
```

**Note** that you also need to have the [`tfastx`](https://github.com/wrpearson/fasta36/) program protein vs dna searches: 
