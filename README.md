# ett_ms
Early terminated trasnscripts and missing proteins reflect artifacts in bacterial proteomes (code for the manuscript)

# Abstract

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
