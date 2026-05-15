
## ett_ms/figures_code

This directory provides all the 'R'-scripts that were used to make the figures for the paper:
Early terminated trasnscripts and missing proteins reflect artifacts in bacterial proteomes.

This file outlines the mapping between the files produced by the
various cluster analysis scripts, and the figures in the papers.  The
scripts used to summarize the cluster data are provided in the `./analysis_code` directory.

## figures_code files:

This directory contains three types of files:

1. The script that downloads the analyzed data files used to produce the plots: `get_results_ins_dir_date_oscode.sh` and the script to produce the plots once the data has been downloaded: `plot_figs_out9_cm0.sh`.

Once the data sumaries have been downloaded, into a `results` directory, running the `../plot_figs_out9_cm0.sh` script will produce the figures and tables.

2. The 'R' scripts used to produce the figures and tables.

3. A set of `yaml` files that specify the files used to produce the plots in the `yaml/` directory.

## 'R' libraries required:

All the plots (figures) and tables were created using R: v. 4.5.1
(2025-06-13) -- "Great Square Root".  In addition to the base 'R' distribution, the following libraries are used:
```
library('dplyr')
library('forcats')
library('getopt')
library('ggExtra')
library('ggplot2')
library('ggtext')
library('optparse')
library('patchwork')
library('purrr')
library('RColorBrewer')
library('scales')
library('stringr')
library('tidyr')
library('xtable')
library('yaml')
```

In addition, many of the scripts use functions from `get_yaml_opts.R`
and `read_omes.R` to read options from the `yaml/` directory and read
the `OSCODE.bad_clust` and `OSCODE.bad_omes.samp.am` files.

In addition, most of the scripts read the `yields.tsv` file, which
provides statistics on the number of proteomes and clusters from each
bacteria, and the `oscode_names.tsv` file, which maps `OSCODE`s to
bacteria names.

## Recreating the figures:

1. The `get_results_ins_dir_date_oscode.sh` script downloads analyzed
data from the EBI cluster.  To produce the massaged/summarized data
files used for plotting, The raw data files on the zenodo repository
must be downloaded and analyzed with the scripts in `../analysis_code`
(a shell script that does the analyses is provided).  Once the
`../analysis_code/` scripts have been run (mostly to do `tfastx36`
comparisons of mode-length proteins against genomes that produce
short-outlier proteins or are annotated to be missing proteins
altogether), the results files can be downloaded to a `results/`
directory by running `../get_results_ins_dir_date_oscode.sh` from that
directory (the scripts assume that the `results/` directory is
contained in the `figures_code/` directory).

2. Once the data summaries have been downloaded, figures can be drawn
using `plot_figs_out9_cm0.sh`.  The plotting script is designed to be
run from a `results/` directory with in the `figures_code/` directory.

The `.R` plotting scripts have can take optional arguments from the
command line, or from a `.yaml` file that is specified by the `-Y
file.yaml` option.  Command line arguments (such as the name of the
`figure.pdf` file, `--pdf figure.pdf`) override options specified in
the `.yaml` file.  All plotting scripts have a `--pub` option (off by
default).  If `--pub` is not specified, plots contain the command used
to produce the plot at the bottom of the page.

