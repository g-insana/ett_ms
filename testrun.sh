#!/usr/bin/env bash
source env.sh
export FA_DIR=/hps/nobackup/martin/uniprot/production/temp/mmseqs/fasta9_incsurv/
$ETT_BIN/bs_bad_clusters_oscode_ins2_cm0.sh TEST
