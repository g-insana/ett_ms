#!/bin/sh

export DATA_DIR=outdir9
export FA_DIR=fasta9_incsurv
export DATA_DATE=2026-05-13

$P_SBATCH_8 bs_bad_clusters_oscode_ins2_cm0.sh $*



