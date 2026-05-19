#!/usr/bin/env bash
source env.sh
mkdir -p outdir
$ETT_BIN/bs_bad_clusters_oscode_ins2_cm0.sh TEST
