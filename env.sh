#!/usr/bin/env bash
set -u
set -e
echo "Setting up environment variables"
export ETT_BIN=$(pwd)/analysis_code/
export DATA_DIR=$(pwd)/data
export BS_ETT_DIR=$(pwd)/outdir
export DATA_DATE=2026-05-18 #replace this with the date of your analysis
