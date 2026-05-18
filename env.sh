#!/usr/bin/env bash
set -u
set -e
echo "Setting up environment variables"
export ETT_BIN=$(pwd)/analysis_code/
export DATA_DIR=$(pwd)/data
export FA_DIR=$(pwd)/fasta
