#!/bin/sh

mode_prot=$1
bad_prot=$2
gca_ome=$3
taxon=$4

export BIN_DIR=/homes/pearson/ett

ssearch36 -E '1e-6 -1' -XG -m8CBl -s MD10 $mode_prot $bad_prot >> ${taxon}_${gca_ome}_short.ss_MD10

