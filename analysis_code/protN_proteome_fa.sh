#!/bin/sh

query=$1
oscode=$2
prot_id=$3
gca_id=$4
suff=$5

fasta36 -E '1e-6 -1' -XG -m8CBl -s MD10 $query $DATA_DIR/$oscode/proteomes/proteome_${prot_id}.fa > ${oscode}_${gca_id}_v_prot.${suff}
