#!/bin/sh

query=$1
s_type=$2
taxon=$3
suff=$4

f=${query%.*}
${ETT_BIN}/down_1gca.sh $f > $f.nt
tfastx36 -E '1e-6 -1' -XG -m8CBl -s MD10 $query $f.nt  > ${taxon}_${f}_${s_type}.${suff}
rm -f $f.nt
