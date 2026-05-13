#!/bin/sh

DIR=$1
shift
OSCODE=$1
shift
acc=$*
## echo $OSCODE $acc

/hps/software/users/martin/uniprot/src/anaconda3/envs/python385/bin/extractor.py -f/hps/nobackup/martin/uniprot/production/temp/mmseqs/$DIR/$OSCODE.FASTA -i /hps/nobackup/martin/uniprot/production/temp/mmseqs/$DIR/$OSCODE.FASTA.idx2 -s $acc 2> down_fasta.err

