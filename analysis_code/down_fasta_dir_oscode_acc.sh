#!/bin/sh

DIR=$1
shift
OSCODE=$1
shift
acc=$*
## echo $OSCODE $acc

extractor.py -f $DIR/$OSCODE.FASTA -i $DIR/$OSCODE.FASTA.idx -s $acc 2> down_fasta.err
