#!/bin/sh

n=$1
curl --silent "https://www.ebi.ac.uk/ena/browser/api/fasta/${n}?download=true&gzip=false"

