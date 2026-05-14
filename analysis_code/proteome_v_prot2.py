#!/usr/bin/env python3

## script to look for full-length proteins in short-protein proteomes

## genome_v_prot.py
## take a file of the form
## >GCA_013423585.1	664	B   <- this is the proteome with short proteins
## 1760960	KOZ39023|560	560 <- mode protein for the cluster that has a short protein
## 1220219	KOZ39495|586	586 <- 
## 589694	KOZ37315|479	479
## 1500335	KOZ38449|329	329
## 687435	KOZ37629|231	231
##
## (1) extract the proteome_acc, count of bad clusters, B/S
## (2) get the list of protein acc's to search
## (3) build a string for down_fasta_ecoli.sh -- send to GCA_nnnnn.fa
## (4) run tfastx with GCA_nnn.fa \!down_1gca.sh+GCA_nnnn > GCA_nnnn_B/S.tfx
##

import fileinput
import sys
import os
import subprocess
import argparse

BIN_DIR="/homes/pearson/ett"

def check_args(test_args=None):
    """
    parse arguments and check for error conditions
    """

    parser = argparse.ArgumentParser(
        description="parse summ2 file, extract cluster queries, run genome searches"
    )

    parser.add_argument(
        "--oscode","-T",
        dest='taxon',
        type=str,
        default='ecoli'
        )

    parser.add_argument(
        "--fa_dir",
        dest='fa_dir',
        type=str,
        default='fasta10_incsurv'
        )

    parser.add_argument(
        "--script",
        dest='down_script',
        type=str,
        default='down_fasta_oscode.sh'
        )

    parser.add_argument(
        "-P",
        "--prot_id_map",
        dest="prot_id_map",
        type=str,
        help="mapping of proteome_ids to ENA_accs",
    )

    parser.add_argument(
        "--suff","-S",
        dest='suff',
        type=str,
        default='tfxg'
        )

    parser.add_argument(
        dest="files",
        type=str,
        nargs='*')

    return parser.parse_args()

## reads proteome mapping/BUSCO/metadata info
def read_proteome_map(map_file):
    ## proteome_map{} dict maps proteome_id's to GCA accessions, but
    ## also provides BUSCO and metadata if available

    proteome_map = {}
    proteome_idx = 0

    ome_fields = "proteome_id up_acc ena_acc".split(" ")

    with open(map_file,'r') as pid_f:
        for line in pid_f:
            proteome_data = line.strip('\n').split('\t')
            proteome_id = proteome_data[proteome_idx]
            proteome_info = dict(zip(ome_fields,proteome_data[:3]))
            proteome_map[proteome_info['ena_acc']] = proteome_id

    return proteome_map

def run_query(gca_acc, prot_ids, fa_dir, taxon, s_type, proteome_map, args):

    # (1) get the sequences
    if (len(prot_ids) > 0):
        query_name = gca_acc

        N_prot_ids = len(prot_ids)

        q_acc_str = "' '".join(prot_ids)

        q_file_name = f'{query_name}.aa'

        if not os.path.exists(q_file_name):
            cmd_str = f'{BIN_DIR}/{args.down_script} {fa_dir} {taxon} \'{q_acc_str}\' > {q_file_name}'
            ## print(cmd_str)
            os.system(cmd_str)

        prot_id = proteome_map[gca_acc]

        subprocess.run([f'{BIN_DIR}/protN_proteome_fa.sh',q_file_name,taxon,prot_id, gca_acc,args.suff])

def main():

    args = check_args()

    proteome_map = {}
    if args.prot_id_map:
        proteome_map = read_proteome_map(args.prot_id_map)

    prot_ids = []
    gca_acc = ''
    for line in fileinput.input(files=args.files):
        if line.startswith('>'):
            ## if have previous set, generate files, do search
            if (gca_acc):
                run_query(gca_acc, prot_ids,args.fa_dir, args.taxon, s_type, proteome_map, args)

            (gca_acc, n_bad_clust, s_type) = line[1:].strip('\n').split('\t')
            prot_ids = []
            continue

        clust_data = line.strip('\n').split('\t')
        if (len(clust_data)==4):
            (clust_id, prot_id, prot_len, mode_prot_id) = clust_data
        else:
            (clust_id, prot_id, prot_len) = clust_data

        prot_ids.append(prot_id)
        
    if (gca_acc):
        run_query(gca_acc, prot_ids, args.fa_dir, args.taxon, s_type, proteome_map, args)

if __name__ == "__main__":

    main()

