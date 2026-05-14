#!/usr/bin/env python3

## 26-Nov-2026
##
## track_missing_prots2.py
## modification of track_missing_prots.py that is aware of clusters with short proteins
##
## ignores clusters with short proteins
##
## using the raw clustering data
## (also need mapping from proteome number to GCA):
##
## (1) read all the clusters
## (2) associate clusters with proteomes
## (3) find proteomes with the smallest number of clusters
## (4) for those proteomes, identify clusters/mode proteins that are missing from those proteomes
## (5) write out a tfx_ex format file with missing mode proteins
##

import sys
import fileinput
import argparse
import statistics as stat
from pathlib import Path
import random
import copy
import re
import math

def check_args(test_args=None):
    """
    parse arguments and check for error conditions
    """

    parser = argparse.ArgumentParser(
        description="identify mode proteins missing from proteomes"
    )

    parser.add_argument(
        "-P",
        "--prot_id_map",
        dest="prot_id_map",
        type=str,
        help="mapping of proteome_ids to ENA_accs",
    )

    parser.add_argument(
        "--miss_fract",
        dest="miss_fract",
        type=float,
        default=0.5,
        help="number of proteins per proteome",
    )

    parser.add_argument(
        "--min_len",
        dest="min_len",
        type=int,
        default=100,
        help="number of proteins per proteome",
    )

    parser.add_argument(
        "--miss_max",
        dest="miss_max",
        type=int,
        default=20,
        help="number of proteins per proteome",
    )

    parser.add_argument(
        "--n_samp",
        dest="n_samp",
        default=50,
        type=str,
        help="number of proteins per proteome",
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
            proteome_info['ome_id'] = proteome_id
            proteome_map[proteome_id] = proteome_info['ena_acc']

    return proteome_map

def main():

    ## format of Proteins_cluster_m50pct.tsv
    f_names = "cluster_id protein_ids proteins_count proteomes_count representative seqlen_range".split(' ')
    ## indexes of fields, to avoid building dictionary of data line
    fname_idx = dict(zip(f_names,range(len(f_names))))

    args = check_args()

    clust_dict = {}
    lclust_dict = {}

    ome_clusts = {}
    ome_lclusts = {}

    proteome_map = {}
    if args.prot_id_map:
        proteome_map = read_proteome_map(args.prot_id_map)

## cluster_file header:
## cluster_id	protein_ids	proteins_count	proteomes_count	representative	seqlen_range	seqlen_mode
##

    if (len(args.files) > 0):

        for this_file in args.files:
            this_oscode = this_file.split('/')[-2]
            in_clust_fields = []

            lclust_cnt = 0

            with open(this_file,'r') as in_clust_fd:
                idx_clust_bad = 0
                for line in in_clust_fd:
                    if (line.startswith('#')):
                        continue
                    in_data = line.strip('\n').split('\t')
                    if (in_data[0] == 'cluster_id'):
                        in_clust_fields = [x for x in in_data]
                        continue

                    clust_data = dict(zip(in_clust_fields, in_data))

                    clust_id = clust_data['cluster_id']

                    sseqlen_mode = re.sub(r'\.\d+$','',clust_data['seqlen_mode'])
                    seqlen_mode = int(sseqlen_mode)

                    if (seqlen_mode > args.min_len):
                        lclust_cnt += 1
                    
                    one_ome_set = set()
                    ## parse cluster info
                    for prot_id in clust_data['protein_ids'].split(' '):

                        (this_proteome, this_prot_info) = prot_id.split(':')

                        if (this_proteome == ''):
                            continue

                        prot_info_fields = this_prot_info.split('|')
                        (this_prot_acc, this_prot_len) = prot_info_fields[0:2]
                        this_prot_len = int(this_prot_len)

                        if (clust_id not in clust_dict and this_prot_len == seqlen_mode):
                            clust_dict[clust_id] = {'info':this_prot_info,'len':this_prot_len}

                        if (seqlen_mode > args.min_len):
                            lclust_dict[clust_id] = {'info':this_prot_info,'len':this_prot_len} 

                        if (this_proteome not in one_ome_set):
                            one_ome_set.add(this_proteome)
                            if (this_proteome not in ome_clusts):
                                ome_clusts[this_proteome] = {'ids':[clust_id]}
                            elif (clust_id not in ome_clusts[this_proteome]):
                                ome_clusts[this_proteome]['ids'].append(clust_id)

                            if (seqlen_mode > args.min_len):
                                if (this_proteome not in ome_lclusts):
                                    ome_lclusts[this_proteome] = {'ids':[clust_id]}
                                elif (clust_id not in ome_lclusts[this_proteome]):
                                    ome_lclusts[this_proteome]['ids'].append(clust_id)

        ## here we have all the clusters (with mode example in
        ## clust_dict[], and all the clusters per proteome in ome_clusts[]

        ome_clust_keys = []
        for ome in ome_clusts.keys():
            ome_clusts[ome]['cnt'] = len(ome_clusts[ome]['ids'])
            ome_clust_keys.append(ome)

        ome_lclust_keys = []
        for ome in ome_lclusts.keys():
            ome_lclusts[ome]['cnt'] = len(ome_lclusts[ome]['ids'])
            ome_lclust_keys.append(ome)

        ome_clust_keys.sort(key=lambda x: ome_clusts[x]['cnt'])

    ## now we need to find the clust_id's that are NOT in the low-cluster proteomes

    total_clusts = set(list(clust_dict.keys()))

    total_lclusts = set(list(lclust_dict.keys()))

## need this format:
## >GCA_001025535.1	151	B_N50L:3.73
## 817	CNQ15512|222	222	KLS42055|164|UPI00064C7AB1
## 6309	CNQ59478|447	447	KLS40861|252|UPI00064CAAC2

    clust_cnt = len(clust_dict.keys())

    ## we are only interested in clusters with mode_len > 100

    lclust_cnt = len(lclust_dict.keys())

    max_miss_thresh = lclust_cnt*args.miss_fract

    missing_omes = []
    for ome in ome_lclust_keys:
        if (ome_lclusts[ome]['cnt'] <= max_miss_thresh):
            missing_omes.append(ome)

    if (len(missing_omes) > args.miss_max):
        missing_omes_s = random.sample(missing_omes,args.miss_max)
    else:
        missing_omes_s = missing_omes

    for ome in missing_omes_s:

        gca_id = proteome_map[ome]

        clust_len = ome_lclusts[ome]['cnt']

        print(f'>{gca_id}\t{clust_len}\tMISS')
        miss_clusts = total_lclusts.difference(ome_lclusts[ome]['ids'])

        miss_clusts_list = list(miss_clusts)
        if (len(miss_clusts_list)>args.n_samp):
            miss_samp = random.sample(miss_clusts_list,50)
        else:
            miss_samp = miss_clusts_list

        for miss in miss_samp:
            mode_prot = clust_dict[miss]['info']
            mode_len =  clust_dict[miss]['len']
            print(f'{miss}\t{mode_prot}\t{mode_len}\t{mode_prot}')

if __name__ == "__main__":

    main()
