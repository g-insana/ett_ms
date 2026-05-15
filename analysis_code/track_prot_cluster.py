#!/usr/bin/env python3

## 26-Nov-2026
##
## track_prot_cluster.py
##
## using three files, bad_clust, bad_omes, and the raw clustering data
## (also need mapping from proteome number to GCA):
##
## (1) identify 10 bad clusters with >60% mode
## (2) read the 2000 sampled proteomes 
## (3) for those 10 clusters, print out the
##     cluster, proteome, protein length, pct_mode for that proteome
##

import sys
import fileinput
import argparse
import statistics as stat
import random
import copy
import re
import math

def check_args(test_args=None):
    """
    parse arguments and check for error conditions
    """

    parser = argparse.ArgumentParser(
        description="summarize proteomes with bad protein lengths"
    )

    parser.add_argument(
        "--mode_pct",
        dest="mode_pct",
        type=float,
        default = 60,
        help="threshold for good cluster mode_pct",
    )

    parser.add_argument(
        "--in_clust_file",
        dest='in_clust_file',
        type=str,
        help="file with bad clusters"
        )

    parser.add_argument(
        "--in_samp_omes",
        dest='in_samp_omes',
        type=str,
        help="file with sampled bad omes"
        )

    parser.add_argument(
        "--bad_file",
        dest="bad_file",
        type=str,
        help="file for bad proteome/protein",
    )

    parser.add_argument(
        "-P",
        "--prot_id_map",
        dest="prot_id_map",
        type=str,
        help="mapping of proteome_ids to ENA_accs",
    )

    parser.add_argument(
        "-X",
        "--exclude",
        dest="exclude",
        type=str,
        default="B",
        help="mapping of proteome_ids to ENA_accs",
    )

    parser.add_argument(
        dest="files",
        type=str,
        nargs='*')

    return parser.parse_args()

def ome_id_gca(ome_id, proteome_map):
    if (ome_id in proteome_map):
        return proteome_map[ome_id]['ena_acc']
    else:
        return ome_id

def gca_ome_id(ena_acc, proteome_map):
    if (ena_acc in proteome_map):
        return proteome_map[ena_acc]['ome_id']
    else:
        return ena_acc

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
            proteome_map[proteome_info['ena_acc']] = proteome_info

    return proteome_map

def main():

    ## format of Proteins_cluster_m50pct.tsv
    f_names = "cluster_id protein_ids proteins_count proteomes_count representative seqlen_range".split(' ')
    ## indexes of fields, to avoid building dictionary of data line
    fname_idx = dict(zip(f_names,range(len(f_names))))

    args = check_args()

    cluster_idx = 0
    rep_idx = fname_idx['representative']
    prots_idx = fname_idx['protein_ids']
    range_idx = fname_idx['seqlen_range']
    proteins_cnt_idx = fname_idx['proteins_count']
    proteomes_cnt_idx = fname_idx['proteomes_count']

    bad_example_cnt = 0

    clust_align_widths = {}
    out_clust_fd = sys.stdout

    bad_cluster_list = []
    bad_cluster_data = {}
    max_clust_bad = 100

    if (args.in_clust_file):
        in_clust_fields = []
        with open(args.in_clust_file,'r') as in_clust_fd:
            idx_clust_bad = 0
            for line in in_clust_fd:
                if (line.startswith('#')):
                    continue
                in_data = line.strip('\n').split('\t')
                if (in_data[0] == 'clust_id'):
                    in_clust_fields = [x for x in in_data]
                    continue

                clust_data = dict(zip(in_clust_fields, in_data))
                if (float(clust_data['pct_mode']) > args.mode_pct and
                    clust_data['bad_SL_flag']=='Short' and
                    int(clust_data['mode_len']) > 100):
                    if (float(clust_data['pct_bad']) < 2.0 and len(bad_cluster_list) > 20):
                        break
                    bad_cluster_list.append(clust_data['clust_id'])
                    bad_cluster_data[clust_data['clust_id']] = clust_data

        if (len(bad_cluster_list) > 10):
            bad_cluster_samp = random.sample(bad_cluster_list,10)
        else:
            bad_cluster_samp = bad_cluster_list

        bad_cluster_samp.sort(key=lambda x: float(bad_cluster_data[x]['pct_bad']),reverse=True)

        bad_samp_cluster_data = {x:bad_cluster_data[x] for x in bad_cluster_samp}

    else:
        sys.stderr.write("no --in_cluster_file")
        exit(1)

    proteome_map = {}
    if args.prot_id_map:
        proteome_map = read_proteome_map(args.prot_id_map)
    else:
        sys.stderr.write("no proteome map")
        exit(1)

    ome_samp_data = {}
    if (args.in_samp_omes):
        in_ome_fields = []
        with open(args.in_samp_omes,'r') as in_clust_fd:
            for line in in_clust_fd:
                if (line.startswith('#')):
                    continue
                in_data = line.strip('\n').split('\t')
                if (in_data[0] == 'proteome_id'):
                    in_ome_fields = [x for x in in_data]
                    continue

                if (in_data[-1]==args.exclude):  ## skip bad omes, focus on sampled
                    continue

                ome_data = dict(zip(in_ome_fields, in_data))
                ome_samp_data[gca_ome_id(ome_data['proteome_id'],proteome_map)]=ome_data

        print("## len(in_samp_omes): %d"%(len(list(ome_samp_data.keys()))),file=sys.stderr)


    else:
        sys.stderr.write("no --in_samp_omes")
        exit(1)


    print("\t".join("clust_id ome_id gca_id prot_len pct_mode".split(" ")))

    for in_file in args.files:

        with open(in_file,'r') as f_in:

            for line in f_in:
                if line.startswith("cluster_id"):
                    continue

                ## grab the cluster line
                clust_data  = line.strip('\n').split('\t')

                cluster_id = clust_data[cluster_idx]

                if (cluster_id not in bad_samp_cluster_data):
                    continue

                ## here we have a targeted bad cluster, now report the protein lengths for the target omes

                print("## "+' '.join([bad_samp_cluster_data[cluster_id][x] for x in ['clust_id','n_omes','p_tot','p_bad','pct_short','pct_mode']]))

                mode_len = int(bad_samp_cluster_data[cluster_id]['mode_len'])

                proteome_set = set()
                for prot_id in clust_data[fname_idx['protein_ids']].split(' '):
                    (this_proteome, this_prot_info) = prot_id.split(':')
                    this_info_split = this_prot_info.split('|')
                    (this_prot_acc, this_prot_len) = this_info_split[:2]

                    if (this_proteome in ome_samp_data and this_proteome not in proteome_set):
                        proteome_set.add(this_proteome)
                        ## have a sampled ome with target cluster
                        print("\t".join([cluster_id, this_proteome,ome_samp_data[this_proteome]['proteome_id'], this_prot_len, "%.1f"%(100.0*int(this_prot_len)/mode_len)]))

if __name__ == "__main__":

    main()
