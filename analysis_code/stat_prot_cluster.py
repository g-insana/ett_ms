#!/usr/bin/env python3

## 31-Jan-2026
## extended to count number of distinct proteins from each proteome (allowing duplicates in a cluster, or in multiple clusters)
## possibly report to stderr duplicates in cluster and duplicates in multiple clusters

## 23-Jan-2026
##
## extended to identify proteomes missing from a cluster, so that the
## mode protein for that cluster, can be matched with a set of
## proteomes that lack the cluster, so see if the mode protein is
## really missing.

## needs to have a way of identifying proteomes that are missing from a cluster

## 6-Jan-2026
##
## modified to ensure that only one protein per proteome is counted.
## may need to check for longest, shortest, to ensure that the correct
## one is saved, or simply randomize.
##

## 26-Nov-2026
##
## stat_prot_cluster.py
##
## the goal is to produce a table that can be used to draw a set of plots similar to 
## the cluster plots from f2ABCD_all_out_good_4p_pub.pdf
## (A) proteomes in cluster / for this plot, clusters in proteome
## (B) outliers per cluster / outliers per proteome
## (C) no outliers (per cluster) / no outliers (per proteome)
## (D) short outliers (per cluster) / short outliers (per proteome)
##
## these stats must include proteomes without outliers -- currently
## the 2000 sampled are sampled from proteomes WITH outliers

## using three files, bad_clust, bad_omes, and the raw clustering data
## (also need mapping from proteome number to GCA):
##
## (1) read the raw clustering data, reporting (for each proteome)
##     (a) how many clusters is the proteome in 
##     (b) how many good clusters is the proteome in
##     (c) how many mode length protein
##     (d) how many outlier clusters (good + outlier)
##     (e) how many short outliers
##     (f) is this a sampled proteome?
## (2) for each bacteria: count the number of proteomes
##     (a) no outliers/ outliers / total
## (3) read the 2000 sampled proteomes
##     (a) 
##

import sys
import fileinput
import argparse
import statistics as stat
import random
import copy
import re
import math
import random
import bisect

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
        "--fract",
        "-f",
        dest="thresh_fract",
        type=float,
        default = 0.75,
        help="threshold canonical length",
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
        help="exclude from sample set (B)",
    )

    parser.add_argument(
        "--miss_samp_ex",
        dest="miss_samp_ex",
        type=int,
        default=20,
        help="number of missing clusters to search")

    parser.add_argument(
        "--miss_samp_ex2",
        dest="miss_samp_ex2",
        type=int,
        default=20,
        help="number of proteomes to search for missing clusters")

    parser.add_argument(
        "--miss_samp_file",
        dest="miss_samp_file",
        type=str,
        help="destination for missing cluster proteomes")

    parser.add_argument(
        dest="files",
        type=str,
        nargs='*')

    return parser.parse_args()

def ome_id_gca(ome_id, proteome_id_map):
    if (ome_id in proteome_id_map):
        return proteome_id_map[ome_id]
    else:
        return ome_id

def gca_ome_id(ena_acc, proteome_gca_map):
    if (ena_acc in proteome_gca_map):
        return proteome_gca_map[ena_acc]
    else:
        return ena_acc

## find the lowest position of value in list using pseudo-binary search

def cluster_bin_search(key_list, target_dict, field, value):

    left = 0
    right = len(key_list) - 1

    while left <= right:
        # Calculate the middle index (using floor division //)
        mid = (left + right) // 2
        
        # If the value is greater than the middle element, ignore the left half
        if target_dict[key_list[mid]][field] < value:
            left = mid + 1
        # If the target is smaller than the middle element, ignore the right half
        else:
            right = mid - 1

    # If the loop finishes without finding the value, it is not in the list
    return left

## reads proteome mapping/BUSCO/metadata info
def read_proteome_map(map_file):
    ## proteome_map{} dict maps proteome_id's to GCA accessions, but
    ## also provides BUSCO and metadata if available

    proteome_id_map = {}
    proteome_gca_map = {}
    proteome_idx = 0

    ome_fields = "proteome_id up_acc ena_acc".split(" ")

    ## ignoring q1_len,med_len,q3_len

    with open(map_file,'r') as pid_f:
        for line in pid_f:
            proteome_data = line.strip('\n').split('\t')
            proteome_id = proteome_data[proteome_idx]

            proteome_info = dict(zip(ome_fields,proteome_data[:3]))
            proteome_info['ome_id'] = proteome_id
            proteome_id_map[proteome_id] = proteome_info['ena_acc']
            proteome_gca_map[proteome_info['ena_acc']] = proteome_id

    return (proteome_id_map, proteome_gca_map)

def main():

    if (len(sys.argv) > 10):
        sys_argv_str = "# " + ' '.join(sys.argv[:9]) + ' ... ' + sys.argv[-1]
    else:
        sys_argv_str = "# "+' '.join(sys.argv)

    print(sys_argv_str)

    ## format of Proteins_cluster_m50pct.tsv
    f_names = "cluster_id protein_ids proteins_count proteomes_count representative seqlen_range".split(' ')

    summ_int_fields = "n_omes p_tot p_bad p_short mode_len p_mode".split(" ")
    summ_float_fields = "pct_bad pct_short pct_mode".split(" ")

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

    all_cluster_list = []
    all_cluster_data = {}
    max_clust_bad = 100

    ## read in .bad_clust summary file, get summary data for all clusters
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

                for f in summ_float_fields:
                    clust_data[f] = float(clust_data[f])
                for f in summ_int_fields:
                    clust_data[f] = int(clust_data[f])

                all_cluster_list.append(clust_data['clust_id'])
                all_cluster_data[clust_data['clust_id']]=clust_data

    else:
        sys.stderr.write("no --in_cluster_file")
        exit(1)

    ## read in proteome_id to GCA mapping
    if args.prot_id_map:
        (proteome_id_map, proteome_gca_map) = read_proteome_map(args.prot_id_map)
    else:
        sys.stderr.write("no proteome map")
        exit(1)

    ## this marks whether a proteome is part of the sample set
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

                ome_data = dict(zip(in_ome_fields, in_data))

                this_ome_id = gca_ome_id(ome_data['proteome_id'],proteome_gca_map)
                if (this_ome_id in ome_samp_data):
                    ome_samp_data[this_ome_id]['s_type'] += ome_data['s_type']
                else:
                    ome_samp_data[this_ome_id] = ome_data

        print("## len(in_samp_omes): %d"%(len(list(ome_samp_data.keys()))),file=sys.stderr)

    else:
        sys.stderr.write("no --in_samp_omes")
        exit(1)


    out_fields="GCA_id ome_id cl_prot_cnt clust_cnt mode_cnt outlier_cnt short_cnt in_sample".split(" ")

    print("\t".join(out_fields))

    ## now get the raw cluster data

    ## args.files have cluster data file(s)
    for in_file in args.files:

        all_cluster_set = set()
        protein_stats = {}
        proteome_stats = {}
        with open(in_file,'r') as f_in:

            for line in f_in:
                if line.startswith("cluster_id"):
                    continue

                ## set up ome_info for this cluster
                this_clust_omes = {}
                ## grab the cluster line
                clust_data  = line.strip('\n').split('\t')
                cluster_id = clust_data[cluster_idx]

                all_cluster_set.add(cluster_id)

                mode_len = int(all_cluster_data[cluster_id]['mode_len'])
                mode_thresh_short = int(mode_len * args.thresh_fract - 0.5)
                mode_thresh_long = int(mode_len/args.thresh_fract + 0.5)

                have_mode_prot = False
                for prot_id in clust_data[fname_idx['protein_ids']].split(' '):
                    (this_proteome, this_prot_info) = prot_id.split(':')
                    (this_prot_acc, this_prot_len) = this_prot_info.split('|')[0:2]
                    this_prot_len = int(this_prot_len)

                    if (not have_mode_prot and this_prot_len == mode_len):
                        have_mode_prot = True
                        all_cluster_data[cluster_id]['mode_protein']=this_prot_info
                        all_cluster_data[cluster_id]['mode_proteome']=this_proteome
                        
                    ## we want: cluster_cnt, outlier, short_outlier, mode_length?
                    if (this_proteome not in this_clust_omes):
                        this_clust_omes[this_proteome] = [this_prot_len]
                    else:
                        this_clust_omes[this_proteome].append(this_prot_len)

                ## done with for (prot_id in clust_data)

                ## need to check if mode threshold is met
                
                mode_cnt = 0
                prot_cnt = 0
                for this_proteome in this_clust_omes.keys():
                    prot_cnt += 1
                    if (this_clust_omes[this_proteome][0] == mode_len):
                        mode_cnt += 1

                if (float(prot_cnt)*args.mode_pct/100.0 >= mode_cnt):
                    continue

            ## now scan through this_clust_omes and update stats

                for this_proteome in this_clust_omes.keys():

                    this_ome_list = this_clust_omes[this_proteome]
                    if (len(this_ome_list) > 1):
                        this_prot_len = random.choice(this_ome_list)
                    else:
                        this_prot_len = this_ome_list[0]

                    if (this_proteome not in proteome_stats):
                        GCA_id = ome_id_gca(this_proteome,proteome_id_map)
                        proteome_stats[this_proteome] = {'GCA_id':GCA_id, 'ome_id':this_proteome, 'cl_prot_cnt':0, 'clust_cnt':0, 'outlier_cnt':0, 'short_cnt':0, 'mode_cnt':0, 'in_sample':False, 'clust_set':set()}

                    ## look up proteome in proteome_id_map, check if in sample
                    if (this_proteome in ome_samp_data):
                        proteome_stats[this_proteome]['in_sample']=ome_samp_data[this_proteome]['s_type']
                    else:
                        proteome_stats[this_proteome]['in_sample']='N'

                    this_data = proteome_stats[this_proteome]
                    this_data['clust_cnt'] += 1
                    this_data['cl_prot_cnt'] += len(this_clust_omes[this_proteome])

                    this_data['clust_set'].add(cluster_id)

                    if (this_prot_len < mode_thresh_short):
                        this_data['outlier_cnt'] += 1
                        this_data['short_cnt'] += 1
                    elif (this_prot_len > mode_thresh_long):
                        this_data['outlier_cnt'] += 1

                    if (this_prot_len == mode_len):
                        this_data['mode_cnt'] += 1

        ## here we have read all of the clusters and looked at all the
        ## proteomes in all the clusters we could sort the proteomes by
        ## the number of clusters they belong to, so that that are missing
        ## dozens (but not hundreds?) of clusters are identified.  we then
        ## need to figure out which clusters they are missing from the
        ## proteomes that are missing the most clusters

        if (args.miss_samp_file and len(args.miss_samp_file)>0):

            miss_samp_fd = open(args.miss_samp_file,'w')

            ## get list of proteomes, from fewest clusters to most clusters
            missed_prot_keys = sorted(proteome_stats.keys(),key=lambda x: proteome_stats[x]['clust_cnt'])

            ## get large-ish number of clusters
            clust90_key = missed_prot_keys[int(0.9*len(missed_prot_keys))]

            ## get a threshold where clusters are probably missing
            good_clust_thresh = int(0.9*proteome_stats[clust90_key]['clust_cnt'])

            ## find the last "good" proteome
            last_good_ix = cluster_bin_search(missed_prot_keys, proteome_stats, 'clust_cnt',good_clust_thresh-1)

            ## this gives us some proteomes that are missing clusters
            ## we now need clusters that are missing those proteomes
            ## for each proteome, we have the clust_set and we have all_cluster_sets, so we can find clusters
            ## missing proteomes

            samp_ix = last_good_ix//2

            ## try sampling around the middle, up and down for
            ## args.miss_samp_ex/2 then, get args.miss_samp_ex2 proteomes
            ## to see if proteins are really missing

            if (samp_ix > args.miss_samp_ex) :
                samp_low = samp_ix - args.miss_samp_ex//2
            else:
                samp_low = max(0,last_good_ix - args.miss_samp_ex)

            samp_high = min(samp_low+args.miss_samp_ex,last_good_ix)

            ## print out the cluster/mode information

            for s_ix in range(samp_low,samp_high):
                this_ome_stat = proteome_stats[missed_prot_keys[s_ix]]
                print('>'+'\t'.join([this_ome_stat['GCA_id'],str(this_ome_stat['clust_cnt']),this_ome_stat['in_sample']]),file=miss_samp_fd)
                missed_cluster_set = all_cluster_set - proteome_stats[missed_prot_keys[s_ix]]['clust_set']

                missed_cluster_list = [ x for x in missed_cluster_set]
                for clust_id in random.sample(missed_cluster_list, min(args.miss_samp_ex2, len(missed_cluster_list))):
                    this_cluster = all_cluster_data[clust_id]
                    if (this_cluster['mode_len'] >= 200):
                        print('\t'.join([clust_id,this_cluster['mode_protein'],str(this_cluster['mode_len'])]),file=miss_samp_fd)

        ## done with result file, write out proteomes

        s_stats_keys = sorted(proteome_stats.keys(),key=lambda x: -proteome_stats[x]['outlier_cnt'])

        for ome_id in s_stats_keys:
            this_ome_data = proteome_stats[ome_id]
            print('\t'.join([str(this_ome_data[x]) for x in out_fields]))

if __name__ == "__main__":

    main()
