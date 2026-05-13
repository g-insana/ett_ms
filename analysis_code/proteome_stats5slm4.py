#!/usr/bin/env python3

## 02-Mar-2026
## print full protein sequence id (possibly acc|len|upi or acc|len) in examples for consistency with down_fasta_oscode_*.sh
##

## 09-Feb-2026 (proteome_stats5slm4.py) modified to accept
## --ome_subset, which provide a set of sampled proteomes so that all
## statistics are based on clusters that only include this subset of proteomes

## 06-Feb-2026
## modified to read and interpret new cluster files, that have a prot_id of
## 12345:P09488|218|UPI_98765  (proteome_id:acc|len|up_acc)
##

## 26-Nov-2025
## two changes:
##  (1) make certain that sampled proteomes have reasonable mode
##  lengths that reflect the overall statistics for short vs long outliers
##  (2) double double check that we are not looking at duplicate proteins per proteome
##

## 10-Aug-2025
## for short sequences, we check to see if the long sequence can be
## found in the short sequence proteome.  So we find 'N' short
## sequences from the same proteome and check to see if there is a
## long (mode) sequence in that proteome.

## for long sequences, we want to see if the long sequence is in the
## mode-sequence proteome, so rather than looking at multiple clusters
## in the same proteome, we need to look at multiple mode-proteomes
## for the same cluster.  Older versions looked at only one
## mode-proteome. This version (proteome_stats4ln.py) reports multiple
## mode proteomes, to see how often the long sequence can be found in
## them.

## 28-Aug-2025
## we are also interested in whether some proteomes produce more long
## sequences, while others produce more short sequences.  So it makes
## some sense to produce a proteome-based set of long_ex as well as a
## cluster set of long_ex.

## 5-Jun-2025
## produce statistics on long (and short) sequences with 'XXXX'?  this
## requires looking at the sequences, and they are not available
## without extraction.  But we could look for non 'XXXX' sequences
## when selecting examples

## 4-Jun-2025
## modify the long cluster examples to avoid clusters with a mode length < 100.
##

## 2-Jun-2025
##
## while short proteins may be a property of the proteome, long
## proteins are probably a property of the cluster, so long protein
## examples (--long_prots) should be selected to (1) have no
## duplicates and (2) sample multiple clusters.  Thus, the logic for
## long protein examples is quite different from that of short proteins.
## Basically "bad_proteomes_long" should be "bad_clusters_long"

## 21-May-2025
##
## modified to extract information on "long bad" sequences.  Need mode
## (reference) sequence and genome GCA, long sequence and genome GCA
## check to see if extended sequence is in reference genome GCA, how
## often is it in other databases?
## add --long_prots --long_tfx_ex, --long_tfx_file
##

## 16-May-2025
## modified to read only two files:  "Protein_clusters_m60000.tsv" and a mapping/BUSCO/metadata file
##

## read a "Protein_clusters_m60000.tsv" file
## fields: cluster_id	protein_ids	proteins_count	proteomes_count	representative	seqlen_range
## count how many proteomes have a length < 0.75 of representative
## ?write out proteome(s) with short ortholog

## modified to check for duplicate proteome_id's -- only want
## short/long proteins with unique proteome_id

## include code to discover the most "median" proteomes.


## 15-May-2025
## modify to mark clusters with > 50% long proteins from short proteins 

## 9-April-2025
## need to make certain that samp_tfx proteomes are in .bad_omes.samp set
## so get samp_tfx proteomes from that set -- this has been confirmed
##
## There are 2 options for the number of bad/sampled proteomes to be displayed (--bad_prots, --samp_prots)
##           2 options for the number of proteomes to be examined with tfx (--bad_tfx_ex, --samp_tfx_ex)
##           1 option for the number of proteins to be examined for bad/samp tfx searches (--tfx_search_cnt)

## 2-Mar-2025
## modified to write cluster stats and proteome stats to files
## (--out_clust_file,--out_omes_file), implement sampling across all proteomes

## 3-Mar-2025

## provide sampled protein examples for worst and randomly selected clusters,
## and for worst and randomly selected proteomes
## add data type for sampled vs worst proteomes

## 21-Mar-2025

## change proteome sampling to sample from different ends of N50
## spectrum IF N50 information available
## Also label 'B' (bad), 'S' (sampled) to include N50 info
##

## 22-Mar-2025 -- filter out short clusters for tfx examples

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
        "--fract",
        "-f",
        dest="thresh_fract",
        type=float,
        default = 0.75,
        help="threshold canonical length",
    )
    parser.add_argument(
        "--bad_low_fract",
        dest="bad_low_fract",
        type=float,
        default = 0.1,
        help="lower fraction for bad display (0.1)",
    )
    parser.add_argument(
        "--bad_prots",
        "-b",
        dest="bad_prots",
        type=int,
        default = 20,
        help="# of worst proteomes to display",
    )
    parser.add_argument(
        "--mode_pct",
        dest="mode_pct",
        type=float,
        default = 60,
        help="threshold for good cluster mode_pct",
    )

    parser.add_argument(
        "--mode_good_only",
        dest="mode_good_only",
        action='store_true',
        default = False,
        help="threshold for good cluster mode_pct",
    )

    parser.add_argument(
        "--samp_prots",
        dest="samp_prots",
        type=int,
        default = 1000,
        help="# sampled proteomes to display",
    )
    parser.add_argument(
        "--bad_tfx_ex",
        dest="bad_tfx_ex",
        type=int,
        default = 20,
        help="# of tfastx samples per bad proteome/cluster",
    )

    parser.add_argument(
        "--long_tfx_ome_file",
        dest="long_tfx_ome_file",
        type=str,
        help="destination for long tfx ome info"
    )

    parser.add_argument(
        "--long_tfx_clust_file",
        dest="long_tfx_clust_file",
        type=str,
        help="destination for long tfx clust info"
    )

    parser.add_argument(
        "--long_tfx_ex",
        dest="long_tfx_ex",
        type=int,
        default=20,
        help="number of mode proteomes for long examples"
    )

    parser.add_argument(
        "--long_prots",
        dest="long_prots",
        type=int,
        default=100,
        help="number of long outlier examples"
    )

    parser.add_argument(
        "--no_excl",
        dest="no_excl",
        action='store_true',
        default = False,
        help="ignore excluded proteomes"
    )

    parser.add_argument(
        "--no_redund",
        dest="no_redun",
        action='store_true',
        default = False,
        help="ignore redundant proteomes"
    )

    parser.add_argument(
        "--show_mode_region",
        dest="mode_region",
        type=int,
        default = 0,
        help="show region around mode",
    )

    parser.add_argument(
        "--samp_N50",
        dest="do_N50",
        action='store_true',
        default = False,
        help="stratify examples by N50"
    )
    parser.add_argument(
        "--samp_tfx_ex",
        dest="samp_tfx_ex",
        type=int,
        default = 20,
        help="# of tfastx samples per sampled proteome/cluster",
    )
    parser.add_argument(
        "--tfx_search_cnt",
        dest="tfx_search_cnt",
        type=int,
        default = 20,
        help="# of tfastx samples per sampled proteome/cluster",
    )
    parser.add_argument(
        "--out_clust_file",
        dest='out_clust_file',
        type=str,
        help="file for bad clusters"
        )

    parser.add_argument(
        "--out_omes_file",
        dest='out_omes_file',
        type=str,
        help="file for bad clusters"
        )

    parser.add_argument(
        "--out_tfx_file",
        dest='out_tfx_file',
        type=str,
        help="file for tfx examples"
    )

    parser.add_argument(
        "--example_cnt",
        "-x",
        dest="example_cnt",
        type=int,
        default = 0,
        help="threshold canonical length",
    )

    parser.add_argument(
        "--target_length",
        dest="target_length",
        type=int,
        default=200,
        help="min length for target acc info",
    )

    parser.add_argument(
        "--bad_file",
        dest="bad_file",
        type=str,
        help="file for bad proteome/protein",
    )

    parser.add_argument(
        "--ome_subset",
        dest="ome_subset",
        type=str,
        help="file with the only proteomes that will be examined",
    )

    parser.add_argument(
        "-P",
        "--prot_id_map",
        dest="prot_id_map",
        type=str,
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

def print_proteomes(this_ome_key, bad_proteomes_cnt, n_clust, proteome_map, have_busco, busco_fields_out, s_type, out_omes_fd):

    if (this_ome_key == ''):
        return

    bad_proteomes_name = ome_id_gca(this_ome_key, proteome_map)

    bad_cnt_total = bad_proteomes_cnt[this_ome_key]['total']
    bad_cnt_short = bad_proteomes_cnt[this_ome_key]['short']

    if (bad_cnt_total  > 0):
        pct_short = 100.0*bad_cnt_short / float(bad_cnt_total)
    else:
        pct_short = 0.0

    print("\t".join((bad_proteomes_name, str(bad_cnt_total),str(bad_cnt_short),"%.3g"%(bad_cnt_total*100.0/n_clust),"%.1f"%(pct_short))),end='',file=out_omes_fd)
    if (have_busco):
        if (this_ome_key in proteome_map):
            if ('B_prot_cnt' in proteome_map[this_ome_key] and int(proteome_map[this_ome_key]['B_prot_cnt'])>0):
                print("\t"+"\t".join([str(proteome_map[this_ome_key][x]) for x in busco_fields_out]),end='',file=out_omes_fd)
            else:
                print("\t"+"\t".join(['0' for x in busco_fields_out]),end='', file=out_omes_fd)
        else:
            sys.stderr.write("***WARNING %s not found in proteome_map\n"%(this_ome_key))

    print("\t"+s_type,file=out_omes_fd)


## print_long_clust_ex(cluster_id,
##              this_bad_cluster_info [{'proteome':,'acc','len','is_short'},....]
##              proteome_map,
##              this_clust_stats
##              n_ex_clust, s_type, args, out_tfx_fd)
## if there is more than one long protein in a cluster, print the longest and the shortest
##

def print_long_clust_ex(cluster_id, this_info, mode_list, proteome_map, mode_stats, s_type, args, out_tfx_fd):

    if (not cluster_id):
        return

    ## this_bad_cluster_info is a list 
    this_info.sort(key=lambda x: -x['len'])
    this_info0 = this_info[0]

    ## list bad cluster, plus some info about it
    long_proteome = this_info0['proteome']

    if (long_proteome in proteome_map):
        long_proteome = proteome_map[long_proteome]['ena_acc']

    print(f">{cluster_id}\t{mode_stats['tot_cnt']}\t{mode_stats['bad_cnt']}\t{mode_stats['bad_cnt_long']}\t{long_proteome}\t{this_info0['seq_id']}",file=out_tfx_fd)

    ## then show the two long proteins (and associated proteomes)
    for mode_prot in mode_list:
        this_ome_key = mode_prot['proteome']

        if (this_ome_key in proteome_map):
            mode_proteome_name = proteome_map[this_ome_key]['ena_acc']
            if (re.search('_N50',s_type) and 'A_n50' in proteome_map[this_ome_key]):
                n50_s_type = "%s:%.2f"%(s_type,math.log10(float(proteome_map[this_ome_key]['A_n50'])))
            else:
                n50_s_type = s_type
        else:
            mode_proteome_name = this_ome_key
            n50_s_type = s_type

        print('\t'.join((mode_proteome_name,s_type,mode_prot['seq_id'])),file=out_tfx_fd)

## modified 4-July-2025 to provide short protein id from short protein genome
def print_tfx_ex(this_ome_key, bad_proteome_info, bad_proteome_cnt, proteome_map, clust_stats, n_ex_clust, s_type, args, out_tfx_fd):

    if (this_ome_key == ''):
        return

    if (this_ome_key in proteome_map):
        bad_proteomes_name = proteome_map[this_ome_key]['ena_acc']
    else:
        bad_proteomes_name = this_ome_key

    ## list bad proteome, plus some info about it -- percent bad clusters
    print(f">{bad_proteomes_name}\t{bad_proteome_cnt}\t{s_type}",file=out_tfx_fd)

    ## only use examples longer than args.target_length
    bad_clust_keys = []
    for bad_prot_k in bad_proteome_info.keys():

        bad_prot_len = bad_proteome_info[bad_prot_k]['len']
        if (bad_prot_len > clust_stats[bad_prot_k]['mode_len']):
            bad_info = bad_proteome_info[bad_prot_k]
            print("[print_tfx_ex] protein %s/%s longer %d than mode %d"%(this_ome_key,bad_info['acc'],bad_info['len'],clust_stats[bad_prot_k]['mode_len']), 
                  file=sys.stderr)
            continue
        
        if (clust_stats[bad_prot_k]['mode_len'] > args.target_length):
            bad_clust_keys.append(bad_prot_k)

    ## if doing N50, sort bad_clust_keys by n50, sample top and bottom
    ## quartiles

    if (len(bad_clust_keys) > n_ex_clust):
        bad_clust_keys = random.sample(bad_clust_keys, n_ex_clust)

    for bad_clust_k in bad_clust_keys:
        bad_stats = clust_stats[bad_clust_k]
        bad_prot_acc = bad_proteome_info[bad_clust_k]['acc']
        bad_prot_len = bad_proteome_info[bad_clust_k]['len']
        # bad_prot_id_str = "%s|%d"%(bad_prot_acc,bad_prot_len)
        bad_prot_id_str = bad_proteome_info[bad_clust_k]['seq_id']
        ## the short protein (fourth column) is 
        print('\t'.join((bad_clust_k, bad_stats['mode_id'], str(bad_stats['mode_len']),bad_prot_id_str)),file=out_tfx_fd)

## modified 4-July-2025 to provide short protein id from short protein genome

def print_short_tfx_ex(this_ome_key, bad_proteome_info, this_bad_ome_cnt, proteome_map, clust_stats, n_ex_clust, s_type, args, out_tfx_fd):

    if (this_ome_key == ''):
        return

    bad_proteome_name = ome_id_gca(this_ome_key, proteome_map)

    ## list bad proteome, plus some info about it -- percent bad clusters
    print(f">{bad_proteome_name}\t{this_bad_ome_cnt}\t{s_type}",file=out_tfx_fd)

    bad_clust_keys = []
    for bad_clust_k in bad_proteome_info.keys():
        if (clust_stats[bad_clust_k]['mode_len'] > args.target_length):
            bad_clust_keys.append(bad_clust_k)

    ## if doing N50, sort bad_clust_keys by n50, sample top and bottom
    ## quartiles

    if (len(bad_clust_keys) > n_ex_clust):
        bad_clust_keys = random.sample(bad_clust_keys, n_ex_clust)

    for bad_clust_k in bad_clust_keys:
        bad_stats = clust_stats[bad_clust_k]
        bad_prot_acc = bad_proteome_info[bad_clust_k]['acc']
        bad_prot_len = bad_proteome_info[bad_clust_k]['len']
##        bad_prot_id_str = "%s|%d"%(bad_prot_acc,bad_prot_len)
        bad_prot_id_str = bad_proteome_info[bad_clust_k]['seq_id']
        print('\t'.join((bad_clust_k, bad_stats['mode_id'], str(bad_stats['mode_len']), bad_prot_id_str)),file=out_tfx_fd)

## modified 28-Aug-2025 for long protein info

def print_long_tfx_ex(this_ome_key, bad_proteome_info, this_bad_ome_cnt, proteome_map, clust_stats, n_ex_clust, s_type, args, out_tfx_fd):

    if (this_ome_key == ''):
        return

    bad_proteome_name = ome_id_gca(this_ome_key, proteome_map)

    bad_clust_keys = []
    for bad_clust_k in bad_proteome_info.keys():
        if (clust_stats[bad_clust_k]['mode_len'] > args.target_length):
            bad_clust_keys.append(bad_clust_k)

    ## if doing N50, sort bad_clust_keys by n50, sample top and bottom
    ## quartiles

    ## list bad proteome, plus some info about it -- percent bad clusters
    print(f">{bad_proteome_name}\t{this_bad_ome_cnt}\t{len(bad_clust_keys)}\t{s_type}",file=out_tfx_fd)

    if (len(bad_clust_keys) > n_ex_clust):
        bad_clust_keys = random.sample(bad_clust_keys, n_ex_clust)

    for bad_clust_k in bad_clust_keys:
        bad_stats = clust_stats[bad_clust_k]
        ## in this loop, bad_stats has:
        ## {'tot_cnt': 1593, 'bad_cnt': 70, 'bad_cnt_short': 0, 'bad_cnt_long': 70, 'bad_fract': 0.04394224733207784,
        ## 'prot_cnt': 1593, 'bad_fract_short': 0.0, 'bad_SL_flag': 'Long', 'ome_cnt': 1592, 'mode_proteome': '1336781',
        ## 'mode_id': 'CPI93547|201', 'mode_len': 201, 'mode_cnt': 1250,
        ## 'width': None, 'q1_len': 201, 'med_len': 201, 'q3_len': 201}
        ## which provides mode_proteome and mode_id

        mode_ome_key = bad_stats['mode_proteome']
        mode_proteome_name = ome_id_gca(mode_ome_key, proteome_map)

        ## bad_proteome_info[bad_clust_k] has:
        ## {'proteome': '276375', 'acc': 'EEC35266', 'len': 302, 'is_short': False}

        ## which provides the long protein and its proteome

        bad_prot_acc = bad_proteome_info[bad_clust_k]['acc']
        bad_prot_len = bad_proteome_info[bad_clust_k]['len']
##        bad_prot_id_str = "%s|%d"%(bad_prot_acc,bad_prot_len)
        bad_prot_id_str = bad_proteome_info[bad_clust_k]['seq_id']

        ## to be similar to the .long_ex (long) file, I need:
        ## v-bad-proteome (with multiple bad clusters)
        ## >GCA_021083485.1	109	B_N50L:4.17
        ## v-cluster v-mode-proteome    mode-protein-id long_protein
        ## 61938     GCA_029310475.1	CPF18766|653	MDF3707975|724

        print('\t'.join((bad_clust_k, mode_proteome_name, str(bad_stats['mode_id']), bad_prot_id_str)),file=out_tfx_fd)

## print a list of _ex entries for bad proteomes
## bad proteome examples can be short or long, so now uses print_tfx_ex_func()
##

def print_tfx_ex_list(ome_list, bad_proteomes, bad_proteome_cnt, b_type, proteome_map,
                      clust_stats, example_cnt, s_type, args, out_tfx_fd, print_tfx_ex_func):

    do_n50_str = re.search('_N50',s_type)

    for tfx_ome in ome_list:
        if (do_n50_str and tfx_ome in proteome_map and 'A_n50' in proteome_map[tfx_ome]):
            n50_s_type = "%s:%.2f"%(s_type,math.log10(float(proteome_map[tfx_ome]['A_n50'])))
        else:
            n50_s_type = s_type

        ## print(ome_id_gca(tfx_ome,proteome_map), b_type, bad_proteome_cnt[tfx_ome],file=sys.stderr)

        print_tfx_ex_func(tfx_ome, bad_proteomes[tfx_ome], bad_proteome_cnt[tfx_ome][b_type], proteome_map,
                     clust_stats, example_cnt, n50_s_type, args, out_tfx_fd)

## reads proteome mapping/BUSCO/metadata info
def read_proteome_map(map_file):
    ## proteome_map{} dict maps proteome_id's to GCA accessions, but
    ## also provides BUSCO and metadata if available

    proteome_map = {}

    ## modify to allow busco scores
    have_busco = False
    busco_fields = 'proteome_id upid ena_acc B_prot_cnt B_assem_lvl B_compl_comb B_compl_single B_compl_dup B_frag B_miss B_compl_descr B_is_rep B_is_ref B_is_excl B_is_redun'.split(' ')

    busco_fields2 = 'upid proteome_id proteome_taxid species_taxid ena_acc biosample proj_acc B_prot_cnt B_assem_lvl B_compl_comb B_compl_single B_compl_dup B_frag B_miss B_compl_descr B_is_rep B_is_ref B_is_excl'.split(' ')

    busco_fields_out = 'B_prot_cnt B_assem_lvl B_compl_comb B_compl_single B_frag B_miss B_is_re_ex B_prot_status'.split(' ')

    assm_fields = 'A_wgs_set A_dt_created A_ft_country A_ft_coll_by A_center_name A_genome_rep A_tot_len A_n50 A_scaf_cnt A_contig_cnt A_contig_n50 A_contig_l50 A_contig_n75 A_contig_n90 A_scaffold_l50 A_scaffold_n75 A_scaffold_n90 A_assem_meth A_coverage A_seq_tech A_annot_date A_annot_pipe A_annot_method A_genes_tot'.split(' ')

    assm_fields_am = 'A_wgs_set A_dt_created A_ft_country A_ft_coll_by A_center_name A_genome_rep A_tot_len A_n50 A_scaf_cnt A_contig_cnt A_contig_n50 A_contig_l50 A_contig_n75 A_contig_n90 A_scaffold_l50 A_scaffold_n75 A_scaffold_n90 A_assem_meth A_coverage A_seq_tech A_annot_date A_annot_pipe A_annot_method A_genes_tot A_am'.split(' ')

    assm_fields_out = assm_fields
    assm_fields_out_am  = assm_fields_am

    busco_assm_fields = busco_fields + assm_fields
    busco_assm_fields_am= busco_fields + assm_fields_am
    busco_assm_fields2 = busco_fields2 + assm_fields
    busco_assm_fields2_am = busco_fields2 + assm_fields_am

    ome_fields = busco_fields[:3]

    proteome_idx = 0

    n1xx_fields = "A_tot_len A_n50 A_scaf_cnt A_coverage A_contig_n50 A_contig_l50 A_contig_n90 A_scaffold_l50 A_scaffold_n75 A_scaffold_n90".split(" ")
    nmax_fields = "A_contig_cnt".split(" ")

    with open(map_file,'r') as pid_f:
        for line in pid_f:
            proteome_data = line.strip('\n').split('\t')
            proteome_id = proteome_data[proteome_idx]
            if (not have_busco):
                if (len(proteome_data) == len(busco_fields)):
                    have_busco = True
                    busco_fields_in = busco_fields
                    proteome_idx = 0
                    continue
                elif (len(proteome_data) == len(busco_assm_fields)):
                    have_busco = True
                    busco_fields_in = busco_assm_fields
                    busco_fields_out = busco_fields_out + assm_fields_out
                    proteome_idx = 0
                    continue
                elif (len(proteome_data) == len(busco_assm_fields_am)):
                    have_busco = True
                    busco_fields_in = busco_assm_fields_am
                    busco_fields_out = busco_fields_out + assm_fields_out_am
                    ome_fields = ['proteome_id','upid','ena_acc']
                    proteome_idx = 0
                    continue
                elif (len(proteome_data) == len(busco_assm_fields2)):
                    have_busco = True
                    busco_fields_in = busco_assm_fields2
                    busco_fields_out = busco_fields_out + assm_fields_out
                    ome_fields = ['proteome_id','upid','ena_acc']
                    proteome_idx = 1
                    continue
                elif (len(proteome_data) == len(busco_assm_fields2_am)):
                    have_busco = True
                    busco_fields_in = busco_assm_fields2_am
                    busco_fields_out = busco_fields_out + assm_fields_out_am
                    ome_fields = ['proteome_id','upid','ena_acc']
                    proteome_idx = 1
                    continue
                else:
                    sys.stdout.write("*** WARNING -- unrecognized busco fields\n")

            if (have_busco):
                proteome_info = dict(zip(busco_fields_in,proteome_data))


                ## check for null fields, put in '1' if no value
                for f in n1xx_fields:
                    if (f in proteome_info and proteome_info[f]==''):
                        proteome_info[f]="1"

                for f in nmax_fields:
                    if (f in proteome_info and proteome_info[f]==''):
                        proteome_info[f]="1000000"

                ## is_rep, is_ref, is_excl, are t/f, is_redund is 1,-1,0
                ## as far as I can tell, you can only be one of these (or included)

                is_re_ex = 'inc'

                prot_status = 1
                if (proteome_info['B_is_rep']=='t'):
                    in_re_ex = 'rep'
                    prot_status = 0
                elif (proteome_info['B_is_ref']=='t'):
                    in_re_ex = 'ref'
                    prot_status = 0
                elif (proteome_info['B_is_excl'] == 't'):
                    in_re_ex = 'excl'
                    prot_status = 3
                elif ('B_is_redun' in proteome_info and proteome_info['B_is_redun'] == '1'):
                    in_re_ex = 'redund'
                    prot_status = 2

                proteome_info['B_is_re_ex'] = is_re_ex
                proteome_info['B_prot_status'] = str(prot_status)
            else:
                proteome_info = dict(zip(ome_fields,proteome_data[:3]))
                proteome_info['B_is_re_ex'] = 'None'
                proteome_info['B_prot_status'] = '1'

            proteome_map[proteome_id] = proteome_info

    return (proteome_map, have_busco, busco_fields_out)

def main():

    ## format of Proteins_cluster_m50pct.tsv
    f_names = "cluster_id protein_ids proteins_count proteomes_count representative seqlen_range seqlen_mode".split(' ')
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
    if (args.out_clust_file):
        out_clust_fd = open(args.out_clust_file,'w')

    out_omes_fd = sys.stdout
    if (args.out_omes_file):
        out_omes_fd = open(args.out_omes_file,'w')

    out_tfx_fd = sys.stdout
    if (args.out_tfx_file):
        out_tfx_fd = open(args.out_tfx_file,'w')

    long_tfx_clust_fd = None
    if (args.long_tfx_clust_file):
        long_tfx_clust_fd = open(args.long_tfx_clust_file,'w')

    long_tfx_ome_fd = None
    if (args.long_tfx_ome_file):
        long_tfx_ome_fd = open(args.long_tfx_ome_file,'w')

    have_ome_subset=False
    ome_subset = set()
    if (args.ome_subset):
        have_ome_subset=True
        with open(args.ome_subset,'r') as fd:
            for line in fd:
                if (line.startswith('GCA_id')):
                    continue
                ome_fields = line.strip('\n').split('\t')
                ome_subset.add(ome_fields[1])

    proteome_map = {}
    if args.prot_id_map:
        (proteome_map, have_busco,busco_fields_out) = read_proteome_map(args.prot_id_map)

    ## start looking for bad sequences
    if (args.thresh_fract > 1.0):
        args.thresh_fract = 1.0/args.thresh_fract

    ## now scan the Protein_m50pct.tsv file(s)

    for in_file in args.files:

        with open(in_file,'r') as f_in:
            bad_proteomes = {}
            bad_proteomes_short = {}
            bad_proteomes_long = {}
            bad_clusters_short = {}
            bad_clusters_long = {}
            bad_clusters_long_modes = {}
            bad_proteome_cnt = {}
            clust_stats = {}
            n_clust = 0
            bad_mode_cnt = 0
            bad_clust = 0
            bad_cluster_set = set()

            for line in f_in:
                if line.startswith("cluster_id"):
                    continue

                ## grab the cluster line
                n_clust += 1
                clust_data  = line.strip('\n').split('\t')

                cluster_id = clust_data[cluster_idx]

                ## get_representative length

                rep_id = clust_data[rep_idx]
                rep_id_list = clust_data[rep_idx].split('|')
                (rep_acc, rep_len) = rep_id_list[0:2]
                rep_len = int(rep_len)

                (clust_min, clust_max) = clust_data[range_idx].split('-')
                clust_min = int(clust_min)
                clust_max = int(clust_max)

                proteins_clust_cnt = int(clust_data[proteins_cnt_idx])

                if (cluster_id in clust_align_widths):
                    clust_width = clust_align_widths[cluster_id]
                else:
                    clust_width = None

                ## we know there is a short/long one
                ## extract list of proteins (with proteome and length)

                bad_cnt = 0
                bad_cnt_short = 0
                bad_cnt_long = 0

                clust_info = []
                uniq_prot_set = set()
                dup_prot_set = set()

                ## scan through cluster proteome:acc|len to check for duplicate proteomes

                ## this scan does not examine protein lengths, so has both good and bad to avoid duplicates

                clust_prot_cnt = 0
                for prot_id in clust_data[fname_idx['protein_ids']].split(' '):
                    (this_proteome, this_prot_info) = prot_id.split(':')

                    if (this_proteome == ''):
                        continue

                    if (have_ome_subset and this_proteome not in ome_subset):
                        continue

                    clust_prot_cnt += 1

                    prot_info_list = this_prot_info.split('|')
                    (this_prot_acc, this_prot_len) = prot_info_list[0:2]
                    this_prot_len = int(this_prot_len)
                    if (this_proteome in uniq_prot_set):
                        dup_prot_set.add(this_proteome)
                        continue
                    else:
                        uniq_prot_set.add(this_proteome)
                        clust_info.append({'proteome':this_proteome,
                                           'acc':this_prot_acc,
                                           'len':this_prot_len,
                                           'seq_id':this_prot_info})

                if (not ome_subset and clust_prot_cnt < proteins_clust_cnt):
                    sys.stderr.write(f"*** cluster_id {cluster_id} proteomes:{clust_data[proteomes_cnt_idx]} : cluster count mismatch {proteins_clust_cnt} != {clust_prot_cnt}\n")

                ## at this point, clust_info has all the distinct proteome accs|lens
                ## sort it to get max, median (possibly mode)

                clust_info.sort(key=lambda x: -x['len'])

                clust_size = len(clust_info);
                clust_midx = int(clust_size/2)
                clust_q1x = int(clust_size/4)
                clust_max = clust_info[0]
                
                clust_med = clust_info[clust_midx]
                clust_q3 = clust_info[clust_q1x]    ## backwards because of sort high to low
                clust_q1 = clust_info[clust_midx + clust_q1x]

                clust_lens = [x['len'] for x in clust_info]

                ## get the mode clust_len
                clust_len_cnt = {}
                for prot in clust_info:
                    if (prot['len'] in clust_len_cnt):
                        clust_len_cnt[prot['len']]['cnt'] += 1
                    else:
                        clust_len_cnt[prot['len']] = {'info':prot, 'cnt':1}

                clust_len_sort = sorted(list(clust_len_cnt.keys()), key=lambda x: -clust_len_cnt[x]['cnt'])
                ## clust_len_sort is largest number (mode)
                mode_cnt = clust_len_cnt[clust_len_sort[0]]['cnt'] 
                mode_info = clust_len_cnt[clust_len_sort[0]]['info']
                mode_len = mode_info['len']
                mode_acc = mode_info['acc']
                mode_proteome = mode_info['proteome']
                # mode_id = f'{mode_acc}|{mode_len}'
                mode_id = mode_info['seq_id']

                ## report median/mode mismatch
                if (clust_med['len'] != mode_len):
                    sys.stderr.write("*** median != mode\n")
                    sys.stderr.write(f'cluster: {cluster_id} max: {clust_max}\n')
                    sys.stderr.write(f'med: {clust_med}\n')
                    sys.stderr.write(f'mode: {clust_len_cnt[clust_len_sort[0]]}\n')

                ## check to see if this is a good (stable) mode
                mode_pct = 100.0*float(mode_cnt)/float(clust_size)
                mode_good = (mode_pct > args.mode_pct)

                if (not mode_good):
                    bad_mode_cnt += 1
                    print(f' *** bad cluster mode {cluster_id} mode_cnt: {mode_cnt} size: {clust_size} {mode_pct:.1f}%',file=sys.stderr)

                    # if (args.mode_good_only):
                    #     print(f' *** skipping {cluster_id}',file=sys.stderr)
                    #     continue

                mode_thresh_short = int(mode_len * args.thresh_fract - 0.5)
                mode_thresh_long = int(mode_len/args.thresh_fract + 0.5)

                mode_prot_ids = []

                ## could have a duplicate here, skip if so
                for this_prot_info in clust_info:
                    if (this_prot_info['proteome'] in dup_prot_set):
                        continue

                    this_proteome = this_prot_info['proteome']
                    this_prot_acc = this_prot_info['acc']
                    this_prot_len = this_prot_info['len']
                    
                    ## get mode proteins
                    if (this_prot_len == mode_len):
                        mode_prot_ids.append(this_prot_info)
                        continue

        	    ## find proteins outside range
                    if (this_prot_len < mode_thresh_short) or (this_prot_len > mode_thresh_long):

			## first example for bad_cluster
                        if (cluster_id not in bad_cluster_set):
                            bad_clust += 1
                            bad_cluster_set.add(cluster_id)

                        bad_cnt += 1
                        if (this_prot_len < mode_thresh_short):
                            is_short = True
                            bad_cnt_short += 1
                        else:
                            is_short = False
                            bad_cnt_long += 1

                        this_prot_info['is_short'] = is_short

                        ## we are now going to save proteome information, but only if mode_good if args.mode_good_only

                        if ((not mode_good) and args.mode_good_only):
                            continue

                        bad_prot_info = this_prot_info.copy()
                        if (this_proteome not in bad_proteomes):
                            bad_proteomes[this_proteome] = { cluster_id : bad_prot_info }
                            bad_proteome_cnt[this_proteome] = {'total':1, 'short':0, 'long':0}
                        else:
                            bad_proteome_cnt[this_proteome]['total'] += 1
                            if (cluster_id not in bad_proteomes[this_proteome]):
                                bad_proteomes[this_proteome][cluster_id] = bad_prot_info

                        if (is_short): 
                            bad_proteome_cnt[this_proteome]['short'] += 1
                            if (this_proteome not in bad_proteomes_short):
                                bad_proteomes_short[this_proteome] = { cluster_id : bad_prot_info }
                            else:
                                if (cluster_id not in bad_proteomes_short[this_proteome]):
                                    bad_proteomes_short[this_proteome][cluster_id] = bad_prot_info
                            if (cluster_id not in bad_clusters_short):
                                bad_clusters_short[cluster_id] = [bad_prot_info]
                            else:
                                bad_clusters_short[cluster_id].append(bad_prot_info)

                        elif (mode_len >= args.target_length):    ## is_long and mode_len >= 100
                            bad_proteome_cnt[this_proteome]['long'] += 1
 
                            ## save bad_long proteome info
                            if (this_proteome not in bad_proteomes_long):
                                bad_proteomes_long[this_proteome] = { cluster_id : bad_prot_info }
                            else:
                                if (cluster_id not in bad_proteomes_long[this_proteome]):
                                    bad_proteomes_long[this_proteome][cluster_id] = bad_prot_info

                            if (len(bad_proteomes_long[this_proteome]) != bad_proteome_cnt[this_proteome]['long']):
                                print(f"*** mismatch bad_cnt {this_proteome} len: {len(bad_proteomes_long[this_proteome])} cnt: {bad_proteome_cnt[this_proteome]['long']}",file=sys.stderr)

                            ## save bad_long cluster info
                            if (cluster_id not in bad_clusters_long):
                                bad_clusters_long[cluster_id] = [ bad_prot_info ]
                            else:
                                bad_clusters_long[cluster_id].append(bad_prot_info)

                    ## done with for loop through clust_info

                ## all done with this list of proteins
                ## tot_cnt = int(clust_data[proteins_cnt_idx])
                tot_cnt = clust_prot_cnt

                ## if we have some bad_clusters_long, then save some mode proteins/proteomes
                ## but this can only be done after the cluster has been fulling examined
                if (long_tfx_clust_fd and cluster_id in bad_clusters_long):
                    if (len(mode_prot_ids) > args.long_tfx_ex):
                        long_mode_prot_ids = random.sample(mode_prot_ids,args.long_tfx_ex)
                    else:
                        long_mode_prot_ids = [x for x in mode_prot_ids]

                    bad_clusters_long_modes[cluster_id] = long_mode_prot_ids

                if (bad_cnt > 0):
                    bad_fract_short = bad_cnt_short/bad_cnt
                    bad_SL_flag = 'Short'
                    if (bad_fract_short < 0.5):
                        bad_SL_flag = 'Long'
                else:
                    bad_fract_short = 0.0
                    bad_SL_flag = 'None'

                clust_stats[cluster_id] = {
                    'tot_cnt': tot_cnt,
                    'bad_cnt': bad_cnt,
                    'bad_cnt_short': bad_cnt_short,
                    'bad_cnt_long': bad_cnt_long,
                    'bad_fract' : float(bad_cnt)/float(tot_cnt),
                    'prot_cnt': clust_prot_cnt,
                    'bad_fract_short': bad_fract_short,
                    'bad_SL_flag' : bad_SL_flag,
                    'ome_cnt': clust_prot_cnt,
                    'mode_proteome': mode_proteome,
                    'mode_id': mode_id,
                    'mode_len': mode_len,
                    'mode_cnt': mode_cnt,
                    'width' : clust_width,
                    'q1_len' : clust_q1['len'],
                    'med_len' : clust_med['len'],
                    'q3_len' : clust_q3['len']
                }

            ## get the distribution of clust_stats
            clust_keys_sort = sorted(clust_stats.keys(), key = lambda x : -clust_stats[x]['bad_fract'])

            print("#! "+' '.join(sys.argv))
            print("## %d clusters [%d total, %d bad_mode] with proteins < %.2f or > %.2f representative length"%(bad_clust, n_clust, bad_mode_cnt, args.thresh_fract,1.0/args.thresh_fract))

            bad95 = 0
            bad10 = 0
            bad95_printed = False
            
            ## print out information on clusters
            for cluster in clust_keys_sort:
                this_cluster = clust_stats[cluster]

                if (this_cluster['bad_fract'] > 0.95):
                    bad95 += 1
                    continue

                if (not bad95_printed):
                    print("## %d clusters with >95%% of proteins outside length range"%(bad95),file=out_clust_fd)
                    print('\t'.join(("clust_id","n_omes","p_tot","p_bad","pct_bad","p_short","pct_short","bad_SL_flag","mode_len","p_mode","pct_mode","q1_len","med_len","q3_len")),file=out_clust_fd)
                    bad95_printed=True

                if (this_cluster['bad_fract'] < args.bad_low_fract):
                    bad10 += 1
                    continue

                if (this_cluster['tot_cnt'] > 0):
                    mode_fract = float(this_cluster['mode_cnt'])/this_cluster['tot_cnt']
                else:
                    mode_fract_ = 0.0

                bad_fract_str =  "%.2f"%(100.0*this_cluster['bad_fract'])
                if (0 < this_cluster['bad_fract'] < 0.01):
                    bad_fract_str =  "%.3g"%(100.0*this_cluster['bad_fract'])

                bad_fract_short_str = "%.2f"%(100.0*this_cluster['bad_fract_short'])
                if (0 < 100.0*this_cluster['bad_fract_short'] < 0.01):
                    bad_fract_short_str =  "%.3g"%(100.0*this_cluster['bad_fract_short'])
                
                print('\t'.join((cluster,str(this_cluster['ome_cnt']),
                                 str(this_cluster['tot_cnt']),
                                 str(this_cluster['bad_cnt']),
                                 bad_fract_str,
                                 str(this_cluster['bad_cnt_short']),
                                 bad_fract_short_str,
                                 this_cluster['bad_SL_flag'],
                                 str(this_cluster['mode_len']),
                                 str(this_cluster['mode_cnt']),
                                 "%.2f"%(100.0*mode_fract),
                                 str(this_cluster['q1_len']),
                                 str(this_cluster['med_len']),
                                 str(this_cluster['q3_len'])
                               )),file=out_clust_fd)

            n_bad_clust = len(bad_proteome_cnt.keys())

            print("## %d clusters with < %.2f%% of proteins outside length range"%(bad10,100.0*args.bad_low_fract),file=out_clust_fd)
            print("")
            ## end of cluster stats

            ## beginning of proteome stats summary
            print("## proteomes with most short/long proteins (%d total) in %d clusters"%(n_bad_clust,n_clust),file=out_omes_fd)
            clust_values = [bad_proteome_cnt[x]['total'] for x in bad_proteome_cnt.keys()]
            print("## cluster min: %d median: %.1f max: %d"%(min(clust_values), stat.median(clust_values), max(clust_values)),file=out_omes_fd)
            print("## cluster quantiles: ",end='',file=out_omes_fd)
            print("## ",end='',file=out_omes_fd)
            print(stat.quantiles(clust_values,n=10),file=out_omes_fd)

            ## beginning of proteome count data
            print("\t".join(('proteome_id','n_clusters','ns_clusters','pct_cluster','pct_short')), end='',file=out_omes_fd)
            if (have_busco):
                print("\t"+"\t".join(busco_fields_out),end='',file=out_omes_fd)
            print("\ts_type",file=out_omes_fd)

            ## restructure this section to
            ## (1) show a large number of the bad proteomes
            ## (1a) write out a sample of bad-proteins from bad proteomes
            ## (2) show a random sample of the remaining proteomes (perhaps 1000)
            ## (2a) write out a sample of bad-proteins from sampled proteomes

            bad_proteomes_keys = list(bad_proteomes.keys())
            bad_proteomes_keys.sort(key = lambda x : -bad_proteome_cnt[x]['total'])

            ## here, it would make sense to do some sampling 
            ## need 100% of the args.bad_prots

            ## print out the bad proteome
            max_bad_omes = min(len(bad_proteomes_keys),args.bad_prots)
            for this_ome_key in bad_proteomes_keys[:max_bad_omes]:
                print_proteomes(this_ome_key, bad_proteome_cnt, n_clust, proteome_map, have_busco, busco_fields_out, 'B', out_omes_fd)

            ## also print out bad_proteome_examples for tfastx

            ## here we print out two types of short protein/proteome data:
            ## (1) the worst proteomes (those with the most short/long
            ##     proteins) then divided by highest/lowest N50
            ## (2) a sample from all bad proteomes

            ## if have_busco, do some n50 samples
            ## to label by n50, need to sort tfx_samp_omes by n50

            bad_proteomes_short_keys = [x for x in bad_proteomes_keys if x in bad_proteomes_short]
            bad_proteomes_short_keys.sort(key = lambda x : -bad_proteome_cnt[x]['short'])

            if (have_busco and args.do_N50):
                bad_omes_n50 = {}
                ## bad_proteomes_short_keys has been sorted by the largest number of short proteins
                for tfx_ome in bad_proteomes_short_keys[:args.bad_prots]:
                    if (tfx_ome in proteome_map):
                        try:
                            bad_omes_n50[tfx_ome] = float(proteome_map[tfx_ome]['A_n50'])
                        except ValueError:
                            print(proteome_map[tfx_ome],file=sys.stderr)
                            proteome_map[tfx_ome]['A_n50'] = 1
                            bad_omes_n50[tfx_ome] = 0.0
                            
                bad_omes_n50_keys = list(bad_omes_n50.keys())
                bad_omes_n50_keys.sort(key=lambda x: bad_omes_n50[x])

                if (len(bad_omes_n50_keys) > args.bad_tfx_ex):
                    ## now take half from top, half from bottom
                    tfx_samp_short_omes = []
                    for ix in range(0,int(args.bad_tfx_ex/2)):
                        tfx_samp_short_omes.append(bad_omes_n50_keys[ix])

                    print_tfx_ex_list(tfx_samp_short_omes,
                                      bad_proteomes_short, bad_proteome_cnt, 'short', proteome_map,
                                      clust_stats, args.tfx_search_cnt,
                                      'B_N50L', args, out_tfx_fd, print_short_tfx_ex)

                    tfx_samp_short_omes = []
                    bad_omes_N = len(bad_omes_n50_keys)
                    for ix in range(0, int(args.bad_tfx_ex/2)):
                        tfx_samp_short_omes.append(bad_omes_n50_keys[bad_omes_N-1-ix])
                    
                    print_tfx_ex_list(tfx_samp_short_omes,
                                      bad_proteomes_short, bad_proteome_cnt, 'short', proteome_map,
                                      clust_stats, args.tfx_search_cnt, 'B_N50H', args, out_tfx_fd,
                                      print_short_tfx_ex)

                else:
                    print_tfx_ex_list(bad_omes_n50_keys,
                                      bad_proteomes_short, bad_proteome_cnt, 'short', proteome_map,
                                      clust_stats, args.tfx_search_cnt, 'B_N50', args, out_tfx_fd,
                                      print_short_tfx_ex)

            else:  ## no do_n50
                tfx_samp_short_omes = [x for x in bad_proteomes_short_keys[:args.bad_tfx_ex]]

                if (len(tfx_samp_short_omes) > args.samp_tfx_ex):
                    tfx_samp_short_omes = random.sample(tfx_samp_short_omes,args.samp_tfx_ex)

                print_tfx_ex_list(tfx_samp_short_omes,bad_proteomes_short, bad_proteome_cnt, 'short', proteome_map,
                                  clust_stats, args.tfx_search_cnt, 'B_N50', args, out_tfx_fd,
                                  print_short_tfx_ex)

            ## get a sample of the worst bad proteomes for tfastx searches
            ## these examples should all be short
            tfx_bad_omes = [x for x in bad_proteomes_short_keys[:args.bad_prots]]
            if (len(tfx_bad_omes) > args.bad_tfx_ex):
                tfx_bad_omes = random.sample(tfx_bad_omes,args.bad_tfx_ex)

            ## now do a sample of everything
            if (args.samp_prots > 1 and 
                args.samp_prots < len(bad_proteomes_keys)):
                ## bad_proteomes_good_keys = [x for   ## get bad proteomes with > thresh proteins ]

                sample_bad_proteomes_keys = random.sample(bad_proteomes_keys,args.samp_prots)
            else:
                sample_bad_proteomes_keys = [ x for x in bad_proteomes_keys]

            sample_bad_proteomes_keys.sort(key = lambda x : -bad_proteome_cnt[x]['total'])

            for this_ome_key in sample_bad_proteomes_keys:
                print_proteomes(this_ome_key, bad_proteome_cnt, n_clust, proteome_map, have_busco, busco_fields_out, 'S', out_omes_fd)


            ## get a sample of the sampled bad proteomes for tfastx searches (now only want short examples)
            tfx_samp_short_omes = [x for x in sample_bad_proteomes_keys if x in bad_proteomes_short]
            if (len(tfx_samp_short_omes) > args.samp_tfx_ex):
                tfx_samp_short_omes = random.sample(tfx_samp_short_omes,args.samp_tfx_ex)

            print_tfx_ex_list(tfx_samp_short_omes,bad_proteomes_short, bad_proteome_cnt, 'short', proteome_map,
                              clust_stats, args.tfx_search_cnt, 'S_N50', args, out_tfx_fd,
                              print_short_tfx_ex)

            ## we have now printed out all the short_tfx_ex examples, also print out long_tfx_ex examples

            ## check match of bad_proteomes_long with bad_proteomes_cnt

            for bp in bad_proteomes_long.keys():
                if (len(bad_proteomes_long[bp]) != bad_proteome_cnt[bp]['long']):
                    print(f"**2 mismatch bad_cnt {this_proteome} len: {len(bad_proteomes_long[this_proteome])} cnt: {bad_proteome_cnt[this_proteome]['long']}",file=sys.stderr)

            if (long_tfx_ome_fd):
                bad_proteomes_long_keys = list(bad_proteomes_long.keys())
                bad_proteomes_long_keys = [ x for x in bad_proteomes_long_keys if bad_proteome_cnt[x]['long'] >= 5]

                bad_proteomes_long_keys.sort(key = lambda x : -bad_proteome_cnt[x]['long'])

                tfx_samp_long_omes = [x for x in bad_proteomes_long_keys[:args.bad_prots]]
                if (len(tfx_samp_long_omes) > args.bad_tfx_ex):
                    tfx_samp_long_omes = random.sample(tfx_samp_long_omes, args.bad_tfx_ex)
                    tfx_samp_long_omes.sort(key = lambda x : -bad_proteome_cnt[x]['long'])

                    for bp in tfx_samp_long_omes:
                        if (len(bad_proteomes_long[bp]) != bad_proteome_cnt[bp]['long']):
                            print(f"**3 mismatch bad_cnt {this_proteome} len: {len(bad_proteomes_long[this_proteome])} cnt: {bad_proteome_cnt[this_proteome]['long']}",file=sys.stderr)

                    print_tfx_ex_list(tfx_samp_long_omes, bad_proteomes_long, bad_proteome_cnt, 'long', proteome_map,
                                      clust_stats, args.tfx_search_cnt, 'L_N50', args, long_tfx_ome_fd,
                                      print_long_tfx_ex)

            ## finally, print information on the long cluster members
            ## get a list of proteomes with long examples
            ## have_bad_clusters_long[cluster_id] = {'ome':bad_proteome, 'prot':bad_prot_info}

            if (long_tfx_clust_fd):
                bad_clusters_long_keys = list(bad_clusters_long.keys())
                bad_clusters_long_keys.sort(key = lambda x : -len(bad_clusters_long[x]))
                nex_bad_long = min(len(bad_clusters_long_keys),args.long_prots)

                n_bad_long = 0
                for bad_long in bad_clusters_long_keys:
                    this_bad_cluster_info = bad_clusters_long[bad_long]
                    if (len(this_bad_cluster_info) <= 5):
                        break

                    this_bad_mode_info = bad_clusters_long_modes[bad_long]
                    if this_bad_mode_info[0]['len'] <= args.target_length:
                        continue

                    this_bad_cluster_info.sort(key=lambda x: -x['len'])

                    print_long_clust_ex(bad_long, this_bad_cluster_info, this_bad_mode_info, proteome_map, clust_stats[bad_long], 'L_N50', args,long_tfx_clust_fd)

                    n_bad_long += 1
                    if (n_bad_long >= nex_bad_long):
                        break

if __name__ == "__main__":

    main()
