#!/usr/bin/env python3

## read a genus_clX.tfxg file

## for each query, identify if full length alignment (with frameshifts)

## if not full length, report >95% matches
## report number of searches, number full-length

## parse_tfxg2s.py sorts results by N50, if asked

import fileinput
import argparse
import sys
import re

def check_args(test_args=None):
    """
    parse arguments and check for error conditions
    """

    parser = argparse.ArgumentParser(
        description="parse summ2 file, extract cluster queries, run genome searches"
    )

    parser.add_argument(
        "-D",
        dest='debug',
        action='store_true',
        default=False
        )

    parser.add_argument(
        "-p","--percid",
        dest='percid_min',
        type=float,
        default=90.0
        )

    parser.add_argument(
        "-S","--sort",
        dest='sort_n50',
        action='store_true',
        default=False
        )

    parser.add_argument(
        "-M","--max",
        dest='max_slen',
        action='store_true',
        default=False
        )

    parser.add_argument(
        "--summ",
        dest='summ_flg',
        action='store_true',
        default=False
        )

    parser.add_argument(
        dest="files",
        type=str,
        nargs='*')

    return parser.parse_args()

## go through a list of hits, identify full-length based on a single alignment (cluster_good_1) and
## full length using overlapping query alignments (cluster_good_N)
##

## now returns cluster1, s_size1, clusterN, s_sizeN(ave)

def process_hit(hit_list, args):

    ## first check for full-length alignment

    cluster_good_1 = 0
    good_1_len = 0

    cluster_good_N = 0
    good_N_len = 0

    q_len_thresh = int(0.95 * hit_list[0]['q_len'])

    for hit in hit_list:
        if hit['alen'] >= q_len_thresh:
            if args.debug:
                print(hit['q_seqid'],hit['s_seqid'],hit['q_len'], hit['s_len'],hit['alen'])
            cluster_good_1 += 1

            term_fs_status = 0
            btop = hit['BTOP']
            if (re.search(r'\*',btop)):
                term_fs_status += 1
            if (re.search(r'\\',btop) or re.search(r'/',btop)):
                term_fs_status += 2

            return ('1', term_fs_status, hit['s_len'])


    ## here if no single alignment that works

    ## sort hits left to right in query
    hit_list.sort(key = lambda x : x['q_start'])

    N_alen = 1
    tot_alen = hit_list[0]['alen']
    tot_slen = hit_list[0]['s_len']
    max_slen = tot_slen

    this_hit = hit_list[0]

    overlap_thresh = int(0.05 * this_hit['q_len'] + 0.5)

##    print(this_hit)
    lap_hits = []
    lap_used = set()
    for ix in range(len(hit_list)-1):
        next_hit = hit_list[ix+1]

##        if (abs(next_hit['q_start'] - this_hit['q_end']) < overlap_thresh):
        if (next_hit['q_start'] - this_hit['q_end'] < overlap_thresh):
##            print(next_hit)
##            print(tot_alen)

            if (ix not in lap_used):
                lap_used.add(ix)
                lap_hits.append(this_hit)

            if (ix+1 not in lap_used):
                lap_used.add(ix+1)
                lap_hits.append(next_hit)

            N_alen += 1
            tot_alen += next_hit['alen']
            tot_slen += next_hit['s_len']
            if (next_hit['s_len'] > max_slen):
                max_slen = next_hit['s_len']
            this_hit = next_hit

    ## check for frameshifts, term codons in lap_hits

    term_fs_status = 0
    for lap_hit in lap_hits:
        btop = lap_hit['BTOP']
        if (re.search(r'\*',btop)):
            term_fs_status += 1
        if (re.search(r'\\',btop) or re.search(r'/',btop)):
            term_fs_status += 2

    len_return = tot_slen/N_alen
    if (args.max_slen):
        len_return = max_slen

    if (tot_alen >= q_len_thresh):
        return('N', term_fs_status, len_return)
    else:
        return('0', term_fs_status, len_return) 

def main():

    args = check_args()

    tfx_fields = 'q_seqid q_len s_seqid s_len percid alen mism gaps q_start q_end s_start s_end expect bits BTOP'.split(' ')
    f_fields = 'percid expect bits'.split(' ')
    i_fields = 'q_len s_len alen mism gaps q_start q_end s_start s_end'.split(' ')

    if (len(sys.argv) > 10):
        sys_argv_str = "# " + ' '.join(sys.argv[:9]) + ' ... ' + sys.argv[-1]
    else:
        sys_argv_str = "# "+' '.join(sys.argv)

    print(sys_argv_str)
    print('\t'.join("file logN50 X1cnt X1fs X1T X1fsT X1len Ncnt Nfs NT NfsT Nlen X0cnt X0len tot_hit tot_cnt pct_good pct_1contig pct_frame".split(' ')))

    tfxg_stats = []

    fs_status_str = ('0','T','fs','Tfs')

    os_code=args.files[0].split('_')[0]

    for tfx_file in args.files:
        with open(tfx_file,'r') as infile:

            in_hit=False
            cluster_cnt = 0
            cluster_good_1tot = 0
            cluster_good_1fs = 0
            cluster_good_1T = 0
            cluster_good_1fsT = 0

            cluster_good_Ntot = 0
            cluster_good_Nfs = 0
            cluster_good_NT = 0
            cluster_good_NfsT = 0

            cluster_bad_Ntot = 0
            cluster_1slen = 0
            cluster_Nslen = 0
            cluster_0slen = 0
            
            for line in infile:

                if line.startswith("# "):
                    if (in_hit):
                        if (len(hit_list) > 0):
                            (align_type, fs_status, s_len) = process_hit(hit_list,args)

                            if (align_type == '1'):
                                cluster_good_1tot += 1
                                if (fs_status & 2):
                                    cluster_good_1fs += 1
                                if (fs_status & 1):
                                    cluster_good_1T += 1
                                if (fs_status == 3):
                                    cluster_good_1fsT += 1

                                cluster_1slen += s_len
                            elif (align_type == 'N'):
                                cluster_good_Ntot += 1

                                if (fs_status & 2):
                                    cluster_good_Nfs += 1
                                if (fs_status & 1):
                                    cluster_good_NT += 1
                                if (fs_status == 3):
                                    cluster_good_NfsT += 1

                                cluster_Nslen += s_len
                            elif (align_type == '0'):
                                cluster_bad_Ntot += 1
                                cluster_0slen += s_len

                    in_hit = False
                    continue
                
                if not in_hit:
                    cluster_cnt += 1
                    hit_list = []

                in_hit = True
                this_hit = dict(zip(tfx_fields, line.strip('\n').split('\t')))
                
                for f in i_fields:
                    this_hit[f] = int(this_hit[f])

                for f in f_fields:
                    this_hit[f] = float(this_hit[f])

                if (this_hit['percid'] >= args.percid_min):
                    hit_list.append(this_hit)

            ## process the last line if it exists
            if (in_hit):
                if (len(hit_list) > 0):
                    (align_type, s_len) = process_hit(hit_list,args)

                    if (align_type == '1'):
                        cluster_good_1tot += 1
                        cluster_1slen += s_len
                    elif (align_type == 'N'):
                        cluster_good_Ntot += 1
                        cluster_Nslen += slen
                    elif (align_type == '0'):
                        cluster_bad_Ntot += 1
                        cluster_0slen += s_len

        if (cluster_cnt > 0): 
            tot_good = cluster_good_1tot + cluster_good_Ntot
            if (tot_good > 0):
                fract_1tot_good = cluster_good_1tot/float(tot_good)
                fract_frame = float(cluster_good_1fs + cluster_good_Nfs)/tot_good
            else:
                fract_1tot_good = 0.0
                fract_frame = 0.0

            if (cluster_good_1tot > 0):
                cluster_slen_1ave = cluster_1slen/cluster_good_1tot
            else:
                cluster_slen_1ave = 0

            if (cluster_good_Ntot > 0):
                cluster_slen_Nave = cluster_Nslen/cluster_good_Ntot
            else:
                cluster_slen_Nave = 0

            if (cluster_bad_Ntot > 0):
                cluster_slen_0ave = cluster_0slen/cluster_bad_Ntot
            else:
                cluster_slen_0ave = 0

            if (len(tfx_file.split(':')) > 1):
                logN50 = re.search(r':(\d+\.\d+).t',tfx_file).group(1)
                logN50 = float(logN50)
            else:
                logN50 = -1.0

            this_result = {'file':tfx_file,
                           'good_1tot':cluster_good_1tot,'good_1fs':cluster_good_1fs,'good_1T':cluster_good_1T,'good_1fsT':cluster_good_1fsT,
                           'good_1len': cluster_slen_1ave, 
                           'good_Ntot':cluster_good_Ntot,'good_Nfs':cluster_good_Nfs,'good_NT':cluster_good_NT,'good_NfsT':cluster_good_NfsT,
                           'good_Nlen': cluster_slen_Nave, 
                           'bad_Ntot':cluster_bad_Ntot, 'bad_0len': cluster_slen_0ave, 
                           'tot_cnt':cluster_cnt, 'tot_good':tot_good,
                           'fract_1tot_good':fract_1tot_good,
                           'fract_frame':fract_frame,
                           'logN50':logN50}

            tfxg_stats.append(this_result)

            ## print('\t'.join([tfx_file, str(cluster_good_1tot),str(cluster_good_Ntot), str(cluster_cnt),'%.1f'%(100.0*tot_good/cluster_cnt),'%.1f'%(100.0*fract_1tot_good)]))

    if (tfxg_stats[0]['logN50'] > 0):
        tfxg_stats.sort(key = lambda x : x['logN50'])

    for r in tfxg_stats:

        if (not args.summ_flg):

            print('\t'.join([r['file'],'%.2f'%(r['logN50']),
                             str(r['good_1tot']),str(r['good_1fs']),str(r['good_1T']),str(r['good_1fsT']),'%.0f'%(r['good_1len']),
                             str(r['good_Ntot']),str(r['good_Nfs']),str(r['good_NT']),str(r['good_NfsT']),'%.0f'%(r['good_Nlen']),
                             str(r['bad_Ntot']),'%.0f'%(r['bad_0len']),
                             str(r['tot_good']),str(r['tot_cnt']),
                             '%.1f'%(100.0*r['tot_good']/r['tot_cnt']),
                             '%.1f'%(100.0*r['fract_1tot_good']),
                             '%.1f'%(100.0*r['fract_frame'])
                             ]))

if __name__ == "__main__":

    main()
