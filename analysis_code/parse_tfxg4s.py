#!/usr/bin/env python3

## (1-July-2025)
## modified parse_tfxg3s.py to incorporate mode vs short pairwise ss

## read an ome_GCA_X.tfxg file
## also read ome_GCA_X_short.ss_md10 file to get overlap info

## for each query, identify if full length alignment (with frameshifts)

## if not full length, report >95% matches
## report number of searches, number full-length

## parse_tfxg2s.py sorts results by N50, if asked

## _1fsT and _NfsT counts are NOT included in _fs and _T counts

import fileinput
import argparse
import copy
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
        "--s_pid",
        dest='s_pid',
        action='store_true',
        default=False
        )

    parser.add_argument(
        "--tfx_ex",
        dest='tfx_ex',
        type=str,
        default=None
        )

    parser.add_argument(
        "--no_X0",
        dest='no_X0',
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

def process_hit(hit_list, args, ss_align_info):

    ## first check for full-length alignment

    cluster_good_1 = 0
    good_1_len = 0

    cluster_good_N = 0
    good_N_len = 0

    q_len_thresh = int(0.95 * hit_list[0]['q_len'])

    for hit in hit_list:
        q_seqid = hit['q_seqid']
        if hit['alen'] >= q_len_thresh:
            cluster_good_1 += 1

            term_fs_status = 0
            btop = hit['BTOP']
            if (re.search(r'\*',btop)):
                term_fs_status += 1
            if (re.search(r'\\',btop) or re.search(r'/',btop)):
                term_fs_status += 2

            if args.debug:
                print(q_seqid,hit['s_seqid'],hit['q_len'], hit['s_len'],hit['alen'], hit['percid'], term_fs_status, end='')
                if (q_seqid in ss_align_info):
                    ss_info = ss_align_info[q_seqid]
                    print("\tss:",ss_info['percid'],ss_info['N_ext_len'],ss_info['N_ext'],ss_info['C_ext_len'],ss_info['C_ext'])
                else:
                    print()

            return (q_seqid,'1', term_fs_status, hit['s_len'])

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
        q_seqid = lap_hit['q_seqid']
        btop = lap_hit['BTOP']
        if (re.search(r'\*',btop)):
            term_fs_status += 1
        if (re.search(r'\\',btop) or re.search(r'/',btop)):
            term_fs_status += 2

        if args.debug:
            print("lap",q_seqid,lap_hit['s_seqid'],lap_hit['q_len'], lap_hit['alen'],lap_hit['percid'],lap_hit['q_start'],lap_hit['q_end'], term_fs_status, end='')
            if (q_seqid in ss_align_info):
                ss_info = ss_align_info[q_seqid]
                print("\tss:",ss_info['percid'],ss_info['N_ext_len'],ss_info['N_ext'],ss_info['C_ext_len'],ss_info['C_ext'])
            else:
                print();

    len_return = tot_slen/N_alen
    if (args.max_slen):
        len_return = max_slen

    if (tot_alen >= q_len_thresh):
        return(q_seqid, 'N', term_fs_status, len_return)
    else:
        if args.debug:
            q_seqid=hit_list[0]['q_seqid']
            print("0H",q_seqid,end=''),
            if (q_seqid in ss_align_info):
                ss_info = ss_align_info[q_seqid]
                print("\tss:",ss_info['percid'],ss_info['N_ext_len'],ss_info['N_ext'],ss_info['C_ext_len'],ss_info['C_ext'])
            else:
                print();

        return(q_seqid, '0', term_fs_status, len_return) 

def process_ss(hit_list, debug):

    ## always only one hit in _short.ss_MD10

    hit = hit_list[0]
    if (len(hit_list)>1):
        sys.stderr.write(f"hit_list > 1 {hit['q_seq_id']} {hit['s_seq_id']}\n")
        return (hit['q_seqid'], -1, -1, -1, 0.0)

    N_term_ext = hit['q_start'] - 1
    C_term_ext =hit['q_len'] - hit['q_end']

    delta = int(float(min(hit['q_len'], hit['s_len']))/50.0 + 0.5)

##    if debug:
##        print("ss",hit['q_seqid'],hit['s_seqid'],N_term_ext, C_term_ext, delta)

    return(hit['q_seqid'], N_term_ext, C_term_ext, delta, hit['percid'])


def read_hit(line, tfx_fields, i_fields, f_fields, percid_min, ss_align_info):

    this_hit = dict(zip(tfx_fields, line.strip('\n').split('\t')))

    for f in i_fields:
        this_hit[f] = int(this_hit[f])

    for f in f_fields:
        this_hit[f] = float(this_hit[f])

    if (this_hit['percid'] >= percid_min):
        return this_hit
##    elif (this_hit['percid'] >= 75.0 and (this_hit['q_seqid'] in ss_align_info and ss_align_info[this_hit['q_seqid']]['percid'] > percid_min)):
##        return this_hit
    else:
        return None

def read_ss_file(tfx_file, tfx_fields, i_fields, f_fields, args):

    ss_align_info = {}

    tfx_dir_name = tfx_file.split('/')
    if (len(tfx_dir_name) > 1):
        tfx_dir_pref = f'{tfx_dir_name[0]}/'
        tfx_file_name = tfx_dir_name[1]
    else:
        tfx_dir_pref = ''
        tfx_file_name = tfx_dir_name[0]

    tfx_file_parts = tfx_file_name.split('_')
    ss_MD10_prefix = tfx_dir_pref+'_'.join(tfx_file_parts[:3])

    with open(ss_MD10_prefix+"_short.ss_MD10",'r') as infile:
        in_hit=False
        hit_list = []

        ## for _short.ss_MD10 files, there will be many queries each with one hit

        clust_C_ext = clust_N_ext = 0
        for line in infile:

            if line.startswith("# "):
                ## if already read a list of results, process and record
                if (in_hit):
                    clust_C_ext = clust_N_ext = 0
                    if (len(hit_list) > 0):
                        (q_seqid, N_ext_len, C_ext_len, delta, percid) = process_ss(hit_list,args.debug)
                        if (N_ext_len > delta):
                            clust_N_ext += 1
                        if (C_ext_len > delta):
                            clust_C_ext += 1

                        ss_align_info[q_seqid] = {'percid':percid, 'N_ext': clust_N_ext, 'N_ext_len':N_ext_len, 'C_ext': clust_C_ext, 'C_ext_len':C_ext_len}

                in_hit = False
                hit_list = []
                continue

            this_hit = read_hit(line, tfx_fields, i_fields, f_fields, args.percid_min,  ss_align_info)
            if this_hit:
                hit_list.append(this_hit)

            ## process the last line if it exists
            clust_C_ext = clust_N_ext = 0
            if (len(hit_list) > 0):
                (q_seqid, N_ext_len, C_ext_len, delta,percid) = process_ss(hit_list,args.debug)
                if (N_ext_len > delta):
                    clust_N_ext += 1
                if (C_ext_len > delta):
                    clust_C_ext += 1
                ss_align_info[q_seqid] = {'percid': percid, 'N_ext': clust_N_ext, 'N_ext_len':N_ext_len, 'C_ext': clust_C_ext, 'C_ext_len':C_ext_len}

    return(ss_align_info)

def main():

    args = check_args()

    tfx_fields = 'q_seqid q_len s_seqid s_len percid alen mism gaps q_start q_end s_start s_end expect bits BTOP'.split(' ')
    f_fields = 'percid expect bits'.split(' ')
    i_fields = 'q_len s_len alen mism gaps q_start q_end s_start s_end'.split(' ')

    if (len(sys.argv) > 10):
        sys_argv_str = "# " + ' '.join(sys.argv[:9]) + ' ... ' + sys.argv[-1]
    else:
        sys_argv_str = "# "+' '.join(sys.argv)

    have_tfx_ex = False
    if (args.tfx_ex):
        tfx_gca_set = set()
        have_tfx_ex = True
        with open(args.tfx_ex,'r') as ex_fd:
            for line in ex_fd:
                if (line.startswith('>')):
                    (gca_acc, cnt, stype) = line[1:].strip('\n').split('\t')
                    tfx_gca_set.add(gca_acc)

    all_hdrs = "file logN50 X1cnt X1fs X1T X1fsT X1len X1Next X1Cext Ncnt Nfs NT NfsT Nlen NNext NCext X0cnt X0len X0Next X0Cext tot_hit tot_cnt pct_good pct_1good pct_fs".split(' ')
    noX0_hdrs = "file logN50 X1cnt X1fs X1T X1fsT X1len X1Next X1Cext Ncnt Nfs NT NfsT Nlen NNext NCext X0cnt tot_hit tot_cnt pct_good pct_1good pct_fs".split(' ')

    if (not args.summ_flg):
        print(sys_argv_str)
        if (not args.no_X0):
            print('\t'.join(all_hdrs))
        else:
            print('\t'.join(noX0_hdrs))

    tfxg_stats = []
    cluster_total = 0

    fs_status_str = ('0','T','fs','Tfs')

    os_code=args.files[0].split('_')[0]

    for tfx_file in args.files:

        if (have_tfx_ex):
            file_parts = tfx_file.split('_')
            gca_part = '_'.join(file_parts[1:3])

            if (gca_part not in tfx_gca_set):
                continue

        ## as originally designed, this section calculates the
        ## statistics for the entire proteome (with multiple
        ## clusters).  That does not work when we need to merge the
        ## tfx results with the _short.ss_MD10 results, because the
        ## details of the alignments have been lost.

        ## need to save all the alignment info, so it can be merged
        ## alternatively, read the ss_MD10 info at the same time, so that it is available to be merged
        ## since it is shorter, why not grab the .ss_MD10 stuff and merge it with this file

        
        ss_align_info = read_ss_file(tfx_file, tfx_fields, i_fields, f_fields, args)

        with open(tfx_file,'r') as infile:

            in_hit=False
            cluster_cnt = 0
            cluster_good_1tot = 0
            cluster_good_1fs = 0
            cluster_good_1T = 0
            cluster_good_1fsT = 0

            no_cluster_cnt = 0

            cluster_good_Ntot = 0
            cluster_good_Nfs = 0
            cluster_good_NT = 0
            cluster_good_NfsT = 0

            cluster_bad_Ntot = 0
            cluster_1slen = 0
            cluster_Nslen = 0
            cluster_0slen = 0
            
            cluster_1Next = 0
            cluster_1Cext = 0

            cluster_NNext = 0
            cluster_NCext = 0

            cluster_0Next = 0
            cluster_0Cext = 0

            for line in infile:

                if line.startswith("# "):
                    if (in_hit):
                        if (len(hit_list) > 0):
                            (q_seqid, align_type, fs_status, s_len) = process_hit(hit_list,args, ss_align_info)

                            if (q_seqid in ss_align_info):
                                N_ext = ss_align_info[q_seqid]['N_ext']
                                C_ext = ss_align_info[q_seqid]['C_ext']
                            else:
                                sys.stderr.write(f"no ss_align_info[{q_seqid}] {tfx_file}\n")
                                N_ext = 0
                                C_ext = 0

                            if (align_type == '1'):
                                cluster_good_1tot += 1
                                if (fs_status & 2):
                                    cluster_good_1fs += 1
                                if (fs_status & 1):
                                    cluster_good_1T += 1
                                if (fs_status == 3):
                                    cluster_good_1fsT += 1

                                cluster_1Next += N_ext
                                cluster_1Cext += C_ext

                                cluster_1slen += s_len

                            elif (align_type == 'N'):
                                cluster_good_Ntot += 1

                                if (fs_status & 1):
                                    cluster_good_Nfs += 1
                                if (fs_status & 2):
                                    cluster_good_NT += 1
                                if (fs_status == 3):
                                    cluster_good_NfsT += 1

                                cluster_NNext += N_ext
                                cluster_NCext += C_ext

                                cluster_Nslen += s_len
                            elif (align_type == '0'):
                                cluster_bad_Ntot += 1
                                cluster_0slen += s_len

                                cluster_0Next += N_ext
                                cluster_0Cext += C_ext

                            else:
                                sys.stderr.write("Uncounted result\m")

                    in_hit = False
                    continue
                
                if not in_hit:
                    cluster_cnt += 1
                    hit_list = []

                in_hit = True
                this_hit = read_hit(line, tfx_fields, i_fields, f_fields, args.percid_min, ss_align_info)
                
                if (this_hit):
                    hit_list.append(this_hit)

            ## process the last line if it exists
            if (in_hit):
                if (len(hit_list) > 0):
                    short_hit = copy.copy(hit_list[0])
                    (q_seqid, align_type, fs_status, s_len) = process_hit(hit_list,args, ss_align_info)

                    if (align_type == '1'):
                        cluster_good_1tot += 1
                        cluster_1slen += s_len
                    elif (align_type == 'N'):
                        cluster_good_Ntot += 1
                        cluster_Nslen += slen
                    elif (align_type == '0'):
                        cluster_bad_Ntot += 1
                        cluster_0slen += s_len
                    else:
                        sys.stderr.write("Uncounted result\m")

        if (cluster_cnt > 0): 
            tot_good = cluster_good_1tot + cluster_good_Ntot
            if (tot_good > 0):
                fract_1tot_good = cluster_good_1tot/float(tot_good)
            else:
                fract_1tot_good = 0.0

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
                logN50 = re.search(r':(\d+\.\d+).tfx',tfx_file).group(1)
                logN50 = float(logN50)
            else:
                logN50 = -1.0

            this_result = {'file':tfx_file,
                           'good_1tot':cluster_good_1tot,'good_1fs':cluster_good_1fs,'good_1T':cluster_good_1T,'good_1fsT':cluster_good_1fsT,
                           'good_1len': cluster_slen_1ave, 'good_1Next':cluster_1Next, 'good_1Cext':cluster_1Cext,
                           'good_Ntot':cluster_good_Ntot,'good_Nfs':cluster_good_Nfs,'good_NT':cluster_good_NT,'good_NfsT':cluster_good_NfsT,
                           'good_Nlen': cluster_slen_Nave, 'good_NNext': cluster_NNext, 'good_NCext' : cluster_NCext,
                           'bad_Ntot':cluster_bad_Ntot, 'bad_0len': cluster_slen_0ave,
                           'bad_0Next' : cluster_0Next, 'bad_0Cext' : cluster_0Cext,
                           'tot_cnt':cluster_cnt, 'tot_good':tot_good,
                           'fract_1tot_good':fract_1tot_good, 'logN50':logN50}

            tfxg_stats.append(this_result)

            ## print('\t'.join([tfx_file, str(cluster_good_1tot),str(cluster_good_Ntot), str(cluster_cnt),'%.1f'%(100.0*tot_good/cluster_cnt),'%.1f'%(100.0*fract_1tot_good)]))

    ## here we have read all the tfxg_l files, now we want to read
    ## all the .ss_MD10 files and merge the results

##    for tfx_result in tfxg_stats:

    if (tfxg_stats[0]['logN50'] > 0):
        tfxg_stats.sort(key = lambda x : x['logN50'])

    tot_bad_Ntot = 0
    tot_good_1tot = 0
    tot_good_Ntot = 0
    tot_good=0

    tot_fs = 0
    tot_Term = 0
    tot_NfsT = 0
    tot_fsT = 0

    tot_N_ext = 0
    tot_C_ext = 0

    for r in tfxg_stats:
        if (not args.summ_flg):
            if (r['tot_good']>0): 
                ## to get total frame-shifts, need _fs + _fsT
                pct_fs = (r['good_1fs'] + r['good_Nfs'] + r['good_1fsT'] +r['good_NfsT'])/r['tot_good']
            else:
                pct_fs = 0.0

            d_file = r['file']
            if (args.s_pid):
                d_file_list = d_file.split('_')
                d_file = '_'.join(d_file_list[1:3])

            if (not args.no_X0):
                print('\t'.join([d_file,'%.2f'%(r['logN50']),
                                 str(r['good_1tot']),str(r['good_1fs']),str(r['good_1T']),str(r['good_1fsT']),'%.0f'%(r['good_1len']),
                                 str(r['good_1Next']), str(r['good_1Cext']),
                                 str(r['good_Ntot']),str(r['good_Nfs']),str(r['good_NT']),str(r['good_NfsT']),'%.0f'%(r['good_Nlen']),
                                 str(r['good_NNext']), str(r['good_NCext']),
                                 str(r['bad_Ntot']),'%.0f'%(r['bad_0len']),
                                 str(r['bad_0Next']), str(r['bad_0Cext']),
                                 str(r['tot_good']),str(r['tot_cnt']),
                                 '%.1f'%(100.0*r['tot_good']/r['tot_cnt']),
                                 '%.1f'%(100.0*r['fract_1tot_good']),
                                 '%.1f'%(100.0*pct_fs)]))
            else:
                print('\t'.join([d_file,'%.2f'%(r['logN50']),
                                 str(r['good_1tot']),str(r['good_1fs']),str(r['good_1T']),str(r['good_1fsT']),'%.0f'%(r['good_1len']),
                                 str(r['good_1Next']), str(r['good_1Cext']),
                                 str(r['good_Ntot']),str(r['good_Nfs']),str(r['good_NT']),str(r['good_NfsT']),'%.0f'%(r['good_Nlen']),
                                 str(r['good_NNext']), str(r['good_NCext']),
                                 str(r['bad_Ntot']),
                                 str(r['tot_good']),str(r['tot_cnt']),
                                 '%.1f'%(100.0*r['tot_good']/r['tot_cnt']),
                                 '%.1f'%(100.0*r['fract_1tot_good']),
                                 '%.1f'%(100.0*pct_fs)]))

        tot_bad_Ntot += r['bad_Ntot']
        tot_good_1tot += r['good_1tot']
        tot_good_Ntot += r['good_Ntot']
        tot_good += r['tot_good']

        tot_N_ext += r['good_NNext'] + r['good_1Next']
        tot_C_ext += r['good_NCext'] + r['good_1Cext']

        tot_fs += r['good_1fs'] + r['good_Nfs']
        tot_Term += r['good_1T'] + r['good_NT']
        tot_fsT += r['good_1fsT'] + r['good_NfsT']

        no_cluster_cnt += r['tot_cnt'] - r['tot_good']

        tot_NfsT = tot_good_1tot + tot_good_Ntot - (tot_fs + tot_Term + 2*tot_fsT)

    if (not args.summ_flg):
        print("# "+"\t".join(("oscode","0cnt","1cnt","Ncnt","fs","Term","fsT","n_fsT","Next","Cext","hits","no_hit","tot_omes")))
        print("# ",end='')

    print(f"{os_code}\t{tot_bad_Ntot}\t{tot_good_1tot}\t{tot_good_Ntot}\t{tot_fs}\t{tot_Term}\t{tot_fsT}\t{tot_NfsT}\t{tot_N_ext}\t{tot_C_ext}\t{tot_good}\t{no_cluster_cnt}\t{len(tfxg_stats)}")

if __name__ == "__main__":

    main()
