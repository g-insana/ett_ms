#!/usr/bin/env python3

## read a *.ok2_MD10 file, report: (1) number of queries, (2) number of hits at various thresholds

import fileinput
import argparse
import subprocess

def check_args(test_args=None):
    """
    parse arguments and check for error conditions
    """

    parser = argparse.ArgumentParser(
        description="parse fasta .ok2 file, extract cluster queries, run genome searches"
    )

    parser.add_argument(
        "--percid",
        dest="pct_thresh",
        type=float,
        default=90.0)

    parser.add_argument(
        "--aln",
        dest="aln_thresh",
        type=float,
        default=0.90)

    parser.add_argument(
        dest="files",
        type=str,
        nargs='*')

    return parser.parse_args()

def main():

    args = check_args()

    tfx_fields = 'q_seqid q_len s_seqid s_len percid alen mism gaps q_start q_end s_start s_end expect bits BTOP'.split(' ')
    f_fields = 'percid expect bits'.split(' ')
    i_fields = 'q_len s_len alen mism gaps q_start q_end s_start s_end'.split(' ')

    print("\t".join("file query_cnt hit_cnt hit_pct good_pct pct_pct good_aln aln_pct aln_pct pct_aln_pct".split(' ')))

    for tfx_file in args.files:
        with open(tfx_file,'r') as infile:

            in_hit=False
            query_cnt = 0
            aln_cnt = 0
            pct_good_cnt = 0
            aln_good_cnt = 0
            pct_aln_good_cnt = 0

            for line in infile:

                if line.startswith("# "):
                    if line.startswith("# Query:"):
                        in_hit = False
                        query_cnt += 1
                    continue

                if not in_hit:
                    aln_cnt += 1
                else:
                    continue

                in_hit = True
                this_hit = dict(zip(tfx_fields, line.strip('\n').split('\t')))
                
                for f in i_fields:
                    this_hit[f] = int(this_hit[f])

                for f in f_fields:
                    this_hit[f] = float(this_hit[f])

                if (this_hit['percid'] >= args.pct_thresh):
                    pct_good_cnt += 1

                a_len_thresh = int(args.aln_thresh* min(this_hit['q_len'],this_hit['s_len']))

                if (this_hit['alen'] >= a_len_thresh):
                    aln_good_cnt += 1
                    if (this_hit['percid'] > args.pct_thresh):
                        pct_aln_good_cnt += 1

            ## done with file

        if (query_cnt > 0 ):
            ppct_aln_cnt = 100.0*aln_cnt/query_cnt
        else:
            ppct_aln_cnt = 0.0

        if (aln_cnt > 0):
            ppct_pct_good_cnt = 100.0*pct_good_cnt/aln_cnt
            ppct_aln_good_cnt =  100.0*aln_good_cnt/aln_cnt
            ppct_pct_aln_good_cnt = 100.0*pct_aln_good_cnt/aln_cnt
        else:
            ppct_pct_good_cnt = 0.0
            ppct_aln_good_cnt = 0.0
            ppct_pct_aln_good_cnt = 0.0

        print('\t'.join([tfx_file, str(query_cnt), str(aln_cnt), "%.1f"%(ppct_aln_cnt),
                                                   str(pct_good_cnt), "%.1f"%(ppct_pct_good_cnt),
                                                   str(aln_good_cnt),'%.1f'%(ppct_aln_good_cnt),
                                                   str(pct_aln_good_cnt),'%.1f'%(ppct_pct_aln_good_cnt)]))

if __name__ == "__main__":

    main()
