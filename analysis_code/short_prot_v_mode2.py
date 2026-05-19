#!/usr/bin/env python3

## modified 18-April-2026 to include fa_dir: fasta10_incsurv
##

## short_prot_v_mode2.py --fa_dir fasta10_incsurv --oscode ECOLA ECOLA_c.....tfx_ex2
## take a file of the form

## a derivative of long_protein_v_mode.py, that works with a slightly different format file:

## proteome_id    n_short prots in proteome s_type
## >GCA_025660355.1	444	B_N50L:3.06
## cluster_id   mode_prot_id    mode_len bad_prot_id
## 1045608	AAN78524|388	388	MCU8648263|234
## 313816	AAN81807|282	282	MCU8651106|169
## 2043279	CTZ65892|237	237	MCU8651645|156
## 63207	AAN79527|248	248	MCU8651910|151
## 2016381	AAN83700|251	251	MCU8647279|186

## we want to run each mode_protein against the bad_prot_id, label output by cluster

import fileinput
import sys
import os
import subprocess
import argparse

BIN_DIR=os.getenv('ETT_BIN')

if BIN_DIR is None:
    print(" ERROR: ETT_BIN not defined. Have you sourced env.sh?", file=sys.stderr)
    sys.exit(1)


def check_args(test_args=None):
    """
    parse arguments and check for error conditions
    """

    parser = argparse.ArgumentParser(
        description="parse short_tfx_ex2 file, extract mode_prot, bad_prot, run ssearch"
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
        "--run",
        dest='run_flag',
        action='store_true',
        default=False
        )

    parser.add_argument(
        "--script",
        dest='down_script',
        type=str,
        default='down_fasta_oscode.sh'
        )

    parser.add_argument(
        dest="files",
        type=str,
        nargs='*')

    return parser.parse_args()

def run_query(mode_prot_list, bad_prot_list, cluster_id, fa_dir, taxon,  run_flag, down_script):

    # (1) get the sequences
    if (len(mode_prot_list) <= 0):
        return

    for ix, mode_prot_id in enumerate(mode_prot_list):

        mode_prot_acc = mode_prot_id.split('|')[0]

        mode_file_name = f'{mode_prot_acc}.aa'
        bad_prot_id = bad_prot_list[ix]
        bad_prot_acc = bad_prot_id.split('|')[0]

        bad_file_name = f'{bad_prot_acc}.aa'

        ## extract the mode protein
        have_mode_file = os.path.exists(mode_file_name) and (os.path.getsize(mode_file_name)>0)
        if (not have_mode_file):
            cmd_str = f'{BIN_DIR}/{down_script} {fa_dir} {taxon} \'{mode_prot_id}\' > {mode_file_name}'
            if (run_flag):
                os.system(cmd_str)
            else:
                print(cmd_str)

        ## extract the bad protein
        have_bad_file = os.path.exists(bad_file_name) and (os.path.getsize(bad_file_name)>0)
        if (not have_bad_file):
            cmd_str = f'{BIN_DIR}/{down_script} {fa_dir} {taxon} \'{bad_prot_id}\' > {bad_file_name}'
            if (run_flag):
                os.system(cmd_str)
            else:
                print(cmd_str)

        cmd_list = [f'{BIN_DIR}/short_prot2_mode_ss.sh',mode_file_name,bad_file_name,cluster_id,taxon]

        if (run_flag):
            subprocess.run(cmd_list)
        else:
            print(' '.join(cmd_list))


def main():

    args = check_args()

    mode_prot_list = []
    bad_prot_list = []

    mode_prot_acc = ''
    gca_ome_id = ''

    for line in fileinput.input(files=args.files):
        if line.startswith('>'):
            ## if have previous set, generate files, do search
            if (gca_ome_id):
                run_query(mode_prot_list, bad_prot_list, gca_ome_id, args.fa_dir, args.taxon, args.run_flag, args.down_script)

            (gca_ome_id, n_clust, s_type) = line[1:].strip('\n').split('\t')

            mode_prot_list = []
            bad_prot_list = []
            continue

        (cluster_id, mode_acc, mode_len, bad_acc) = line.strip('\n').split('\t')

        cluster_str = f'cl_{cluster_id}'
        mode_prot_list.append(mode_acc)
        bad_prot_list.append(bad_acc)
        
    if (gca_ome_id):
        run_query(mode_prot_list, bad_prot_list, gca_ome_id, args.fa_dir, args.taxon, args.run_flag, args.down_script)

if __name__ == "__main__":

    main()

