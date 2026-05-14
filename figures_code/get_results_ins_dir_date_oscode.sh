#!/bin/sh


if [[ $# < 2 ]]; then
    echo $0 outdir9_incsurv_cm0/x 2026-05-13 OSCODES1 OSCODES2 ...
    exit 1
fi

if [[ -z $1 || "$1" == "x" ]]; then
    data_dir='outdir9_incsurv_cm0'
else
    data_dir=$1
fi

shift

date_str=$1
shift

r_dir=/homes/pearson/bs_tmp/ett_ms/$data_dir
i_dir=/homes/insana/ett2/$data_dir

scp -p pearson@codon_m:${i_dir}/yields.tsv .

for os in $* ; do

    os_dir="${r_dir}/${os}_${date_str}_incsurv_cm0"
    ios_dir="${i_dir}/${os}"

    echo "$os : $os_dir : $ios_dir"

    ## get Insana files:
    scp -p pearson@codon_m:${ios_dir}/${os}_c0.5_p0.90_cm0.bad_clust .
    scp -p pearson@codon_m:${ios_dir}/${os}_c0.5_p0.90_cm0.bad_omes.samp.am .

    ## get Pearson files:
    scp -p pearson@codon_m:${os_dir}/${os}_proteome.stats .
    ## tfx short protein files
    scp -p pearson@codon_m:${os_dir}/${os}_bad_omes_2000.tfxg_stats_\* .
    ## missing prot tfx files
    scp -p pearson@codon_m:${os_dir}/${os}_miss_tfxg.summ_r3 .

done
