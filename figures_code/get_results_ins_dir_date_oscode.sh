#!/bin/sh

if [[ $# < 2 ]]; then
    echo "Usage: $0 DATA_DATE OSCODE1 [OSCODE2...]"
    echo "E.g.: $0 2026-05-18 TEST NEIGO"
    echo "Note: Replace 2026-05-18 with the date you defined in ../env.sh"
    exit 1
fi

date_str=$1

shift

r_dir=../outdir
i_dir=../data


ln -sf ${i_dir}/yields.tsv .

for os in $* ; do

    os_dir="${r_dir}/${os}_${date_str}"
    ios_dir="${i_dir}/${os}"

    if [ ! -d "${os_dir}" ]; then
      echo "ERROR: no such dir '${os_dir}'"
      exit 2
    fi

    echo "$os : $os_dir : $ios_dir"

    ## get input files:
    ln -sf ${ios_dir}/${os}_c0.5_p0.90_cm0.bad_clust .
    ln -sf ${ios_dir}/${os}_c0.5_p0.90_cm0.bad_omes.samp.am .

    ## get analysis files:
    ln -sf ${os_dir}/${os}_proteome.stats .
    ## cluster length files
    ln -sf ${os_dir}/${os}_clust_proteome_lens.samp_2B .
    ## tfx short protein files
    ln -sf ${os_dir}/${os}_bad_omes_2000.tfxg_stats_* .
    ## missing prot tfx files
    ln -s ${os_dir}/${os}_miss_tfxg.summ_r3 .

done
