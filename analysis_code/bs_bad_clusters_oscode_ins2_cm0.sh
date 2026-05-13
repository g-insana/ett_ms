#!/bin/sh

## version of analysis script that only runs, and then parses, the
## .short_tfx_ex accesssions vs the appropriate genomic DNA (genome_v_prot2.py)
## also compares short protein to mode protein with ssearch (short_prot_v_mode2.py)

if [[ -z "$DATA_DIR" ]]; then
    data_dir=outdir9
else
    data_dir=$DATA_DIR
fi

if [[ -z "$DATA_EXT" ]]; then
    data_ext=incsurv_cm0
else
    data_ext=$DATA_EXT
fi

if [[ -z "$DATA_DATE" ]]; then
    data_date=`date +%Y-%m-%d`
else
    data_date=$DATA_DATE
fi

if [[ -z "$BS_ETT_DIR" ]]; then
    bs_ett_dir="bs_tmp/ett_ms"
else
    bs_ett_dir=$BS_ETT_DIR
fi

if [[ -z "$FA_DIR" ]]; then
    fa_dir=fasta9_incsurv
else
    fa_dir=$FA_DIR
fi

echo "using: bs_ett_dir: $bs_ett_dir data_dir: $data_dir data_ext: $data_ext data_date: $data_date fa_dir: $fa_dir"

export INS_PROT=/homes/insana/ett2/${data_dir}_${data_ext}
export ETT_BIN=/homes/pearson/ett/ett_ms/analysis_code
export TMP_DIR=/homes/pearson/$bs_ett_dir/${data_dir}_${data_ext}

if [[ ! -d $TMP_DIR ]]; then
    mkdir $TMP_DIR
fi

proteome_stats_bin=proteome_stats5slm4.py

param_name="c0.5_p0.90_cm0"

for OSCODE in $*; do

    tfx_dir=$TMP_DIR/${OSCODE}_${data_date}_${data_ext}

    tfx_ex_file=${INS_PROT}/${OSCODE}/${OSCODE}_${param_name}.short_tfx_ex

    if [[ ! -d $tfx_dir ]]; then
	mkdir $tfx_dir
    fi

    echo $tfx_dir

    tfx_name=${OSCODE}_${param_name}

    pushd $tfx_dir

    ## this version works with the files produced in the initial analysis in ${INS_PROT}/${OSCODE}/${OSCODE}_${param_name}.short_tfx_ex
    ## so it does not run $proteome_stats_bin

    ##    echo "program: " $proteome_stats_bin

    ## $ETT_BIN/$proteome_stats_bin --mode_good_only --bad_low_fract 0.0 --target_length=200 -P $INS_PROT/$OSCODE/${OSCODE}_proteomeids_busco_assemblymetadata.tsv --out_clust_file=${OSCODE}_${param_name}.bad_clust --out_omes_file=${OSCODE}_${param_name}.bad_omes.samp.am --out_tfx_file=$tfx_name.short_tfx_ex --samp_N50 --samp_prots=2000 --bad_tfx_ex=100 --bad_prots=2000 --samp_tfx_ex=100 --long_prots=100 --long_tfx_clust_file=$tfx_name.long_ex_clust --long_tfx_ome_file=$tfx_name.long_ex_ome $INS_PROT/$OSCODE/Protein_clusters_m50pct.tsv
    ## echo "$OSCODE proteome_stats done" `date`

    ## produce OSCODE_proteome.stats files (Suppl. Fig. 1)
    ${ETT_BIN}/stat_prot_cluster.py --in_clust_file $data_prefix.bad_clust --in_samp_omes $data_prefix.bad_omes.samp.am -P $INS_PROT/${OSCODE}/${OSCODE}.proteomes.tsv $INS_PROT/${OSCODE}/Protein_clusters_m50pct.tsv --miss_samp_file ${OSCODE}_miss.tab > ${OSCODE}_proteome.stats

    ## do the .tfx searches
    rm -f ${OSCODE}_bad_omes_2000.tfxg_stats_BHL3ss ${OSCODE}_bad_omes_2000.tfxg_stats_S1K3ss
    ${ETT_BIN}/genome_v_prot2.py --script down_fasta_dir_oscode_acc.sh --fa_dir $fa_dir --oscode $OSCODE $tfx_ex_file
    ## genome_v_prot2.py needs down_fasta_dir_oscode_acc.sh, 

    echo "$OSCODE genome_v_prot2.py short done" `date`

    ## now do the short vs mode ssearch
    rm -f ${OSCODE}_cl_*_short.ss_MD10
    ${ETT_BIN}/short_prot_v_mode2.py --script down_fasta_dir_oscode_acc.sh --fa_dir $fa_dir --oscode $OSCODE --run $tfx_ex_file

    echo "$OSCODE short_prot_v_mode2.py short done" `date`

    ## combine protein overlap info with genomic overlap info
    ${ETT_BIN}/parse_tfxg4s.py -S --tfx_ex $tfx_ex_file ${OSCODE}_GCA*_B_N50*.tfxg > ${OSCODE}_bad_omes_2000.tfxg_stats_BHL3ss
    ${ETT_BIN}/parse_tfxg4s.py -S --tfx_ex $tfx_ex_file ${OSCODE}_GCA*_S_N50*.tfxg > ${OSCODE}_bad_omes_2000.tfxg_stats_S1K3ss

    echo "$OSCODE tfxg files parsed, tfxg_stats files produced"

    popd

done
