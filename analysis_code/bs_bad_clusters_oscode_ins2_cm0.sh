#!/bin/sh

## version of analysis script that only runs, and then parses, the
## .short_tfx_ex accesssions vs the appropriate genomic DNA (genome_v_prot2.py)
## also compares short protein to mode protein with ssearch (short_prot_v_mode2.py)

do_stats=0
do_tfastx=1
do_miss=1

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

echo "using: out_dir: $bs_ett_dir data_dir: $DATA_DIR data_date: $data_date"

export OUT_DIR=$bs_ett_dir

if [[ ! -d $OUT_DIR ]]; then
    mkdir -p $OUT_DIR
fi

proteome_stats_bin=proteome_stats5slm4.py

param_name="c0.5_p0.90_cm0"

for OSCODE in $*; do

    tfx_dir=$OUT_DIR/${OSCODE}_${data_date}

    tfx_ex_file=${DATA_DIR}/${OSCODE}/${OSCODE}_${param_name}.short_tfx_ex
    
    clusters_file=${DATA_DIR}/${OSCODE}/${OSCODE}_Protein_clusters_m50pct.tsv

    if [[ ! -d $tfx_dir ]]; then
	mkdir $tfx_dir
    fi

    echo $tfx_dir

    tfx_name=${OSCODE}_${param_name}

    pushd $tfx_dir

    data_prefix=${OSCODE}_${param_name}
    ins_data_prefix=$DATA_DIR/$OSCODE/$data_prefix

    map_file=$DATA_DIR/$OSCODE/${OSCODE}_proteomeids_busco.tsv 
    if [[ ! -f $map_file ]]; then
	map_file=$DATA_DIR/$OSCODE/${OSCODE}.proteomes.tsv
	if [[ ! -f $map_file ]]; then
	    echo "Mapping file: $map_file not found"
	    exit 1
	fi
    fi

    echo "working on $OSCODE/$data_prefix"

    if [[ $do_stats == 1 ]]; then
    ## this version works with the files produced in the initial analysis in ${DATA_DIR}/${OSCODE}/${OSCODE}_${param_name}.short_tfx_ex
    ## so it does not run $proteome_stats_bin

        echo "program: " $proteome_stats_bin

	$ETT_BIN/$proteome_stats_bin --mode_good_only --bad_low_fract 0.0 --target_length=200 -P $map_file --out_clust_file=${OSCODE}_${param_name}.bad_clust --out_omes_file=${OSCODE}_${param_name}.bad_omes.samp.am --out_tfx_file=$tfx_name.short_tfx_ex --samp_N50 --samp_prots=2000 --bad_tfx_ex=100 --bad_prots=2000 --samp_tfx_ex=100 --long_prots=100 --long_tfx_clust_file=$tfx_name.long_ex_clust --long_tfx_ome_file=$tfx_name.long_ex_ome ${clusters_file}
	echo "$OSCODE proteome_stats done" `date`
    fi

    ## produce OSCODE_proteome.stats files (Suppl. Fig. 1)
    echo "produce ${OSCODE}_proteome.stats"
    ${ETT_BIN}/stat_prot_cluster.py --in_clust_file $ins_data_prefix.bad_clust --in_samp_omes $ins_data_prefix.bad_omes.samp.am -P $map_file --miss_samp_file ${OSCODE}_miss.tab ${clusters_file} > ${OSCODE}_proteome.stats
    
    echo "produce ${OSCODE}_clust_proteome_lens.samp_2B"
    ${ETT_BIN}/track_prot_cluster.py -X='B' --mode_pct 80.0 --in_clust_file $ins_data_prefix.bad_clust --in_samp_omes $ins_data_prefix.bad_omes.samp.am -P $map_file ${clusters_file} > ${OSCODE}_clust_proteome_lens.samp_2B

    if [[ $do_tfastx == 1 ]]; then
	## do the .tfx searches
	rm -f ${OSCODE}_bad_omes_2000.tfxg_stats_BHL3ss ${OSCODE}_bad_omes_2000.tfxg_stats_S1K3ss
	${ETT_BIN}/genome_v_prot2.py --script down_fasta_acc.sh --oscode $OSCODE $tfx_ex_file
	## genome_v_prot2.py needs down_fasta_acc.sh, 

	echo "$OSCODE genome_v_prot2.py short done" `date`

	## now do the short vs mode ssearch
	rm -f ${OSCODE}_cl_*_short.ss_MD10
  echo "short vs mode"
	echo "${ETT_BIN}/short_prot_v_mode2.py --script down_fasta_acc.sh --oscode $OSCODE --run $tfx_ex_file"
	${ETT_BIN}/short_prot_v_mode2.py --script down_fasta_acc.sh --oscode $OSCODE --run $tfx_ex_file

	echo "$OSCODE short_prot_v_mode2.py short done" `date`

	## combine protein overlap info with genomic overlap info
	${ETT_BIN}/parse_tfxg4s.py -S --tfx_ex $tfx_ex_file ${OSCODE}_GCA*_B_N50*.tfxg > ${OSCODE}_bad_omes_2000.tfxg_stats_BHL3ss
	${ETT_BIN}/parse_tfxg4s.py -S --tfx_ex $tfx_ex_file ${OSCODE}_GCA*_S_N50*.tfxg > ${OSCODE}_bad_omes_2000.tfxg_stats_S1K3ss

	echo "$OSCODE tfxg files parsed, tfxg_stats files produced"

    fi

    if [[ $do_miss == 1 ]]; then

	## look for proteins missing from proteomes
	
	echo "${ETT_BIN}/track_missing_prots2.py --miss_fract=0.67 -P $map_file ${clusters_file} > ${OSCODE}_miss_prots.miss_tfx_ex_r3"

	${ETT_BIN}/track_missing_prots2.py --miss_fract=0.67 -P $map_file ${clusters_file} > ${OSCODE}_miss_prots.miss_tfx_ex_r3

	rm -f *.miss_tfx_r3
	${ETT_BIN}/genome_v_prot2.py --script down_fasta_acc.sh --oscode $OSCODE --suff miss_tfx_r3 ${OSCODE}_miss_prots.miss_tfx_ex_r3

	${ETT_BIN}/parse_tfxg3s.py *.miss_tfx_r3 > ${OSCODE}_miss_tfxg.summ_r3

	echo "$OSCODE missing protein tfx search done"

	${ETT_BIN}/proteome_v_prot2.py -P $map_file --script down_fasta_acc.sh --oscode $OSCODE --suff miss_ok2_r3 ${OSCODE}_miss_prots.miss_tfx_ex_r3

	${ETT_BIN}/parse_fa36.py *.miss_ok2_r3 > ${OSCODE}_pmiss_ok2.summ_r3

	echo "$OSCODE missing protein fasta search done"
    fi

    popd
done
