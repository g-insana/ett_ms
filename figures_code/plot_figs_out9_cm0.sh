#!/bin/sh

## plot_figs_cm0.sh (run in results)

## modified 7-April-2026 to get proteome/cluster statistics from outdir8_incsurv_cm0.tsv
## 

BIN_DIR=../

PUB_STR="--pub"
export YAML=../yaml

## figure 1 -- cluster qualities across 20 taxa
echo "fig1:: $BIN_DIR/clust_tot_qual_all5r.R -Y $YAML/clust_lens_20_cm0.yaml --pdf f1_clust_stat_cm0.pdf"
$BIN_DIR/fig1_clust_qual_5p.R -Y $YAML/clust_lens_20_cm0.yaml --stats yields.tsv $PUB_STR --pdf f1_clust_stats_cm0_pub.pdf
##ln -s -f `pwd`/f1_clust_stats_cm0_pub.pdf `pwd`/suppl_f1_clust_stats_cm0_pub.pdf ../figs

## figure 2 -- use current version
## echo "$BIN_DIR/clust_one_dist.R -Y clust_one_dist_4bad.yaml --stats yields.tsv --pub --pdf f2_clust_one_4bad_pub.pdf"
## $BIN_DIR/fig2_clust_one_dist.R -Y $YAML/f2_clust_one_dist_4bad.yaml --pub --pdf f2_clust_one_4bad_pub.pdf

## figure 3 -- cluster quality for 4 distant taxa
echo "fig3:: $BIN_DIR/fig3_clust_dist2m.R -Y clust_dist_4bad_cm0.yaml --stats yields.tsv --pdf f3_clust_dist_4far.pdf"
$BIN_DIR/fig3_clust_dist2p.R -Y $YAML/clust_dist_4bad_cm0.yaml --stats yields.tsv $PUB_STR --pdf f3_clust_dist_4bad_cm0_pub.pdf
##ln -s -f `pwd`/f3_clust_dist_4bad_cm0_pub.pdf ../figs

## suppl_figure 1 -- proteome qualities across 20 taxa
echo "suppl_fig1:: $BIN_DIR/suppl_f1_omes_stats_qual6p.R -Y $YAML/clust_omes_stats_20_cm0.yaml --stats yields.tsv --pdf suppl_f64_omes_stats.pdf"
$BIN_DIR/suppl_f1_omes_stats_qual6p.R -Y $YAML/clust_omes_stats_20_cm0.yaml --stats yields.tsv $PUB_STR --pdf suppl_f1_omes_stats6_cm0_pub.pdf
##ln -s -f `pwd`/suppl_f4_omes_stats6_cm0_pub.pdf ../figs

echo "suppl_fig2 (miss) $BIN_DIR/suppl_f2_miss_stats.R -Y ../yaml/miss_omes.yaml $PUB_STR --pdf suppl_f2_miss_summ_r3_pub.pdf"
$BIN_DIR/suppl_f2_miss_stats.R -Y ../yaml/miss_omes.yaml $PUB_STR --pdf suppl_f2_miss_summ_r3_pub.pdf
##ln -s -f `pwd`/suppl_f2_miss_summ_r3_pub.pdf ../figs

## suppl figure 3 -- distribution of outliers across proteomes
echo "suppl_fig5:: $BIN_DIR/suppl_f3_bad_omes_dist_clust.R -Y $YAML/bad_omes_3bad_clust2B_cm0.yaml --stats yields.tsv --pdf suppl_f3_omes_dist_lens_cm0.pdf"
$BIN_DIR/suppl_f3_bad_omes_dist_clust.R -Y $YAML/bad_omes_3bad_clust2B_cm0.yaml --stats yields.tsv $PUB_STR --pdf suppl_f3_omes_dist_lens_cm0_pub.pdf
##ln -s -f `pwd`/suppl_f3_omes_dist_lens_cm0_pub.pdf ../figs

## figure 4 -- distribution of proteome cluster quality
echo "fig4:: $BIN_DIR/fig4_bad_omes_dist_busco2p.R -Y $YAML/bad_omes_3bad_cm0.yaml --stats yields.tsv --pdf f4_omes_dist_busco_2far.pdf"
$BIN_DIR/fig4_bad_omes_dist_busco2p.R -Y $YAML/bad_omes_3bad_cm0.yaml --stats yields.tsv $PUB_STR --pdf f4_omes_dist_busco_3bad_cm0_pub.pdf
##ln -s -f `pwd`/f4_omes_dist_busco_3bad_cm0_pub.pdf ../figs

echo "suppl_f4:: suppl_f4_bad_omes_busco_all.R -Y $YAML/bad_omes_3bad_cm0.yaml --stats yields.tsv $PUB_STR --pdf suppl_f4_omes_all_busco_3bad_cm0.pdf"
$BIN_DIR/suppl_f4_bad_omes_busco_all.R -Y $YAML/bad_omes_3bad_cm0.yaml --stats yields.tsv $PUB_STR --pdf suppl_f4_omes_all_busco_3bad_cm0_pub.pdf
##ln -s -f `pwd`/suppl_f4_omes_all_busco_3bad_cm0_pub.pdf ../figs

## fig 5 -- alignment figure fig5_bglr_puuc_sized.pdf

## figure 6 -- lots of measures of proteome quality
echo "fig6:: $BIN_DIR/fig6_bad_omes_meta_tfx.R -Y $YAML/bad_omes_ECOLX_cm0.yaml --stats yields.tsv --pdf f7_omes_dist_ECOLX.pdf"
$BIN_DIR/fig6_bad_omes_meta_tfx.R -Y $YAML/bad_omes_ECOLX_cm0.yaml --stats yields.tsv $PUB_STR --pdf f6_omes_dist_ECOLX_cm0_pub.pdf

echo "suppl_fig6a :: $BIN_DIR/suppl_f6_bad_omes_meta_tfx.R -Y $YAML/bad_omes_dist3m6_cm0.yaml --stats yields.tsv --pdf f7_omes_dist_ECOLX.pdf"
$BIN_DIR/fig6_bad_omes_meta_tfx.R -Y $YAML/bad_omes_dist3m6_cm0.yaml --stats yields.tsv $PUB_STR --pdf suppl_f6_omes_dist_max_cm0_pub.pdf
##ln -s -f `pwd`/suppl_f6_omes_dist_max_cm0_pub.pdf ../figs

## figure 7 -- distribution of short/long alignment errors/corrections
echo "suppl_fig7:: $BIN_DIR/suppl_f7_align_props_box5p.R -Y $YAML/align_20_short/long.yaml --stats yields.tsv  --pdf f7_align20_short/long_box_cm0.pdf"
$BIN_DIR/suppl_f7_align_props_box5p.R -Y $YAML/align_20_short.yaml --stats yields.tsv  $PUB_STR --pdf suppl_f7_align20_short_box_cm0_pub.pdf
##ln -s -f `pwd`/suppl_f7_align20_short_box_cm0_pub.pdf ../figs

echo ""
echo "table1:: ../clust_tables2.R -Y ../yaml/clust_lens_20_cm0.yaml --stats yields.tsv --table1 clust_info.tb1.tex --table2 clust_info.tb2.tex"
../clust_tables2.R -Y ../yaml/clust_lens_20_cm0.yaml --stats yields.tsv --table1 clust_info.tb1.tex

echo ""
echo "suppl table1:: ../clust_table_suppl.R --file=yields.tsv --supp_table1 ett_suppl_tbl1_outdir9_cm0.tex"
../clust_table_suppl.R --file=yields.tsv --supp_table1 ett_suppl_tbl1_outdir9_cm0.tex
