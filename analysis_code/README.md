
## ett_ms/analysis_code

This directory provides all of the scripts that were used to analyze the data for the paper:
Early terminated transcripts and missing proteins reflect artifacts in bacterial proteomes.

This file outlines the mapping between the files produced by the
various cluster analysis scripts, and the figures in the papers.  The
'R' code used to produce the figures is provided in the `../figures_code` directory

## Overall analysis strategy.

1. The [clustering and labelling pipeline](https://github.com/g-insana/ProteomeCluster) produces a clusters file:
`Protein_clusters_m50pct.tsv`, with tab-delimited fields starting with a header line:
```
cluster_id protein_ids proteins_count proteomes_count representative seqlen_range seqlen_mode
```

Most of the analysis of clusters is done by parsing the `protein_ids` field, which consists of a space delimited string of the form:
```
1405376:CQR69793|367 1406373:KPA48333|367 1410261:CRF37971|367 1913545:KRS08592|367|UPI00028D7880 ...
```
where each entry consists of a `proteome_id:protein_accession|protein_length` and some entries include a third field, `|upi_accession` (this accession can be used to retrieve the sequence from UniParc).

2. The `proteome_stats5slm3.py` script reads the `Proteins_m50pct.tsv`
file and produces two summaries of outlier frequencies:
`OSCODE...bad_clust`, which summaries outlier numbers and frequences
for each individual cluster, and `OSCODE....bad_omes.samp.am`, which
samples up to 2000 proteomes and reports outlier frequencies for the
clusters those proteomes contribute to.

- Fig. 1 is produced from the `OSCODE...bad_clust` files, which is produced by `proteome_stats5slm4.py`

- Suppl. Fig. 1 is produced from the `OSCODE_proteome.stats` files, which are produced by `stat_prot_cluster.py`.

- Suppl. Fig 2 is produced from the `OSCODE_miss_tfxg.summ_r3` files, which are produced by `track_missing_prots2.py`, `genome_v_prot2.py`, and `parse_tfxg3s.py`.

- Fig. 2 is produced from `OSCODE_cl_XXXXXXX.clust` files (XXXXX is
  a cluster number for that OSCODE), which are produced by
  `one_cluster.sh Protein_clusters_m50pct.tsv XXXXXX` where `XXXXXX`
  is a cluster id from OSCODE....bad_clust that has > 60% mode proteins and > 200 aa mode length, and plotted by `fig2_clust_one_dist.R`

- Fig. 3 is produced from `OSCODE....bad_clust` files (`ECOLX,KLEPN,SALER,NEIGO`) plotted by `fig3_clust_dist2p.R`

3. the `proteome_stats5slm3.py` script produces an
`OSCODE....short_tfx_ex` file that contains up to 2000 examples of
both bad (proteomes with bad clusters with many outliers) or sampled
proteomes.  The `.short_tfx_ex` files include ENA genome IDs
(`GCA_1235567.1`) and a set of short outlier proteins from that
proteome/genome, and a set of mode-length proteins from the same
cluster.  The mode-length proteins are compared to the genome using
`tfastx36`, to look for mode-length alignments, with
`genome_v_prot2.py`, and the short proteins are compared to the
mode-length proteins with `short_prot_v_mode2.py` to check for
N-terminal or C-terminal extensions by the mode protein. The results
from these searches are parsed and combined by `parse_tfxg4s.py`, and
the results written to files: `OSCODE_bad_omes_2000.tfxg_stats_BHL3ss`
(bad proteomes) and `OSCODE_bad_omes_2000.tfxg_stats_S1K3ss`.

- Fig. 4 and Suppl. Fig. 4 are produced by the script
  `fig4_bad_omes_dist_busco2p.R` using the
  `OSCODE_bad_omes_2000.tfxg_stat_BHL3ss` and `_S1K3ss` together with the `OSCODE.....bad_omes.samp.am` files.

- Fig. 6 is produced from the same set of files from `ECOLX`
  (`ECOLX_c0.5_p0.90_cm0.bad_omes.samp.am`,
  `ECOLX_bad_omes_2000.tfxg_stats_BHL3ss`, and
  `ECOLX_bad_omes_2000.tfxg_stats_S1K3ss`) using the `fig6_bad_omes_meta_tfx.R` script.

- Suppl. Fig. 6 is produced with the same files and script as Fig. 6, but displaying data from 5 additional bacteria.

- Suppl. Fig. 7 summaries the properties of corrected short alignments
  using the sampled data from all 20 bacteria
  (`OSCODE_bad_omes_2000.tfxg_stats_S1K3ss`) and the `suppl_f7_align_props_box5p.R` script.


## How to run

The `bs_bad_clusters_oscode_ins2_cm0.sh` script runs all the necessary scripts.

Suggested operation:
```
source env.sh
mkdir -p outdir
$ETT_BIN/bs_bad_clusters_oscode_ins2_cm0.sh OSCODE
```

(`env.sh` sets up the necessary environment variables so that `bs_bad_clusters_oscode_ins2_cm0.sh` can run properly)

`$ETT_BIN/bs_bad_clusters_oscode_ins2_cm0.sh OSCODE` runs the analysis for one bacteria.

`$ETT_BIN/bs_bad_clusters_oscode_ins2_cm0.sh $(cat analysis_code/oscodes20)` runs the analysis for all 20 bacteria.

Note that cluster information needs to be downloaded and unpacked. a script in the main directory (`fetch_datasets.sh`) can be used for that purpose.

Proteome files also need to be present, to run the searches, these can be prepared and downloaded to `proteomes/` named subfolders under each `OSCODE/` directory, using either upid or gca_set_acc information present in the `OSCODE/OSCODE.proteomes.tsv` files.
