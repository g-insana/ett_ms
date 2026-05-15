#!/usr/bin/env Rscript --vanilla

################
## clust_tables.R ecoli...bad_clust
## 
################
## reads in *.bad_clust, --stats file --names file
## produces table 1:
## bacteria OSCODE orig_clust good_clust %good_clust med_good_clust
## 
## produces table 2:
## OSCODE N_good_clust N_out_50 N_out_75 N_out_95 %_out_50 %_out_75 %_out_95
##

library('dplyr')
library('purrr')
library('stringr')
library('optparse')
library('yaml')
library('xtable')

get_yaml_file = 'get_yaml_opts.R'

if (file.exists(get_yaml_file)) {
   source(get_yaml_file)
} else if (file.exists(paste0("../",get_yaml_file))) {
   source(paste0("../",get_yaml_file))
} else {
  cat(paste("cannot open", get_yaml_file))
  sys.exit(1)
}

## .bad_clust
## 0 clusters with >95% of proteins outside length range
## clust_id	n_omes	p_tot	p_bad	pct_bad	p_short	pct_short	bad_SL_flag	mode_len	p_mode	pct_mode	q1_len	med_len	q3_len
## 8159 	1345	1346	980	72.81	0	0.00	Long	354	297	22.07	354	587	590
## 38126	974	974	643	66.02	0	0.00	Long	90	331	33.98	90	135	160

## clust_qual_all3.R 

p.name<-'clust_tables'
p.full_name <- paste0(p.name,'.R')

args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.full_name,"\n",args),collapse=' ', sep=' ')

opt_list = list(
    make_option(c("--files"),type='character',action='store',help='comma separated *.bad_clust, files REQUIRED',default=NA),
    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file', default=NA),
    make_option(c("--table1"),type='character',action='store',help='table1 file name', default='table1.tex'),
    make_option(c("--table2"),type='character',action='store',help='table1 file name', default=NA),
    make_option(c("-S","--stats"),type='character',action='store',help='os_stats_file', default=NA),
    make_option(c("-N","--names"),type='character',action='store',help='os_names_file', default=NA),
    make_option(c("-F","--fa_stats"),type='character',action='store',help='fa_stats_file', default=NA),
    make_option(c("-M","--mode_good"),action='store_true',help='only good modes', default=TRUE),
    make_option(c("--mode_pct_thresh"),type='double', action='store',help='mode pct threshold', default=60.0),
    make_option(c("--pub"),action='store_true',help='publication plot', default=FALSE)
)

opt <- get_yaml_opts(args, opt_list, paste0(p.name,'.yaml'))

if (length(opt$files)==1) {
    file_list <- strsplit(opt$files,',')[[1]]
} else {
    file_list <- opt$files
}

have_fa_stats <- FALSE

os_stats <- NULL
if (!is.na(opt$stats) && length(opt$stats)>0) {
   os_stats <- read.table(opt$stats, header=TRUE, sep='\t')
   ## print("os_stats")
   ## print(head(os_stats))

   os_stats_names = names(os_stats)

   if ('OSCODE' %in% os_stats_names && 'proteomes_total' %in% os_stats_names) {
      ## change the names in os_stats_name_edits

      os_stats_name_edits = c('OSCODE'='oscode','proteomes_total'='count','clusters_50pct'='clust_cnt',
                              'ps_min'='fa_min','ps_q1'='fa_1q','ps_median'='fa_med','ps_q3'='fa_3q','ps_max'='fa_max','ps_mean'='fa_mean')

      new_os_names <- str_replace_all(os_stats_names, os_stats_name_edits)   
      names(os_stats) <- new_os_names
      have_fa_stats <- TRUE
      os_stats$clust_fn <- os_stats$clust_cnt/os_stats$fa_med
   }

   os_stats <- os_stats[order(-os_stats$count),]
   ## print(os_stats)
}

os_names <- NULL
if (!is.na(opt$names) && length(opt$names)>0) {
   os_names <- read.table(opt$names, header=TRUE, sep='\t')
   ## print("os_stats")
   ## print(head(os_stats))
   os_stats <- merge(os_stats,os_names,by='oscode')
   print(os_stats)
}

fa_stats <- NULL
if (!have_fa_stats && !is.na(opt$fa_stats) && length(opt$fa_stats)>0) {
   fa_stats <- read.table(opt$fa_stats, header=TRUE, sep='\t')
   print(paste("fa_stats",opt$fa_stats))
   print(head(os_stats))
   fa_stats <- fa_stats[order(-fa_stats$ome_cnt),]
   print(fa_stats)
} else {
   fa_stats <- os_stats[,c('oscode','count','fa_min','fa_1q','fa_med','fa_3q','fa_max')]
   colnames(fa_stats)[2]<-'ome_cnt'
}


clust_data <- NULL
leg.labels.oscode = c()
leg.labels.nrow = c()
taxon.names = c()

## as we read each cluster file, we want the number of clusters with different levels of errors
## we need a place to store those numbers
## taxon thresh N_bad pct_bad N_short pct_short N_long pct_long
##

max_cluster = 0

mode_good_fn = c()
mode_good_df = NULL

## table for all taxa
clust_data = NULL
clust_thresh = NULL

for (clust_file in file_list) {

    ## print(clust_file)
    
    if (file.exists(clust_file)) {

        this_taxon = strsplit(clust_file,'_')[[1]][1]
    	taxon.names <- append(taxon.names, this_taxon)

	this_os_info <- os_stats[os_stats$oscode==this_taxon,]
##	print(this_os_info)

	this_data <- read.table(clust_file,header=TRUE,sep='\t')

    	## print(head(this_data))

    	this_data$taxon = this_taxon 

    	old_nrow = nrow(this_data)
    	new_nrow <- old_nrow

        old_nrow = nrow(this_data)
        this_data <- this_data[this_data$pct_mode > opt$mode_pct_thresh,]
        new_nrow = nrow(this_data)
	
	mode_str = sprintf("%d\t%d\t%.1f",old_nrow, new_nrow,100.0*new_nrow/old_nrow)
	mode_str <- setNames(mode_str, this_taxon)

	med_Nomes <- median(this_data$n_omes)
	q1_Nomes <- quantile(this_data$n_omes, probs=c(0.25))

	med_Nmode <- median(this_data$p_mode)
	q1_Nmode <- quantile(this_data$p_mode, probs=c(0.25))

	print(this_os_info)
	print(new_nrow)
	print(quantile(this_data$p_mode, probs=c(0.1,0.25,0.5))/this_os_info$count)
	
	med_pct_omes <- 100.0*med_Nomes/this_os_info$count
	q1_pct_omes <- 100.0*q1_Nomes/this_os_info$count

	med_pct_mode <- 100.0*med_Nmode/this_os_info$count
	q1_pct_mode <- 100.0*q1_Nmode/this_os_info$count

##	print(c(med_Nomes, med_pct_omes))

	mode_good_fn <- append(mode_good_fn, mode_str)

	print("this_os_info:")
	print(head(this_os_info))

	this_good_df <- data.frame('s_name'=this_os_info$short_name, 'taxon'=this_taxon,
	                           'N_omes'=this_os_info$count, 'N_ori_clust'=old_nrow, 'N_good_clust'=new_nrow,
	                           'pct_good_clust_n'=round(100*new_nrow/old_nrow,1),
	                           'pct_good_clust'=sprintf("%.1f",100*new_nrow/old_nrow),
				   'N50_omes'=sprintf("%.0f",med_Nomes), 'N25_omes'=sprintf("%.0f",q1_Nomes),
				   'pct50_omes'=sprintf("%.1f",med_pct_omes), 'pct25_omes'=sprintf("%.1f",q1_pct_omes),
				   'pct50_mode'=med_pct_mode, 'pct25_mode'=q1_pct_mode
				   )
							   
	mode_good_df <- rbind(mode_good_df, this_good_df)
    }

    clust_data <- rbind(clust_data, this_data)
}



if (nrow(fa_stats) > 0) {
   mode_good_df <- merge(mode_good_df, fa_stats[,c('oscode','ome_cnt','fa_med')], by.x='taxon', by.y='oscode')
   mode_good_df$pct_fa_good <- round(100.0 * mode_good_df$N_good_clust/mode_good_df$fa_med,1)
   print("fa_stats merge mode_good_df")
   print(mode_good_df)
  ome_stat_cols=c('N_omes','fa_med','N_good_clust','pct_fa_good','pct_good_clust_n')
  float_fields = c(6:7)+1
  table1_fields = c('s_name','taxon','N_omes','fa_med', 'N_good_clust','pct_fa_good','pct_good_clust_n')
} else {
  ome_stat_cols=c('N_omes','N_good_clust','pct_good_clust_n')
  float_fields = c(5)+1
  table1_fields = c('s_name','taxon','N_omes', 'N_good_clust','pct_good_clust')
}

## italisize s_name
mode_good_df$s_name = paste0('\\textit{',mode_good_df$s_name,'}')

mode_good_medians <- sapply(mode_good_df[,ome_stat_cols],median)
median_df <- data.frame('s_name'='median','taxon'=NA, round(t(mode_good_medians),1))

print("head(mode_good_df)")
print(head(mode_good_df[,table1_fields]))
print("")
print("*** summary(mode_good_df[,'N_good_clust']")
print(summary(mode_good_df[,'N_good_clust']))
print("median_df")
print(median_df[,table1_fields])

mode_good_df_tbl <- rbind(mode_good_df[order(-mode_good_df$N_omes),table1_fields],median_df[,table1_fields])

## print out opt$table1 -- ori_clust, good_clust, pct_good, pct_prot_coverage

out_fields_digits=rep(0,length(names(mode_good_df_tbl))+1)
print("out_fields_digits:")
print(out_fields_digits)
print("float_fields:")
print(float_fields)

out_fields_digits[c(1:3)] = NA
out_fields_digits[float_fields]=1
print("out_fields_digits:")
print(out_fields_digits)

sink(opt$table1)
print(xtable(mode_good_df_tbl, digits=out_fields_digits), include.rownames=FALSE, NA.string='', format.args=list(big.mark=' '),sanitize.text.function=identity)
sink()

this_quantile <- function(x) {
   quants <- quantile(x, probs=c(0.5, 0.75, 0.95))
   q_names = names(quants)
   names(quants) <- q_names
   quants
}

## calculate the columns for table 1 In addition to ori_clusters,
## good_clusters, and %good_clusters, I need the median number of
## proteomes for each cluster/total_proteomes for taxa

taxa_medians <- aggregate(pct_bad~taxon, data=clust_data,FUN=this_quantile)

taxa_medians_ag <- data.frame(taxon=taxa_medians$taxon, 
		'pct50_bad'=taxa_medians$pct_bad[,1],
		'pct75_bad'=taxa_medians$pct_bad[,2],
		'pct95_bad'=taxa_medians$pct_bad[,3])

taxa_medians_ag <- merge(taxa_medians_ag, mode_good_df[,c('taxon','N_good_clust','pct50_mode','pct25_mode')],by='taxon',all.x=TRUE)

names(taxa_medians_ag) = c('taxon','pct50_bad','pct75_bad','pct95_bad','N_good','pct50_mode','pct25_mode')


taxa_medians_ag <- within(taxa_medians_ag, {N95_bad <- pct95_bad*N_good/100
		                          N75_bad <- pct75_bad*N_good/100
					  N50_bad <- pct50_bad*N_good/100 })


taxa_median_medians <- sapply(taxa_medians_ag[,c('pct50_bad','pct75_bad','pct95_bad','pct50_mode','N50_bad','N75_bad','N95_bad')],median)

print("taxa_median_medians")
print(taxa_median_medians)

taxa_median_df <- data.frame('taxon'='median',t(taxa_median_medians))
print(taxa_median_df)

taxa_medians_p <- within(taxa_medians_ag, {pct50_bad<-sprintf("%0.2f",pct50_bad)
	                                   pct95_bad<-sprintf("%0.2f",pct95_bad)
					   pct50_mode <- sprintf("%0.1f",pct50_mode)
					   pct25_mode <- sprintf("%0.1f",pct25_mode)
	                                   N50_bad<-sprintf("%.0f",N50_bad)
	                                   N75_bad<-sprintf("%.0f",N75_bad)
	                                   N95_bad<-sprintf("%.0f",N95_bad)
					   }
					   )

table2_fields = c("taxon","pct50_mode","N50_bad","N75_bad","N95_bad","pct50_bad","pct75_bad","pct95_bad")

## taxa_medians_p <- rbind(taxa_medians_p[-order(taxa_medians_p$taxon),table2_fields], taxa_median_df)
## taxa_medians_p <- rbind(taxa_medians_p[,table2_fields], taxa_median_df)

if (!is.na(opt$table2) && length(opt$table2)>0) {
  sink(opt$table2)
  print(xtable(taxa_medians_p[,table2_fields]),include.rownames=FALSE,format.args=list(big.mark=' '))
  sink()
}

short_stats <- clust_data %>% group_by(taxon) %>%
    summarize(
        count_S = sum(bad_SL_flag=='Short'),
        count_L = sum(bad_SL_flag=='Long'),
	count_Z_fct = sum(p_bad < 1)/n(),
	short_fct = ifelse(count_S + count_L > 0, count_S/(count_S + count_L), NA)

    ) %>%
    arrange (desc(short_fct)) %>%
    ungroup()

print("short_stats")
print(short_stats)
print("summary(short_stats)")
summary(short_stats$short_fct)

