#!/usr/bin/env Rscript --vanilla

################
## clust_tables.R ecoli...bad_clust
## 
################
## reads in 
## produces suppl table 1:  cluster yields
##
## 

library('ggplot2')
library(dplyr)
library(purrr)
## library(cowplot)
library(RColorBrewer)
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


p.name<-'clust_tables'
p.full_name <- paste0(p.name,'.R')

args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.full_name,"\n",args),collapse=' ', sep=' ')

opt_list = list(
    make_option(c("--file"),type='character',action='store',help='comma separated *.bad_clust, files REQUIRED',default=NA),
    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file', default=NA),
    make_option(c("--supp_table1"),type='character',action='store',help='supp_table1 file name', default='supp_table1.tex'),
    make_option(c("--pub"),action='store_true',help='publication plot', default=FALSE)
)

opt <- get_yaml_opts(args, opt_list, paste0(p.name,'.yaml'))

if (length(opt$file)==1) {
    file_list <- strsplit(opt$file,',')[[1]]
} else {
    file_list <- opt$file
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

yield_names=c('OSCODE','N_omes','N_clust','N_seq','N_clust_ns','N_seq_ns','N_clust_10pct','N_seq_10pct','N_clust_50pct','N_seq_50pct')

yield_df = NULL

for (clust_file in file_list) {

    print(clust_file)
    
    if (file.exists(clust_file)) {

	this_data <- read.table(clust_file,header=TRUE,sep='\t')
##	print(names(this_data))
##	print(yield_names)

	os_stats_names = names(this_data)
        if ('OSCODE' %in% os_stats_names && 'proteomes_total' %in% os_stats_names) {
           ## change the names in yield_names
	   fa_names <- c('ps_min'='fa_min','ps_q1'='fa_q1','ps_median'='fa_med','ps_q3'='fa_q3','ps_max'='fa_max','ps_mean'='fa_mean')
	   yield_names <- c(yield_names, fa_names)
        } 

        names(this_data) = yield_names

	## calculate percentages

	this_data$pct_seq_ns = round(100.0*this_data$N_seq_ns/this_data$N_seq,1)
	this_data$pct_clust_ns = round(100.0*this_data$N_clust_ns/this_data$N_clust,1)
	this_data$pct_seq_10pct = round(100.0*this_data$N_seq_10pct/this_data$N_seq,1)
	this_data$pct_clust_10pct = round(100.0*this_data$N_clust_10pct/this_data$N_clust,1)
	this_data$pct_seq_50pct = round(100.0*this_data$N_seq_50pct/this_data$N_seq,1)
	this_data$pct_clust_50pct = round(100.0*this_data$N_clust_50pct/this_data$N_clust,1)

    	## print(head(this_data))
	
	yield_df <- rbind(yield_df, this_data)
    }
}

## print(yield_df)

yield_medians <- sapply(yield_df,median)

print("yield medians")
print(yield_medians)

median_df <- data.frame(round(t(yield_medians),1))
median_df$OSCODE='median'
median_df[,c(2:4)]=NA

## median_df_names = names(median_df)
## print(median_df)

## print("yield_df:")
## print(yield_df[-order(yield_df$N_omes),])

yield_df_tbl <- rbind(yield_df[order(-yield_df$N_omes),],median_df)
print("yield_df_tbl:")
yield_df_tbl

## print out opt$table1 -- ori_clust, good_clust, pct_good, pct_prot_coverage

## out_fields = c('OSCODE', 'N_omes','N_clust','N_seq','N_clust_ns','N_clust_10pct','N_clust_50pct','pct_clust_ns','pct_clust_10pct','pct_clust_50pct','N_seq_ns','N_seq_10pct','N_seq_50pct','pct_seq_ns','pct_seq_10pct','pct_seq_50pct')
## float_fields = c(8,9,10,14,15,16)+1

out_fields = c('OSCODE', 'N_omes','N_clust','N_seq','pct_clust_ns','N_clust_10pct','pct_clust_10pct','N_clust_50pct','pct_clust_50pct','pct_seq_ns','pct_seq_10pct','pct_seq_50pct')
float_fields = c(5,7,9,10:12)+1

out_fields_digits=rep(0,length(out_fields)+1)
print("out_fields:")
print(out_fields)
print("float_fields:")
print(float_fields)

out_fields_digits[float_fields]=1
print("out_fields_digits:")
print(out_fields_digits)

sink(opt$supp_table1)
print(xtable(yield_df_tbl[,out_fields], digits=out_fields_digits),include.rownames=FALSE,NA.string='',format.args=list(big.mark=' '))
sink()

