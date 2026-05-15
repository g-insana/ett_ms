#!/usr/bin/env Rscript --vanilla

################
## clust_all_qual_all.R ecoli...bad_clust
## 
################
## reads in *.bad_clust
## extracts percent of clusters with 0.0, 0.1, 1.0, and 10.0% outliers
## (later include long/short statistics?)
##

library('ggplot2')
library('scales')
library('dplyr')
library('tidyr')
library('stringr')
library('purrr')
library('scales')
library('RColorBrewer')
library('patchwork')
library('optparse')
library('yaml')

options(tibble.width=Inf)

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

p.name<-'clust_mode_qual_all5'
p.full_name <- paste0(p.name,'.R')

args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.full_name,"\n",args),collapse=' ', sep=' ')

opt_list = list(
    make_option(c("--files"),type='character',action='store',help='comma separated *.bad_clust, files REQUIRED',default=NA),
    make_option(c("-D","--debug"),action='store_true',help='debug', default=FALSE),
    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file', default=NA),
    make_option(c("-P","--pdf"),type='character',action='store',help='PDF file name', default=NA),

    make_option(c("-S","--stats"),type='character',action='store',help='os_stats_file', default=NA),
    make_option(c("-F","--fa_stats"),type='character',action='store',help='fa_stats_file', default=NA),
    make_option(c("-M","--mode_good"),action='store_true',help='only good modes', default=FALSE),
    make_option(c("-L","--log"), action='store_true', help='log panel A',default=FALSE),
    make_option(c("--mode_pct_thresh"),type='double', action='store',help='mode pct threshold', default=60.0),
    make_option(c("--pub"),action='store_true',help='publication plot', default=FALSE)
)

opt <- get_yaml_opts(args, opt_list, paste0(p.name,'.yaml'))

## print(opt)

if (length(opt$files)==1) {
    file_list <- strsplit(opt$files,',')[[1]]
} else {
    file_list <- opt$files
}


## this section of code is more complicated because there are old stats files with "oscode count clust_cnt"  (e.g. oscode_ome_clust_cnt.tsv)
## but also the more complete stats file that also contains fa_stats  (outdir8_incsurv_cm0.tsv)

os_stats <- NULL
have_fa_stats <- FALSE

if (!is.na(opt$stats) && length(opt$stats)>0) {
   os_stats <- read.table(opt$stats, header=TRUE, sep='\t')
   if (opt$debug) {
      print("os_stats")
      print(head(os_stats))
   }

   os_stats_names = names(os_stats)

   if ('OSCODE' %in% os_stats_names && 'proteomes_total' %in% os_stats_names) {
      ## change the names in os_stats_name_edits

      os_stats_name_edits = c('OSCODE'='oscode','proteomes_total'='count','clusters_50pct'='clust_cnt',
                              'ps_min'='fa_min','ps_q1'='fa_q1','ps_median'='fa_med','ps_q3'='fa_q3','ps_max'='fa_max','ps_mean'='fa_mean')

      new_os_names <- str_replace_all(os_stats_names, os_stats_name_edits)   
      names(os_stats) <- new_os_names
      have_fa_stats <- TRUE
      os_stats$clust_fn <- os_stats$clust_cnt/os_stats$fa_med
   }

   os_stats <- os_stats[order(-os_stats$count),]
   if (opt$debug) {
      print(head(os_stats))
   }
}

if (!have_fa_stats && !is.na(opt$fa_stats) && length(opt$fa_stats)>0) {
   fa_stats <- read.table(opt$fa_stats, header=TRUE, sep='\t')
   have_fa_stats <- TRUE

   if (opt$debug) {
      print("fa_stats")
      print(fa_stats)
   }

   if (have_fa_stats && nrow(os_stats) > 0 & nrow(fa_stats)) {
      os_stats <- merge(os_stats, fa_stats, by='oscode')
   }
   os_stats$clust_fn <- os_stats$clust_cnt/os_stats$fa_med
   print("*** os_stats/summary")
   ##print(os_stats[,c('oscode','count','ome_cnt','clust_cnt','clust_fn','fa_med','clust_fn')])
   ## print(summary(os_stats[,c('oscode','count','ome_cnt','clust_cnt','clust_fn','fa_med','clust_fn')]))
   print(os_stats)
   print(summary(os_stats))
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

## need to set up combinations of colors and symbols to label up to 20 taxa
## 5 colors
## 4 shapes (square=0, circle=1, triangle=2, diamond=5)

## this is brewer.pal(5, 'Set1)
## red, blue, green, purple, orange
ds_colors = rep(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"),4)

ds_shapes = rep(c(0,1,2,5),each=5)
lds_shapes = rep(c(15, 16, 17, 18), each=5)
##lds_sizes = 1.2*c(rep(2,15),rep(2.8,5))
lds_sizes = 1.2*c(rep(2,15),rep(2.8,5))

mode_good_fn = c()
mode_good_cnt = NULL

## table for all taxa
clust_data = NULL
clust_thresh = NULL

for (clust_file in file_list) {

    if (opt$debug) {
        print(clust_file)
    }
    
    if (file.exists(clust_file)) {

    this_taxon = strsplit(clust_file,'_')[[1]][1]
    taxon.names <- append(taxon.names, this_taxon)

    this_os_stat = os_stats[os_stats$oscode==this_taxon,]
##    print(this_os_stat[,c('oscode','count','clust_cnt','short_name','ome_cnt','fa_med')])

    this_data <- read.table(clust_file,header=TRUE,sep='\t')

##    print(head(this_data))

    this_data$taxon = this_taxon 

##    if (this_taxon == 'SALER') {
##       print(this_data[this_data$clust_id==95583,c('clust_id','n_omes','p_tot','p_bad')])
##    }

    new_nrow = nrow(this_data)

    if (opt$mode_good) {
        old_nrow = nrow(this_data)
        this_data <- this_data[this_data$pct_mode > opt$mode_pct_thresh,]
        new_nrow = nrow(this_data)
	
##	print("after mode_good")
##	print(head(this_data))

	mode_str = sprintf("%d\t%d\t%.1f",old_nrow, new_nrow,100.0*new_nrow/old_nrow)
	mode_str <- setNames(mode_str, this_taxon)

	mode_good_fn <- append(mode_good_fn, mode_str)
	if (opt$debug) {
	  print(sprintf("%s opt$mode_good filter: old: %d new: %d: %.1f%%\n",this_taxon, old_nrow, new_nrow,100.0*new_nrow/old_nrow))
	}
	mode_good_cnt <- rbind(mode_good_cnt, data.frame('taxon'=this_taxon, 'N_good_clust'=new_nrow))
        }
    }

    this_ome_cnt <- this_os_stat$count
##    print(paste("this_ome_cnt",this_ome_cnt))

    this_pct_mode_good <- new_nrow/old_nrow
##    print(paste("this_pct_mode_good",this_pct_mode_good))

    this_data$ome_pct_all <-  100.0 * this_data$n_omes / this_ome_cnt

    if (sum(this_data$ome_pct_all > 100.0)>0) {
        print(paste(this_taxon,"over 100: ",this_ome_cnt))
	print(head(this_data[this_data$ome_pct_all>100.,c('clust_id','n_omes','p_tot','ome_pct_all')]))
    }

    this_data$ome_pct_mode <- 100.0 * this_data$p_mode / this_ome_cnt

    clust_data <- rbind(clust_data, this_data)
}

## print("summary(clust_data)")
##  clust_data |> group_by(taxon) |>
##  	   summarize(across(c(ome_pct_all), 
## 	                    .fns=list(Min=min, Med=median, Max=max),na.rm=TRUE),
##	   .groups='drop')


if (opt$debug) {
   print("mode_good cluster count")
   print(mode_good_cnt)
   print(summary(mode_good_cnt))
}

clust_data$pct_bad_nz <- ifelse(clust_data$pct_bad > 0, clust_data$pct_bad, 0.001)

## taxa_medians <- aggregate(pct_bad~taxon, data=clust_data,FUN=median)
## print("taxa medians")
## print(taxa_medians)

clust_taxa_quants <- clust_data %>%
   group_by(taxon) %>%
   summarize(
     quants = list(quantile(pct_mode, probs=seq(0,1,0.25),na.rm=TRUE,))
   ) %>%
   unnest_wider(quants)

if (opt$debug) {
   print(clust_taxa_quants)
   print("summary(clust_taxa_quants$`50%`)")
   print(summary(clust_taxa_quants$`50%`))
}

if (opt$debug) {
   print(paste("clust 25% > 80",nrow(clust_taxa_quants[clust_taxa_quants$`25%`>80.0,])))
## print(clust_taxa_quants[clust_taxa_quants$`25%`>80.0,])
   print(paste("clust 25% > 90",nrow(clust_taxa_quants[clust_taxa_quants$`25%`>90.0,])))
## print(clust_taxa_quants[clust_taxa_quants$`25%`>90.0,])
   print(paste("clust 25% > 95",nrow(clust_taxa_quants[clust_taxa_quants$`25%`>95.0,])))
## print(clust_taxa_quants[clust_taxa_quants$`25%`>95.0,])
   print(paste("clust 50% > 95",nrow(clust_taxa_quants[clust_taxa_quants$`50%`>95.0,])))
## print(clust_taxa_quants[clust_taxa_quants$`50%`>95.0,])
   print(paste("clust 50% > 98",nrow(clust_taxa_quants[clust_taxa_quants$`50%`>98.0,])))
## print(clust_taxa_quants[taxa_quants$`50%`>98.0,])
}

all_taxa_quants <- clust_data %>%
   group_by(taxon) %>%
   summarize(
     quants = list(quantile(ome_pct_all, probs=seq(0,1,0.25),na.rm=TRUE,))
   ) %>%
   unnest_wider(quants)

if (opt$debug) {
   print(all_taxa_quants)
   print("summary(all_taxa_quants$`50%`)")
   print(summary(all_taxa_quants$`50%`))

   print(paste("all 25% > 80",nrow(all_taxa_quants[all_taxa_quants$`25%`>80.0,])))
## print(all_taxa_quants[all_taxa_quants$`25%`>80.0,])
   print(paste("all 25% > 90",nrow(all_taxa_quants[all_taxa_quants$`25%`>90.0,])))
## print(all_taxa_quants[all_taxa_quants$`25%`>90.0,])
   print(paste("all 50% > 95",nrow(all_taxa_quants[all_taxa_quants$`50%`>95.0,])))
## print(all_taxa_quants[all_taxa_quants$`50%`>95.0,])
   print(paste("all 50% > 98",nrow(all_taxa_quants[all_taxa_quants$`50%`>98.0,])))
## print(all_taxa_quants[all_taxa_quants$`50%`>98.0,])
}

mode_taxa_quants <- clust_data %>%
   group_by(taxon) %>%
   summarize(
     quants = list(quantile(ome_pct_mode, probs=seq(0,1,0.25),na.rm=TRUE,))
   ) %>%
   unnest_wider(quants)

if (opt$debug) {
   print(mode_taxa_quants)
   print("summary(mode_taxa_quants$`50%`)")
   print(summary(mode_taxa_quants$`50%`))

print(paste("ome 25% > 80",nrow(mode_taxa_quants[mode_taxa_quants$`25%`>80.0,])))
## print(mode_taxa_quants[mode_taxa_quants$`25%`>80.0,])
print(paste("ome 25% > 90",nrow(mode_taxa_quants[mode_taxa_quants$`25%`>90.0,])))
## print(mode_taxa_quants[mode_taxa_quants$`25%`>90.0,])
print(paste("ome 50% > 95",nrow(mode_taxa_quants[mode_taxa_quants$`50%`>95.0,])))
## print(mode_taxa_quants[mode_taxa_quants$`50%`>95.0,])
print(paste("ome 50% > 98",nrow(mode_taxa_quants[mode_taxa_quants$`50%`>98.0,])))
## print(mode_taxa_quants[mode_taxa_quants$`50%`>98.0,])
}

prot_missed <- clust_data %>% group_by(taxon) %>%
    summarize(
        N_miss_all60 = sum(ome_pct_all<60.0),
        pct_miss_all60 = 100.0*N_miss_all60/n(),
        N_miss_all75 = sum(ome_pct_all<75.0),
        pct_miss_all75 = 100.0*N_miss_all75/n(),
        N_miss_mode75 = sum(pct_mode<75.0),
        pct_miss_mode75 = 100.0*N_miss_mode75/n(),
        N_miss_all80 = sum(ome_pct_all<80.0),
        pct_miss_all80 = 100.0*N_miss_all80/n(),
        N_miss_mode80 = sum(pct_mode<80.0),
        pct_miss_mode80 = 100.0*N_miss_mode80/n(),
    ) %>% ungroup()

print("prot_missed")
print(prot_missed)
colSums(prot_missed[,c(2,4,6,8,10)])
summary(prot_missed)

short_stats <- clust_data %>% group_by(taxon) %>%
    summarize(
        count_S = sum(bad_SL_flag=='Short'),
        count_L = sum(bad_SL_flag=='Long'),
	count_Z_fct = sum(p_bad < 1)/n(),
	bad_fct_med = median(pct_bad),
	short_fct = ifelse(count_S + count_L > 0, count_S/(count_S + count_L), 0.0)
    ) %>%
    arrange (desc(short_fct)) %>%
    ungroup()

   print("short_stats")
   print(short_stats)
   print("summary(fraction of outliers short)")
   summary(short_stats$short_fct)

## taxa.ordered = taxa_medians[order(c(taxa_medians$pct_bad)),]$taxon
## order taxa by decreasing proteome count
taxa.ordered = os_stats[order(-os_stats$count),]$oscode
## taxa.ordered

short_stats$taxon <- factor(short_stats$taxon, levels=taxa.ordered, ordered=TRUE)

clust_data$taxon <- factor(clust_data$taxon, levels=taxa.ordered, ordered=TRUE)
clust_thresh$taxon <- factor(clust_thresh$taxon, levels=taxa.ordered, ordered=TRUE)

clust_thresh$thresh <- factor(clust_thresh$thresh, levels=c('pct0','pct0.1','pct1.0','pct10.0','pct100'), ordered=TRUE)

ds_colors = ds_colors[1:length(taxa.ordered)]
names(ds_colors) = taxa.ordered
##print(ds_colors)

ds_shapes = ds_shapes[1:length(taxa.ordered)]
names(ds_shapes) = taxa.ordered
##print(ds_shapes)

lds_shapes = lds_shapes[1:length(taxa.ordered)]
lds_sizes = lds_sizes[1:length(taxa.ordered)]
names(lds_shapes) = taxa.ordered
##print(lds_shapes)

## print(taxa.ordered)
## print(ds_colors)

q_color = scale_color_manual(values=ds_colors)
q_shape = scale_shape_manual(values=ds_shapes,labels=taxa.ordered)
lq_shape = scale_shape_manual(values=lds_shapes,labels=taxa.ordered)
lq_size = scale_size_manual(values=lds_sizes,labels=taxa.ordered)

theme_set(theme_linedraw(base_size=14))
theme.leg_id <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='longdash'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='longdash'),
  plot.title=element_text(face='plain', hjust=0,size=12),
#  plot.title=element_blank(),
  plot.title.position='plot',
  plot.margin=margin(t=2,b=2,r=10),
  plot.caption=element_text(size=9,hjust=0),
  axis.text.x = element_text(size=11),
  axis.title.x = element_text(size=12),
  axis.text.y = element_text(size=11),
  axis.title.y = element_text(size=12),
  legend.key=element_blank(),
  legend.position.inside=c(0.55,0.4),
  legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
  legend.text=element_text(size=9,hjust=0),
  legend.justification=c(0,1),
  legend.title=element_blank())

cs_guide <- guide_legend(nrow=10, override.aes=list(shape=lds_shapes, size=lds_sizes, alpha=1.0))


top_delta = 0.5
top100_delta = 100.0-top_delta
log50_delta = log10(top_delta) - log10(50)

xtop_scale = scale_x_discrete(position='top')

ylog_scale_bad = scale_y_log10('percent',limits=c(0.001, 100.0), breaks=10^seq(-3,2.0,1.0),labels=c('0.00','0.01','0.1','1.0','10.0','100'))
ypct_scale_bad = scale_y_continuous(limits=c(0.0, 100.0), breaks=seq(0,100,20))

yprot_scale_bad60 = scale_y_continuous(limits=c(60, 100.0))

yprot_scale_bad60_log = scale_y_continuous(limits=c(50, 100.0), breaks=c(60,75,90,95,98,99,100),labels=c('60','75','90','95','98','99','100'))
yprot_scale_bad60_log = scale_y_continuous(limits=c(50, 100.0), breaks=c(60,75,90,95,98,99,100),labels=c('60','75','90','95','98','99','100'))

ypct_scale_short = scale_y_continuous(limits=c(75.0, 100.0), breaks=seq(75.0,100,5))

## str(clust_data)

sp_plot_margin = c(0.1, 0.1, 0.25, 0.1)*1.00


taxon_medians <- clust_data %>% group_by(taxon) %>%
                 summarize(
	            all_median = median(ome_pct_all),
	            all_q1 = fivenum(ome_pct_all)[2],
		    ome_mode_median = median(ome_pct_mode),
		    ome_mode_q1 = fivenum(ome_pct_mode)[2],
		    clust_mode_median = median(pct_mode),
		    clust_mode_q1 = fivenum(ome_pct_mode)[2],
		    clust_bad_median = median(pct_bad)
		 )

print("taxon_medians")
print(taxon_medians)

taxon_med_meds <- taxon_medians %>% summarize(across(where(is.numeric), \(x) median(x,na.rm=TRUE)))
print(taxon_med_meds)

## panel A. -- pct of omes in clusters

p_clust_all_box <- ggplot(clust_data, aes(x=taxon, y=ome_pct_all, color=taxon)) + geom_boxplot(outlier.shape=NA) + q_color + yprot_scale_bad60 + 
             geom_hline(yintercept=taxon_med_meds$all_median,linetype='longdash') +
	     theme.leg_id +
	     ggtitle("A. proteomes per cluster ") + ## xtop_scale +
	     labs(y='percent') + 
	     theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1.0, vjust=1.0, color=ds_colors), legend.position='None') +
     	     theme(plot.margin=unit(sp_plot_margin,'cm'))

print(paste("med_omes/clust",taxon_med_meds$all_median))


## panel B -- pct of proteins at mode / cluster

p_clust_mode_box <- ggplot(clust_data, aes(x=taxon, y=pct_mode, color=taxon)) + 
		 geom_boxplot(outlier.shape=NA) + q_color + yprot_scale_bad60 + 
		 geom_hline(yintercept=taxon_med_meds$clust_mode_median,linetype='longdash') +
		 theme.leg_id + 
   	         ggtitle("B. mode proteins per cluster") +
	         labs(y='percent') + xtop_scale +
	         theme(axis.title.x=element_blank(),
		       axis.text.x=element_blank(),
		       axis.ticks.x=element_blank(),
		       legend.position='None') +
	         theme(plot.margin=unit(sp_plot_margin,'cm'))

print(paste("med_at_mode",taxon_med_meds$clust_mode_median))

## panel C -- pct of outliers / cluster

p_bad_box <- ggplot(clust_data, aes(x=taxon, y=pct_bad_nz, color=taxon)) + 
	  geom_boxplot(outlier.shape=5,outlier.alpha=0.2,outlier.size=2.0) + q_color + ylog_scale_bad + 
	  theme.leg_id + 
	  ggtitle("C. outliers per cluster") + labs(y='percent') +
          geom_hline(yintercept=taxon_med_meds$clust_bad_median,linetype='longdash') +
	  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position='None') +
	  theme(plot.margin=unit(sp_plot_margin,'cm'))

print(paste("med_outliers:",taxon_med_meds$clust_bad_median))

med_not_bad <- 100.0*median(short_stats$count_Z_fct)
print(paste("med_not_bad:",med_not_bad))

## panel D -- no outlier clusters

p_not_bad <- ggplot(short_stats, aes(x=taxon, y=100.0*count_Z_fct, color=taxon)) + geom_point(aes(shape=taxon,size=taxon)) + ## geom_point(shape=3,size=3) + 
	     geom_hline(yintercept=med_not_bad, linetype="longdash") +
	     theme.leg_id + 
	  q_color + ypct_scale_bad + lq_shape + lq_size +
          ggtitle("D. no outlier clusters") +
	  labs(y='percent') + 
	  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position='None') +
          theme(plot.margin=unit(sp_plot_margin,'cm'))

med_short_pct <- 100.0*median(short_stats$short_fct)
print(paste("med_short_pct:",med_short_pct))

p_bad_short <- ggplot(short_stats, aes(x=taxon, y=100.0*short_fct, color=taxon)) + geom_point(aes(shape=taxon,size=taxon)) + ## geom_point(shape=3,size=3) +
	     geom_hline(yintercept=med_short_pct, linetype="longdash") +
	     theme.leg_id + 
  	     q_color + ypct_scale_short + lq_shape + lq_size +
             ggtitle("E. short outliers per cluster") +
	     labs(y='percent') + 
	     theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1.0, vjust=1.0, color=ds_colors), legend.position='None') +
	     theme(plot.margin=unit(sp_plot_margin,'cm'))

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (!is.na(opt$pdf)) {
   pdf_file = opt$pdf
} else {
   pdf_file = paste0('bfig_',taxa.ordered[1],"_len_dist.pdf")
}

p_not_bad_v <- p_not_bad

  n_panel <- 5
  p_pct_cnt1 <- ( p_clust_all_box / p_clust_mode_box / p_bad_box / p_not_bad_v / p_bad_short) 

rel_heights=c(5,4,5,4,4)


if (! opt$pub) {

   p_pct_cnt1 <- (p_pct_cnt1 / doc_panel) + plot_layout(ncol=1,heights=c(rel_heights,0.1))

   ggsave(file=pdf_file, plot=p_pct_cnt1, width=5.0, height= 1.2 +(n_panel * 1.8))

} else {

   p_pct_cnt1 <- p_pct_cnt1 + plot_layout(ncol=1,heights=rel_heights)

   if (!grepl('_pub',pdf_file)) {
      pdf_file_pref <- strsplit(pdf_file,'\\.')[[1]][1]
      pdf_file <- paste0(pdf_file_pref,"_pub.pdf")
   }

   ggsave(file=pdf_file, plot=p_pct_cnt1, width=5.0, height= 1.0 +(n_panel * 1.8))
}

## warnings()
