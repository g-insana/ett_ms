#!/usr/bin/env Rscript --vanilla

################
## align_props_box_short.R all_*.summ_4s all_*.summ_l2
## derived from comp_summ2.R, but takes per-cluster (proteome) data for each taxa and builds box plot

## align_props_box_short.R --files *.tfxg_stats_S1K3ss,*_long_clust.tfxg_stats
##

################
## align_props_short.R takes files of the form:
## *.tfxg_stats_S1K3ss
## file	                        logN50	X1cnt	X1fs	X1T	X1fsT	X1len	X1Next	X1Cext	Ncnt	Nfs	NT	NfsT	Nlen	NNext	NCext	X0cnt	X0len	X0Next	X0Cext	tot_hit	tot_cnt	pct_good	pct_1good	pct_fs
## 9ENTR_GCA_030358845.1	4.14	1	1	0	0	10168	1	0	4	0	0	0	5367	3	1	15	3390	7	10	5	20	25.0	20.0	20.0
## 9ENTR_GCA_030358725.1	4.46	0	0	0	0	0	0	0	6	0	0	0	11742	6	0	5	4275	4	2	6	11	54.5	0.0	0.0
## 9ENTR_GCA_001518515.1	4.79	1	0	0	0	74896	1	0	11	0	0	0	38704	6	6	1	50085	1	0	12	13	92.3	8.3	0.0
## 9ENTR_GCA_016428285.1	4.83	1	1	0	0	8239	1	0	0	0	0	0	0	0	0	2	6680	2	0	1	3	33.3	100.0	100.0

## *_long_clust.tfxg_stats
## file         	logN50	X1cnt	X1fs	X1T	X1fsT	X1len	X1Next	X1Cext	Ncnt	Nfs	NT	NfsT	Nlen	NNext	NCext	X0cnt	X0len	X0Next	X0Cext	tot_hit	tot_cnt	pct_good	pct_1good	pct_fs
## 9ENTR_cl_24258	-1.00	20	6	20	6	137656	0	20	0	0	0	0	0	0	0	0	0	0	0	20	20	100.0	100.0	30.0
## 9ENTR_cl_44223	-1.00	20	1	20	1	110665	20	0	0	0	0	0	0	0	0	0	0	0	0	20	20	100.0	100.0	5.0
## 9ENTR_cl_64383	-1.00	20	0	0	0	139274	20	0	0	0	0	0	0	0	0	0	0	0	0	20	20	100.0	100.0	0.0

## and generates a set of box plots for the 20 proteomes looking at:
##
## hit fraction (recovery of full length for short, recovery of long in short proteome)
## 1cnt vs Ncnt
## Term vs fs
## Next vs Cext

## it might make sense to plot these with the proteome (ordered by fraction short) on the x-axis, with 4 panels for the 4 measures, left short, right long
## (need ordering based on 1/2 transition), which can be found in oscode_f50.tab

library('ggplot2')
library('dplyr')

library('stringr')
library('patchwork')
library('optparse')
library('yaml')

get_yaml_file = 'get_yaml_opts.R'

if (file.exists(get_yaml_file)) {
   source(get_yaml_file)
} else if (file.exists(paste0("../",get_yaml_file))) {
   source(paste0("../",get_yaml_file))
} else {
  cat(paste("cannot open", get_yaml_file))
  sys.exit(1)
}

p.name<-'align_props_box_short'
p.full_name <- paste0(p.name,'.R')

args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.full_name,"\n",args),collapse=' ', sep=' ')

opt_list = list(
    make_option(c("--files"),type='character',action='store',help='comma separated *.bad_clust, files REQUIRED',default=NA),
    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file', default=NA),
    make_option(c("-P","--pdf"),type='character',action='store',help='PDF file name', default=NA),
    make_option(c("-S","--stats"),type='character',action='store',help='os_stats_file', default=NA),
    make_option(c("--pub"),action='store_true',help='publication plot', default=FALSE)
)

opt <- get_yaml_opts(args, opt_list, paste0(p.name,'.yaml'))

if (length(opt$files)==1) {
    file_list <- strsplit(opt$files,',')[[1]]
} else {
    file_list <- opt$files
}

## print(opt)
## print(file_list)

os_stats <- NULL
if (!is.na(opt$stats) && length(opt$stats)>0) {
   os_stats <- read.table(opt$stats, header=TRUE, sep='\t')
   print("head(os_stats)")

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
   print(head(os_stats))
}

summ_data <- NULL
leg.labels.oscode = c()
leg.labels.nrow = c()

## need to set up combinations of colors and symbols to label up to 20 taxa
## 5 colors
## 4 shapes (square=0, circle=1, triangle=2, diamond=5)

## for short alignments
ds_colors = rep(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"),4)
ds_shapes = c(rep(c(0,1,2,5),each=5))
ds_sizes = rep(2,20)

## for long alignments
lds_colors = rep(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"),4)
lds_shapes = rep(c(15, 16, 17, 18), each=5)
lds_sizes = c(rep(2,15),rep(2.8,5))

## use list of summ df's for analysis, to start, just read one
##

stats_df = NULL

for (stats_file in file_list) {

    print(stats_file)
    if (file.exists(stats_file)) {

       ts_df <- read.table(stats_file,header=TRUE,sep='\t')

       if (sum(str_detect(stats_file,'long_clust'))>0){
           this_s_type = 'long'
       } else {
           this_s_type = 'short'
       }

       this_taxon = strsplit(stats_file, '_')[[1]][1]

       print(head(ts_df))

       if (nrow(ts_df) > 0) {
           ts_df$s_type = this_s_type
           ts_df$taxon = this_taxon
	   sum_exts <- with(ts_df,X1Next+NNext+X1Cext+NCext)

           ts_df <- within(ts_df, {
	                   pct_Next <- ifelse(sum_exts > 0, 100.0*(X1Next+NNext)/sum_exts, 0.0);
	                   pct_err  <- ifelse(tot_hit>0, 100.0*(X1fs + X1T + Nfs + NT - X1fsT - NfsT)/tot_hit, 0.0)
			   })

	} else {
	   print(sprintf("*** %s *** no rows"))
	   ts_df <- data.frame('taxon'=this_taxon, 's_type'=this_s_type, 'pct_good'=0.0, 'pct_1good'=0.0, 'pct_fs'=0.0, 'pct_err'=0.0, 'pct_Next'=0.0)
        }	   

     print("head(ts_df)")
     print(head(ts_df[,c('taxon','s_type','pct_good','pct_1good','pct_fs','pct_Next', 'pct_err')]))

     stats_df <- rbind(stats_df,ts_df)
    }
}

if (!is.null(stats_df) && nrow(stats_df) > 0) {
  stats_df <- within(stats_df, {taxon_sl<-paste0(taxon,'_',s_type)})
}

summary(stats_df)

taxa.ordered = os_stats$oscode

## print(taxa.ordered)

stats_df$taxon <- factor(stats_df$taxon, levels=taxa.ordered, ordered=TRUE)
stats_df$s_type <- factor(stats_df$s_type, levels=c('short','long'),ordered=TRUE)

tax_sl_levels = c(paste0(taxa.ordered,"_short"),paste0(taxa.ordered,"_long"))
##print(tax_sl_levels)

stats_df$taxon_sl <- factor(stats_df$taxon_sl, levels=tax_sl_levels, ordered=TRUE)

## print("stats_df")
## print(stats_df, digits=4)
## print(summary(stats_df[as.character(stats_df$s_type)=='short',]))
## print(summary(stats_df[as.character(stats_df$s_type)=='long',]))

##names(ds_colors) = paste0(taxa.ordered,'_short')
names(ds_colors) = taxa.ordered
names(lds_colors) = paste0(taxa.ordered,'_long')

## ds_shapes = ds_shapes[1:length(taxon.names)]
## names(ds_shapes) = paste0(taxa.ordered,'_short')
names(ds_shapes) = taxa.ordered
##names(ds_sizes) = paste0(taxa.ordered,'_short')
names(ds_sizes) = taxa.ordered

names(lds_shapes) = paste0(taxa.ordered,'_long')
names(lds_sizes) = paste0(taxa.ordered,'_long')

## a_ds_colors <- c(ds_colors, lds_colors)
a_ds_colors <- ds_colors
## a_ds_shapes <- c(ds_shapes, lds_shapes)
a_ds_shapes <- ds_shapes
## a_ds_sizes <- c(ds_sizes, lds_sizes)
a_ds_sizes <- ds_sizes

## print(a_ds_colors)
## print(a_ds_shapes)
## print(a_ds_sizes)

q_color = scale_color_manual(values=a_ds_colors,labels=taxa.ordered)
q_shape = scale_shape_manual(values=a_ds_shapes,labels=taxa.ordered)
q_size = scale_size_manual(values=a_ds_sizes,labels=taxa.ordered)

theme_set(theme_linedraw(base_size=14))
theme.leg_id <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_text(face='plain', hjust=0,size=12),
  ## plot.margin=margin(t=2,r=8,b=2),
  plot.caption=element_text(size=9,hjust=0),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.text.y = element_text(size=11),
  axis.title.y = element_text(size=12),
  legend.key=element_blank(),
  legend.text=element_text(size=9,hjust=0),
  legend.title=element_blank())

# cs_guide <- guide_legend(nrow=10, override.aes=list(shape=lds_shapes, size=lds_sizes, alpha=1.0))
cs_guide <- guide_legend()

cs_guide <- guide_legend(nrow=5,override.aes=list(alpha=1.0,size=1.5))
this_cs_guide = guides(color=cs_guide)

y_scale = scale_y_continuous(breaks=100.0 * seq(0,1.0, 0.2),limits= 100.0 *c(0.0,1.0))
y_scale_r = scale_y_continuous(breaks=100.0 * seq(0,1.0,0.2),limits=100.0 * c(0.0,1.0),position='right')
## y_scale_r = y_scale

x_scale = scale_x_continuous('proteome fraction long',breaks=seq(0,1.0,0.2))

x_var = 'taxon'

## stats_df_short <- stats_df[stats_df$s_type == 'short',]
## stats_df_long <- stats_df[stats_df$s_type == 'long',]

if (is.null(stats_df) | nrow(stats_df)==0) {
  print("no stats_df")
  quit()
}

taxa_medians <- stats_df %>% group_by(taxon) %>%
               	  summarize(
		     hit_mean = mean(pct_good),
		     hit1_mean = mean(pct_1good),
		     fs_mean = mean(pct_fs),
		     err_mean = mean(pct_err),
		     Next_mean = mean(pct_Next),

		     hit_median = median(pct_good),
		     hit1_median = median(pct_1good),
		     fs_median = median(pct_fs),
		     err_median = median(pct_err),
		     Next_median = median(pct_Next))
print(taxa_medians)

taxa_med_meds <- taxa_medians %>% 
                   summarize(
		      hit_mean_med = median(hit_mean),
		      hit1_mean_med = median(hit1_mean),
		      fs_mean_med = median(fs_mean),
		      err_mean_med = median(err_mean),
		      Next_mean_med = median(Next_mean),

		      hit_med_med = median(hit_median),
		      hit1_med_med = median(hit1_median),
		      fs_med_med = median(fs_median),
		      err_med_med = median(err_median),
		      Next_med_med = median(Next_median))
print(taxa_med_meds[,1:5])
print(taxa_med_meds[,6:10])

ave_hit_fn <- mean(stats_df$pct_good)
med_hit_fn <- median(stats_df$pct_good)
print(paste("ave hit:",ave_hit_fn, med_hit_fn))

ave_hit_fn=unlist(taxa_med_meds[1,1])
med_hit_fn=unlist(taxa_med_meds[1,6])
print(paste("ave hit:",ave_hit_fn, med_hit_fn))

print(paste("nrow(stats_df):",nrow(stats_df)))

nrow_cnt <- stats_df %>% group_by(taxon) %>% summarize(nrows=n())
print(nrow_cnt)
print(summary(nrow_cnt))

symbol_alpha <- 0.25
if (nrow(stats_df) < 500) {
   print(paste("changing symbol_alpha",nrow(stats_df)))
   symbol_alpha=1.0
}

p_hit_fn <- ggplot(stats_df,aes(x=.data[[x_var]], y=pct_good, color=taxon, shape=taxon)) +
  ggtitle("A. long recovered") +
  theme.leg_id + 
  this_cs_guide +
  geom_boxplot(outlier.shape=NA,alpha=0.8,linewidth=0.5, median.linewidth=1) + 
##  geom_jitter(alpha=symbol_alpha,size=1, width=0.33,height=0.0) + 
##  stat_summary(fun=mean, geom='point', aes(size=taxon)) +
##  stat_summary(fun.data=mean_se, geom='errorbar',width=0.2) +
  geom_hline(yintercept=ave_hit_fn,linetype="dashed") + 
  geom_hline(yintercept=med_hit_fn,linetype="longdash") + 
  y_scale + q_color + q_shape + q_size + 
  ylab("percent")


ave_X1_fn <- mean(stats_df$pct_1good)
med_X1_fn <- median(stats_df$pct_1good)
print(paste("ave X1:",ave_X1_fn, med_X1_fn))

ave_X1_fn=unlist(taxa_med_meds[1,2])
med_X1_fn=unlist(taxa_med_meds[1,7])
print(paste("ave X1:",ave_X1_fn, med_X1_fn))


p_X1_fn <- ggplot(stats_df,aes(x=.data[[x_var]], y=pct_1good, color=taxon, shape=taxon)) +
  ggtitle("B. one contig") +
  theme.leg_id + 
  this_cs_guide +
  geom_boxplot(outlier.shape=NA,alpha=0.8,linewidth=0.5, median.linewidth=1) + 
##  geom_jitter(alpha=symbol_alpha,size=1, width=0.33,height=0.0) + 
##  stat_summary(fun=mean, geom='point', aes(size=taxon)) +
##  stat_summary(fun.data=mean_se, geom='errorbar',width=0.2) +

  geom_hline(yintercept=ave_X1_fn,linetype="dashed") + 
  geom_hline(yintercept=med_X1_fn,linetype="longdash") + 
  y_scale + q_color + q_shape + q_size + 
  ylab("percent")

ave_fs_fn <- mean(stats_df$pct_fs)
med_fs_fn <- median(stats_df$pct_fs)
print(paste("ave fs:",ave_fs_fn,med_fs_fn))

ave_fs_fn=unlist(taxa_med_meds[1,3])
med_fs_fn=unlist(taxa_med_meds[1,8])
print(paste("ave fs:",ave_fs_fn,med_fs_fn))

p_fs_fn <- ggplot(stats_df,aes(x=.data[[x_var]], y=pct_fs, color=taxon, shape=taxon)) +
  ggtitle("C. frame-shifts") +
  theme.leg_id + 
  this_cs_guide +
  geom_boxplot(outlier.shape=NA,alpha=0.8,linewidth=0.5, median.linewidth=1) + 
##  geom_jitter(alpha=symbol_alpha,size=1, width=0.33,height=0.0) + 
##  stat_summary(fun=mean, geom='point', aes(size=taxon)) +
##  stat_summary(fun.data=mean_se, geom='errorbar',width=0.2) +
  geom_hline(yintercept=ave_fs_fn,linetype="dashed") + 
  geom_hline(yintercept=med_fs_fn,linetype="longdash") + 
  y_scale + q_color + q_shape + q_size + 
  ylab("percent")

ave_err_fn <- mean(stats_df$pct_err)
med_err_fn <- median(stats_df$pct_err)
print(paste("ave error:",ave_err_fn, med_err_fn))

ave_err_fn=unlist(taxa_med_meds[1,4])
med_err_fn=unlist(taxa_med_meds[1,9])
print(paste("ave error:",ave_err_fn, med_err_fn))

p_err_fn <- ggplot(stats_df,aes(x=.data[[x_var]], y=pct_err, color=taxon, shape=taxon)) +
  ggtitle("D. errors (fs+Term)" ) +
  theme.leg_id + 
  this_cs_guide +
  geom_boxplot(outlier.shape=NA,alpha=0.8,linewidth=0.5, median.linewidth=1) + 
##  geom_jitter(alpha=symbol_alpha,size=1, width=0.33,height=0.0) + 
##  stat_summary(fun=mean, geom='point', aes(size=taxon)) +
##  stat_summary(fun.data=mean_se, geom='errorbar',width=0.2) +
  geom_hline(yintercept=ave_err_fn,linetype="dashed") + 
  geom_hline(yintercept=med_err_fn,linetype="longdash") + 
  y_scale + q_color + q_shape + q_size + 
  ylab("percent")

ave_Next_fn <- mean(stats_df$pct_Next)
med_Next_fn <- median(stats_df$pct_Next)
print(paste("ave Next:",ave_Next_fn,med_Next_fn))

ave_Next_fn=unlist(taxa_med_meds[1,5])
med_Next_fn=unlist(taxa_med_meds[1,10])
print(paste("ave Next:",ave_Next_fn,med_Next_fn))

p_Next_fn <- ggplot(stats_df,aes(x=.data[[x_var]], y=pct_Next, color=taxon, shape=taxon)) +
  ggtitle("E. N-extensions") +
  theme.leg_id + 
  this_cs_guide +
  geom_boxplot(outlier.shape=NA,alpha=0.8,linewidth=0.5, median.linewidth=1) + 
##  geom_jitter(alpha=symbol_alpha,size=1, width=0.33,height=0.0) + 
##  stat_summary(fun=mean, geom='point', aes(size=taxon)) +
##  stat_summary(fun.data=mean_se, geom='errorbar',width=0.2) +
  geom_hline(yintercept=ave_Next_fn,linetype="dashed") + 
  geom_hline(yintercept=med_Next_fn,linetype="longdash") + 
  y_scale + q_color + q_shape + q_size + 
  ylab("percent")

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (!is.na(opt$pdf)) {
   pdf_file = opt$pdf
} else {
   pdf_file = paste0('lfig_',taxon.names[1],"_len_dist.pdf")
}

p_pct_cnt1 <- p_hit_fn + p_X1_fn + p_fs_fn + p_err_fn + p_Next_fn + guide_area() + plot_layout(guides='collect', ncol=2)

if (! opt$pub) {
   p_pct_cnt1 <-  p_pct_cnt1 / doc_panel + plot_layout(heights=c(5,0.1))
   ggsave(file=pdf_file, plot=p_pct_cnt1, width=9.0, height=7.2)
} else {
   ggsave(file=pdf_file, plot=p_pct_cnt1, width=9.0, height=7)
}

warnings()
