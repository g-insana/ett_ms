#!/usr/bin/env Rscript --vanilla

################
##
## miss_stats.R *_miss_tfx.summ
##
## calculate statistics of properties of missing proteins
##
##

library('ggplot2')
library('dplyr')
library('stringr')
library('patchwork')
library('getopt')
library('yaml')
library('optparse')
library('ggExtra')

get_yaml_file = 'get_yaml_opts.R'

if (file.exists(get_yaml_file)) {
   source(get_yaml_file)
} else if (file.exists(paste0("../",get_yaml_file))) {
   source(paste0("../",get_yaml_file))
} else {
  cat(paste("cannot open", get_yaml_file))
  sys.exit(1)
}

## read files of the form:

p.name<-'miss_stats'
p.full_name<-paste0(p.name,'.R')

cmd_args<- commandArgs(trailingOnly=TRUE)

## print(cmd_args)

plabel=paste(c(p.full_name,"\n",cmd_args),collapse=' ', sep=' ')

option_list = list(
  	    make_option(c("-D","--debug"),action='store_true',help='debug',default=FALSE),
  	    make_option(c("-f","--files"),type='character',action='store',help='comma separated *.summ, files REQUIRED'),
	    make_option(c("-P","--pdf"),type='character',action='store',help='PDF file name', default=NA),
	    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file', default=NA),
	    make_option(c("--pub"),action='store_true',help='publication plot', default=FALSE),
	    make_option(c("-S","--stats"),type='character',action='store',help='os_stats_file', default=NA),
	    make_option(c("-J","--jitter"),action='store_true',help='jitter', default=FALSE),
	    make_option(c("--a_label"),action='store',help='panel a label', default=NA),
	    make_option(c("--b_label"),action='store',help='panel b label', default=NA)
	    )

opt <- get_yaml_opts(args, option_list, paste0(p.name,'.yaml'))

f_names=c('file','X1cnt','Ncnt','tot_cnt','pct_good','pct_frame')
q_fields=c('X1cnt','Ncnt','tot_cnt','pct_good','pct_frame')


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

miss_stats_df = NULL

taxon.names=c()

for (tfx_file in opt$files) {

    this_taxon = strsplit(tfx_file,'_')[[1]][1]
    taxon.names <- append(taxon.names, this_taxon)

##    print(paste(tfx_file, this_taxon))

    this_stats_df = read.table(tfx_file,sep='\t',header=TRUE)

##    print(paste(this_taxon, nrow(this_stats_df)))
    
    if (nrow(this_stats_df)) {
        this_stats_df$taxon = this_taxon
	this_stats_names = names(this_stats_df)
    } else {
##       this_stats_df = data.frame('file'=this_taxon,'X1cnt'=0, 'X1fs'=0, 'X1T'=0, 'X1len'=0, 
##                                  'Ncnt'=0, 'Nfs'=0, 'NT'=0, 'NfsT'=0, 'Nlen'=0'tot_cnt'=0, 'pct_good'=100.0, 'pct_frame'=0.0, 'taxon'=this_taxon, 
         
	 this_stat_vec = rep(0, length(this_stats_names))
	 names(this_stat_vec) = this_stats_names

	 this_stats_df <- data.frame(t(this_stat_vec))
	 this_stats_df[,c('file','taxon')] <- this_taxon
##	 print(this_stats_df)
    }
    miss_stats_df = rbind(miss_stats_df,this_stats_df)
}

taxa.ordered = os_stats[order(-os_stats$count),]$oscode
miss_stats_df$taxon <- factor(miss_stats_df$taxon, levels=taxa.ordered, ordered=TRUE)

## print(tapply(miss_stats_df[,q_fields], miss_stats_df$taxon, summary),width=200)

## get medians across the various measures

taxon_medians <- miss_stats_df %>% group_by(taxon) %>%
                 summarize(
		    n_samp = n(), 
		    q1_good = fivenum(pct_good)[2],
	            pct_good = median(pct_good),
		    q3_good = fivenum(pct_good)[4],
	            pct_fs = median(pct_frame)
		 )

print("taxon_medians")
print(taxon_medians)

taxon_med_meds <- taxon_medians %>% summarize(across(where(is.numeric), \(x) median(x,na.rm=TRUE)))
print(taxon_med_meds)

################
## do some box plots

## this is brewer.pal(5, 'Set1)
## red, blue, green, purple, orange
ds_colors = rep(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"),4)

ds_shapes = rep(c(0,1,2,5),each=5)
lds_shapes = rep(c(15, 16, 17, 18), each=5)
##lds_sizes = 1.2*c(rep(2,15),rep(2.8,5))
lds_sizes = 1.2*c(rep(2,15),rep(2.8,5))

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

xtop_scale = scale_x_discrete(position='top')

ypct_scale = scale_y_continuous(limits=c(0.0, 100.0), breaks=seq(0,100,20))

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

if (is.na(opt$a_label)){
  a_label="A. missing proteins recovered "
} else {
  a_label=opt$a_label
}

if (is.na(opt$b_label)){
  b_label="B. missing proteins caused by frame-shifts "
} else {
  b_label=opt$b_label
}


p_pct_miss <- ggplot(miss_stats_df, aes(x=taxon, y=pct_good, color=taxon)) + geom_boxplot(outlier.shape=NA) + q_color + ypct_scale + 
             geom_hline(yintercept=taxon_med_meds$pct_good,linetype='longdash') +
	     theme.leg_id +
	     ggtitle(a_label) + ## xtop_scale +
	     labs(y='percent') + 
	     theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1.0, vjust=1.0, color=ds_colors), legend.position='None')

if (opt$jitter) {
    p_pct_miss <- p_pct_miss + geom_jitter(alpha=0.5,shape=1,size=1)
}

p_pct_fs <- ggplot(miss_stats_df, aes(x=taxon, y=pct_frame, color=taxon)) + geom_boxplot(outlier.shape=NA) + q_color + ypct_scale + 
             geom_hline(yintercept=taxon_med_meds$pct_fs,linetype='longdash') +
	     theme.leg_id +
	     ggtitle(b_label) + ## xtop_scale +
	     labs(y='percent') + 
	     theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45, hjust=1.0, vjust=1.0, color=ds_colors), legend.position='None')
if (opt$jitter) {
    p_pct_fs <- p_pct_fs + geom_jitter(alpha=0.5,shape=1, size=1)
}

## print(paste("med_pct_good",taxon_med_meds$all_median))

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

n_panel = 2

if (!is.na(opt$pdf)) {
   pdf_file = opt$pdf
} else {
   pdf_file = paste0('fig_',taxa.ordered[1],"_missed_dist.pdf")
}

p_pct_all <- p_pct_miss / p_pct_fs

rel_heights = c(5,5)


if (! opt$pub) {

   p_pct_all <- (p_pct_all / doc_panel) + plot_layout(ncol=1,heights=c(rel_heights,0.1))

   ggsave(file=pdf_file, plot=p_pct_all, width=5.0, height= 1.2 +(n_panel * 1.8))

} else {

   p_pct_all <- p_pct_all + plot_layout(ncol=1,heights=rel_heights)

   if (!grepl('_pub',pdf_file)) {
      pdf_file_pref <- strsplit(pdf_file,'\\.')[[1]][1]
      pdf_file <- paste0(pdf_file_pref,"_pub.pdf")
   }

   ggsave(file=pdf_file, plot=p_pct_all, width=5.0, height= 1.0 +(n_panel * 1.8))
}
