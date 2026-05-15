#!/usr/bin/env Rscript --vanilla

################
## clust_one_dist.R --pdf --files=9ENTR_c0.5_p0.90_cm2.bad_clust,ACIBA_c0.5_p0.90_cm2.bad_clust,...
##

## 6-Sept-2025 -- reads a file of protein lengths

################
## reads in *.clust

## 1205927:KOZ37682	472
## 1207527:KOZ89088	472
## 1207753:KPH34246	472
## 1207994:KOZ65201	472
## 1208965:KOZ04643	472
## 1209275:KOZ57974	472

################
## add plot to examine cluster quality vs cluster length

library('ggplot2')
library('ggtext')
library('optparse')
library('yaml')
library('dplyr')
library('purrr')
library('RColorBrewer')
library('patchwork')

get_yaml_file = 'get_yaml_opts.R'

if (file.exists(get_yaml_file)) {
   source(get_yaml_file)
} else if (file.exists(paste0("../",get_yaml_file))) {
   source(paste0("../",get_yaml_file))
} else {
  cat(paste("cannot open", get_yaml_file))
  sys.exit(1)
}

## clust_one_dist.R 

p.name<-'clust_one_dist'
p.full_name <- paste0(p.name,'.R')

args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.full_name,"\n",args),collapse=' ', sep=' ')

option_list = list(
	    make_option(c("--files"),type='character',action='store',help='comma separated *.bad_clust, files REQUIRED'),
	    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file', default=NA),
	    make_option(c("-P","--pdf"),type='character',action='store',help='PDF file name', default=NA),
	    make_option(c("--pub"),action='store_true',help='publication plot', default=FALSE),
	    make_option(c("-S","--stats"),type='character',action='store',help='os_stats_file', default=NA),
	    make_option(c("--mode_pct_thresh"),type='double', action='store',help='mode pct threshold', default=60.0),
	    make_option(c("-M","--mode_good"),action='store_true', default=FALSE)
	    )

opt <- get_yaml_opts(args, option_list, paste0(p.name,'.yaml'))

if (length(opt$files)==1) {
    file_list <- strsplit(opt$files,',')[[1]]
} else {
    file_list <- opt$files
}

if (length(opt$files)==1) {
  file_list <- strsplit(opt$files,',')[[1]]
} else {
  file_list <- opt$files
}

print(opt)

os_stats <- NULL
have_full_name = FALSE
have_short_name = FALSE
if (!is.na(opt$stats) && length(opt$stats)>0) {
   os_stats <- read.table(opt$stats, header=TRUE, sep='\t')
##   print("head(os_stats)")
##   print(head(os_stats))
}

os_names <- NULL
if (!is.na(opt$names) && length(opt$names)>0) {
   os_names <- read.table(opt$names, header=TRUE, sep='\t')
   ## print("os_stats")
   ## print(head(os_stats))
   os_stats <- merge(os_stats,os_names,by='oscode')
   print(os_stats)
}

have_full_name <- ('full_name' %in% colnames(os_stats))
have_short_name <- ('short_name' %in% colnames(os_stats))
clust_data = NULL
leg.labels.oscode = c()
leg.labels.nrow = c()
taxon.names = c()

print(file_list)

mode_good_fn = c()
mode_good_cnt = NULL

report_fields=c('n_omes','p_tot','p_bad','p_short','pct_short','mode_len','p_mode','pct_mode','mode_pct_omecnt')

mmode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x,ux))
  ux[tab == max(tab)]
}

for (file in file_list) {

    print(file)

    taxon = strsplit(file,'_')[[1]][1]

    taxon.names <- append(taxon.names, taxon)

    this_df <- read.table(file, sep='\t', col.names=c('acc','len'),header=FALSE)

    this_df$taxon <- taxon

    this_label = taxon
    if (taxon == 'NEIGO') {
        this_label = 'NEIGO (5X)'
    }

    leg.labels.oscode <- append(leg.labels.oscode, this_label)

    this_mode = mmode(this_df$len)
    print(paste(taxon, "mode:", this_mode))
    this_df$mode_pct <- 100.0 * this_df$len/this_mode

    if (taxon=='NEIGO') {
        long_df <- this_df[rep(seq_len(nrow(this_df)),each=5),]
	print(nrow(long_df))
        this_df <- rbind(this_df, long_df )
    }

    print(summary(this_df))

    clust_data <- rbind(clust_data, this_df)
}

clust_data$taxon <- factor(clust_data$taxon, levels=taxon.names, ordered=TRUE)

## summary(clust_data)
## print(head(clust_data))

theme_set(theme_linedraw(base_size=14))
theme.leg_id <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_text(face='plain', hjust=0,size=12),
  plot.title.position='panel',
  plot.caption=element_text(size=9,hjust=0),
  axis.text.x = element_text(size=11),
  axis.title.x = element_text(size=12),
  axis.text.y = element_text(size=11),
  axis.title.y = element_text(size=12),
  legend.key=element_blank(),
  legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
  legend.justification=c(0,1),
  legend.title=element_blank())

y_scale_pct = scale_y_log10("percent bad proteins",breaks=c(0.0001, 0.0010, 0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.000", "0.001", "0.01", "0.1", "1.0", "10.0", "100.0"),limits=c(0.00010,50.0))

x_scale_pct = scale_x_log10("percent bad proteins",breaks=c(0.0001, 0.0010, 0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.000", "0.001", "0.01", "0.1", "1.0", "10.0", "100.0"),limits=c(0.00010,50.0))

y_scale_cnt = scale_y_log10("number of  bad proteins",breaks=c(0.1, 1, 10, 100,1000,10000,100000),labels=c("0", "1", "10", "100", "1,000", "10,000", "100,000"),position='right')

x_scale = scale_x_continuous('fraction of clusters',breaks=seq(0,1.0,0.2))

if (is.null(opt$colors)) {
    q_color_l = scale_color_brewer(palette = "Set1",labels=leg.labels.oscode)
    q_color = scale_color_brewer(palette = "Set1", labels=leg.labels.oscode)
} else {
    print("have colors")

    print(unlist(opt$colors))

    v_colors = unlist(opt$colors)

    d_colors = v_colors[seq(2,length(v_colors),2)]
    d_colors = setNames(d_colors,v_colors[seq(1,length(v_colors)-1,2)])

    print(d_colors)

    q_color_l = scale_color_manual(values=d_colors,labels=leg.labels.oscode)
    q_color = scale_color_manual(values=d_colors,labels=leg.labels.oscode)
}

color_guide= guide_legend(position='inside')
color_pos = c(0.60,0.975)

print(leg.labels.oscode)

p_mode_pct <- ggplot(clust_data) + 
	   geom_freqpoly(aes(x=mode_pct, color=taxon),fill='transparent',binwidth=2) + 
##	   geom_histogram(aes(x=mode_pct, color=taxon),fill='transparent',position='identity',binwidth=2) +
	   q_color_l + theme.leg_id + xlab("percent of mode length") + ylab("number of proteins") +
	   scale_x_continuous(limits=c(40,160),breaks=seq(40,160,20)) +
	   scale_y_continuous() +
	   geom_vline(xintercept=75,linetype='dashed') +
	   geom_vline(xintercept=133,linetype='longdash') +
           guides(color=color_guide) + theme(legend.position.inside=color_pos,legend.text=element_markdown(size=8),legend.key.spacing.y=unit(-4,'pt'))

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (!is.na(opt$pdf)) {
   pdf_file = opt$pdf
} else {
   pdf_file = clust_data[1,]$taxon
}

if (! opt$pub) {
   p_pct_mode <- ( p_mode_pct) / doc_panel + 
	          plot_layout(heights=c(10, 0.1))

   ggsave(file=pdf_file, plot=p_pct_mode, width=4.0, height=3.8)
} else {

  ## automatically add _pub to file name if necessary

   if (!grepl('_pub',pdf_file)) {
      pdf_file_pref <- strsplit(pdf_file,'\\.')[[1]][1]
      pdf_file <- paste0(pdf_file_pref,"_pub.pdf")
   }

   p_pct_mode <- ( p_mode_pct) 
   ggsave(file=pdf_file, plot=p_mode_pct, width=4.0, height=3.5)
}


warnings()
