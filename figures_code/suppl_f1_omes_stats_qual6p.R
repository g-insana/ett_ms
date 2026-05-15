#!/usr/bin/env Rscript --vanilla

## modified 13-Mar-2026 to produce single 6-panel plot for supplemental figure

################
## omes_stats_mode_qual_all.R ecoli...bad_clust
## 
## derivative of clust_mode_qual_all5.R that does statistics on omes, rather than clusters
##
################
## reads in *_proteome.stats files
##
## 
## modified to read --stats=oscode_ome_clust_cnt.tsv, no names --taxon_order
##

library('ggplot2')
library('stringr')
library('scales')
library('dplyr')
library('tidyr')
library('purrr')
library('RColorBrewer')
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

## *proteome.stats
## GCA_id	ome_id	clust_cnt	mode_cnt	outlier_cnt	short_cnt	in_sample
## GCA_001258975.1	1346192	1690	1127	205	197	BS
## GCA_001254255.1	1342561	1971	1388	198	193	BS

p.name<-'omes_stats_mode_qual_all6'
p.full_name <- paste0(p.name,'.R')

args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.full_name,"\n",args),collapse=' ', sep=' ')

opt_list = list(
    make_option(c("--files"),type='character',action='store',help='comma separated *.bad_clust, files REQUIRED',default=NA),
    make_option(c("-D","--debug"),action='store_true',help='debug', default=FALSE),
    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file', default=NA),
    make_option(c("-P","--pdf"),type='character',action='store',help='PDF file name', default=NA),
    make_option(c("-U","--sup_pdf"),type='character',action='store',help='PDF file name', default=NA),
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

have_fa_stats <- FALSE
os_stats <- NULL
if (!is.na(opt$stats) && length(opt$stats)>0) {
   os_stats <- read.table(opt$stats, header=TRUE, sep='\t')
   if (opt$debug) {
      print("os_stats")
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

##   print(head(os_stats))
   os_stats <- os_stats[order(-os_stats$count),]
   if (opt$debug) {
      print(os_stats)
   }
}

if (!have_fa_stats && !is.na(opt$fa_stats) && length(opt$fa_stats)>0) {
   fa_stats <- read.table(opt$fa_stats, header=TRUE, sep='\t')

   if (opt$debug) {
      print("fa_stats")
##   print(head(fa_stats))
      print(fa_stats)
   }

   if (nrow(os_stats) > 0 & nrow(fa_stats)) {
      os_stats <- merge(os_stats, fa_stats, by='oscode')
      os_stats$ome_fn <- os_stats$ome_cnt/os_stats$fa_med
      if (opt$debug) {
         print("os_stats")
         print(os_stats[,c('oscode','count','ome_cnt','ome_cnt','fa_med','ome_fn')])
         print("summary(os_stats$ome_fn)")
         print(summary(os_stats$ome_fn))
      }
   }
}

leg.labels.oscode = c()
leg.labels.nrow = c()
taxon.names = c()

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

## table for all taxa
ome_data <- NULL
ome_thresh = NULL

if (opt$debug) {
  print("file_list")
  print(file_list)
}

for (ome_file in file_list) {

    if (opt$debug) {
      print("ome_file:")
      print(ome_file)
    }
    
    if (file.exists(ome_file)) {

    this_taxon = strsplit(ome_file,'_')[[1]][1]
    taxon.names <- append(taxon.names, this_taxon)

    this_os_stat = os_stats[os_stats$oscode==this_taxon,]
##    print(this_os_stat[,c('oscode','count','ome_cnt','short_name','ome_cnt','fa_med')])

    this_data <- read.table(ome_file,header=TRUE,sep='\t')

##    print("head(this_data)")
##    print(head(this_data))

    this_data$taxon = this_taxon 

    new_nrow = nrow(this_data)

    this_ome_cnt <- this_os_stat$count

    if (opt$debug) {
      print(paste("this_ome_cnt",this_ome_cnt))
      print(paste("pre-nrow:", nrow(this_data)))
    }

    this_data <- this_data[this_data$clust_cnt > this_os_stat$clust_cnt/5,]

    if (opt$debug) {
      print(paste("post-nrow:", nrow(this_data)))
    }

    this_max_clust <- max(this_data$clust_cnt)

    if (opt$debug) {
      print(paste("taxon: ",this_taxon,"max_clust:",this_max_clust,"os_clust_cnt:",this_os_stat$clust_cnt))
    }
    this_data$pct_clust <- 100 * this_data$clust_cnt / this_os_stat$clust_cnt

    this_data$pct_mode <- 100.0 * this_data$mode_cnt / this_data$clust_cnt
    this_data$pct_bad <- 100.0 * this_data$outlier_cnt / this_data$clust_cnt
    this_data$pct_short <- 100.0 * this_data$short_cnt / this_data$outlier_cnt

    ## print(paste(this_taxon,"cluster summary"))
    ## print(summary(this_data))

##    print(paste("Nrows:",nrow(ome_data),nrow(this_data)))

    ome_data <- rbind(ome_data, this_data)
    }
}

if (opt$debug) {
  print("head(ome_data)")
  print(head(ome_data))
}

ome_data$pct_bad_nz <- ifelse(ome_data$pct_bad > 0, ome_data$pct_bad, 0.001)

short_stats <- ome_data %>% group_by(taxon) %>%
     summarize(
	count_Z_fct = sum(outlier_cnt < 1)/n(),
 	short_fct = median(ifelse(outlier_cnt > 0, short_cnt/outlier_cnt, 0.0))
     ) %>%
     arrange (desc(short_fct)) %>%
     ungroup()

if (opt$debug) {
  ## print(short_stats)
  print("summary(fraction of outliers short)")
  summary(short_stats$short_fct)
  #head(short_stats)
}

## taxa.ordered = taxa_medians[order(c(taxa_medians$pct_bad)),]$taxon
if (!is.na(opt$taxa_order) && length(opt$taxa_order)>0) {
   taxa_order <- read.table(opt$taxa_order, header=TRUE, sep='\t')
   taxa.ordered = taxa_order$oscode
   if (opt$debug) {
      print("taxon_order")
      print(taxa.ordered)
   }
} else {
## order taxa by decreasing proteome count
  taxa.ordered = os_stats[order(-os_stats$count),]$oscode
}
## taxa.ordered

short_stats$taxon <- factor(short_stats$taxon, levels=taxa.ordered, ordered=TRUE)

ome_data$taxon <- factor(ome_data$taxon, levels=taxa.ordered, ordered=TRUE)
ome_thresh$taxon <- factor(ome_thresh$taxon, levels=taxa.ordered, ordered=TRUE)

ome_thresh$thresh <- factor(ome_thresh$thresh, levels=c('pct0','pct0.1','pct1.0','pct10.0','pct100'), ordered=TRUE)

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
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
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

y_scale<-function(x) {
  x <- ifelse(x >= top100_delta,top100_delta, x)

  xs <- ifelse (x <= 50.0, x, 
                ifelse(x>top100_delta,100.0,50.0*(log10(100.0-x)-log10(50))/log50_delta + 50.0))
  xs   
}

y_scale_inv<-function(x) {


  xs <- ifelse(x<=50.0, x, 
               ifelse(x > top100_delta, 100.0,
	              100 - 50*(10**(x-50.0) * (-log50_delta)/50.0)))
  xs
}

sq100_trans <- function() {

  trans_new("sq100", function(x) y_scale(x), function(x) y_scale_inv(x))
}

xtop_scale = scale_x_discrete(position='top')

ylog_scale_bad = scale_y_log10('percent',limits=c(0.001, 100.0), breaks=10^seq(-3,2.0,1.0),labels=c('0.00','0.01','0.1','1.0','10.0','100'))

ylog_scale_bad_r = scale_y_log10('percent',limits=c(0.001, 100.0), breaks=10^seq(-3,2.0,1.0),labels=c('0.00','0.01','0.1','1.0','10.0','100'),position='right')

ypct_scale_bad = scale_y_continuous(limits=c(0.0, 100.0), breaks=seq(0,100,20))

yprot_scale_bad60 = scale_y_continuous(limits=c(60, 100.0))

yprot_scale_bad60_log = scale_y_continuous(limits=c(50, 100.0), breaks=c(60,75,90,95,98,99,100),labels=c('60','75','90','95','98','99','100'))

yprot_scale_bad60_log_r = scale_y_continuous(limits=c(50, 100.0), breaks=c(60,75,90,95,98,99,100),labels=c('60','75','90','95','98','99','100'), position='right')

ypct_scale_short = scale_y_continuous(limits=c(40.0, 100.0), breaks=seq(40.,100.,10))
ypct_scale_short_r = scale_y_continuous(limits=c(40.0, 100.0), breaks=seq(40.,100.,10), position='right')

## str(ome_data)

sp_plot_margin = c(0.1, 0.1, 0.25, 0.1)*1.00

taxon_medians <- ome_data %>% group_by(taxon) %>%
	           summarize(
	             clust_median = median(pct_clust),
	             clust_q1 = fivenum(pct_clust)[2],
	             mode_median = median(pct_mode),
	             mode_q1 = fivenum(pct_mode)[2],
	             bad_median = median(pct_bad),
	             short_median = median(pct_short,na.rm=TRUE)
		   )

print("taxon_medians")
print(taxon_medians)

taxon_med_meds <- taxon_medians %>% summarize(across(where(is.numeric), \(x) median(x,na.rm=TRUE)))
print(taxon_med_meds)


p_ome_clust_box <- ggplot(ome_data, aes(x=taxon, y=pct_clust, color=taxon)) + geom_boxplot(outlier.shape=NA) + q_color + yprot_scale_bad60 + 
	     geom_hline(yintercept=taxon_med_meds$clust_median, linetype="longdash") +
	     theme.leg_id + 
	     ggtitle("A. clusters per proteome ") + ## xtop_scale +
	     labs(y='percent') + 
	     theme(axis.title.x=element_blank(),
	           axis.text.x=element_text(angle=45, hjust=1.0, vjust=1.0, color=ds_colors),
	           legend.position='None') +
     	     theme(plot.margin=unit(sp_plot_margin,'cm'))

ome_medians <- ome_data %>% group_by(taxon) %>%
               summarize(ome_median = median(pct_mode))

## ome_med_median <- median(ome_medians$ome_median)
## print(paste("ome_med_median:",ome_med_median))

p_ome_mode_box <- ggplot(ome_data, aes(x=taxon, y=pct_mode, color=taxon)) + geom_boxplot(outlier.shape=NA) + q_color + yprot_scale_bad60 + 
	     theme.leg_id + 
	     geom_hline(yintercept=taxon_med_meds$mode_median, linetype="longdash") +
	     ggtitle("B. mode proteins per proteome ") + ## xtop_scale +
	     labs(y='percent') + 
	     theme(axis.title.x=element_blank(),
	         axis.text.x=element_blank(),
		 legend.position='None') +
     	     theme(plot.margin=unit(sp_plot_margin,'cm'))

p_bad_box <- ggplot(ome_data, aes(x=taxon, y=pct_bad_nz, color=taxon)) + geom_boxplot(outlier.shape=5,outlier.alpha=0.2,outlier.size=2.0) + q_color + ylog_scale_bad + 
	     geom_hline(yintercept=taxon_med_meds$bad_median, linetype="longdash") +
	     theme.leg_id + 
	     ggtitle("C. outliers per proteome") + labs(y='percent') +
	     theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position='None') +
	     theme(plot.margin=unit(sp_plot_margin,'cm'))


med_not_bad <- 100.0*median(short_stats$count_Z_fct)
print(paste("med_not_bad:",med_not_bad))

p_not_bad <- ggplot(short_stats, aes(x=taxon, y=100.0*count_Z_fct, color=taxon)) + geom_point(aes(shape=taxon,size=taxon)) + ## geom_point(shape=3,size=3) + 
	     geom_hline(yintercept=med_not_bad, linetype="longdash") +
	     theme.leg_id + 
	  q_color + ypct_scale_bad + lq_shape + lq_size +
          ggtitle("D. no outlier proteomes") +
	  labs(y='percent') + 
	  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position='None') +
          theme(plot.margin=unit(sp_plot_margin,'cm'))

## median(ome_data$pct_short,na.rm=TRUE)
## med_short_pct <- median(ome_data$pct_short,na.rm=TRUE)
## print(head(ome_data$pct_short))
## print(paste("med_short_pct:",med_short_pct))

## med_short_pct <- 100.0*median(short_stats$short_fct)
## print(paste("med_short_pct2:",med_short_pct))

p_bad_short <- ggplot(short_stats, aes(x=taxon, y=100*short_fct, color=taxon)) + 
	    ## geom_boxplot(outlier.shape=5,outlier.alpha=0.2,outlier.size=2) + 
	     geom_point(aes(shape=taxon,size=taxon)) +
	     geom_hline(yintercept=taxon_med_meds$short_median, linetype="longdash") +
	     theme.leg_id + 
 	     q_color + ypct_scale_short + lq_shape + lq_size +
             ggtitle("E. short outliers per proteome") +
	     labs(y='percent') + 
	     theme(axis.title.x=element_blank(),
	           axis.text.x=element_text(angle=45, hjust=1.0, vjust=1.0, color=ds_colors),
                   legend.position='None') +
	     theme(plot.margin=unit(sp_plot_margin,'cm'))

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (!is.na(opt$pdf)) {
   pdf_file = opt$pdf
} else {
   pdf_file = paste0('bfig_',taxa.ordered[1],"_len_dist.pdf")
}

  p_ome_clust_box_m = p_ome_clust_box + theme(plot.margin=unit(c(0,10,0,0),'pt'))
  p_not_bad_m = p_not_bad + theme(plot.margin=unit(c(0,10,0,0),'pt'))
  p_ome_mode_box_m = p_ome_mode_box + theme(plot.margin=unit(c(0,10,0,0),'pt'))

  n_panel <- 5
  widths_cnt1 <- c(10,1,10)
  p_pct_cnt1 <- (p_ome_clust_box / p_ome_mode_box  / p_bad_box / p_not_bad / p_bad_short )
  rel_heights <- c(5,4,5,4,4)

## p_pct_cnt2 <- ( p_ome_mode_box / p_ome_mode_box_ls )


if (! opt$pub) {

   p_pct_cnt1 <- (p_pct_cnt1 / doc_panel) + plot_layout(heights=c(rel_heights,0.1))

   ggsave(file=pdf_file, plot=p_pct_cnt1, width=5.0, height= 1.2 +(n_panel * 1.8))

} else {

   p_pct_cnt1 <- p_pct_cnt1 + plot_layout(heights=rel_heights, widths=widths_cnt1)

   if (!grepl('_pub',pdf_file)) {
      pdf_file_pref <- strsplit(pdf_file,'\\.')[[1]][1]
      pdf_file <- paste0(pdf_file_pref,"_pub.pdf")
   }

   ggsave(file=pdf_file, plot=p_pct_cnt1, width=5, height= 1.0 +(n_panel * 1.8))
}

## warnings()
