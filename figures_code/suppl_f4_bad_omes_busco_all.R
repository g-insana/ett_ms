#!/usr/bin/env Rscript --vanilla

################
## bad_omes_dist4.R --omes ecoli....bad_omes.samp --tfx_b=ecoli_bad_omes200.tfxg_stats_BHL --tfx_s=ecoli_bad_omes200.tfxg_stats_S1K
## 
## this script is only used for looking at reference vs other/redundant/excluded proteomes

################
## reads in *.bad_clust.samp.am
##
## 13-May-2025 modified to read field names from .samp.am file
##

library('ggplot2')
library('ggtext')
library('stringr')
library('forcats')  ## for fct_relevel
library('RColorBrewer')
library('patchwork')
library('getopt')
library('optparse')
library('yaml')
library('stringr')

get_yaml_file = 'get_yaml_opts.R'
read_omes_file= 'read_omes.R'

if (file.exists(get_yaml_file)) {
   source(get_yaml_file)
   source(read_omes_file)
} else if (file.exists(paste0("../",get_yaml_file))) {
   source(paste0("../",get_yaml_file))
   source(paste0("../",read_omes_file))
} else {
  cat(paste("cannot open", get_yaml_file))
  sys.exit(1)
}

## .bad_omes
## proteomes with most short/long proteins (119905 total) in 3871 clusters
## cluster min: 1 median: 5.0 max: 675
## cluster quantiles: ## [2.0, 3.0, 4.0, 4.0, 5.0, 6.0, 8.0, 9.0, 12.0]
## proteome_id	n_clusters	%_cluster	B_prot_cnt	B_compl_comb	B_compl_single	B_frag	B_miss
## GCA_013423585.1	675	17.4	11884	60.2	58.9	30.9	8.9
## GCA_030296335.1	576	14.9	7517	51.1	50.9	29.8	19.1
## GCA_030316365.1	557	14.4	7583	50.2	50.0	30.9	18.9
## GCA_025660355.1	478	12.3	9185	28.0	26.8	26.8	45.2
## GCA_010725135.1	474	12.2	9094	29.8	29.5	41.1	29.1

## also could have genome assembly statistics

p.name<-'bad_omes_dist_busco2m'
p.full_name <- paste0(p.name,".R")

cmd_args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.full_name,"\n",cmd_args),collapse=' ', sep=' ')

option_list = list(
	    make_option(c("--omes"),type='character',action='store',help='comma separated *.summ, files REQUIRED'),
	    make_option(c("--tfx_bad"),type='character', action='store',help='tfx .stats_BHL file', default=NA),
	    make_option(c("--tfx_samp"),type='character', action='store',help='tfx .stats_S1K file', default=NA),
	    make_option(c("-P","--pdf"),type='character',action='store',help='PDF file name', default=NA),
	    make_option(c("-T","--title"),type='character',action='store',help='plot title'),
	    make_option(c("--stats"),type='character',action='store',help='proteome stats',default=NA),
	    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file'),
	    make_option(c("--pub"),action='store_true',help='yaml file',default=FALSE)
	    )

opt <- get_yaml_opts(cmd_args, option_list, paste0(p.name,".yaml"))

if (length(opt$files)==1) {
    file_list <- strsplit(opt$omes,',')[[1]]
} else {
    file_list <- opt$omes
}

if (length(opt$cc_cnt) == 0) {    
   opt$cc_cnt = 10
}

B_compl_thresh = 60.0


## print(opt)

os_stats <- NULL
have_full_name=FALSE
have_short_name=FALSE
if (length(opt$stats)>0) {
   print(paste("read.table", opt$stats))
   os_stats <- read.table(opt$stats, header=TRUE, sep='\t')

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

ome_fields <- c('proteome_id','n_clusters','pct_cluster','pct_cluster_nz','B_prot_cnt','B_assem_lvl','B_compl_comb','B_compl_single','B_frag','B_miss')

tfx_fields <- c('file','logN50','x1cnt','Ncnt','tot_cnt','pct_good','pct_frame')
tfx_fields3 <- c('file','logN50','x1cnt','x1len', 'Ncnt','Nlen','x0cnt','x0len', 'tot_good','tot_cnt','pct_good','pct_frame')

ome_fields_save = c(ome_fields,'B_prot_status', 'A_n50', 'A_contig_cnt', 'A_annot_date', 'A_coverage', 's_type','A_am','B_prot_status_f')

taxon.names = c()
taxon_labels = c()

## b_val_parms = c('B_frag', 'B_compl_comb','B_compl_single','B_miss')
b_val_parms = c('B_frag', 'B_compl_comb','B_miss')

## b_val_strings = c('fragment score','completeness (combined)','completeness (single)','missing score')
b_val_strings = c('fragment score','completeness (combined)','missing score')

n_omes <- length(file_list)
n_parm <- length(b_val_parms)
n_plots <- n_omes * n_parm

if (length(opt$pub) == 0 || ! opt$pub) {
   opt$pub = FALSE
   plot_list = vector("list",n_plots+1)
} else {
   plot_list = vector("list",n_plots)
}

## universal plot initiation

theme_set(theme_linedraw(base_size=14))
theme.leg_id <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_markdown(size=14, hjust=0),
  plot.caption=element_text(size=9,hjust=0),
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  legend.position.inside=c(0.95,0.05), legend.justification.inside=c(1,0),
  legend.text=element_text(size=9, hjust=0),
  legend.key=element_blank(),
  legend.background=element_rect(fill="transparent", color='black',linetype='solid',linewidth=0.4),
  legend.box.background=element_rect(fill="transparent"),
  legend.justification=c(0,1),
  legend.title=element_blank())

theme.leg_id_NL <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_text(face='plain', size=14, hjust=0),
  plot.caption=element_text(size=9,hjust=0),
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  legend.position='none',
  legend.title=element_blank())


ly_scale_pct = scale_y_log10("bad clusters per proteome (%)",breaks=c(0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.00", "0.1", "1.0", "10.0", "100.0"),limits=c(0.01,30))

ly_scale_pct_r = scale_y_log10("bad clusters per proteome (%)",breaks=c(0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.00", "0.1", "1.0", "10.0", "100.0"),limits=c(0.01,30),position='right')

y_scale_pct = scale_y_continuous("bad clusters per proteome (%)",limits=c(0.0,30))
y_scale_pct_r = scale_y_continuous("bad clusters per proteome (%)",limits=c(0.0,30),position='right')

y_scale_cnt = scale_y_log10("bad clusters per proteome (N)",breaks=c(0.1, 1, 10, 100,1000),labels=c("0", "1", "10", "100", "1,000"))

x_scale_bus = scale_x_continuous("BUSCO fragment score",limits=c(-0.2,35),breaks=seq(0,30,10))

lx_scale_bus = scale_x_log10(breaks=c(0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.00", "0.1", "1.0", "10.0", "100.0"),limits=c(0.01,100))

x_scale_compl = scale_x_continuous(limits=c(B_compl_thresh,100))

x_scale = scale_x_continuous('sampled proteomes')
x_scale_pct = scale_x_continuous('fraction of sampled proteomes', limits=c(0,1),breaks=seq(0,1.0,0.2))

colors2 <- c('S'= '#00BFC4', 'B'= '#F8766D' )
colors4 <- brewer.pal(4, 'Set1')
names(colors4) = c('S0','S1','S2','S3')

stat_labels= c('S0'='reference','S1'='other','S2'='redundant','S3'='excluded')

s_color_sb <- scale_color_manual(values=colors2, labels=c('S'='total','B'='worst'))

drop_flag=TRUE
s_color_stat <- scale_color_manual(values=colors4,labels=stat_labels,drop=drop_flag)
s_alpha_stat <- scale_alpha_manual(values=rev(c(1.5,0.8,1.5,0.8)),labels=stat_labels,drop=drop_flag)
s_size_stat  <-  scale_size_manual(values=rev(c(2.0,1.0,2.0,1.0)),labels=stat_labels,drop=drop_flag)
s_shape_stat  <-  scale_shape_manual(values=rev(c(16,1,16,1)),labels=stat_labels,drop=drop_flag)

s_shape_sb <- scale_shape_manual(values=c('S'=1,'B'=0), labels=c('S'='all','B'='worst'))

title_idx = 1

taxon_titles = rep("",n_omes)

for (ix in 1:n_omes) {

    ome_file <- file_list[ix]

    print(paste(c(ix,ome_file,file_list[ix])))

    taxon = strsplit(ome_file,'_')[[1]][1]
    if (length(opt$title)!=0) {
       taxon_name = opt$title
    } else {
       taxon_name = taxon
    }      

    taxon.names <- append(taxon.names, taxon)

    this_os_info <- os_stats[os_stats$oscode==taxon,]

    ome_data <- read_ome_file(ome_file, taxon, ome_fields_save, dedup=FALSE)
    
    ome_data <- sample_dedup(ome_data, this_os_info$count, 0.2, explain=TRUE)

    taxon_label = taxon

    taxon_cnt = ''
    if (taxon %in% os_stats$oscode) {
       this_os_info <- os_stats[os_stats$oscode==taxon,]
       taxon_cnt = this_os_info$count
       taxon_clust_cnt = this_os_info$clust_cnt
       if (have_short_name) {
           taxon_label = sprintf("%d proteomes\n%d clusters",taxon_cnt,taxon_clust_cnt)
	   taxon_title = sprintf("*%s* (%s)",this_os_info$short_name,taxon)
       } else {
           taxon_label = sprintf("%d proteomes\n%d clusters",taxon_cnt,taxon_clust_cnt)
	   taxon_title = sprintf("%s",taxon)
       }
       taxon_cnt_str = sprintf("(%d)",taxon_cnt)
       taxon_labels <- append(taxon_labels, taxon_label)
       taxon_titles[[ix]] <- taxon_title
    }

    ome_data$pct_cluster_nz = ifelse(ome_data$pct_cluster < 0.01, 0.01, ome_data$pct_cluster)
    ome_data$s_type = factor(ome_data$s_type,levels=c('B','S'),ordered=TRUE)

    low_prot_cnt = nrow(ome_data[ome_data$B_prot_cnt < 1,])
    if (low_prot_cnt > 1) {
        print(paste("taxon:", taxon_name, " rows with B_prot_cnt < 1", low_prot_cnt))
        print(ome_data[ome_data$B_prot_cnt < 1,c('pct_cluster_nz','B_frag','B_prot_status','B_prot_status_f')])
    }

    print(paste("taxon_name:",taxon_name))
    d_cols = c(c('s_type','pct_cluster','pct_cluster_nz'),b_val_parms)
    print(d_cols)
    print(head(ome_data[,d_cols]))
    print(summary(ome_data[,d_cols]))

    ## all done with ome_data, now get tfx_data

    N_omes = nrow(ome_data)
    NS_omes = nrow(ome_data[ome_data$s_type=='S',])
    NB_omes = nrow(ome_data[ome_data$s_type=='B',])

    ## print("all, NS, NB rows")
    ## print(paste(N_omes, NS_omes, NB_omes))

    ome_data_samp <- ome_data[ome_data$s_type=='S',]
    ome_data_samp <- ome_data_samp[order(ome_data_samp$pct_cluster_nz),]
    ome_data_samp$xrank = 1:nrow(ome_data_samp)/nrow(ome_data_samp)

    ome_data_bad <- ome_data[ome_data$s_type=='B',]
    ome_data_bad <- ome_data_bad[order(ome_data_bad$pct_cluster_nz),]
    ome_data_bad$xrank = 1:nrow(ome_data_bad)/nrow(ome_data_bad)

    ## we would like the "bad" omes to sample no more than the worst 5% of the sample
    ## so if os_stats$count == 120,000 (E. coli), then 2,000 is fine
    ## but if os_stats$count == 20,000, then we only want 1000

    bad_fract = as.integer(min(taxon_cnt/20, 2000))
    bad_fract_pos = max(nrow(ome_data_bad)-bad_fract + 1,1)

    ## print(paste("bad_fract_pos", taxon_name, bad_fract_pos))
    ## print(ome_data_bad[bad_fract_pos,c('proteome_id','pct_cluster_nz','xrank','s_type')])


    ome_data_s <- ome_data[ome_data$B_prot_status_f=='S0',]
    print(paste("nrows S0:",nrow(ome_data_s)))
    ome_data_s <- rbind(ome_data_s,ome_data[sample(nrow(ome_data),nrow(ome_data)/2),])
    ome_data_s <- ome_data_s[!is.na(ome_data_s$B_prot_status_f),]

    ## plots against other busco scores:
    ## B_compl_comb, B_compl_single, B_miss

    legend_rev_flag=TRUE
    ## legend_rev_flag=FALSE

    ome_data_s$B_prot_status_f <- fct_relevel(ome_data_s$B_prot_status_f, 'S3','S2','S1','S0')

    g_busco1 <- guide_legend(position='inside',reverse=legend_rev_flag)

    b_idx = 1
    for (b_val in b_val_parms ) {

        ## calculate correlations
        clust_pct_corr <- cor(ome_data$pct_cluster_nz,ome_data[[b_val]])
	if (is.na(clust_pct_corr)) {
	    print("*** using sampled data ***")
	    ome_data_samp_t <- ome_data[sample(nrow(ome_data),nrow(ome_data)/10),]
	    clust_pct_corr <- cor(ome_data_samp_t$pct_cluster_nz,ome_data_samp_t[[b_val]])
	}

	if (is.na(clust_pct_corr)) {
	   print("*** using worst data ***")
	   ome_data_bad <- ome_data[ome_data$s_type == 'B',]
	   clust_pct_corr <- cor(ome_data_bad$pct_cluster_nz,ome_data_bad[[b_val]])
	}

	print(paste("clust_pct_corr",b_val))
	print(clust_pct_corr)

	b_val_nz = paste0(b_val,"_nz")

	ome_data_s[,b_val_nz] = ifelse(ome_data_s[,b_val] < 0.01, 0.01, ome_data_s[,b_val])

	this_scale_x <- lx_scale_bus
	this_scale_y <- ly_scale_pct
	ann_pos_x = 0.015
	if (grepl('B_compl',b_val)) {
	    print(paste("bad B_compl",b_val))
	    if (nrow(ome_data_s[ome_data_s[,b_val_nz] < B_compl_thresh,]) > 1) {
	        print(paste("re-scaled",b_val,nrow(ome_data_s[ome_data_s[,b_val_nz] < B_compl_thresh,])))
	    }
	    ome_data_s[,b_val_nz] <- ifelse( ome_data_s[,b_val_nz] < B_compl_thresh, B_compl_thresh, ome_data_s[,b_val_nz])
	    this_scale_x <- x_scale_compl
	    ## this_scale_y <- y_scale_pct
	    ann_pos_x = 62.0
        }
	this_guide <- guides(color=g_busco1, shape=g_busco1, size=g_busco1, alpha=g_busco1)
	if (b_idx == 1) {
	    this_theme <- theme.leg_id + theme(legend.position.inside=c(0.95,0.05), legend.justification.inside=c(1,0))
	    this_title = ggtitle(sprintf("%s. %s",LETTERS[title_idx],taxon_titles[ix]))
	} else {
	    this_theme <- theme.leg_id_NL
	    this_title = ggtitle(sprintf("%s.",LETTERS[title_idx]))
	}

        p_this_bval <- ggplot(ome_data_s,aes(y=pct_cluster_nz,x=.data[[b_val_nz]], color=B_prot_status_f,alpha=B_prot_status_f,size=B_prot_status_f,shape=B_prot_status_f)) +
	this_title +
        geom_point(data=filter(ome_data_s,B_prot_status_f != 'S0'), position=position_dodge(0.01)) +
        geom_point(data=filter(ome_data_s,B_prot_status_f == 'S0'), position=position_dodge(0.01)) +
        this_scale_y + this_scale_x + s_color_stat + s_alpha_stat + s_size_stat + s_shape_stat +
	xlab(sprintf("BUSCO %s", b_val_strings[b_idx])) +
        annotate('text',label=paste("r^2 ==",sprintf("%.3f",clust_pct_corr)),y=0.015, x=ann_pos_x,parse=TRUE,hjust=0) +
        this_theme + this_guide

	this_idx = (ix - 1)*n_parm + b_idx
	plot_list[[this_idx]] = p_this_bval

	b_idx = b_idx + 1
	title_idx = title_idx + 1
    }
}

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

big_plot = Reduce('/', plot_list) + plot_layout(ncol=n_parm)

if (is.na(opt$pdf)) {
   ## print(opt$omes)
   print(strsplit(opt$omes,'_'))
   fig3_file=paste0("fig3_omes_busco_",strsplit(opt$omes,"\\.")[[1]][1],".pdf")
} else {
  fig3_file = opt$pdf
}

if (! opt$pub) {
   big_plot <- big_plot / doc_panel + plot_layout(heights=c(10.,0.1))
   ggsave(file=fig3_file, plot=big_plot, width=4.0*n_parm, height=(4.0*n_omes + 0.2))
} else {

   if (!grepl('_pub',fig3_file)) {
      fig3_file_pref <- strsplit(fig3_file,'\\.')[[1]][1]
      fig3_file <- paste0(fig3_file_pref,"_pub.pdf")
   }

   ggsave(file=fig3_file, plot=big_plot, width=4.0*n_parm, height=(4.0*n_omes))
}


warnings()
