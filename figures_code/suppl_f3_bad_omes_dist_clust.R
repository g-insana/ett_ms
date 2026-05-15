#!/usr/bin/env Rscript --vanilla

################
## bad_omes_dist_clust.R --omes ecoli....bad_omes.samp --tfx_b=ecoli_bad_omes200.tfxg_stats_BHL --tfx_s=ecoli_bad_omes200.tfxg_stats_S1K
## 
## this script is only used for looking at reference vs other/redundant/excluded proteomes

################
## reads in *.bad_clust.samp.am
##
## 13-May-2025 modified to read field names from .samp.am file
##

library('ggplot2')
library('ggtext')
library('dplyr')
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

p.name<-'bad_omes_dist_clust'
p.full_name <- paste0(p.name,".R")

cmd_args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.full_name,"\n",cmd_args),collapse=' ', sep=' ')

option_list = list(
	    make_option(c("--omes"),type='character',action='store',help='comma separated *.summ, files REQUIRED'),
	    make_option(c("-C","--bad_clust"),type='character',action='store',help='*.bad_clust, files REQUIRED'),
	    make_option(c("-P","--pdf"),type='character',action='store',help='PDF file name', default=NA),
	    make_option(c("-T","--title"),type='character',action='store',help='plot title'),
	    make_option(c("--stats"),type='character',action='store',help='proteome stats',default=NA),
	    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file'),
	    make_option(c("--cl_lens"),type='character',action='store',help='clust_proteome_len.samp files'),
	    make_option(c("--pub"),action='store_true',help='yaml file',default=FALSE)
	    )

opt <- get_yaml_opts(cmd_args, option_list, paste0(p.name,".yaml"))

if (length(opt$files)==1) {
    file_list <- strsplit(opt$omes,',')[[1]]
} else {
    file_list <- opt$omes
}

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

ome_fields <- c('proteome_id','n_clusters','ns_clusters','pct_cluster','pct_short','pct_cluster_nz','B_prot_cnt','B_assem_lvl','B_compl_comb','B_compl_single','B_frag','B_miss')

tfx_fields <- c('file','logN50','x1cnt','Ncnt','tot_cnt','pct_good','pct_frame')
tfx_fields3 <- c('file','logN50','x1cnt','x1len', 'Ncnt','Nlen','x0cnt','x0len', 'tot_good','tot_cnt','pct_good','pct_frame')

ome_fields_save = c(ome_fields,'B_prot_status', 'A_n50', 'A_contig_cnt', 'A_annot_date', 'A_coverage', 's_type','A_am','B_prot_status_f')

taxon.names = c()

n_omes <- length(file_list)

if (length(opt$pub) == 0 || ! opt$pub) {
   opt$pub = FALSE
   plot_list = vector("list",n_omes+1)
} else {
   plot_list = vector("list",n_omes)
}

## universal plot initiation

theme_set(theme_linedraw(base_size=14))
theme.leg_id <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_markdown(face='plain', size=14, hjust=0),
  plot.caption=element_text(size=9,hjust=0),
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  legend.position='None',
  legend.text=element_blank(),
  legend.key=element_blank(),
  legend.background=element_rect(fill="transparent", color='black',linetype='solid',linewidth=0.4),
  legend.box.background=element_rect(fill="transparent"),
  legend.justification=c(0,1),
  legend.title=element_blank())

## if (opt$pub) {
##   theme.leg_id_NL <- theme.leg_id + theme(legend.position='None', strip.text=element_blank())
## } else {
   theme.leg_id_NL <- theme.leg_id + theme(legend.position='None')
## }

ly_scale_pct = scale_y_log10("bad clusters per proteome (%)",breaks=c(0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.00", "0.1", "1.0", "10.0", "100.0"),limits=c(0.01,30))

ly_scale_pct_r = scale_y_log10("bad clusters per proteome (%)",breaks=c(0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.00", "0.1", "1.0", "10.0", "100.0"),limits=c(0.01,30),position='right')

y_scale_pct_r = scale_y_continuous("bad clusters per proteome (%)",limits=c(0.0,30),position='right')

y_scale_cnt = scale_y_log10("bad clusters per proteome (N)",breaks=c(0.1, 1, 10, 100,1000),labels=c("0", "1", "10", "100", "1,000"))

x_scale_bus = scale_x_continuous("BUSCO fragment score",limits=c(-0.2,35),breaks=seq(0,30,10))
lx_scale_bus = scale_x_log10("BUSCO fragment score", breaks=c(0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.00", "0.1", "1.0", "10.0", "100.0"),limits=c(0.01,30))


x_scale = scale_x_continuous('sampled proteomes')
x_scale_pct = scale_x_continuous('fraction of sampled proteomes', limits=c(0,1),breaks=seq(0,1.0,0.2))

colors2 <- c('S'= '#00BFC4', 'B'= '#F8766D' )
colors4 <- brewer.pal(4, 'Set1')
names(colors4) = c('S0','S1','S2','S3')

colors3 <- c('mode'='lightgray','not-mode'='#00BFC4','outlier'='#F8766D')
alpha3 <- c('mode'=0.75,'not-mode'=0.9,'outlier'=1.0)

stat_labels= c('S0'='reference','S1'='other','S2'='redundant','S3'='excluded')

s_color_sb <- scale_color_manual(values=colors2, labels=c('S'='total','B'='worst'))

s_color_stat <- scale_color_manual(values=colors4,labels=stat_labels,drop=FALSE)
s_alpha_stat <- scale_alpha_manual(values=rev(c(1.5,0.8,0.5,0.5)),labels=stat_labels,drop=FALSE)
s_size_stat  <-  scale_size_manual(values=rev(c(2.0,1.0,1.0,1.0)),labels=stat_labels,drop=FALSE)
s_shape_stat  <-  scale_shape_manual(values=rev(c(16,1,1,1)),labels=stat_labels,drop=FALSE)

s_color_len <- scale_color_manual(values=colors3)
s_alpha_len <- scale_alpha_manual(values=alpha3)

s_shape_sb <- scale_shape_manual(values=c('S'=1,'B'=0), labels=c('S'='all','B'='worst'))

title_idx = 1

for (ix in 1:n_omes) {
    ome_file <- file_list[ix]
    clust_file <- opt$bad_clust[ix]
    summ_file <- opt$cl_lens[ix]
    
    print("read files:")
    print(paste(c(ix,ome_file,file_list[ix],clust_file,summ_file)))

    taxon = strsplit(ome_file,'_')[[1]][1]
    if (length(opt$title)!=0) {
       taxon_name = opt$title
    } else {
       taxon_name = taxon
    }      

    taxon.names <- append(taxon.names, taxon)

    this_os_info <- os_stats[os_stats$oscode==taxon,]

    print(paste("read_ome_file",ome_file))
    ome_data <- read_ome_file(ome_file, taxon, ome_fields_save, dedup=FALSE)
    print(nrow(ome_data))
    
    print(paste("read_clust_file",clust_file))
    clust_data <- read_cluster_file(clust_file, taxon)
    print(nrow(clust_data))

    print(paste("len_data",summ_file))
    len_data <- read.table(summ_file,sep='\t',header=TRUE)
    print(nrow(len_data))

    ome_data_samp <- ome_data[ome_data$s_type=='S',]
    N_samp = nrow(ome_data_samp)
    ome_data <- ome_data_samp

    print(paste(taxon_name, nrow(ome_data_samp),nrow(len_data)))

    print(summary(ome_data_samp[,c('n_clusters','ns_clusters','pct_cluster','pct_short')]))

    p_quants = quantile(ome_data_samp$pct_cluster,p=seq(0,1.0,0.1))
    print(p_quants)
    print(sprintf("%0.1f",(N_samp*p_quants/100)))

    print(paste("taxon_name:",taxon_name))
    ## print(head(ome_data))

    taxon_label = taxon_name

    taxon_cnt = ''
    if (taxon %in% os_stats$oscode) {
       this_os_info <- os_stats[os_stats$oscode==taxon,]
       print("this_os_info:")
       print(this_os_info)

       taxon_cnt = this_os_info$count
       taxon_clust_cnt = this_os_info$clust_cnt
       if (have_short_name) {
           taxon_label = sprintf("%d proteomes\n%d clusters",taxon_cnt,taxon_clust_cnt)
	   ## "*" produces italics with element_markdown 
	   taxon_title = sprintf("*%s* (%s)",this_os_info$short_name,taxon_name)
       } else {
           taxon_label = sprintf("%d proteomes\n%d clusters",taxon_cnt,taxon_clust_cnt)
	   taxon_title = sprintf("*%s*",taxon_name)
       }
       taxon_cnt_str = sprintf("(%d)",taxon_cnt)
    }

    ome_data$pct_cluster_nz = ifelse(ome_data$pct_cluster < 0.01, 0.01, ome_data$pct_cluster)
    ome_data$B_frag_nz = ifelse(ome_data$B_frag < 0.01, 0.01, ome_data$B_frag)

    ome_data$s_type = factor(ome_data$s_type,levels=c('B','S'),ordered=TRUE)

    low_prot_cnt = nrow(ome_data[ome_data$B_prot_cnt < 1,])
    if (low_prot_cnt > 1) {
        print(paste("taxon:", taxon_name, " rows with B_prot_cnt < 1", low_prot_cnt))
        print(ome_data[ome_data$B_prot_cnt < 1,c('pct_cluster_nz','B_frag','B_prot_status','B_prot_status_f')])
    }

    N_omes = nrow(ome_data)
    NS_omes = nrow(ome_data[ome_data$s_type=='S',])
    NB_omes = nrow(ome_data[ome_data$s_type=='B',])

    print("all, NS, NB rows")
    print(paste(N_omes, NS_omes, NB_omes))

    ome_data_samp <- ome_data_samp[order(ome_data_samp$pct_cluster_nz),]
    ome_data_samp$xrank = 1:nrow(ome_data_samp)/nrow(ome_data_samp)

    ome_data_c2 = ome_data_samp

    ## summary(ome_data_c2[,c('xrank','pct_cluster_nz','s_type')])
    print("ggplot p_bad_pct nrows:")
    print(paste("S:",nrow(ome_data_c2[ome_data_c2$s_type=='S',])))

    p_bad_pct <- ggplot(ome_data_c2,aes(x=xrank,y=pct_cluster_nz,color=s_type)) +
      theme.leg_id + geom_point(shape=1) + ly_scale_pct + x_scale_pct + 
      ggtitle(sprintf("%s. %s",LETTERS[title_idx],taxon_title)) + 
      s_color_sb+ guides(color=guide_legend(position='inside')) +
      theme(legend.position.inside=c(0.95,0.05), legend.justification.inside=c(1,0)) +
      annotate('text',label=taxon_label,x=0.0, y=30.0,hjust=0,vjust=1,size=5)

    title_idx = title_idx + 1

    this_plot = p_bad_pct 

    ## for the lengths plot, merge in the xranks, then facet on clust_id

    u_clust = data.frame('clust_id'=unique(len_data$clust_id),'xpos'=0.05,'ypos'=30)
    u_clust <- merge(u_clust, clust_data[,c('clust_id','p_bad','pct_bad','p_short','pct_short','pct_mode')],by='clust_id')

    len_data <- merge(len_data, ome_data_c2[ome_data_c2$s_type=='S',c('proteome_id', 'xrank')],by.x='gca_id',by.y='proteome_id')

    print("len_data")
    print(head(len_data))

    u_order <- order(-u_clust$pct_bad)
    f_rank <- rank(-u_clust$pct_bad)

    names(f_rank) = u_clust$clust_id
    clust_rank = sprintf("%d",f_rank)
    names(clust_rank) = names(f_rank)
    print("f_rank")
    print(f_rank)

    print("clust_rank")
    print(clust_rank)
    


    u_clust$clust_id <- factor(u_clust$clust_id,levels=u_clust[u_order,]$clust_id, ordered=TRUE)

    print("head(u_clust)")
    print(head(u_clust))

    len_data$clust_id <- factor(len_data$clust_id,levels=u_clust[u_order,]$clust_id, ordered=TRUE)
    
    ## print(len_data)

    len_data$pct_bad = ifelse(((len_data$pct_mode < 75) | (len_data$pct_mode > 133)), 
      'outlier',
      ifelse(len_data$pct_mode==100,'mode','not-mode'))

    len_data$pct_bad = factor(len_data$pct_bad, levels=c('mode', 'not-mode','outlier'),ordered=TRUE)

    tot_counts <- len_data %>% group_by(clust_id) %>% summarize(tot_row=sum(n()))

    out_counts <- len_data %>% group_by(clust_id, pct_bad) %>% summarize(out_row=sum(n())) %>% filter(pct_bad=='outlier')

    print("out_counts")
    print(out_counts,n=Inf)

    u_clust <- u_clust %>% left_join(out_counts, by='clust_id') %>% left_join(tot_counts, by='clust_id')

##    u_clust$clust_label = sprintf("%d (%.1f%%; m=%0.1f%%) s=%d/%d (%.1f%%)",
##      u_clust$p_bad,u_clust$pct_bad.x,u_clust$pct_mode,u_clust$out_row,u_clust$tot_row,(100.0*u_clust$out_row/u_clust$tot_row))
    u_clust$clust_label = sprintf("N~bad~: %d (%.1f%%) mode: %0.1f%%",  u_clust$p_bad,u_clust$pct_bad.x,u_clust$pct_mode)

    print("u_clust")
    print(u_clust)

    n_clust = length(unique(len_data$clust_id))

    p_clust_lens <- ggplot(len_data,aes(x=xrank, y=pct_mode,color=pct_bad,alpha=pct_bad)) +
      theme.leg_id_NL +
      geom_point(shape=3,size=1) + 
      scale_y_continuous(limits=c(20,110),breaks=seq(25,100,25)) +
      x_scale_pct +
      s_color_len + s_alpha_len +            
      ylab("protein length (% mode)") +
      xlab("fraction of sampled proteomes") +
      geom_richtext(data=u_clust, aes(x=xpos,y=ypos-10,label=clust_label),color='black',alpha=1.0,size=4,hjust=0,vjust=0,fill=NA,label.color=NA) +
##      facet_grid(clust_id ~ .)
      facet_wrap(vars(clust_id),strip.position='right',ncol=1, labeller=as_labeller(clust_rank))
 
    plot_list[[ix]] = p_bad_pct / p_clust_lens + plot_layout(height=c(2,6))
}

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (! opt$pub) {
   big_plot <- (Reduce('|', plot_list)) / doc_panel + plot_layout(heights=c(10,0.1))
##   str(big_plot)
} else {
   big_plot = (Reduce('|', plot_list))
   p_heights = c(4, 6)
}

if (is.na(opt$pdf)) {
   ## print(opt$omes)
   print(strsplit(opt$omes,'_'))
   fig3_file=paste0("fig3_omes_len",strsplit(opt$omes,"\\.")[[1]][1],".pdf")
} else {
  fig3_file = opt$pdf
}

if (! opt$pub) {
   ggsave(file=fig3_file, plot=big_plot, width=5*n_omes, height=(1.25*n_clust + 4.2))
} else {
   if (!grepl('_pub',fig3_file)) {
      fig3_file_pref <- strsplit(fig3_file,'\\.')[[1]][1]
      fig3_file <- paste0(fig3_file_pref,"_pub.pdf")
   }
   ggsave(file=fig3_file, plot=big_plot, width=5*n_omes, height=(1.25*n_clust + 4))
}


warnings()
