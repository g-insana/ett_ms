#!/usr/bin/env Rscript --vanilla

################
## bad_omes_dist5m.R --yaml bad_omes_dist3m.yaml
## 
## differs from bad_omes_dist4m.R by not deduplicating in read_ome_file(), but waiting until
##   after tfx_file's are read
##
## (1) determine how many bad_omes to plot -- should be the worst 10% of proteomes
## (2) also add tfx searches (separately -- how many??)
## (3) likewise, sample 10% of proteomes
## (4) should be no need to de-duplicate

## much more complex version of bad_omes_dist3r2.R that seeks to do the analysis for multiple datasets simultaneously
##
## to start, yaml file will have lists of sample.am, *.tfxg_stats_S1K, and *.tfxg_stats_BHL
##
################

library('ggplot2')
library('stringr')
## library(tidyr)
## library(cowplot)
library('patchwork')
library('yaml')
library('getopt')
library('optparse')
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

p.name<-'bad_omes_dist5m'
p.full_name <- paste0(p.name,'.R')

cmd_args<- commandArgs(trailingOnly=TRUE)

print(paste(cmd_args))

option_list = list(
	    make_option(c("-D","--debug"),action='store_true',help='debug flag',default=FALSE),
	    make_option(c("--omes"),type='character',action='store',help='comma separated *.summ, files REQUIRED'),
	    make_option(c("--tfx_bad"),type='character', action='store',help='tfx .stats_BHL file', default=NA),
	    make_option(c("--tfx_samp"),type='character', action='store',help='tfx .stats_S1K file', default=NA),
	    make_option(c("-P","--pdf"),type='character',action='store',help='PDF file name', default=NA),
	    make_option(c("-Y","--yaml"),type='character',action='store',help='yaml file', default=NA),
	    make_option(c("-S","--stats"),type='character',action='store',help='os_stats_file', default=NA),
	    make_option(c("--busco"),action='store_true',help='include busco panel', default=FALSE),
	    make_option(c("--pub"),action='store_true',help='pdb for publication', default=FALSE),
	    make_option(c("--frame_thresh"),type='double',help='threshold for 1_seg', default=50.0),
	    make_option(c("--cc_cnt"),type='double',action='store',help='contig count threshold', default=10.0)
	    )

opt <- get_yaml_opts(cmd_args, option_list,"bad_omes_dist3m5b.yaml" )

if (length(opt$files)==1) {
    file_list <- strsplit(opt$files,',')[[1]]
} else {
    file_list <- opt$files
}

if (length(opt$cc_cnt) == 0) {    
   opt$cc_cnt = 10
}

## print(opt)

plabel=paste(c(p.name,"\n",cmd_args),collapse=' ', sep=' ')

## need a function to read each .bad_ome.samp.am file

taxon.names = c()

n_omes = length(opt$omes)
if (length(opt$pub) == 0 || ! opt$pub) {
   opt$pub = FALSE
   plot_list = vector("list",n_omes+1)
} else {
  plot_list = vector("list",n_omes)
}
    
## set up stuff once for plots:
    theme_set(theme_linedraw(base_size=14))
    theme.leg_id <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
    panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
    panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
    plot.title=element_text(face='plain', hjust=0, size=14),
    plot.caption=element_text(size=9,hjust=0),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.position.inside=c(0.05,0.90),
    legend.text=element_text(size=9, hjust=0),
    legend.key=element_blank(),
    legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
    legend.justification=c(0,1),
    legend.title=element_blank())

    theme.leg_id_NL <- theme.leg_id + theme(legend.position='None')

    y_scale_pct = scale_y_log10("bad clusters in proteome (%)",limits=c(0.03,30), breaks=c(0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.00", "0.1", "1.0", "10.0", "100.0"))

    y_scale_pct_r = scale_y_log10("bad clusters in proteome (%)",limits=c(0.03,30), breaks=c(0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.00", "0.1", "1.0", "10.0", "100.0"),position='right')
    
    y_scale_cnt = scale_y_log10("bad clusters in proteome (N)",breaks=c(0.1, 1, 10, 100,1000),labels=c("0", "1", "10", "100", "1,000"))

    x_scale = scale_x_continuous('sampled proteomes')

    colors2 <- c('S'= '#00BFC4', 'B'= '#F8766D' )
    ## colors2 <- c('S'= '#1F78B4', 'B'= '#E31A1C' )

    s_color <- scale_color_manual(values=colors2, labels=c('S'='total','B'='worst'))

    alpha_sshape_labels=c('ND',sprintf('<%.0f 1 segment',opt$frame_thresh),sprintf('>%.0f 1 segment',opt$frame_thresh))

    s_shape <- scale_shape_manual(values=c('no_frame'=1,'H_frame'=2,'L_frame'=6),labels=c('ND','>80% frame','<80% frame'))
    ss_shape <- scale_shape_manual(values=c('no_seg'=0,'1_seg'=1,'N_seg'=3),labels=alpha_sshape_labels)

    s_alpha <- scale_alpha_manual(values=c('no_frame'=0.4,'H_frame'=1.0,'L_frame'=1.0),labels=alpha_sshape_labels)
    ss_alpha <- scale_alpha_manual(values=c('no_seg'=0.2,'1_seg'=1.0,'N_seg'=1.0),labels=alpha_sshape_labels)

## done with plot setup

os_stats <- NULL
if (length(opt$stats)>0) {
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

   if (opt$debug) {
     print("head(os_stats)")
     print(head(os_stats))
   }
}
    
ome_fields <- c('proteome_id','n_clusters','pct_cluster','B_prot_cnt','B_assem_lvl','B_compl_comb','B_compl_single','B_frag','B_miss')
ome_fields_save = c('taxon', ome_fields, 'A_n50', 'A_contig_cnt', 'A_annot_date', 'A_coverage', 's_type','A_am')

title_idx = 1

for (ix in seq_along(opt$omes)) {    

    ome_file = opt$omes[ix]

    print(paste("Reading ",ome_file))

    taxon = strsplit(ome_file,'_')[[1]][1]
    taxon.names <- append(taxon.names, taxon)
    ## get ome_data for taxon

    this_os_info <- os_stats[os_stats$oscode==taxon,]

    ome_data <- read_ome_file(ome_file, taxon, ome_fields_save, dedup=FALSE)

    ## ome_data is a mixture of bad omes and sampled omes
    ## (1) isolate the bad omes (vs the sampled omes)

    ome_data <- sample_dedup(ome_data, this_os_info$count, 0.2, explain=TRUE)

    ## get tfx_bad
    tfx_bad_file <- opt$tfx_bad[ix]
    tfx_bad <- read_tfx_file(tfx_bad_file, 'BHL')

    ## get tfx_samp
    tfx_samp_file <- opt$tfx_samp[ix]
    tfx_samp <- read_tfx_file(tfx_samp_file,'S1K')

    ## extract the GCA acc from the file name
    ## label s_type

    tfx_all <- rbind(tfx_bad, tfx_samp)

    tfx_all$pct_fsTtot <- 0.0
    # tfx_all$pct_fsTtot <- ifelse(tfx_all$tot_hits > 0, 100.0*tfx_all$X1cnt/tfx_all$tot_hit, 0.0)
    tfx_all$pct_fsTtot <- ifelse(tfx_all$tot_hit > 0, 100.0*tfx_all$X1cnt/tfx_all$tot_hit, 0.0)


##    tfx_all$frame_fact = ifelse(tfx_all$pct_fs >= opt$frame_thresh, 'H_frame','L_frame')
    tfx_all$frame_fact = ifelse(tfx_all$pct_fsTtot >= opt$frame_thresh, '1_seg','N_seg')
    tfx_all$frame_fact = factor(tfx_all$frame_fact, levels=c('no_seg','1_seg','N_seg'))

    print("aggregate(tfx_all$pct_fsTtot)")
    print(aggregate(pct_fsTtot ~ frame_fact + s_type, tfx_all,NROW))

    ome_data <- merge(ome_data, tfx_all[,c('proteome_id','pct_fs','pct_fsTtot')],all.x=TRUE)

    ## add factor reflecting fraction 1_seg
    ome_data$frame_fact <- 'no_seg'
    ome_data$frame_fact <- ifelse(!is.na(ome_data$pct_fsTtot), ifelse(ome_data$pct_fsTtot > opt$frame_thresh, '1_seg','N_seg'),'no_seg')
    ome_data$frame_fact <- factor(ome_data$frame_fact, levels=c('no_seg','1_seg','N_seg'))

    clust_comb_corr <- with(ome_data[ome_data$s_type=='S',],cor(pct_cluster,B_frag))
    
    ## str(clust_comb_corr)
    ## summary(clust_comb_corr)

    ## ome_data_s <- ome_data[as.character(ome_data$pct_fsTtot)!='no_seg',]
    ome_data_s <- ome_data

    ## print("ome_data_s nrow(ND, 1_seg, N_seg)")
    ## print(nrow(ome_data[as.character(ome_data$frame_fact)=='no_seg',]))
    ## print(nrow(ome_data[as.character(ome_data$frame_fact)=='1_seg',]))
    ## print(nrow(ome_data[as.character(ome_data$frame_fact)=='N_seg',]))

    ## print("head/tail ome_data_s 1_seg N_seg")
    ## ome_1_seg = ome_data_s[as.character(ome_data_s$frame_fact)=='1_seg',]
    ## ome_N_seg = ome_data_s[as.character(ome_data_s$frame_fact)=='N_seg',]
    ## print(colnames(ome_1_seg))
    ## print(head(ome_1_seg[,c('proteome_id','logN50','pct_cluster','B_frag','A_contig_cnt','A_coverage','frame_fact')]))
    ## print(tail(ome_1_seg[,c('proteome_id','logN50','pct_cluster','B_frag','A_contig_cnt','A_coverage','frame_fact')]))
    ## print(head(ome_N_seg[,c('proteome_id','logN50','pct_cluster','B_frag','A_contig_cnt','A_coverage','frame_fact')]))
    ## print(tail(ome_N_seg[,c('proteome_id','logN50','pct_cluster','B_frag','A_contig_cnt','A_coverage','frame_fact')]))

    taxon_label_omes = sprintf("%s (%d proteomes)",taxon.names[ix],os_stats[os_stats$oscode==taxon.names[ix],]$count)

    clust_pct_corr <- with(ome_data,cor(pct_cluster,B_frag))
    if (is.na(clust_pct_corr)) {
        print("*** using sampled data ***")
        clust_pct_corr <- with(ome_data[sample(nrow(ome_data),nrow(ome_data)/10),],cor(pct_cluster,B_frag))
    }

    if (is.na(clust_pct_corr)) {
        print("*** using worst data ***")
	clust_pct_corr <- with(ome_data[ome_data$s_type == 'B',],cor(pct_cluster,B_frag))
    }

    print("clust_pct_corr")
    nrow(ome_data[ome_data$B_prot_status != 'S2',])
    print(clust_pct_corr)

    ome_data_samp <- ome_data[ome_data$s_type=='S',]

    print(paste("clust_assm_corr sum:",sum(is.na(log10(ome_data_samp$pct_cluster))),sum(is.na(log10(ome_data_samp$A_n50)))))
    print(nrow(ome_data_samp))
    print(sum(is.na(log10(ome_data_samp$pct_cluster & is.na(log10(ome_data_samp$A_n50))))))
    print(head(ome_data_samp[,c('pct_cluster','A_n50')]))

    clust_assm_corr <- with(ome_data[ome_data$s_type=='B',],cor(log10(pct_cluster),log10(A_n50),use="complete.obs"))
    print("clust_assm_corr B")
    print(clust_assm_corr)

    clust_assm_corr <- with(ome_data[ome_data$s_type=='S',],cor(log10(pct_cluster),log10(A_n50),use="complete.obs"))
    print("clust_assm_corr S")
    print(clust_assm_corr)
    
    clust_assm_corr_short <- with(ome_data[ome_data$logN50 < 6,],cor(log10(pct_cluster),log10(A_n50),use="complete.obs"))
    print("clust_assm_corr_short")
    print(clust_assm_corr_short)

    p_bad_logN50 <- ggplot(ome_data_s,aes(y=pct_cluster,x=logN50, color=s_type, shape=frame_fact, alpha=frame_fact)) +
    ggtitle(sprintf("%s. %s",LETTERS[title_idx],taxon_label_omes)) +
    theme.leg_id_NL + 
    geom_point(size=1,position=position_dodge(0.01)) + labs(x="log10(N50)") + s_color + ss_shape + ss_alpha + y_scale_pct + scale_x_continuous(limits=c(3,7)) +
    annotate('text',label=paste("r^2==",sprintf("%.3f",clust_assm_corr)),x=3.0, y=0.08,parse=TRUE, hjust=0) +
    annotate('text',label=paste("r^2==",sprintf("%.3f",clust_assm_corr_short),"~(log(N50)<6)"),x=3.0, y=0.04,parse=TRUE, hjust=0)

    title_idx = title_idx + 1

##    print("head(ome_data)")
##    print(head(ome_data[ome_data$A_contig_cnt > opt$cc_cnt,c('taxon','pct_cluster','A_contig_cnt')]))

    clust_afrag_corr <- with(ome_data[ome_data$A_contig_cnt>0,],cor(log10(pct_cluster),log10(A_contig_cnt),use="complete.obs"))
    
    clust_afrag_corr_cc10 <- with(ome_data[ome_data$A_contig_cnt > opt$cc_cnt,],cor(log10(pct_cluster),log10(A_contig_cnt),use="complete.obs"))
    print("clust_afrag_corr (log10)")
    print(clust_afrag_corr)
    print(clust_afrag_corr_cc10)

    p_bad_contig <- ggplot(ome_data_s[ome_data_s$A_contig_cnt>0,],aes(y=pct_cluster,x=A_contig_cnt, color=s_type, shape=frame_fact, alpha=frame_fact)) +
    ggtitle(sprintf("%s.",LETTERS[title_idx])) +
    theme.leg_id_NL +geom_point(size=1,position=position_dodge(0.01)) +
    y_scale_pct +theme(axis.text.y=element_blank(),axis.title.y=element_blank()) + scale_x_log10("contig count",limits=c(1,5000)) + s_color + ss_shape + ss_alpha +
    ##  guides(color=guide_legend(position='inside'))+theme(legend.position.inside=c(0.2,0.4)) +
    annotate('text',label=paste("r^2 ==",sprintf("%.3f",clust_afrag_corr)),y=0.08, x=1,parse=TRUE,hjust=0) +
    annotate('text',label=paste("r^2 ==",sprintf("%.3f~~(cc>%.0f)",clust_afrag_corr_cc10,opt$cc_cnt)),y=0.04, x=1.0,parse=TRUE,hjust=0)

    title_idx = title_idx + 1

    print("head(ome_data)")
    print(head(ome_data))
    summary(ome_data)

##    clust_gcov_corr <- with(ome_data[ome_data$s_type=='S' & ome_data$A_coverage>1,],cor(pct_cluster,log10(A_coverage),use="complete.obs"))
    print("head(ome_data pct_cluster A_coverage)")
    print(head(ome_data[ome_data$A_coverage>1,c('pct_cluster','A_coverage')]))

    clust_gcov_corr <- with(ome_data[ome_data$A_coverage>1,],cor(pct_cluster,log10(A_coverage),use="complete.obs"))
    print("clust_gcov_corr")
    print(clust_gcov_corr)

    p_bad_gcov <- ggplot(ome_data_s[ome_data_s$A_coverage > 1,],aes(y=pct_cluster,x=A_coverage, color=s_type, shape=frame_fact, alpha=frame_fact)) +
    ggtitle(sprintf("%s.",LETTERS[title_idx])) +
    theme.leg_id_NL +geom_point(size=1,position=position_dodge(0.01)) +
    y_scale_pct +theme(axis.text.y=element_blank(), axis.title.y=element_blank()) + scale_x_log10("genome coverage", limits=c(1.0,1000)) + s_color + ss_shape + ss_alpha +
    ##  guides(color=guide_legend(position='inside'))+theme(legend.position.inside=c(0.7,0.7)) +
    annotate('text',label=paste("r^2 ==",sprintf("%.3f~~(cov > 1)",clust_gcov_corr)),y=0.04, x=1.0 ,parse=TRUE,hjust=0)

    title_idx = title_idx + 1

    p_bad_busco <- ggplot(ome_data_s,aes(y=pct_cluster,x=B_frag, color=s_type, shape=frame_fact, alpha=frame_fact)) +
    ggtitle(sprintf("%s.",LETTERS[title_idx])) +
      theme.leg_id_NL + geom_point(size=1,position=position_dodge(0.01)) +
      y_scale_pct_r + scale_x_continuous("BUSCO fragment score",limits=c(0,42)) + s_color + ss_shape + ss_alpha +
      annotate('text',label=paste("r^2 ==",sprintf("%.3f",clust_pct_corr)),y=0.08, x=2, hjust=0, parse=TRUE)

    if (opt$busco) {
        title_idx = title_idx + 1
    }	

    ## doc_panel = ggplot() + theme_void() + labs(caption=plabel)

    if (!opt$busco) {
        p_col = 3
        p_qual_pct <-  p_bad_logN50 + p_bad_contig + p_bad_gcov + plot_layout(ncol=p_col)
#	p_col = 1
#	p_qual_pct <-  p_bad_gcov + plot_layout(ncol=p_col)

    } else {
       p_col = 4
       p_qual_pct <-  p_bad_logN50 + p_bad_contig + p_bad_gcov + p_bad_busco + plot_layout(ncol=p_col)
    }

    plot_list[[ix]] = p_qual_pct
}

doc_panel = ggplot() + theme_void() + labs(caption=plabel)


## now use wrap_plots plot_list

if (is.na(opt$pdf)) {
   print(opt$omes)
   print(strsplit(opt$omes,'_'))
   fig1i_file=paste0("fig1r_omes_",strsplit(opt$omes,"\\.")[[1]][1],".pdf")
   fig1i0_file=paste0("fig1r0_omes_",strsplit(opt$omes,"\\.")[[1]][1],".pdf")
   fig1i00_file=paste0("fig1r00_omes_",strsplit(opt$omes,"\\.")[[1]][1],".pdf")
} else {
  fig1i_file = opt$pdf
}

if (! opt$pub) {
   plot_list[[n_omes+1]] = doc_panel   
   p_heights = c(rep(10.0,n_omes),0.1)
} else {
   if (!grepl('_pub',fig1i_file)) {
      fig1i_file_pref <- strsplit(fig1i_file,'\\.')[[1]][1]
      fig1i_file <- paste0(fig1i_file_pref,"_pub.pdf")
   }
   p_heights = c(rep(10.0,n_omes))
}

## print(fig1i_file)
## print(length(plot_list))

big_plot = Reduce('/', plot_list) + plot_layout(heights=p_heights)
    
ggsave(file=fig1i_file, plot=big_plot, width=p_col*3.5, height=4.0*length(opt$omes))

warnings()
