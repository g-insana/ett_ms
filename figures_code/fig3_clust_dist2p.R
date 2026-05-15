#!/usr/bin/env Rscript --vanilla

################
## clust_dist2m.R --pdf --files=9ENTR_c0.5_p0.90_cm2.bad_clust,ACIBA_c0.5_p0.90_cm2.bad_clust,...
##

## 22-May-2025 -- modified to read field names from .bad_clust file

################
## reads in *.bad_clust

################
## add plot to examine cluster quality vs cluster length

library('ggplot2')
library('stringr')
library('ggtext')
library('optparse')
library('yaml')
library('dplyr')
library('purrr')
library('RColorBrewer')
library('patchwork')

get_yaml_file = 'get_yaml_opts.R'
read_omes_file = 'read_omes.R'

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

## showtext_auto()

## .bad_clust
## ## 0 clusters with >95% of proteins outside length range
## clust_id	n_omes	p_tot	p_bad	pct_bad	p_short	pct_short	bad_SL_flag	mode_len	p_mode	pct_mode	q1_len	med_len	q3_len
## 151465	15183	16434	4864	29.60	159	3.27	Long	212	6749	41.07	212	212	306
## 204569	15447	15454	4314	27.92	4314	100.00	Short	217	9004	58.26	130	217	217
## 205	10431	10445	2752	26.35	2742	99.64	Short	164	7330	70.18	111	164	164

## clust_dist2m.R 

p.name<-'clust_dist2m'
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

## print(opt)

have_full_name <- FALSE
have_short_name <- FALSE


have_fa_stats <- FALSE

os_stats <- NULL
if (!is.na(opt$stats) && length(opt$stats)>0) {
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

   ## print("head(os_stats)")
   ## print(head(os_stats))
}

os_names <- NULL
if (!is.na(opt$names) && length(opt$names)>0) {
   os_names <- read.table(opt$names, header=TRUE, sep='\t')
   ## print("os_stats")
   ## print(head(os_stats))
   os_stats <- merge(os_stats,os_names,by='oscode')
##   print(os_stats)
}

have_full_name <- ('full_name' %in% colnames(os_stats))
have_short_name <- ('short_name' %in% colnames(os_stats))


clust_data = NULL
leg.labels.oscode = c()
leg.labels.nrow = c()
taxon.names = c()

## print(file_list)

mode_good_fn = c()
mode_good_cnt = NULL

report_fields=c('n_omes','p_tot','p_bad','p_short','pct_short','mode_len','p_mode','pct_mode','mode_pct_omecnt')

for (file in file_list) {

    print(file)

    taxon = strsplit(file,'_')[[1]][1]

    taxon.names <- append(taxon.names, taxon)

    this_df <- read_cluster_file(file, taxon)

    ## also calculated modified row fraction that accounts for zero-bad clusters

    leg.labels.oscode <- append(leg.labels.oscode, taxon)
    this_df$mode_pct_omecnt <- 100.0 * this_df$p_mode/os_stats[os_stats$oscode==taxon,]$count

##    print(summary(this_df[report_fields]))

    if (opt$mode_good) {
        old_nrow = nrow(this_df)
        this_df <- this_df[this_df$pct_mode > opt$mode_pct_thresh,]
        new_nrow = nrow(this_df)
	
	mode_str = sprintf("%s\t%d\t%d\t%.1f\n",taxon, old_nrow, new_nrow,100.0*new_nrow/old_nrow)
	mode_str <- setNames(mode_str, taxon)

	mode_good_fn <- append(mode_good_fn, mode_str)
        ## print(sprintf("opt$mode_good filter: old: %d new: %d: %.1f%%\n",old_nrow, new_nrow,100.0*new_nrow/old_nrow))

    	## print(paste(taxon,old_nrow, new_nrow, sprintf("%.1f",100.0*new_nrow/old_nrow)))
##        print(summary(this_df[report_fields]))

	mode_good_cnt <- rbind(mode_good_cnt, data.frame('taxon'=taxon, 'N_orig_clust'=old_nrow,'N_good_clust'=new_nrow, 'pct_good'=100.0*new_nrow/old_nrow))
    }

    this_df$p_long = this_df$p_bad - this_df$p_short
    this_df$pct_long = 100.0 * ifelse(this_df$p_bad > 0, this_df$p_long/this_df$p_bad, 0.0)

    if (have_short_name) {
        this_os_info = os_stats[os_stats$oscode==taxon,]
	this_row_label <- sprintf("*%s* (%s, %d)",this_os_info$short_name, taxon, nrow(this_df))
    } else {
	this_row_label <- sprintf("%s (%d)", taxon, nrow(this_df))
    }
    leg.labels.nrow <- append(leg.labels.nrow, this_row_label)

    clust_data <- rbind(clust_data, this_df)
}

clust_data$taxon <- factor(clust_data$taxon, levels=taxon.names, ordered=TRUE)
clust_data$bad_SL_flag <- factor(clust_data$bad_SL_flag, levels=c('None','Short','Long'), ordered=TRUE)

## summary(clust_data)
## print(head(clust_data))

clust_data$pct_bad_nz <- ifelse(clust_data$pct_bad > 0.0, clust_data$pct_bad, 1e-4)
clust_data$pn_bad_nz <- ifelse(clust_data$p_bad > 0, clust_data$p_bad, 0.1)

clust_data_nz <- clust_data[clust_data$p_bad > 0,]
clust_data_bq1q3 <- clust_data_nz[clust_data_nz$q1q3_fn>0.0,]

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

y_scale_pct = scale_y_log10("percent outlier proteins",breaks=c(0.0001, 0.0010, 0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.000", "0.001", "0.01", "0.1", "1.0", "10.0", "100.0"),limits=c(0.00010,50.0))

x_scale_pct = scale_x_log10("percent outlier proteins",breaks=c(0.0001, 0.0010, 0.01, 0.1, 1.0, 10.0, 100.0),labels=c("0.000", "0.001", "0.01", "0.1", "1.0", "10.0", "100.0"),limits=c(0.00010,50.0))

y_scale_cnt = scale_y_log10("number of outlier proteins",breaks=c(0.1, 1, 10, 100,1000,10000,100000),labels=c("0", "1", "10", "100", "1,000", "10,000", "100,000"),position='right')

x_scale = scale_x_continuous('fraction of clusters',breaks=seq(0,1.0,0.2))

if (is.null(opt$colors)) {
    if (length(leg.labels.nrow) <= 5) {
        q_color_l = scale_color_brewer(palette = "Set1",labels=leg.labels.nrow)
	q_color = scale_color_brewer(palette = "Set1", labels=leg.labels.oscode)
    } else {
        ds_colors = rep(c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"),4)
	q_color_l = scale_color_manual(values=ds_colors[1:length(leg.labels.nrow)],labels=leg.labels.nrow)
	q_color = scale_color_manual(values=ds_colors[1:length(leg.labels.oscode)],labels=leg.labels.oscode)
    }
} else {
    print("have colors")

    print(unlist(opt$colors))

    v_colors = unlist(opt$colors)

    d_colors = v_colors[seq(2,length(v_colors),2)]
    d_colors = setNames(d_colors,v_colors[seq(1,length(v_colors)-1,2)])

    print(d_colors)

    q_color_l = scale_color_manual(values=d_colors,labels=leg.labels.nrow)
    q_color = scale_color_manual(values=d_colors,labels=leg.labels.oscode)
}

ds_shapes = rep(c(0,1,2,5),each=5)
names(ds_shapes[1:length(leg.labels.nrow)]) = leg.labels.nrow
q_shape = scale_shape_manual(values=ds_shapes[1:length(leg.labels.oscode)],labels=leg.labels.oscode)


clust_data_samp20= clust_data[sample(nrow(clust_data),nrow(clust_data)/20),]
clust_data_samp10 = clust_data[sample(nrow(clust_data),nrow(clust_data)/10),]

## how many clusters? how to sample?

clust_data_samp <- clust_data_samp10

color_guide= guide_legend(position='inside',override.aes=list(alpha=1))

## print(leg.labels.nrow)

color_pos = c(0.025,0.975)
theme_pos <- theme(legend.position.inside=color_pos,legend.text=element_markdown(size=8),legend.key.spacing.y=unit(-4,'pt'))
guide_pos <- guides(color=color_guide, shape=color_guide)
if (length(leg.labels.nrow) > 5) {
   print(length(leg.labels.nrow))
   theme_pos <- theme(legend.text=element_markdown(size=8))
   guide_pos <- guides()
}

p_bad_pct <- ggplot(clust_data,aes(x=f_row_pct,y=pct_bad_nz)) +
  theme.leg_id + geom_line(aes(color=taxon)) + y_scale_pct + x_scale + q_color_l + ggtitle("A. cluster quality (percent)") +
  guide_pos + theme_pos

p_bad_cnt <- ggplot(clust_data_samp,aes(x=f_row_pct,y=pn_bad_nz)) +
  theme.leg_id + geom_point(aes(color=taxon, shape=taxon),alpha=0.25,size=2.0) + q_shape + y_scale_cnt + x_scale + q_color + ggtitle("B. cluster quality (number)") + 
  guide_pos + theme_pos

bad_len_cor = with(clust_data, cor(pct_bad_nz, mode_len))
print("bad_len_cor")
print(bad_len_cor)

options(dplyr.width=Inf)
print(mode_good_cnt)

clust_data_stats <- clust_data %>%    
   group_by(taxon) %>%
     summarize(
     med_pct_omecnt = median(mode_pct_omecnt),
     mn_pct_mode = min(pct_mode),
     q1_pct_mode = quantile(pct_mode,p=0.25),
     med_pct_mode = median(pct_mode),
     q3_pct_mode = quantile(pct_mode,p=0.75),
     mx_pct_mode = max(pct_mode),
     med_p_bad = median(p_bad),
     q3_p_bad = quantile(p_bad,prob=0.75),
     med_pct_bad = median(pct_bad),
     q3_pct_bad = quantile(pct_bad,prob=0.75),
     av_pct_short = mean(pct_short),
     med_pct_short = median(pct_short),
     av_pct_long = mean(pct_long),
     med_pct_long = median(pct_long),
     q3_pct_short = quantile(pct_short,prob=0.75)
     )

clust_data_stats <- inner_join(mode_good_cnt, clust_data_stats, by=c("taxon"))

print("median stats")
print(clust_data_stats)
summary(clust_data_stats)

## calculate cummulative p_bad, normmalize to summ=100%, then quantiles of p_bad

cumm_df <- clust_data %>%
     group_by(taxon)  %>%
     mutate(row_num = row_number(),
            fn_row_num = 100.0*row_number()/n(),
	    cum_pct_bad = cumsum(pct_bad),
            total_pct_bad = sum(pct_bad),
	    freq_cum_pct_bad = cum_pct_bad/total_pct_bad) %>%
     ungroup() %>%
     select(c('taxon','clust_id','pct_bad','row_num','fn_row_num','cum_pct_bad','freq_cum_pct_bad'))

print("cumm_df 1%,5%,10%")
cumm_slice <- cumm_df %>% group_by(taxon) %>% slice(c(as.integer(n()/100),as.integer(n()/20),as.integer(n()/10,n())))
print(cumm_slice,n=Inf)

print("cumm_thresh50")
cumm_thresh50 <- cumm_df %>% group_by(taxon) %>% filter(freq_cum_pct_bad > 0.5) %>% slice(1)
print(cumm_thresh50)
print("cumm_thresh25")
cumm_thresh25 <- cumm_df %>% group_by(taxon) %>% filter(freq_cum_pct_bad > 0.25) %>% slice(1)
print(cumm_thresh25)

p_mode_pct <- ggplot(clust_data) + 
	   geom_freqpoly(aes(x=mode_pct_omecnt, color=taxon),binwidth=2) + 
	   q_color_l + theme.leg_id + xlab("cluster proteomes / total proteomes (%)") + ylab("number of clusters") +
	   scale_x_continuous(limits=c(50,100)) +
           guide_pos + theme_pos


doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (!is.na(opt$pdf)) {
   pdf_file = opt$pdf
} else {
   pdf_file = clust_data[1,]$taxon
}

## print(pdf_file)
## print(strsplit(pdf_file,'.',fixed=TRUE))

## print(pdf_file)
pdf_prefix <- strsplit(pdf_file,'.', fixed=TRUE)[[1]]

## print(paste("pdf_prefix:",pdf_prefix[1]))
pdf_mode = paste0(pdf_prefix[1],"_pmode.pdf")

p_width1=7.8
p_width2=4.0
p_height = 4.2
if (length(leg.labels.nrow)>5) {
   p_width1 <- p_width1 + 6.0
   p_width2 <- p_width2 + 4.0
   p_height <- p_height + 2.0
}

if (! opt$pub) {
   p_pct_cnt1 <- (p_bad_pct + p_bad_cnt ) / doc_panel + 
	          plot_layout(heights=c(10, 0.1))

   p_pct_mode <- ( p_mode_pct) / doc_panel + 
	          plot_layout(heights=c(10, 0.1))

   ggsave(file=pdf_file, plot=p_pct_cnt1, width=p_width1, height=p_height+0.2)
##   ggsave(file=pdf_mode, plot=p_pct_mode, width=p_width2, height=p_height+0.2)
} else {

   if (!grepl('_pub',pdf_file)) {
      pdf_file_pref <- strsplit(pdf_file,'\\.')[[1]][1]
      pdf_file <- paste0(pdf_file_pref,"_pub.pdf")
   }

   if (!grepl('_pub',pdf_mode)) {
      pdf_mode_pref <- strsplit(pdf_mode,'\\.')[[1]][1]
      pdf_mode <- paste0(pdf_mode_pref,"_pub.pdf")
   }

   p_pct_cnt1 <- (p_bad_pct + p_bad_cnt)
   p_pct_mode <- ( p_mode_pct) 
   ggsave(file=pdf_file, plot=p_pct_cnt1, width=p_width1, height=p_height)
##   ggsave(file=pdf_mode, plot=p_mode_pct, width=p_width2, height=p_height)
}

warnings()
