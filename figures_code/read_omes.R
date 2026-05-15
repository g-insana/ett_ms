##
## 16-July-2025 -- revised read_omes.R to add check_column_values() to
##                 consistently check conditions and correct or remove
##

## check column values applies a cond_vector to a df, and either
## removes rows with cond_vector==TRUE (effect='del'), or assigns the
## effect value to df[cond_vector,col_name], reporting changes with
## explain==TRUE
##

library('dplyr')

check_column_values <- function(df, col_name, cond_vector, effect, description, explain) {

    if (missing(explain)) {
        explain = FALSE
    }

    if (sum(cond_vector) > 0) {
       if (explain) {
       	  cat(sprintf(description,sum(cond_vector)))
       }

       if (effect=='del') {
           df <- df[!cond_vector,]
       } else {
           df[cond_vector, col_name] <- effect
       }
    }

    df
}

read_cluster_file <- function(cluster_file, taxon) {

    this_data <- read.table(cluster_file,header=TRUE,sep='\t')

    this_data$taxon <- taxon
    this_data$row_pct <- nrow(this_data):1
    this_data$f_row_pct <- this_data$row_pct/nrow(this_data)
    this_data$row_pn  <- rank(this_data$p_bad)
    this_data$f_row_pn  <- this_data$row_pn/nrow(this_data)
    this_data$q1q3_fn <- with(this_data,(q3_len - q1_len)/med_len)

    if (have_short_name) {
        this_os_info = os_stats[os_stats$oscode==taxon,]
##	print(paste(taxon, this_os_info, nrow(this_data)))
        this_row_label <- sprintf("*%s* (%s, %d)",this_os_info$short_name, taxon, nrow(this_data))
##	print(this_row_label)
    } else {
        this_row_label <- sprintf("%s (%d)",taxon, nrow(this_data))
    }

    this_row_label <- setNames(this_row_label, taxon)

    clust_return = this_data

##    clust_return = list('df'=this_data, 'row_label'=this_row_label)
    clust_return
}

read_ome_file <- function(ome_file, taxon, ome_fields_save, dedup, explain) {

    if (missing(explain)) {
        explain = FALSE
    }

    if (missing(dedup)) {
        dedup = FALSE
    }

    print(paste("read_ome_file():",ome_file,"\n",str(ome_file)))

    ome_data <- read.table(ome_file,sep='\t',quote='',header=TRUE)

    if (dedup) {
        print(paste0("nrows before dedup:",nrow(ome_data)))
        print(nrow(ome_data[ome_data$s_type=='B',]))
	print(nrow(ome_data[ome_data$s_type=='S',]))
        ome_data <- ome_data[!duplicated(ome_data[,"proteome_id"]),]
        print(paste0("after dedup:",nrow(ome_data)))
        print(nrow(ome_data[ome_data$s_type=='B',]))
        print(nrow(ome_data[ome_data$s_type=='S',]))
    }   

    ome_data$taxon <- taxon

    ## fix 'A_am'
    colnames(ome_data)[colnames(ome_data)=='A_assem_meth']<-'A_am'
    ome_data$A_am <- factor(ome_data$A_am)

    if ('B_prot_status' %in% colnames(ome_data)) {
        ome_data$B_prot_status_f <- sprintf('S%d',ome_data$B_prot_status)

	ome_data$B_prot_status_f <- factor(ome_data$B_prot_status_f,levels=rev(c('S0','S1','S2','S3')),ordered=TRUE)

##	print(paste("testing S2:",taxon))
##	print(summary(ome_data[,c('taxon','B_prot_status','B_prot_status_f')]))

##	ome_data$B_prot_status_f <- factor(ome_data$B_prot_status_f,levels=c('S0','S1','S2','S3'),ordered=TRUE)

	## the taxon=='SALER' coding of B_prot_status_f is wrong.  Need to change
	## S1 ('other') to S2 ('redundant')
##	if (taxon == 'SALER') {
##	   ome_data[ome_data$B_prot_status_f=='S1',]$B_prot_status_f='S2'
##        }
    }

    if ('pct_cluster_nz' %in% ome_fields_save) {
        ome_data$pct_cluster_nz = ifelse(ome_data$pct_cluster < 0.01, 0.01, ome_data$pct_cluster)
    }

    ome_data_cols = colnames(ome_data)
    for (save_col in ome_fields_save) {
        if (!(save_col %in% ome_data_cols)) {
	   print(paste("missing column",save_col))
	}
    }

    ## print(paste("Nrow ome_data (S):",nrow(ome_data[ome_data$s_type=='S',])))
    ## print(head(ome_data[ome_data$s_type=='S',c('pct_cluster','A_n50','A_contig_cnt','A_coverage','B_frag', 's_type')]))

    ## print(paste("Nrow ome_data (B):",nrow(ome_data[ome_data$s_type=='B',])))
    ## print(head(ome_data[ome_data$s_type=='B',c('pct_cluster','A_n50','A_contig_cnt','A_coverage','B_frag','s_type')]))

    ome_data <- check_column_values(ome_data, 'A_coverage',
     ome_data$A_coverage=='',"1x",
     "correcting %d A_coverage.'_'\n", explain)

    ome_data$A_coverage = as.numeric(gsub('x$','',ome_data$A_coverage))

    ome_data <- check_column_values(ome_data, 'A_coverage',
        is.na(ome_data$A_coverage),1,
	"correcting %d A_coverage.NAs after gsub\n", explain)

    ome_data <- check_column_values(ome_data, 'A_coverage',
        ome_data$A_coverage<1, 1,
        "correcting %d A_coverage < 1\n", explain)

    ome_data <- check_column_values(ome_data, 'A_n50',
        is.na(ome_data$A_n50), "del",
	"removing %d A_n50.NAs\n", explain)

    ome_data <- check_column_values(ome_data, A_n50,
        ome_data$A_n50 < 1, 'del',
	"removing %d A_n50<1\n", explain)

    ome_data <- check_column_values(ome_data, A_n50,
        is.na(ome_data$A_contig_cnt),'del',
	"removing %d A_contig_cnt.NAs\n",explain)

    ome_data <- check_column_values(ome_data, A_contig_cnt,
        ome_data$contig_cnt < 1,'del',
	"removing %d A_contig_cnt<1\n", explain)

##    print(colnames(ome_data))
##    print(summary(ome_data[,c('A_coverage','A_n50','A_contig_cnt')]))

    ome_data <- ome_data[,ome_fields_save]
    ## print("read_omes ome_data[ome_fields_save]")
    ## print(head(ome_data))

    ome_data$s_type = factor(ome_data$s_type,levels=c('B','S'),ordered=TRUE)
    ome_data$B_prot_cnt <- ifelse(ome_data$B_prot_cnt>0, ome_data$B_prot_cnt, 1)

    ome_data$A_annot_date <- ifelse(is.na(ome_data$A_annot_date),
        '07/01/1970',ome_data$A_annot_date)
    ome_data$A_annot_date <- ifelse(ome_data$A_annot_date == '',
        '07/01/1970',ome_data$A_annot_date)

    ome_data$A_annot_date <- as.Date(ome_data$A_annot_date,format="%m/%d/%y")

    ## print(summary(ome_data[,c('A_annot_date','A_coverage')]))

    ## unique(ome_data[!is.numeric(ome_data$A_coverage),]$A_coverage)

    roundUp <- function(x) ifelse(!is.na(x),2^ceiling(log2(x)),8)
    ome_data$genome_cov_rnd <- roundUp(ome_data$A_coverage)
    ome_data$A_n50 <- ifelse(ome_data$A_n50 < 1, 1, ome_data$A_n50)
    ome_data$logN50 <- log10(ome_data$A_n50)

    ## print("read_omes summary(ome_data)")
    ## print(summary(ome_data))
    
    ome_data
}

## all done with ome_data, now get tfx_data
## tfx_fields are required because the .tfxg_stats_* files do NOT have a header

read_tfx_file <- function(tfx_file, type) {

    print(paste("In read_tfx_file():",tfx_file))

    taxon = strsplit(tfx_file,'_')[[1]][1]

    tfx_data <- read.table(tfx_file,header=TRUE,sep='\t',quote='')
    
##    print(head(tfx_data))

    tfx_data$taxon <- taxon
    tfx_data$proteome_id <- sub("^.*(GCA_.*)_[BS]_N50.*$","\\1",tfx_data$file)
    tfx_data$s_type <- type
    
##    print(head(tfx_data[,c('file','proteome_id')]))

    tfx_data
}

sample_dedup <- function(ome_data, proteome_count, samp_fract, explain) {

    if (missing(explain)) {
        explain=FALSE
    }

    n_bad = proteome_count * samp_fract
    n_bad <- min(n_bad, nrow(ome_data)/2)

    ome_data_bad <- ome_data[ome_data$s_type=='B',]
    ome_data_bad <- ome_data_bad[1:n_bad,]

    ome_data_samp <- ome_data[ome_data$s_type=='S',]
    ome_data_samp <- ome_data_samp[sample(n_bad),]

    if (explain) {
        print(sprintf("Nrows bad: %d, samp: %d",nrow(ome_data_bad),nrow(ome_data_samp)))
##        print(head(ome_data_bad[,c('proteome_id','pct_cluster','B_frag','A_n50','A_contig_cnt','s_type')]))
##        print(head(ome_data_samp[,c('proteome_id','pct_cluster','B_frag','A_n50','A_contig_cnt','s_type')]))
    }
    ## ome_data_s <- rbind(ome_data_samp, ome_data_bad)

    ome_data <- rbind(ome_data_bad, ome_data_samp)

    cat(sprintf("nrows(ome_data) before sample: %d\n",nrow(ome_data)))

    ome_data <- ome_data[sample(nrow(ome_data)),]

    if (explain) { print(sprintf("Before de-dup %d",nrow(ome_data)))}
    ome_data <- ome_data[!duplicated(ome_data[,"proteome_id"]),]

    if (explain) { print(sprintf("After de-dup %d",nrow(ome_data))) }

    ome_data
}
