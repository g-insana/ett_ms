
get_yaml_opts <- function(args, opt_list, yaml_file) {
      
   if (length(args)==0) {
      opt<-NULL

      if (!is.null(yaml_file) && file.exists(yaml_file)) {
          opt$yaml <- yaml_file
      } else {
	  cat("Cannot get program options\n")
	  sys.exit(1)
      }
   } else {
	opt_parser = OptionParser(option_list=opt_list);
	opt = parse_args(opt_parser);
   }

   if (!is.na(opt$yaml) && length(opt$yaml)>0) {
      y.fd <- file(opt$yaml,'r')
      yaml_opt <- read_yaml(y.fd)
      close(y.fd)

##      print("opt")
##      print(opt)

      ## transfer yaml options ONLY if not on command line
      for (y_col in names(yaml_opt)) {
##          print(paste("yaml_opt",y_col,opt[[y_col]]))
          if ((is.null(opt[[y_col]]) || is.na(opt[[y_col]]) ||
	     (is.logical(opt[[y_col]])) && (!opt[[y_col]]) || ! length(opt[[y_col]])>0 )) {
##	     print(paste("yaml logic", is.null(opt[[y_col]]),is.na(opt[[y_col]]),is.logical(opt[[y_col]]),opt[[y_col]],length(opt[[y_col]])))
##	     print(paste("yaml changing ",y_col,"to",yaml_opt[[y_col]]))
	     opt[y_col] = yaml_opt[y_col]
	  }
      }
   }

#   print(opt)

   opt
}
