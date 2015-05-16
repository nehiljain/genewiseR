library(Hmisc)

#' if the out file path is not given then it returns the datatable else writes a tsv on the new path
get_significant_snps <- function(df, threshold, column_name, out_file_path = NULL) {
  signifant_df <- df[get(column_name) < threshold]
  if (is.null(out_file_path)) {
    return(signifant_df)
  } else {
    write.table(x = signifant_df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
}




#' if the out file path is not given then it returns the datatable else writes a tsv on the new path
get_nlp <- function(df, column_name, out_file_path = NULL) {
  nlp_column <- paste0(column_name, ".nlp")
  df[, (nlp_column) := -1 * log10(get(column_name))]
  if (is.null(out_file_path)) {
    return(df)
  } else {
    write.table(x = df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
}


#' calculates snp count, max and mean on given column name and groups all the counts by chromosome
#' if the out file path is not given then it returns the datatable else writes a tsv on the new path
get_max_and_mean <- function(df, column_name, out_file_path = NULL) { 
  max_column <- paste0( "chr_max_",column_name)
  mean_column <- paste0( "chr_mean_",column_name)
  
  df [, (max_column) := max(get(column_name)), by="ensemble_gene_id"]
  df [, (mean_column) := mean(get(column_name)),  by="ensemble_gene_id"]
  df [, snp_count := .N,  by="ensemble_gene_id"]
  #tricky step
  df <- unique(df[,.(ensemble_gene_id, maxT = get(max_column), meanT = get(mean_column), snp_count)])
  
  if (is.null(out_file_path)) {
    return(df)
  } else {
    write.table(x = df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
}




get_quartile <- function(df, column_name) {
  df <- as.data.frame(df)
  if (dim(df)[1] == 1) {
    
    df$quartile <- df[,column_name]
    
  } else{
    cat("Before Quartile", str(df))
    df$quartile <- with(df, cut(get(column_name),
                                breaks=unique(quantile(get(column_name), probs=seq(0,1, by=0.25))), 
                                include.lowest=TRUE))
    df$quartile <- as.factor(df$quartile)
    last_quartile <- tail(levels(df$quartile), n=1)
    cat("\n\nlast quartile\n\n",last_quartile)
    cat("\n\nlevels of quartiles\n\n",levels(df$quartile))
    df$quartile_ch <- as.character(df$quartile)
    df <- filter(df, as.character(quartile_ch) == as.character(last_quartile))
    cat("\n\ndf after fildering\n\n",str(df))
  }
  df <- as.data.table(df)
#   cat("\n\nbefore mean of all\n\n",str(df))
  df[, topQ_nlp := mean(get(column_name), na.rm = TRUE)]
  df <- df[, .(ensemble_gene_id, topQ_nlp)]
  return(df)
}


get_topQ <- function(df, column_name, out_file_path = NULL) {
  result_sign_snp_topq_df <- ddply(df, "ensemble_gene_id", function(df) {
    #   print(str(df))
    return(get_quartile(df, column_name))
  })
  
  result_sign_snp_topq_df <- unique(result_sign_snp_topq_df)
  str(result_sign_snp_topq_df)
  if (is.null(out_file_path)) {
    return(result_sign_snp_topq_df)
  } else {
    write.table(x = result_sign_snp_topq_df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
}

