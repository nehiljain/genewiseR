library(Hmisc)
library(assertthat)


#' if the out file path is not given then it returns the datatable else writes a tsv on the new path
get_significant_snps <- function(df, threshold, column_name, out_file_path = NULL) {
  
  see_if(is.number(threshold))
#   assert_that( column_name %in% names(df), info = "The column names are not present in the datatable", label = NULL)
  
  signifant_df <- df[get(column_name) < threshold]
  
  if (!is.null(out_file_path)) {
    write.table(x = signifant_df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(signifant_df)
}




#' if the out file path is not given then it returns the datatable else writes a tsv on the new path
get_nlp <- function(df, column_name, out_file_path = NULL) {
  
  expect_true( column_name %in% names(df), info = "The column names are not present in the datatable", label = NULL)
  
  df[is.na(df[,column_name]), column_name] <- 0
  nlp_column <- paste0(column_name, ".nlp")
  df[, (nlp_column) := -1 * log10(get(column_name))]
  if (!is.null(out_file_path)) {
    write.table(x = df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(df)
}


#' calculates snp count, max and mean on given column name and groups all the counts by chromosome
#' if the out file path is not given then it returns the datatable else writes a tsv on the new path
get_max_and_mean <- function(df, column_name, out_file_path = NULL) { 
  
  expect_true( column_name %in% names(df), info = "The column names are not present in the datatable", label = NULL)
  df[is.na(df[,column_name]), column_name] <- 0
  max_column <- paste0( "chr_max_",column_name)
  mean_column <- paste0( "chr_mean_",column_name)
  
  df [, (max_column) := max(get(column_name),na.rm = T), by="ensemble_gene_id"]
  df [, (mean_column) := mean(get(column_name),na.rm = T),  by="ensemble_gene_id"]
  df [, snp_count := .N,  by="ensemble_gene_id"]
  
  #tricky step
  df <- unique(df[,.(ensemble_gene_id, maxT = get(max_column), meanT = get(mean_column), snp_count)])
  
  if (!is.null(out_file_path)) {
    write.table(x = df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(df)
}




get_quartile <- function(df, column_name, quartile = 25) {
  
  assert_that(sum(column_name %in% names(df)) == length(column_name))
  
  df <- as.data.frame(df)
  df[is.na(df[,column_name]), column_name] <- 0
  
  
  
  
  if (dim(df)[1] == 1 | length(unique(df[,column_name])) == 1) {
    
    df$quartile <- df[,column_name]
    
  } else{
    df$quartile <- with(df, cut(get(column_name),
                                breaks=unique(quantile(get(column_name), probs=seq(0,1, by=(quartile/100)))), 
                                include.lowest=TRUE))

    df$quartile <- as.factor(df$quartile)
    last_quartile <- tail(levels(df$quartile), n=1)
    df$quartile_ch <- as.character(df$quartile)

    df <- filter(df, as.character(quartile_ch) == as.character(last_quartile))

  }
  
  df <- as.data.table(df)
  topq_column_name <- paste0("topQ_",quartile,"_nlp")
  df[, (topq_column_name) := mean(get(column_name), na.rm = TRUE)]
  df <- df[, .(ensemble_gene_id, get(topq_column_name))]  
  setnames(df, "V2", topq_column_name)
  
  return(df)
}


get_topX_sample <- function(df, column_name, quartile = 25) {
  
  
  assert_that(sum(column_name %in% names(df)) == length(column_name))
  
  df <- as.data.table(df)
  df <- df[,.(ensemble_gene_id, nlp = get(column_name))]
  df[is.na(nlp), nlp := 0]
  df[order(-nlp)]
  
  top_subset <- round(quartile/100 * dim(df)[1])
 
  topq_column_name <- paste0("topQ_",quartile,"_nlp")
  df[1:top_subset, (topq_column_name) := mean(nlp, na.rm = TRUE)]
  df <- df[, .(ensemble_gene_id, get(topq_column_name))]  
  setnames(df, "V2", topq_column_name)
  
  return(df)
}







get_topX_subset <- function(df, column_name, percent = 25) {
  
  
  assert_that(sum(column_name %in% names(df)) == length(column_name))
  
  df <- as.data.table(df)
  df <- df[,.(ensemble_gene_id, nlp = get(column_name))]
  df[is.na(nlp), nlp := 0]
  df[order(-nlp)]
  
  threshold_nlp <- percent/100 * max(df[,nlp])
  
  
  topq_column_name <- paste0("topQ_",quartile,"_nlp")
  df[nlp > threshold_nlp, (topq_column_name) := mean(nlp, na.rm = TRUE)]
  df <- df[, .(ensemble_gene_id, get(topq_column_name))]  
  setnames(df, "V2", topq_column_name)
  
  return(df)
}








#' Top Q statistics.
#' Uses Hmisc cut2 to get mean of nlp(negative log p-value) of snps in the top quartile of each gene
#' 
#'  @param df data.table with the snps maps to genes
#'  @param column_name the name of the column which have negative log p-values
#'  @param quartile the quartile you want to use for topQ stats
#'  @param out_file_path The file where  output will be writen to
#' Note: Internally calls another function get_quartile

get_topQ <- function(df, column_name, threshold = 25, out_file_path = NULL) {
  assert_that(is.data.table(df))
  expect_true( column_name %in% names(df), info = "The column names are not present in the datatable", label = NULL)
#   print(str(df))
  result_sign_snp_topq_df <- ddply(df, "ensemble_gene_id", function(df) {
#       print(str(df))
    return(get_quartile(df, column_name, threshold))
  })
  
  result_sign_snp_topq_df <- unique(result_sign_snp_topq_df)
#   str(result_sign_snp_topq_df)
  if (!is.null(out_file_path)) {
    write.table(x = result_sign_snp_topq_df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(result_sign_snp_topq_df)
}

