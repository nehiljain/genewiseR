

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
  df[, (nlp_column) := -1 * log(get(column_name))]
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
  new_column_names <- c(,max_column, mean_column)
  df [, (max_column) := max(get(column_name)), by="ensembl_gene_id"]
  df [, (mean_column) := mean(get(column_name)),  by="ensembl_gene_id"]
  df [, snp_count := .N,  by="ensembl_gene_id"]
  if (is.null(out_file_path)) {
    return(df)
  } else {
    write.table(x = df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
}