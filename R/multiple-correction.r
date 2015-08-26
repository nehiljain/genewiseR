
#' This function add a new column to the datatable created by bonferroni correction
#' should have a column 'chr' which is integers only
#' 
#' the column name is p_adj_genome.
#'  
#' @param in_un_adj_p_val_snps_data_file_path is tab separated filepath. This is the file with 
#' unadjusted p valu columns.    
#' @param col_names Vector of columns to be acted upon
#' @param out_file_path Optional parameter. The file path to store the output.
#' @return Datable with additional columns for genomewide correction of each column vector

p_adjustment_genomewide <- function (in_un_adj_p_val_snps_data_file_path, col_names , out_file_path=NULL ) {
  

  
  assert_that(is.readable(in_un_adj_p_val_snps_data_file_path))
  snp_stats_dt <- fread(in_un_adj_p_val_snps_data_file_path, sep="\t", sep2="auto", header=T, na.strings="NA",
                          stringsAsFactors = FALSE, verbose =T)
  
  if (!is.integer(snp_stats_dt[,chr])) {
    stop(paste0("The columns 'chr' should be integer, change 'X', 'Y' to numbers and remove any other characters in the column values"))
  }
  
  expect_true( col_names %in% names(snp_stats_dt), info = "The column names are not present in the datatable", label = NULL)
  length_dt <- nrow(snp_stats_dt)
  for (ch in col_names) {
    adj_name <- paste0(ch,".adj_genome_wide")
    flog.debug(sprintf("Col being adjusted @ genomewide scale %s and Output Col - %s", ch, adj_name))
    snp_stats_dt[, (adj_name) := p.adjust(get(ch), "fdr", length_dt)]
  }
  
  setnames(snp_stats_dt, names(snp_stats_dt), norm_var_names(names(snp_stats_dt)))

  if (!is.null(out_file_path)) {
    assert_that(is.writeable(out_file_path))
    flog.debug(sprintf("Output File path %s", out_file_path))
    write.table(x = snp_stats_dt, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(snp_stats_dt)
 
}



#' column names in order chromosome number, snp id, snp position or base pair, pvalue and adjusted pvalue
#' 

p_adjustment_summary <- function(df, chr_name, snp_id_name, snp_pos_name, p_val_name, report_file_name=NULL) {
  
  df <- data.table(df)
  assert_that(is.data.table(df))
  assert_that(is.numeric(df[, get(chr_name)]))
  assert_that(is.character(df[, get(snp_id_name)]))
  assert_that(is.numeric(df[, get(snp_pos_name)]))
  assert_that(is.numeric(df[, get(p_val_name)]))
  
  setnames(df, c(chr_name,snp_id_name,snp_pos_name,p_val_name), c("CHR","SNP","BP","P"))
  
  
  # plot("p-val-ajustment-summary.pdf")
  manhattan(df)
  qq(df$P)
  qplot(df$P) 
  dev.copy(pdf,"p-val-ajustment-summary.pdf", width=4, height=4)
  dev.off()
  
}






#' This function add a new column to the datatable created by bonferroni correction
#' the column name is p_adj_genome. COLUMN required cmh_p_val  
#' @param SNP p value stats data table
#' @return Datable with additional columns for genomewide correction of each column vector

p_adjustment_chrwide <- function (in_un_adj_p_val_snps_data_file_path, out_file_path=NULL, col_names = c("cmh_p_val"), p_adj_method = "fdr") {
  assert_that(is.readable(in_un_adj_p_val_snps_data_file_path))
  assert_that(is.writeable(out_file_path))
  
  snp_stats_dt <- fread(in_un_adj_p_val_snps_data_file_path, sep="\t", sep2="auto", header=T, na.strings="NA",
                        stringsAsFactors = FALSE, verbose =T)
  
  expect_true( col_names %in% names(snp_stats_dt), info = "The column names are not present in the datatable", label = NULL)
  length_dt <- nrow(snp_stats_dt)
  for (ch in col_names) {
    adj_name <- paste0(ch,".p_adj_chr_wide")
    flog.debug(sprintf("Col being adjusted @ chromosomewide scale %s and Output Col - %s", ch, adj_name))
    snp_stats_dt[, (adj_name) := p.adjust(get(ch), p_adj_method, length_dt), by = chr_no]
  }
  if (!is.null(out_file_path)) {
    flog.debug(sprintf("Output File path %s", out_file_path))
    write.table(x = snp_stats_dt, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(snp_stats_dt)
}


