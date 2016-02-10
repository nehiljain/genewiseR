
#' This function add a new column to the datatable created by bonferroni correction
#' should have a column 'chr' which is integers only
#' 
#' the column name is p_adj_genome.
#'  
#' @param in_un_adj_p_val_snps_data_file_path is tab separated filepath. This is the file with 
#' unadjusted p valu columns.    
#' @param col_names Vector of columns to be acted upon
#' @param out_file_path Optional parameter. The file path to store the output.
#' @return DataTable with additional columns for genomewide correction of each column vector

p_adjustment <- function (study_df, 
                                     col_names,
                                     level = "genome",
                                     test = "fdr",
                                     out_file_path=NULL ) {
  
  if (!(col_names %in% names(study_df))) {
    stop(paste0("The column 'chr_no' does not exist in the dataframe.
                Please refer to documentation"))
  }
  
  if(!(level %in% c("genome","chromosome") )) {
    stop(paste0("The value of level is not correct.
                Please refer to documentation"))
  }
  
  if(!data.table::is.data.table(study_df)) {
    study_df <- data.table::data.table(study_df)
  }

  if (!is.integer(study_df[,chr_no])) {
    max_int_chr_no <- study_df[!is.na(as.numeric(chr_no)), chr_no]
    max_int_chr_no <- max(as.numeric(max_int_chr_no))
    study_df[tolower(chr_no) == tolower('X'), chr_no := (max_int_chr_no + 1)]
    study_df[tolower(chr_no) == tolower('Y'), chr_no := (max_int_chr_no + 2)]
    study_df[, chr_no := as.integer(chr_no)]
    if (!is.integer(study_df[,chr_no])) {
      stop(paste0("The columns 'chr_no' should be integer, change 'X', 'Y' to numbers and remove any other characters in the column values"))
    }
  }
  
  expect_true( col_names %in% names(study_df), 
               info = "The column names are not present in the datatable", 
               label = NULL)
  length_dt <- nrow(study_df)
  for (ch in col_names) {
    adj_name <- paste0(ch,".adj_",level,"_wide")
    flog.debug(sprintf("Col being adjusted @ genomewide scale %s and Output Col - %s", ch, adj_name))
    
    if (level == "genome") {
      study_df[, (adj_name) := p.adjust(get(ch), test, length_dt)]  
    } else if (level == "chromosome") {
      study_df[, (adj_name) := p.adjust(get(ch), test, length_dt), by = chr_no]
    }
  }
  
  setnames(study_df, names(study_df), norm_var_names(names(study_df)))

  if (!is.null(out_file_path)) {
    assert_that(is.writeable(out_file_path))
    flog.debug(sprintf("Output File path %s", out_file_path))
    write.table(x = study_df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(study_df)
 
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


