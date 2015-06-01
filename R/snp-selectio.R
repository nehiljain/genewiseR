

#' if the out file path is not given then it returns the datatable else writes a tsv on the new path
get_significant_snps <- function(df, threshold, column_name, out_file_path = NULL) {
  
  see_if(is.number(threshold))
  expect_true( column_name %in% names(df), info = "The column names are not present in the datatable", label = NULL)
  
  signifant_df <- df[get(column_name) < threshold]
  
  if (!is.null(out_file_path)) {
    write.table(x = signifant_df, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(signifant_df)
}


chr_list <- as.character(c(1:29,"X"))


ld_df <- read.table("/home/data/tmp/chr1i/blocks.blocks.det", header=T)
snp_df <- fread("~/sample_nlp_all_snps_in_genes.tsv", sep="\t", sep2="auto", header=T, na.strings="NA", stringsAsFactors = FALSE, verbose =T)
