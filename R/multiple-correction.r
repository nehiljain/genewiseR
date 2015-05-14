
library(data.table)

window_size <- 1000



#' This function add a new column to the datatable created by bonferroni correction
#' the column name is p_adj_genome. COLUMN required cmh_p_val  
#' @param SNP p value stats data table
#' @return Datable with additional columns for genomewide correction of each column vector

p_adjustment_genomewide <- function (in_un_adj_p_val_snps_data_file_path, out_genome_p_adj_file_path, col_names) {
  snp_stats_dt <- fread(in_un_adj_p_val_snps_data_file_path, sep="\t", sep2="auto", header=T, na.strings="NA",
                          stringsAsFactors = FALSE, verbose =T)
  length_dt <- dim(snp_stats_dt)[1]
  for (ch in col_names) {
    print(ch)
    adj_name <- paste0(ch,".p_adj_genome_wide")
    snp_stats_dt[, (adj_name) := p.adjust(get(ch), "fdr", length_dt)]
  }
  write.table(x = snp_stats_dt, file=out_genome_p_adj_file_path, quote = F, sep = "\t", row.names = F)
}


#' This function add a new column to the datatable created by bonferroni correction
#' the column name is p_adj_genome. COLUMN required cmh_p_val  
#' @param SNP p value stats data table
#' @return Datable with additional columns for genomewide correction of each column vector

p_adjustment_chrwide <- function (in_un_adj_p_val_snps_data_file_path, out_genome_p_adj_file_path, col_names) {
  snp_stats_dt <- fread(in_un_adj_p_val_snps_data_file_path, sep="\t", sep2="auto", header=T, na.strings="NA",
                        stringsAsFactors = FALSE, verbose =T)
  length_dt <- dim(snp_stats_dt)[1]
  for (ch in col_names) {
    print(ch)
    adj_name <- paste0(ch,".p_adj_chr_wide")
    snp_stats_dt[, (adj_name) := p.adjust(get(ch), "fdr", length_dt), by = chr_no]
  }
  write.table(x = snp_stats_dt, file=out_genome_p_adj_file_path, quote = F, sep = "\t", row.names = F)
}



# #' This function is the main driver of the other functions. 
execute_script <- function () {
  p_adjustment_chrwide("/home/data/nehil_snp_annotated_study_all_snp_ids.tsv", 
                                   "/home/data/nehil_genome_p_adj_snp_annotated_study.tsv", c("cc_geno","cc_trend"))
}

