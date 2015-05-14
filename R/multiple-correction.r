
library(data.table)

window_size <- 1000



#' This function add a new column to the datatable created by bonferroni correction
#' the column name is p_adj_genome. COLUMN required cmh_p_val  
#' @param SNP p value stats data table
#' @return Datable with additional columns for genomewide correction of each column vector

bonferroni_correction_genomewide <- function (in_un_adj_p_val_snps_data_file_path, out_genome_p_adj_file_path) {
  snp_stats_dt <- fread(in_un_adj_p_val_snps_data_file_path, sep="\t", sep2="auto", header=T, na.strings="NA",
                          stringsAsFactors = FALSE, verbose =T, nrow=100000)
  length_dt <- dim(snp_stats_dt)[1]
  for (ch in p_val_character) {
    print(ch)
    adj_name <- paste0(ch,".p_adjusted")
    snp_stats_dt[, (adj_name) := p.adjust(get(b), "fdr", length_dt)]
  }
  write.table(x = snp_stats_dt, file=out_genome_p_adj_file_path, quote = F, sep = "\t", row.names = F)
}




# #' This function is the main driver of the other functions. 
# execute_script <- function () {
#   bonferroni_correction_genomewide("/home/data/nehil_snp_annotated_study_all_snp_ids.tsv", 
#                                    "/home/data/nehil_genome_p_adj_snp_annotated_study.tsv", c("cc_geno","cc_trend"))
# }

