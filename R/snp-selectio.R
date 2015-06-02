

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


ld_df <- read.table("/home/data/tmp/chr1i/blocks.blocks.det", header=T)
snp_df <- fread("~/sample_nlp_all_snps_in_genes.tsv", sep="\t", sep2="auto", header=T, na.strings="NA", stringsAsFactors = FALSE, verbose =T)
snp_df <- snp_df[,.(chr_no,ensemble_gene_id,snp_pos,cmh_p_val,cmh_p_val.p_adj_genome_wide)]

library(plyr)
library(dplyr)

chr_list <- as.character(c(1:29,"X"))

setnames(ld_df, 
         names(ld_df), 
         c("chr_no","gene_start", "gene_end","KB" ,"NSNPS","SNPS"))
print(names(ld_df))
#       ld_df <- ld_df[chromosome_name %in% chr_list]

result_dt <- map_snps_to_gene(snp_df, ld_df, window_size=0)
result_dt[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL]
write.table(result_dt, out_snps_in_genes_file_path, sep=sep3, row.names=F, quote = F)




snp_select_df <- ddply()

