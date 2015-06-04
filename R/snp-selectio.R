

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
ld_dt <- as.data.table(ld_df)
snp_df <- fread("~/sample_nlp_all_snps_in_genes.tsv", sep="\t", sep2="auto", header=T, na.strings="NA", stringsAsFactors = FALSE, verbose =T)
snp_df <- snp_df[,.(chr_no, snp_pos, cmh_p_val, cmh_p_val.p_adj_genome_wide.nlp, ensemble_gene_id)]

library(plyr)
library(dplyr)
library(logging)
basicConfig()

addHandler(writeToFile, logger="genewise", file="~/coderepo/genewise/genewise_logs.log")
setLevel(10, getHandler('writeToFile', logger='genewise'))

chr_list <- as.character(c(1:29,"X"))

setnames(ld_dt, 
         names(ld_dt), 
         c("chr_no","gene_start", "gene_end","KB" ,"NSNPS","SNPS"))
ld_dt$chr_no <- as.character(ld_dt$chr_no)

length_ids <- dim(ld_dt)[1]
custom_names <- rep("pgi_dal_ld_", length_ids)
custom_names <- paste0(custom_names, seq(1,length_ids))
ld_dt[, ld_id := custom_names]
print(names(ld_dt))


result_dt <- map_snps_to_gene(genome_dt = snp_df[,.(chr_no, snp_pos)], ref_dt = ld_dt, window_size=0)
result_dt[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL]
test <-  merge(x = snp_df, y= result_dt, by = c("chr_no","snp_pos"), all.x=T)
test[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL]




newTest <- test %>% 
          group_by(chr_no, ensemble_gene_id) %>%
          arrange(desc(cmh_p_val.p_adj_genome_wide.nlp), cmh_p_val)

newTest <- newTest[chr_no == "1"]
selected_snps <- newTest[1]

system.time(selected_snps <- ddply(newTest, .(snp_pos), function(df) {
  loginfo(sprintf("ld - %s", df$ld_id), logger="genewise.snp-selction")
  if (length(df$ld_id) > 1) {
    logerror("lenght of the row in d_ply loop is not 1", logger="genewise.snp-selction")
  } 
  if (!is.na(df$ld_id)) {
    loginfo(sprintf("ld found - %s", df$ld_id), logger="genewise.snp-selction")
    if (!df$ld_id %in% selected_snps[ensemble_gene_id == df$ensemble_gene_id, ld_id]) {
      logdebug(sprintf("case : New LD, highest significance snp distance not considered %s", df$ld_id), 
               logger="genewise.snp-selction")
      selected_snps <- rbindlist(list(selected_snps, df))
    }
  } else {
      logdebug(sprintf("case : LD is NA for snp at pos %i", df$snp_pos), logger="genewise.snp-selction")
      if (!any(abs(selected_snps[ensemble_gene_id == df$ensemble_gene_id, snp_pos] - df$snp_pos) <= 1000)) {
        logdebug(sprintf("case : SNp at distance greater than 1000 in gene %s for snp at pos %i", df$ensemble_gene_id, df$snp_pos), 
                 logger="genewise.snp-selction")
        selected_snps <- rbindlist(list(selected_snps, df))
      }
  }
  return(selected_snps)
}))

selected_snps <- unique(selected_snps)
selected_snps <- as.data.table(selected_snps)