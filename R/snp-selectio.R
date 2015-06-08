library(data.table)
library(assertthat)
library(plyr)
library(dplyr)
library(logging)


snp_df <- fread("/home/tmp/nlp_all_snps_in_genes.tsv", sep="\t", sep2="auto", header=T, na.strings="NA", stringsAsFactors = FALSE, verbose =T)
snp_df <- snp_df[,.(chr_no, snp_pos, cmh_p_val, cmh_p_val.p_adj_genome_wide.nlp, ensemble_gene_id)]
snp_df <- snp_df[chr_no == "1"]

basicConfig()

addHandler(writeToFile, logger="genewise", file="~/coderepo/genewise/genewise_logs.log")
setLevel(50, getHandler('writeToFile', logger='genewise'))


#` file path to data from plink for one chromosome

snp_selection <- function(snps_data, ld_blocks_file_path) {
  
  
  if (!is.data.table(snps_data)) {
    snp_df <- as.data.table(snp_df)
  }
  
  
  ld_df <- read.table(ld_blocks_file_path, header=T)
  ld_dt <- as.data.table(ld_df)
  setnames(ld_dt, 
           names(ld_dt), 
           c("chr_no","gene_start", "gene_end","KB" ,"NSNPS","SNPS"))
  
  ld_dt$chr_no <- as.character(ld_dt$chr_no)
  snp_df$chr_no <- as.character(snp_df$chr_no)
  
  length_ids <- dim(ld_dt)[1]
  custom_names <- rep("pgi_dal_ld_", length_ids)
  custom_names <- paste0(custom_names, seq(1,length_ids))
  ld_dt[, ld_id := custom_names]
  
  result_dt <- map_snps_to_gene(genome_dt = snp_df[,.(chr_no, snp_pos)], ref_dt = ld_dt, window_size=0)
  result_dt[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL]
  result_dt <-  merge(x = snp_df, y= result_dt, by = c("chr_no","snp_pos"), all.x=T)
  result_dt[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL]
  
  result_dt <- result_dt %>% 
    group_by(chr_no, ensemble_gene_id) %>%
    arrange(desc(cmh_p_val.p_adj_genome_wide.nlp), cmh_p_val)
  result_dt <- result_dt[chr_no == "1"]
  result_dt <- unique(result_dt)
  selected_snps <- result_dt[1]
  
  d_ply(newTest, .(ensemble_gene_id, snp_pos), function(df) {
    
    loginfo(sprintf("Processing snp at pos %d ld - %s", df$snp_pos, df$ld_id), logger="genewise.snp-selction")
    if (length(df$ld_id) > 1) {
      logerror("lenght of the row in d_ply loop is not 1", logger="genewise.snp-selction")
    } 
    if (!is.na(df$ld_id) & length(df$ld_id) <= 1) {
      loginfo(sprintf("\nld found - %s with geneid %s and number of enteries in selected snp are %i", df$ld_id, df$ensemble_gene_id, length(selected_snps[ensemble_gene_id == df$ensemble_gene_id, ld_id])), logger="genewise.snp-selction")
      
      if (!(df$ld_id %in% selected_snps[ensemble_gene_id == df$ensemble_gene_id, ld_id])) {
        logdebug(sprintf("\ncase : New LD  %s, highest significance snp, distance not considered", df$ld_id), 
                 logger="genewise.snp-selction")
        cat(!(df$ld_id %in% selected_snps[ensemble_gene_id == df$ensemble_gene_id, ld_id]), str(selected_snps[ensemble_gene_id == df$ensemble_gene_id, ld_id]))
        selected_snps <<- rbindlist(list(selected_snps, df))
      }
    } else {
      logdebug(sprintf("\ncase : LD is NA for snp at pos %i", df$snp_pos), logger="genewise.snp-selction")
      if (!any(abs(selected_snps[ensemble_gene_id == df$ensemble_gene_id, snp_pos] - df$snp_pos) <= 1000)) {
        logdebug(sprintf("\ncase : SNp at distance greater than 1000 in gene %s for snp at pos %i", df$ensemble_gene_id, df$snp_pos), 
                 logger="genewise.snp-selction")
        selected_snps <<- rbindlist(list(selected_snps, df))
      }
    }
  })
  
  
}




# 
# snp_df$chr_no <- as.character(snp_df$chr_no)
# length_ids <- dim(ld_dt)[1]
# custom_names <- rep("pgi_dal_ld_", length_ids)
# custom_names <- paste0(custom_names, seq(1,length_ids))
# 
# print(names(ld_dt))


# result_dt <- map_snps_to_gene(genome_dt = snp_df[,.(chr_no, snp_pos)], ref_dt = ld_dt, window_size=0)
# result_dt[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL]
# test <-  merge(x = snp_df, y= result_dt, by = c("chr_no","snp_pos"), all.x=T)
# test[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL]
# 
# 
# 
# 
# newTest <- test %>% 
#           group_by(chr_no, ensemble_gene_id) %>%
#           arrange(desc(cmh_p_val.p_adj_genome_wide.nlp), cmh_p_val)
# 
# newTest <- newTest[chr_no == "1"]
# newTest <- unique(newTest)
# selected_snps <- newTest[1]
# 
# 
# system.time()
# 
# selected_snps <- as.data.table(selected_snps)
