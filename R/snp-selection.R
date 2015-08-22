
#` @examples
#` snp_df <- fread("/home/tmp/nlp_all_snps_in_genes.tsv", sep="\t", sep2="auto", header=T, na.strings="NA", stringsAsFactors = FALSE, verbose =T)
#` snp_df <- snp_df[,.(chr_no, snp_pos, cmh_p_val, cmh_p_val.p_adj_genome_wide.nlp, ensemble_gene_id)]
#` snp_df <- snp_df[chr_no == "1"]
#` np_selection(snps_data = snp_df, ld_blocks_file_path = "/home/data/tmp/chr1i/blocks.blocks.det")
#` 
#` @param ld_blocks_file_path file path to data from plink for one chromosome

snp_selection <- function(snps_data, ld_blocks_file_path, significance_threshold = -1, chr_no_i, p_val_col_name = "cmh_p_val", out_file_path = NULL) {

  
  if (!is.data.table(snps_data)) {
    snps_data <- as.data.table(snps_data)
  }
  
  ld_df <- read.table(ld_blocks_file_path, header=T)
  ld_df <- as.data.table(ld_df)
  setnames(ld_df, 
           names(ld_df), 
           c("chr_no","gene_start", "gene_end","KB" ,"NSNPS","SNPS"))
  
  assert_that(is.numeric(chr_no_i) && !is.null(chr_no_i))
  assert_that(is.data.table(ld_df))
  assert_that(is.data.table(snps_data))
  
  ld_df[,chr_no := as.character(chr_no)]
  snps_data[,chr_no := as.character(chr_no)]
  
  if (significance_threshold > 0) {
    flog.debug(sprintf("\nsignificant snps only\n Threshold: %i", p_val_col_name))
    snps_data <- get_significant_snps(snps_data, significance_threshold, p_val_col_name)
  }
  
  length_ids <- dim(ld_df)[1]
  custom_names <- rep("pgi_dal_ld_", length_ids)
  custom_names <- paste0(custom_names, seq(1,length_ids))

  ld_df[,ld_id := custom_names]
  
  
  snps_data <- as.data.table(snps_data)
  temp <- snps_data[,.(chr_no, snp_pos)]

  result_dt <- map_snps_to_gene(genome_dt = temp , ref_dt = ld_df, window_size=0)
  result_dt[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL, with=FALSE]
  
  ld_snp_merge_dt <-  merge(x = snps_data, y= result_dt, by = c("chr_no","snp_pos"), all.x=T, )
  ld_snp_merge_dt[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL, with=FALSE]

  ld_snp_merge_dt <- ld_snp_merge_dt %>% 
                      group_by(chr_no, ensemble_gene_id) %>%
                      arrange(desc(cmh_p_val.p_adj_genome_wide.nlp), cmh_p_val)
  ld_snp_merge_dt <- ld_snp_merge_dt[chr_no == chr_no_i]
  ld_snp_merge_dt <- unique(ld_snp_merge_dt)
  ld_snp_merge_dt[,gene_start.y := NULL]
  ld_snp_merge_dt[,gene_end.y := NULL]
  selected_snps <- ld_snp_merge_dt[1]
 
  
  d_ply(ld_snp_merge_dt[!is.na(cmh_p_val)], .(ensemble_gene_id, snp_pos), function(df) {
    
    flog.info(sprintf("Processing snp at pos %d ld - %s", df$snp_pos, df$ld_id))
    if (length(df$ld_id) > 1) {
      flog.error("lenght of the row in d_ply loop is not 1")
    } 
    if (!is.na(df$ld_id) & length(df$ld_id) <= 1) {
      flog.debug(sprintf("\nld found - %s with geneid %s and number of enteries in selected snp are %i", df$ld_id, df$ensemble_gene_id, length(selected_snps[ensemble_gene_id == df$ensemble_gene_id, ld_id])))
      
      if (!(df$ld_id %in% selected_snps[ensemble_gene_id == df$ensemble_gene_id, ld_id])) {
        flog.debug(sprintf("\ncase : New LD  %s, highest significance snp, distance not considered", df$ld_id), 
                 "")
        
        selected_snps <<- rbindlist(list(selected_snps, df))
      }
    } else {
      flog.debug(sprintf("\ncase : LD is NA for snp at pos %i", df$snp_pos))
      if (!any(abs(selected_snps[ensemble_gene_id == df$ensemble_gene_id, snp_pos] - df$snp_pos) <= 1000)) {
        flog.debug(sprintf("\ncase : SNp at distance greater than 1000 in gene %s for snp at pos %i", df$ensemble_gene_id, df$snp_pos))
        selected_snps <<- rbindlist(list(selected_snps, df))
      }
    }
  })
  
  selected_snps[
                order(-cmh_p_val.p_adj_genome_wide.nlp, cmh_p_val), 
                snp_ranking := 1:.N, 
                by=.(chr_no, ensemble_gene_id)
              ]
  
  selected_snps <- selected_snps %>% 
                      group_by(chr_no, ensemble_gene_id) %>%
                      arrange(snp_ranking)
  
  #output and returning
  if (!is.null(out_file_path)) {
    write.table(x = selected_snps, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(selected_snps)
  
}






