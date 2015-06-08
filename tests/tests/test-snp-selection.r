context("Selecting Snps")

test_that("length of each ld _id", {
  
  snp_df <- fread("/home/tmp/nlp_all_snps_in_genes.tsv", sep="\t", sep2="auto", header=T, na.strings="NA", stringsAsFactors = FALSE, nrows = 100)
  snp_df <- snp_df[,.(chr_no, snp_pos, cmh_p_val, cmh_p_val.p_adj_genome_wide.nlp, ensemble_gene_id)]
  snp_df <- snp_df[chr_no == "1"]
  result <- snp_selection(snps_data = snp_df, ld_blocks_file_path = "/home/data/tmp/chr1i/blocks.blocks.det")
  
  ld_id_freq_df <- ddply(result, .(ensemble_gene_id), function(df) {
    ld_id_freq <- as.data.frame(table(df$ld_id))
    ld_id_freq
  })
  
  expect_true(all(ld_id_freq_df$Freq == 1))
  
  
  #cleanup after yourself
  rm(snp_df, result, ld_id_freq_df)
  
})

