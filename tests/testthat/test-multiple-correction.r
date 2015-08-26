context('Multiple corrections')

test_that('Genomewise correction', {
  flog.info(getwd())
  adj_genome <- p_adjustment_genomewide(in_un_adj_p_val_snps_data_file_path = "mock_data/genome.gwas", col_names = "p_val")
  expect_false(is.null(adj_genome$p_val_adj_genome_wide))
  expect_true(is.numeric(adj_genome$p_val_adj_genome_wide))
  rm(adj_genome)
})



test_that('Genomewise correction Report', {
  flog.info(getwd())
  adj_genome <- p_adjustment_genomewide(in_un_adj_p_val_snps_data_file_path = "mock_data/genome.gwas", col_names = "p_val")
  p_adjustment_summary(adj_genome, chr_name = "chr", snp_id_name = "snp_id", snp_pos_name="bp", p_val_name="p_val")
  rm(adj_genome)
})

