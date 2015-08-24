context('Multiple corrections')

test_that('Genomewise correctin', {
  flog.info(getwd())
  adj_genome <- p_adjustment_genomewide(in_un_adj_p_val_snps_data_file_path = "mock_data/genome.gwas", col_names = "p_val")
  expect_false(is.null(adj_genome$p_val_adj_genome_wide))
  expect_true(is.numeric(adj_genome$p_val_adj_genome_wide))
  rm(adj_genome)
})

