context("Gene Statistics")

test_that("topQ for each gene", {
  flog.info(getwd())
  gene_df <- fread("../../mock_data/genes_nlp.tsv", sep="\t", sep2="auto", header=T, na.strings="NA", stringsAsFactors = FALSE)
  
  top1 <- get_topQ(df = gene_df, threshold = 1, column_name = "cmh_p_val.p_adj_genome_wide.nlp")
  top50 <- get_topQ(df = gene_df, threshold = 50, column_name = "cmh_p_val.p_adj_genome_wide.nlp")
  
  expect_true(length(unique(gene_df$ensemble_gene_id)) == length(unique(top50$ensemble_gene_id)))

  #cleanup after yourself
  rm(gene_df, top1, top50)
  
})


