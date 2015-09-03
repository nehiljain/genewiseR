context('Snp ids')

test_that('get new snp ids', {
  flog.info(getwd())
  snp_data <- read_csv('mock_data/missing_snp_id.csv')
  snp_data <- generate_new_ids(df = snp_data, prefix = "test", col_name = "id")
  rm(snp_data)
})

