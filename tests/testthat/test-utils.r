context('File Utils')

test_that('rbind files', {
  flog.info(getwd())
  rbind_data <- dir_rbind('mock_data/dir_rbind',  header=T)
  test2 <- read_tsv('mock_data/dir_rbind/test2.tsv')
  test3 <- read_tsv('mock_data/dir_rbind/test3.tsv')
  test4 <- read_tsv('mock_data/dir_rbind/test4.tsv')
  expect_true(nrow(rbind_data) == nrow(test2) + nrow(test3) + nrow(test4))
  expect_true(ncol(rbind_data) == ncol(test2))
  expect_true(ncol(rbind_data) == ncol(test3))
  rm(rbind_data, test2, test3, test4)
})


context('Utils')

test_that('normalise variable names', {
  test_str <- c("asdf Bads", "a-s.d{fa", "df$pct.spent")
  test_str <- norm_var_names(test_str)
  expect_equivalent(test_str, c("asdf_bads","a_s_d_fa","df_pct_spent"))
  rm(test_str)
})