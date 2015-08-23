context("Selecting Snps")

# test_that("length of each ld _id", {
#   snp_df <- fread("/home/tmp/nlp_all_snps_in_genes.tsv", sep="\t", sep2="auto", header=T, na.strings="NA", stringsAsFactors = FALSE, nrows = 1000)
#   snp_df <- snp_df[,.(chr_no, snp_pos, cmh_p_val, cmh_p_val.p_adj_genome_wide.nlp, ensemble_gene_id)]
#   snp_df <- snp_df[chr_no == "1"]
#   result <- snp_selection(snps_data = snp_df, chr_no_i = 1, ld_blocks_file_path = "/home/data/tmp/chr1i/blocks.blocks.det")
#   result[ , snp_pos_diff := snp_pos - shift(snp_pos, 1L, type="lead")]
#   ld_id_freq_df <- ddply(result, .(ensemble_gene_id), function(df) {
#     ld_id_freq <- as.data.frame(table(df$ld_id))
#     ld_id_freq
#   })
#   
#   expect_true(all(ld_id_freq_df$Freq == 1))
#   
#   distance_condition <- dlply(result, .(ensemble_gene_id), function(df) {
#     df <- as.data.table(df)
#     df[ , snp_pos_diff := snp_pos - shift(snp_pos, 1L, type="lead")]
#     bool_result <- na.omit(df[ , !(abs(snp_pos_diff) < 1000 & is.na(shift(ld_id, 1L, type="lead")))])
#     false_list <- which(bool_result == F)
#     difference_in_ld_if_condition_met <- l_ply(false_list, function(index) {
#       c((is.na(df[index, ld_id]) != is.na(df[index+1, ld_id])),
#         (is.na(df[index, ld_id]) != is.na(df[index-1, ld_id]))
#       )
#     })
#     if(all(difference_in_ld_if_condition_met)) {
#       TRUE
#     } else {
#       all(difference_in_ld_if_condition_met)
#     }
#   })
#   
#   expect_true(all(distance_condition))
#   
#          
#   #cleanup after yourself
#   rm(snp_df, result, ld_id_freq_df)
#   
# })

