

#' To get the snp ides of the snps found in the study. Using columns chr_no, snp_pos, ref_allele, in alt_allele 
#' using these conditions rows are joined between the two dataframes provided to this function
#' @param  df1 A csv file (absolute)path  with data about snps and stats from different studies
#' @param  ref_df A csv refernce file (absolute)path 
#' @return returns a data table with all the ref. snp ids joined to each row

get_snp_ids <- function(df1, ref_df) {
  assert_that(is.data.table(df1))
  assert_that(is.data.table(ref_df))
  df1 <- unique(df1, by=c("chr_no", "pos", "ref"))
  ref_df <- unique(ref_df, by=c("chr_no", "pos", "ref"))

  setkey(df1, chr_no, pos)
  setkey(ref_df, chr_no, pos)
  
  result_df <- merge(x = df1, y = ref_df, all.x = T,
                     by = c("chr_no" , "pos" , "ref"), suffixes=c(".study", ".ref"))

  return(result_df)
}






#' The function takes in the data frame after attaching snp ids from the ref. 
#' using get_snp_ids(). The function computes the length of missing snp_name
#' @param the dataframe df
#' @return a dataframe with all missing snp names changed to pgi_dal_snp1,...

generate_new_ids <- function(df) {
  df <- as.data.table(df)
  length_missing_ids <- dim(df[is.na(alt.ref) & is.na(snp_name)])[1]
  custom_names <- rep("pgi_dal_snp", length_missing_ids)
  custom_names <- paste0(custom_names, seq(1,length_missing_ids))
  df[is.na(alt.ref) & is.na(snp_name), snp_name := custom_names]
  return(df)
}

#' execute_script() the driver function for the script
#' @param  csv_file_path A csv file (absolute)path  with data about snps and stats from 
#' different studies
#' @param  ref_tsv_file_path A tsv refernce file (absolute)path 
execute_script <- function(in_csv_file_path, 
                           in_ref_tsv_file_path, 
                           out_snp_name_annotated_study_snps_file_path,
                           in_ref_snps_dir,
                           out_combine_ref_snp_tsv_file_path) {
  study_data <- fread(in_csv_file_path, sep=",", sep2="auto", header=T, na.strings="NA",
                      stringsAsFactors = FALSE, verbose =T)
  ref_data <- fread(in_ref_tsv_file_path, sep="\t", header=T, na.strings="NA",
                    stringsAsFactors = FALSE, verbose =T)
  
  result_data <- get_snp_ids(study_data, ref_data)
  result_data <- unique(result_data)
  result_data <- generate_new_ids(result_data)
  write.table(x = result_data, file=out_snp_name_annotated_study_snps_file_path, quote = F, sep = "\t", row.names = F)
  
#   combined_ref_data <- combine_files_in_dir(in_ref_snps_dir, col_names = c("chr_no","pos","snp_name","ref","alt")  )
#   write.table(x = combined_ref_data, file=out_combine_ref_snp_tsv_file_path,quote = F, sep = "\t", row.names = F)
}


# 
# execute_script(in_csv_file_path = "/home/data/nehil_combine_data/combine_gwas_vcf.csv", 
#      in_ref_tsv_file_path = "/home/data/nehil_combine_all_chromosome_ref_snps.tsv",
#      out_snp_name_annotated_study_snps_file_path = "/home/data/nehil_snp_annotated_study_all_snp_ids.tsv",
#      in_ref_snps_dir = "/home/data/kacper_ref_snp_list/tsv",
#      out_combine_ref_snp_tsv_file_path = "/home/data/nehil_combine_all_chromosome_ref_snps.tsv"
#   )
