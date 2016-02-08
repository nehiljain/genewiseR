

#' Get new snp ids for snps not found in ref db.
#' 
#' 
#' The function takes in the data frame after attaching snp ids from the ref. 
#' using get_snp_ids(). The function computes the length of missing snp_name
#' @param the dataframe df
#' @param prefix, a chacracter value used as prefix to numerical index of the snps, default is 'dal_snp'
#' @param col_name, a character value used as column name of the snp ids
#' @return a datatable with all missing snp names changed to dal_snp1, dal_snp2,...

generate_new_snp_ids <- function(df, 
                             prefix_str = "dal_snp", 
                             col_name = "snp_id") {
  dt <- as.data.table(df)
  empty_columns_names <- dt[,!nzchar(get(col_name))]
  flog.debug(paste0("number of missing ids found", sum(empty_columns_names)))
  length_missing_ids <- sum(empty_columns_names)
  flog.debug(paste0("number of missing ids found", length_missing_ids))
  if (length_missing_ids > 0) {
    custom_names <- rep(prefix_str,length_missing_ids)
    custom_names <- paste0(custom_names, seq(1,length_missing_ids))
    flog.debug(paste0("Missing ids replaced with custom names", length_missing_ids, custom_names[0]))
    dt[empty_columns_names, (col_name) := custom_names]
  }
  setnames(dt, names(dt), norm_var_names(names(dt)))
  return(dt)
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



