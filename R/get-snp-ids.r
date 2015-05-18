reset = T

if (reset) {
  rm(list = ls())  
}

library(data.table)
library(plyr)
library(dplyr)
library(stringr)


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
  print(sum(duplicated(df1[, c("chr_no", "pos", "ref"), with=FALSE])))
  print(sum(duplicated(ref_df[, c("chr_no", "pos", "ref"), with=FALSE])))

  setkey(df1, chr_no, pos)
  setkey(ref_df, chr_no, pos)
  
  result_df <- merge(x = df1, y = ref_df, all.x = T,
                     by = c("chr_no" , "pos" , "ref"), suffixes=c(".study", ".ref"))

  return(result_df)
}

#' This function normalises the input string vector
#' _ and small case output, no spaces, no . etc
#' @param a character vector of names
#' @return a normalised character vector of names
norm_var_names <- function(vars, sep="_") {
  if (sep == ".") sep <- "\\."
  
  # Replace all _ and . and ' ' with the nominated separator.
  
  pat  <- '_|\\.| |,'
  rep  <- sep
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # Replace any all capitals words with Initial capitals
  
  pat  <- stringr::regex('(?<!\\p{Lu})(\\p{Lu})(\\p{Lu}*)')
  rep  <- '\\1\\L\\2'
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # Replace any capitals not at the beginning of the string with _ 
  # and then the lowercase letter.
  
  pat  <- stringr::regex('(?<!^)(\\p{Lu})')
  rep  <- paste0(sep, '\\L\\1')
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # WHY DO THIS? Replace any number sequences not preceded by an
  # underscore, with it preceded by an underscore. The (?<!...) is a
  # lookbehind operator.
  
  pat  <- stringr::regex(paste0('(?<![', sep, '\\p{N}])(\\p{N}+)'))
  rep  <- paste0(sep, '\\1')
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # Remove any resulting initial or trailing underscore or multiples:
  #
  # _2level -> 2level
  
  vars <- stringr::str_replace(vars, "^_+", "")
  vars <- stringr::str_replace(vars, "_+$", "")
  vars <- stringr::str_replace(vars, "__+", "_")
  
  # Convert to lowercase
  
  vars <- tolower(vars)
  
  # Remove repeated separators.
  
  pat  <- paste0(sep, "+")
  rep  <- sep
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  return(vars)
}


#' This function combines the files in the geiven directory
#' The assumption is that it does not have a header (default)
#' 
combine_files_in_dir <- function(dir_path, header = F, col_names = NULL) {
  filename_list <- list.files(dir_path, full.names = T)
  
  cat("Reading Init",)
  combine_data <- fread(filename_list[[1]], sep="\t", header = header, na.strings="NA",
                             stringsAsFactors = FALSE, verbose =T)
  if (is.null(header) & is.null(col_names)) {
    warning("Header and Col Names are both NULL")
    col_names <- seq(1:dim(combine_data)[1], by=1)
  }
  for(l in filename_list) {
    cat("reading:",l)
    data <- fread(l, sep="\t", header = header, na.strings="NA",
                  stringsAsFactors = FALSE, verbose =T)
    combine_data <- rbind(combine_data, data)
  }
  cat("reading over",l)
  combine_data <- unique(combine_data)
  setnames(combine_data, names(combine_data), norm_var_names(col_names))
  return(combine_data)
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
