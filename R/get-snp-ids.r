reset = T

if (reset) {
  rm(list = ls())  
}

library(data.table)
library(plyr)
library(dplyr)
library(stringr)


#' main() the driver function for the script
#' @param  csv_file_path A csv file (absolute)path  with data about snps and stats from different studies
#' @param  ref_csv_file_path A csv refernce file (absolute)path 
main <- function(csv_file_path, ref_csv_file_path) {
  study_data <- fread("/home/data/nehil_combine_data/combine_gwas_vcf.csv", sep=",", sep2="auto", header=T, na.strings="NA",
                         stringsAsFactors = FALSE, verbose =T)
  ref_data <- fread("/home/data/nehil_combine_all_chromosome_ref_snps.tsv", sep="\t", header=T, na.strings="NA",
                        stringsAsFactors = FALSE, verbose =T)
  result_data <- get_snp_ids(study_data, ref_data)
  result_data <- unique(result_data)
  
  write.table(x = result_data, file="/home/data/nehil_combine_study_snp_ids.tsv",quote = F, sep = "\t", row.names = F)
  
  combined_ref_data <- combine_files_in_dir("/home/data/kacper_ref_snp_list/tsv", col_names = c("chr_no","pos","snp_name","ref","alt")  )
  write.table(x = combined_ref_data, file="/home/data/nehil_combine_all_chromosome_ref_snps.tsv",quote = F, sep = "\t", row.names = F)
  
}

#' To get the snp ides of the snps found in the study. Using columns chr_no, snp_pos, ref_allele, in alt_allele 
#' using these conditions rows are joined between the two dataframes provided to this function
#' @param  df1 A csv file (absolute)path  with data about snps and stats from different studies
#' @param  ref_df A csv refernce file (absolute)path 

get_snp_ids <- function(df1, ref_df) {
  str(df1)
  str(ref_df)
#   setnames(df1 , "alt", "study_alt")
  join1_data <- df1[1:dim(df1)[1]/2] %>%
    inner_join(ref_df, by = c("chr_no" = "chr_no", "pos" = "pos", "ref" = "ref"))
  join2_data <- df1[(dim(df1)[1]/2) +1 : dim(df1)[1]] %>%
    inner_join(ref_df, by = c("chr_no" = "chr_no", "pos" = "pos", "ref" = "ref"))
  return(rbind(join1_data,join2_data))
}

#' This function normalises the input string vector
#' _ and small case output, no spaces, no . etc
norm_var_names <- function(vars, sep="_") {
  if (sep == ".") sep <- "\\."
  
  # Replace all _ and . and ' ' with the nominated separator.
  
  pat  <- '_|\\.| |,'
  rep  <- sep
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # Replace any all capitals words with Initial capitals
  
  pat  <- stringr::perl('(?<!\\p{Lu})(\\p{Lu})(\\p{Lu}*)')
  rep  <- '\\1\\L\\2'
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # Replace any capitals not at the beginning of the string with _ 
  # and then the lowercase letter.
  
  pat  <- stringr::perl('(?<!^)(\\p{Lu})')
  rep  <- paste0(sep, '\\L\\1')
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  # WHY DO THIS? Replace any number sequences not preceded by an
  # underscore, with it preceded by an underscore. The (?<!...) is a
  # lookbehind operator.
  
  pat  <- stringr::perl(paste0('(?<![', sep, '\\p{N}])(\\p{N}+)'))
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
combine_files_in_dir <- function(dir_path, header = F, col_names = NULL) {
  filename_list <- list.files(dir_path, full.names = T)
  

  combine_data <- fread(filename_list[[1]], sep="\t", header = header, na.strings="NA",
                             stringsAsFactors = FALSE, verbose =T)
  if (is.null(header) & is.null(col_names)) {
    warning("Header and Col Names are both NULL")
    col_names <- seq(1:dim(combine_data)[1], by=1)
  }
  for(l in filename_list) {
    data <- fread(l, sep="\t", header = header, na.strings="NA",
                  stringsAsFactors = FALSE, verbose =T)
    combine_data <- rbind(combine_data, data)
  }
  combine_data <- unique(combine_data)
  setnames(combine_data, names(combine_data), norm_var_names(col_names))
  return(combine_data)
}



