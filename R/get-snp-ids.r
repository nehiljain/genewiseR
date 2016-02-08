

#' To get the Snp IDs of the snps found in the study. Using columns chr_no, snp_pos, ref_allele, in alt_allele 
#' using these conditions rows are joined between the two dataframes provided to this function
#' @param  data_path A csv file path (absolute) with data about snps and stats from different studies
#' @param  ref_df A dataframe object representing vcf file or txt file with annotation columns
#' @return returns a data table with all the ref. snp ids joined to each row
#' @examples 
#' 
#' combine_gwas_df <- dir_rbind("/Users/nehiljain/code/personal/genewiseR_data/raw_data/", 
#' header = F,col_names = c("chr_no","snp_pos","allele","p_value"))
#' 
#' ref_df <- read_tsv("~/code/personal/genewiseR_data/ref/indels.Bos_taurus.vcf", comment = "##", 
#' progress = T, trim_ws = T)
#' result_df <- get_snp_ids(combine_gwas_df, ref_df )
#' 

get_snp_ids <- function(study_df, 
                        ref_df, 
                        out_file_path=NULL) {
  
  if(!data.table::is.data.table(study_df)) {
    study_df <- data.table::data.table(study_df)
  }
  
  if(!data.table::is.data.table(ref_df)) {
    ref_df <- data.table::data.table(ref_df)
  }
  
  assert_that(data.table::is.data.table(study_df))
  assert_that(data.table::is.data.table(ref_df))

  
  #mapping the chromosome column to chr_no
  #NEHIL: THIS NEEDS TO BE CONSTANTLY UPDATED WITH NEW FILE FORMAT DISCOVERIES
  possbile_chr_names_list <- tolower(c("Chromosome", "Chromsome Name", "#CHROM", "Chr_no", "CHR"))
  
  chr_col_index <- fastmatch::fmatch(possbile_chr_names_list, tolower(names(ref_df)))
  chr_col_index <- na.omit(chr_col_index)[1]
  setnames(ref_df, chr_col_index, norm_var_names("chr_no"))
  
  chr_col_index <- fastmatch::fmatch(possbile_chr_names_list, tolower(names(study_df)))
  # print(names(study_df), chr_col_index)
  chr_col_index <- na.omit(chr_col_index)[1]
  print(chr_col_index)
  setnames(study_df, chr_col_index, norm_var_names("chr_no"))
  
  #mapping the snp position column to snp_pos
  possbile_pos_names_list <- tolower(c("POS", "Snp_pos"))
  
  chr_col_index <- fastmatch::fmatch(possbile_pos_names_list, tolower(names(ref_df)))
  chr_col_index <- na.omit(chr_col_index)[1]
  data.table::setnames(ref_df, chr_col_index, norm_var_names("snp_pos"))
  
  chr_col_index <- fastmatch::fmatch(possbile_pos_names_list, tolower(names(study_df)))
  chr_col_index <- na.omit(chr_col_index)[1]
  data.table::setnames(study_df, chr_col_index, norm_var_names("snp_pos"))
  
  #converting chr_no column to character
  study_df <- study_df[, chr_no:=as.character(chr_no)]
  ref_df <- ref_df[, chr_no:=as.character(chr_no)]
  
  study_df <- unique(study_df, by=c("chr_no", "snp_pos"))
  ref_df <- unique(ref_df, by=c("chr_no", "snp_pos"))
  
  data.table::setkey(study_df, chr_no, snp_pos)
  data.table::setkey(ref_df, chr_no, snp_pos)
  result_df <- merge(x = study_df, y = ref_df, all.x = T,
                     by = c("chr_no", "snp_pos"), suffixes=c(".study", ".ref"))
  data.table::setnames(result_df, names(result_df), norm_var_names(names(result_df)))
  
  if (!is.null(out_file_path)) {
    # assert_that(is.writeable(out_file_path))
    flog.debug(sprintf("Output File path %s", out_file_path))
    write_tsv(x = result_df, file=out_file_path)
  }
  return(result_df)
}

