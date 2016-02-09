

#' Get new snp ids for snps not found in ref db.
#' 
#' 
#' The function takes in the data frame after attaching snp ids from the ref. 
#' using get_snp_ids(). The function computes the length of missing snp_name. It assumes
#' @param study_dfthe study dataframe/data.table 
#' @param prefix_str a chacracter value used as prefix to numerical index of the snps, default is 'dal_snp'
#' 
#' @return a datatable with all missing snp ids changed to dal_snp1, dal_snp2,.... 
#' The column with the ids is snp_id
#' 
#' @examples 
#' 
#' combine_gwas_df <- dir_rbind("../genewiseR_data/raw_data/",
#' header = F,col_names = c("chr","pos","allele","p_value"))
#' 
#' ref_df <- read_tsv("../genewiseR_data/ref/indels.Bos_taurus.vcf", 
#' comment = "##", progress = T, trim_ws = T, col_types = "cicccnnc")
#' 
#' result_df <- get_snp_ids(combine_gwas_df, ref_df, out_file_path = "../genewiseR_data/tmp.tsv" )

generate_new_snp_ids <- function(study_df, 
                             prefix_str = "dal_snp") {
  if(!data.table::is.data.table(study_df)) {
    study_df <- data.table::data.table(study_df)
  }
  
  #mapping the snp position column to snp_pos
  possbile_id_names_list <- tolower(c("id", "Snp_id"))
  
  chr_col_index <- fastmatch::fmatch(possbile_id_names_list, tolower(names(study_df)))
  chr_col_index <- na.omit(chr_col_index)[1]
  data.table::setnames(study_df, chr_col_index, norm_var_names("snp_id"))
  
  study_df <- study_df[, snp_id:=as.character(snp_id)]
  empty_snp_id <- study_df[is.na(snp_id)]
  flog.debug(paste0("number of missing ids found = ", nrow(empty_snp_id)))
  length_missing_ids <- nrow(empty_snp_id)
  
  if (length_missing_ids > 0) {
    custom_names <- rep(prefix_str,length_missing_ids)
    custom_names <- paste0(custom_names, seq(1,length_missing_ids))
    flog.debug(paste0("Missing id replaced with custom names = ", length_missing_ids, " : ",head(custom_names,5)))
    flog.debug(paste0("\n...\n...\n...\n"))
    flog.debug(paste0("Missing id replaced with custom names = ", length_missing_ids, " : ",tail(custom_names,5)))
    study_df[empty_snp_id, snp_id := custom_names]
  }
  setnames(study_df, names(study_df), norm_var_names(names(study_df)))
  return(study_df)
}

