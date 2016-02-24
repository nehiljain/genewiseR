#' It finds all the snps in genome that are in  gene +/- window_size
#' 
#' 
#' @param study_df : a datatable with snp information identified in the study. Should have column names - chr_no, snp_pos
#' @param ref_df : is the Bos taraus data which has gene_start, gene_end and chr_no column names
#' @param window_size : gene start and end +/- window_size
#' @return a dt with all snps in gene windows

map_snps_to_gene <- function(study_df, ref_df, window_size=1000, out_file_path=NULL) {
  
  if(!data.table::is.data.table(study_df)) {
    study_df <- data.table::data.table(study_df)
  }
  
  if(!data.table::is.data.table(ref_df)) {
    ref_df <- data.table::data.table(ref_df)
  }
  
  #mapping the chromosome column to chr_no
  #NEHIL: THIS NEEDS TO BE CONSTANTLY UPDATED WITH NEW FILE FORMAT DISCOVERIES
  possbile_chr_names_list <- tolower(c("Chromosome", "Chromsome Name", "#CHROM", "Chr_no", "CHR"))
  
  chr_col_index <- fastmatch::fmatch(possbile_chr_names_list, tolower(names(ref_df)))
  chr_col_index <- na.omit(chr_col_index)[1]
  setnames(ref_df, chr_col_index, norm_var_names("chr_no"))
  chr_list <- as.character(c(1:29,"X","Y"))
  setnames(ref_df, 
           names(ref_df), 
           norm_var_names(names(ref_df)))
  setnames(ref_gene_id_data, c("gene_start_(bp)", "gene_end_(bp)"), c("gene_start", "gene_end"))
  
  assert_that(all(c("snp_pos","chr_no") %in% names(study_df)))
  assert_that(all(c("gene_start", "gene_start", "chr_no") %in% names(ref_df)))  
  
  study_df[, fake_gene_start := snp_pos]
  study_df[, fake_gene_end := snp_pos]
  ref_df[, fake_gene_start := (gene_start - window_size)]
  ref_df[, fake_gene_end := (gene_end + window_size)]
  setkey(ref_df, chr_no, fake_gene_start, fake_gene_end)
  result_dt <- foverlaps(study_df, ref_df, type="within", nomatch = 0L)
  result_dt[,c("fake_gene_start","fake_gene_end","i.fake_gene_start","i.fake_gene_end"):= NULL]
  
  if (!is.null(out_file_path)) {
    write_tsv(x = result_dt, file=out_file_path)
  }
  return(result_dt)
} 




