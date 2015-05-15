
library(data.table)
library(stringr)

#NOTE: NEHIL NEEDS TO CHANGE NORM_VAR_NAMES FUNCTION


execute <- function (in_genome_data_file_path, sep1, 
                     in_ref_gene_id_file_path, sep2, 
                     out_snps_in_genes_file_path, sep3, 
                     window_size = 1000) {
  chr_list <- as.character(c(1:29,"X"))

  genome_data <- fread(in_genome_data_file_path, sep=sep1, sep2="auto", header=T, na.strings="NA",
                       stringsAsFactors=FALSE, verbose=TRUE)  
#   setnames(genome_data,names(genome_data),norm_var_names(names(genome_data)))
  col_names <- c("chr_no","source","feature_type","gene_start_bp","gene_end_bp")
  genome_data[,chr_no := as.factor(chr_no)]
  
  ref_gene_id_data <- fread(in_ref_gene_id_file_path, sep=sep2, sep2="auto", na.strings="NA",
                          stringsAsFactors=FALSE, verbose=TRUE)
#   setnames(ref_gene_id_data,names(ref_gene_id_data), norm_var_names(names(ref_gene_id_data)))
  setnames(ref_gene_id_data, 
           names(ref_gene_id_data), 
           c("ensemble_gene_id","description", "chromosome_name", "gene_start_(bp)", "gene_end_(bp)", "strand","band"
             ,"associated_gene_name", "associated_gene_source","gene_biotype","source","status","version"))
#   print(names(ref_gene_id_data))
  ref_gene_id_data <- ref_gene_id_data[chromosome_name %in% chr_list] 
  ref_gene_id_data[,chromosome_name := as.factor(chromosome_name)]

  setnames(ref_gene_id_data, c("chromosome_name","gene_start_(bp)", "gene_end_(bp)"), c("chr_no","gene_start", "gene_end"))
  result_dt <- map_snps_to_gene(genome_data, ref_gene_id_data, window_size)
  write.table(result_dt, out_snps_in_genes_file_path, sep=sep3, row.names=F, quote = F)
}

#' The function gets genome data.table and reference data table
#' It finds all the snps in genome that are in  gene +/- window_size
#' @param genome_dt this is a datatable after bonferroni correction. Should have column names, chr_no, pos
#' @param ref_dt this is the Bos taraus dt which has gene_start, gene_end and chr_no column names
#' @param window_size
#' @return a dt with all snps in genes
map_snps_to_gene <- function(genome_dt, ref_dt, window_size) {
  genome_dt[, fake_gene_start := snp_pos]
  genome_dt[, fake_gene_end := snp_pos]
  ref_dt[, fake_gene_start := (gene_start - window_size)]
  ref_dt[, fake_gene_end := (gene_start + window_size)]
  setkey(ref_dt, chr_no, fake_gene_start, fake_gene_end)
  result_dt <- foverlaps(genome_dt, ref_dt, type="within", nomatch = 0L)
  return(result_dt)
} 


# 
# execute( in_genome_data_file_path = "/home/data/nehil_genome_p_adj_snp_annotated_study.tsv",
# in_ref_gene_id_file_path = "/home/data/reference/77/cattle_gene_list-UMD3.1.77.txt",
# out_snps_in_genes_file_path = "/home/data/cattle_gene_level.tsv")
