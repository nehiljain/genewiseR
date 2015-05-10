
library(data.table)

find_snps_in_genes_files <- function (in_genome_data_file_path, in_ref_gene_id_file_path, out_snps_in_genes_file_path) {
  genome_data <- fread(in_genome_data_file_path, sep="\t", sep2="auto", header=T, na.strings="NA",
                       stringsAsFactors=FALSE, verbose=TRUE)
  ref_gene_id_data <- fread(in_ref_gene_id_file_path, sep="\t", sep2="auto", na.strings="NA",
                            stringsAsFactors=FALSE, verbose=TRUE)
  setnames(genome_data,names(genome_data),norm_var_names(names(genome_data)))
  col_names <- c("chr_no","source","feature_type","gene_start_bp","gene_end_bp","")
  setnames(ref_gene_id_data,names(ref_gene_id_data), norm_var_names(names(ref_gene_id_data)))
}

#' The function gets genome data table which has columns pos, chr_no
#' ref_dt has gene_start, gene_end, chr_no
#' @param genome_dt this is a datatable after bonferroni correction. Should have column names, chr_no, pos
#' @param ref_dt this is the Bos taraus dt which has gene_start, gene_end and chr_no column names
#' @return a dt with all snps in genes
find_snps_in_genes <- function(genome_dt, ref_dt, window_size) {
  genome_dt$fake_gene_start <- genome_dt$pos
  genome_dt$fake_gene_end <- genome_dt$pos
  ref_gene_id_data$fake_gene_start <- ref_gene_id_data$gene_start - window_size
  ref_gene_id_data$fake_gene_end <- ref_gene_id_data$gene_end + window_size
  setkey(ref_gene_id_data, chr_no, fake_gene_start, fake_gene_end)
  result_dt <- foverlaps(genome_data, ref_gene_id_data, type="within", nomatch = 0L)
  return(result_dt)
} 






in_genome_data_file_path = "/home/data/nehil_genome_p_adj_snp_annotated_study.tsv"
in_ref_gene_id_file_path = "/home/data/reference/77/Bos_taurus.UMD3.1.77.gtf"
out_snps_in_genes_file_path = "/home/data/cattle_gene_level.tsv"
