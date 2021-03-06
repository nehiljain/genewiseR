% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_new_snp_ids.R
\name{generate_new_snp_ids}
\alias{generate_new_snp_ids}
\title{Get new snp ids for snps not found in ref db.}
\usage{
generate_new_snp_ids(study_df, prefix_str = "dal_snp")
}
\arguments{
\item{prefix_str}{a chacracter value used as prefix to numerical index of the snps, default is 'dal_snp'}

\item{study_dfthe}{study dataframe/data.table}
}
\value{
a datatable with all missing snp ids changed to dal_snp1, dal_snp2,.... 
The column with the ids is snp_id
}
\description{
The function takes in the data frame after attaching snp ids from the ref. 
using get_snp_ids(). The function computes the length of missing snp_name. It assumes
}
\examples{

# combine_gwas_df <- dir_rbind("../genewiseR_data/raw_data/",
#header = F,col_names = c("chr","pos","allele","p_value"))

#ref_df <- read_tsv("../genewiseR_data/ref/indels.Bos_taurus.vcf", 
#comment = "##", progress = T, trim_ws = T, col_types = "cicccnnc")

#result_df <- get_snp_ids(combine_gwas_df, ref_df, out_file_path = "../genewiseR_data/tmp.tsv" )

#with_ids_dt <- generate_new_snp_ids(result_df)
}

