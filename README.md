# genewise


<SOme text to describe genewiseR>

# Setup

1. library(devtools)
2. install_github('nehiljain/genewiseR')




# List of functions exported by the package
1. get_snp_ids - To get the snp ides of the snps found in the study. Using columns chr_no, snp_pos, ref_allele, in alt_allele
2. generate_new_ids - Get new snp ids for snps not found in ref db.
3. p_adjustment_genomewide - genomewide multiple correction [fdr hard coded]
4. p_adjustment_chrwide - chromosomewide multiple correction any method (bon, fdr, etc)
5. p_adjustment_summary - summary plots of comparison between padjusted and raw values genomewide and chromosomeewide
6. get_significant_snps - filter significant snps
7. get_nlp - add column with negative log p value
8. get_max_and_mean - calculates snp count, max and mean on given column name and groups all the counts by chromosome
9. get_topX_sample - Get mean of nlp(negative log p-value) of snps in the top x quartile of each gene 
10. get_topQ - Get mean of nlp(negative log p-value) of snps in the top quartile of each gene
11. explore_topQ - explore topq for 1,5,10,20,25,50
12. snp_selection - snp-selection based on algorithm
13. map_snps_to_gene - It finds all the snps in genome that are in  gene +/- window_size
14. dir_rbind - Rowise combine all the files in a directory on a distributed cluster
15. dir_merge - Combines all the files in a directory using a Full Outer Join `merge(.., all=T)`
16. norm_var_names - Converts character vector to sanitised varirable names

## order of execution:

- 14 and 15? same function?
- 1
- 2
- 3
- 4 it is additional option? not always required?
- 5
- 13
- 6 is this function was used after Hein correction of topQ? is it only statistic?
- 7
- 8
- 9 where x is? 1,5,10,20,25,50 as default?
- 10 ................ 9 and 10 same function (its 9 for top25%)
- 12

# Example process based on indel dataset

combine_gwas_df <- dir_rbind("/Users/nehiljain/code/personal/genewiseR_data/raw_data/", 
header = F,col_names = c("chr_no","snp_pos","allele","p_value"))
 
ref_df <- read_tsv("~/code/personal/genewiseR_data/ref/indels.Bos_taurus.vcf", comment = "##", 
progress = T, trim_ws = T)

result_df <- get_snp_ids(combine_gwas_df, ref_df )
