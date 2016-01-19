# genewise

[![Build Status](https://travis-ci.org/nehiljain/genewise.svg?branch=master)](https://travis-ci.org/nehiljain/genewiseR)


# Setup

1. library(devtools)
2. install_github('nehiljain/genewiseR')
3. 

get_snp_ids - To get the snp ides of the snps found in the study. Using columns chr_no, snp_pos, ref_allele, in alt_allele 
generate_new_ids - Get new snp ids for snps not found in ref db.
p_adjustment_genomewide - genomewide multiple correction [fdr hard coded]
p_adjustment_chrwide - chromosomewide multiple correction any method (bon, fdr, etc)
p_adjustment_summary - summary plots of comparison between padjusted and raw values genomewide and chromosomeewide
get_significant_snps - filter significant snps
get_nlp - add column with negative log p value
get_max_and_mean - calculates snp count, max and mean on given column name and groups all the counts by chromosome
get_topX_sample - Get mean of nlp(negative log p-value) of snps in the top x quartile of each gene 
get_topQ - Get mean of nlp(negative log p-value) of snps in the top quartile of each gene
explore_topQ - explore topq for 1,5,10,20,25,50
snp_selection - snp-selection based on algorithm
map_snps_to_gene - It finds all the snps in genome that are in  gene +/- window_size
dir_rbind - Rowise combine all the files in a directory on a distributed cluster
dir_merge - Combines all the files in a directory using a Full Outer Join `merge(.., all=T)`
norm_var_names - Converts character vector to sanitised varirable names

