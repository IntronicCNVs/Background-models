

# Background-models

This repository contains the scripts and data that were used to generate the results presented in Rigau et al. (Intronic CNVs cause gene expression variation in human population, https://doi.org/10.1101/171165). 

We include code to generate different types of randomised sets of CNVs, that will be used in subsequent scripts to assess the relationship between CNVs and different genomic features, focusing on introns. These are the descriptions of the different R scripst and data files: 

Random_intervals_global_randomization.R

This script generates random sets of CNVs by relocating them all over the genome except for low-mappability regions.
The random sets are stored as matrices with the chromosomes and the start coordinates (CNV size will be the same as in the original map). 


Random_intervals_local_randomization.R

This script generates random sets of CNVs by relocating them all within a 10Mb window, except for low-mappability regions.
The random sets are stored as matrices with the start coordinates (chromosomes and CNV sizes will be the same as in the original map). 


enrichment_general.R

Here the observed number of intronic CNVs is counted in the real CNV maps and in the randomized maps. Then, the enrichment is calculated using the function in calculate_pvalue_randomizations.R 


enrichment_by_intron_size_or_age.R	

The number of CNVs overlapping introns is calculated by groups of introns classified by: (1) evolutionary ages or (2 sizes. A P-value is calculate for each subgroup of introns. 


enrichment_coding_general.R

Here the observed number of CNVs that overlaps with exons is counted in the real CNV maps and in the randomized maps. Then, the enrichment is calculated using the function in calculate_pvalue_randomizations.R 


enrichment_coding_by_age.R

Here the observed number of CNVs that overlaps with exons (and the number of genes affected by exon-overlapping CNVs) is counted in the real CNV maps and in the randomized maps, separating by gene evolutionary ages. Then, the enrichment is calculated using the function in calculate_pvalue_randomizations.R 

calculate_pvalue_randomizations.R 
Function used to calculate the P-value of the enrichment/depletion. This function is extracted from the function permTest() from R package RegioneR (version 1.6.2)

deletions_introns_genes_ages.RData	

RData object containing coordinates of deletions, intronic ranges, exons ranges and gene evolutionary ages, necessary for running all analyses.

Table_S7_list_of_CNVs.xls 

Supplementary table S7 from the article, containing all CNVs used taken into account in the study and their impact on genes. 
