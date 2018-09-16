# Background-models

Random_intervals_global_randomization.R
This script generates random sets of CNVs by relocating them all over the genome except for low-mappability regions.
The random sets are stored as matrices with the chromosomes and the start coordinates (CNV size will be the same as in the original map). 

Random_intervals_local_randomization.R
This script generates random sets of CNVs by relocating them all within a 10Mb window, except for low-mappability regions.
The random sets are stored as matrices with the start coordinates (chromosomes and CNV sizes will be the same as in the original map). 

enrichment_general.R
Here the observed number of CNVs completely within introns (or all overlapping introns) is counted for the real CNV maps and the random maps. Then, the enrichment is calculated using the function in calculate_pvalue_randomizations.R 

enrichment_by_intron_size_or_age.R	
The number of CNVs overlapping introns is calculated by groups of introns classified by evolutionary ages or by sizes. A P-value is calculate for each subgroup of introns. 

deletions_introns_genes_ages.RData	
