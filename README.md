# Background-models

Random_intervals_global_randomization.R
This script generates random sets of CNVs by relocating them all over the genome except for low-mappability regions.
The random sets are stored as matrices with the chromosomes and the start coordinates (CNV size will be the same as in the original map). 

Random_intervals_local_randomization.R
This script generates random sets of CNVs by relocating them all within a 10Mb window, except for low-mappability regions.
The random sets are stored as matrices with the start coordinates (chromosomes and CNV sizes will be the same as in the original map). 

