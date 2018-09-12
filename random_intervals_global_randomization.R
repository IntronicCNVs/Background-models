#!/apps/R/3.4.0/bin/Rscript


# args          = commandArgs(trailingOnly=TRUE)
# input_dir     = args[1]
# output_dir    = args[2]
# map_num       = as.integer(args[3])
# total_permuts = as.integer(args[4])
# size_block    = as.integer(args[5])

input_dir     <- "/home/mrigau/_projects/A_FINAL_VERSION_CNVs/Revision_plos_genetics/extension_randomizations/global_randomization/"
output_dir    <- "/home/mrigau/_projects/A_FINAL_VERSION_CNVs/Revision_plos_genetics/extension_randomizations/global_randomization/output/"
total_permuts <- 10000

library( "regioneR" )
library( "BSgenome.Hsapiens.UCSC.hg19.masked" )
library( "gtools" )
library( "seqbias")


objects <- load("/home/mrigau/_projects/A_FINAL_VERSION_CNVs/Revision_plos_genetics/data_for_randomizations_revision2.RData")

# # Datasets
all_CNV_sets <- c("Phase3_DELS",
                  "Handsaker_DELS",
                  "Zarrei_DELS",
                  "Abyzov_RANGES",
                  "Sud15_DELS")




## Retrieving genome and mask
# Get genome and mask
human.genome    <- getGenomeAndMask(genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$genome
human.autosomal <- filterChromosomes(human.genome, organism="hg", chr.type="autosomal")

# Filter genome and mask, only autosomes
human.mask <- getGenomeAndMask(genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$mask
human.mask.autosomal <- filterChromosomes(human.mask, organism="hg", chr.type="autosomal")

names(human.autosomal) <- seqnames(human.autosomal)

### CNV datasets preparation
# Filtering CNV datasets
for(CNV_set in all_CNV_sets) {
  CNV_ranges <- get(CNV_set)
  # Pass all ranges to UCSC seqname style
  seqlevelsStyle(CNV_ranges)      <- "UCSC"
  # Only CNVs in autosomes
  CNV_ranges <- filterChromosomes(CNV_ranges, chr.type = "autosomal")
  CNV_ranges <- CNV_ranges[!overlapsAny(CNV_ranges, human.mask.autosomal)]
  assign(CNV_set, CNV_ranges)
}

# CREATE GRanges WITH AVAILABLE REGIONS  FOR RANDOMIZATION
# unmasked genome
masked_coverage <- coverage(append(human.autosomal, reduce(human.mask.autosomal)))

coverages_RANGES <- GRanges()
for(i in paste0("chr", 1:22)) {
  coverage_Rle <- masked_coverage[[i]]
  chr_RANGES <- GRanges(seqnames = as.character(i), ranges= ranges(coverage_Rle), coverage = coverage_Rle@values)
  coverages_RANGES <- append(coverages_RANGES, chr_RANGES)
}

# Coverage 2 when there is mask 
# Coverage 1 when no mask : available regions for  randomization
unmasked_genome <- coverages_RANGES[coverages_RANGES$coverage == 1]

unmasked_genome_for_randomization <- unmasked_genome
unmasked_genome_for_randomization$original_chr <- seqnames(unmasked_genome)

# Transform chromosomes from chr- to number
seqlevelsStyle(unmasked_genome_for_randomization) <- "NCBI"
chrs_in_unmasked_genome <- factor(as.character(seqnames(unmasked_genome_for_randomization)), 
                                  levels = 1:length(unmasked_genome_for_randomization))


seqlevels(unmasked_genome_for_randomization) <- as.character(1:length(unmasked_genome_for_randomization))
seqnames(unmasked_genome_for_randomization)  <- 1:length(unmasked_genome_for_randomization)
names(unmasked_genome_for_randomization) <- as.character(seqnames(unmasked_genome_for_randomization))


    
## Randomization function
global_randomization <- function(CNV_map, seed) {
  
  set.seed(seed)
  random_int0 <- random.intervals(unmasked_genome_for_randomization,
                                  n=length(CNV_map),
                                  ms=width(CNV_map)-1)
  
  random_int <- GRanges(seqnames = unmasked_genome_for_randomization[as.character(seqnames(random_int0))]$original_chr,
                        ranges = IRanges(start = start(random_int0) + start(unmasked_genome_for_randomization[seqnames(random_int0)]),
                                         end   = end(random_int0) + start(unmasked_genome_for_randomization[seqnames(random_int0)])))
  
  random_int
}

  
for(map_num in c(1,3,4,5,2)) {
  CNV_set <- all_CNV_sets[map_num]
  CNV_ranges <- get(CNV_set)
  
  all_random_starts <- c()
  all_random_chrs   <- c()
  
  for(i in 1:total_permuts) {
    random_ranges <- global_randomization(CNV_map = CNV_ranges, seed = map_num*1000000 + i)
    all_random_starts <- cbind(all_random_starts, start(random_ranges))
    all_random_chrs   <- cbind(all_random_chrs,   as.character(seqnames(random_ranges)))
    if(i%%500 == 0) print(i)
  }
  
  write.table(all_random_starts, quote = F, row.names = F, col.names = F,
       file = paste0(output_dir, "randomintervals/global_random_intervals_random_starts_", CNV_set, "_",
                     total_permuts, ".txt"))
  rm(all_random_starts)  

  write.table(all_random_chrs, quote = F, row.names = F, col.names = F,
              file = paste0(output_dir, "randomintervals/global_random_intervals_random_chrs_", CNV_set, "_",
                            total_permuts, ".txt"))
  rm(all_random_chrs)
  
  save(CNV_ranges, 
       file = paste0(output_dir, "randomintervals/global_random_intervals_original_", CNV_set, "_GRanges.RData"))
  
  
  print(CNV_set)
}

