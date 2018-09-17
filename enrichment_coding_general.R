library( "regioneR" )
library( "BSgenome.Hsapiens.UCSC.hg19.masked" )
library( "gtools")

# #### Data loading and preparation
# Set working directory to Background-models
setwd("Background-models/")

if (!file.exists("enrichment_coding")) dir.create(file.path(".", "enrichment_coding"))
for(subfolder_name in c("global_random_intervals", "local_random_intervals")) {
  if (!file.exists(paste0("enrichment_coding/", subfolder_name))) dir.create(file.path("enrichment_coding", subfolder_name))
}



# #### Data loading and preparation
objects <- load("deletions_introns_genes_ages.RData")

all_coding_RANGES <- append(ppal_exons_RANGES, alt_exons_RANGES)
all_coding_RANGES <- all_coding_RANGES[seqnames(all_coding_RANGES) %in% 1:22]
all_coding_RANGES <- all_coding_RANGES[!duplicated(all_coding_RANGES)]

seqlevelsStyle(all_coding_RANGES) <- "UCSC"

source("calculate_pvalue_randomizations.R")

# Get genome and mask
human.genome    <- getGenomeAndMask(genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$genome
human.autosomal <- filterChromosomes(human.genome, organism="hg", chr.type="autosomal")

# Filter genome and mask, only autosomes
human.mask <- getGenomeAndMask(genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$mask
human.mask.autosomal <- filterChromosomes(human.mask, organism="hg", chr.type="autosomal")

# Mask with protein-coding genes (for randomizations excluding overlaps with genes)
# protein_coding_RANGES <- filterChromosomes(protein_coding_RANGES, organism="hg", chr.type="autosomal")
# masked.genes.and.gaps <- reduce(append(human.mask.autosomal, protein_coding_RANGES))

# Number of permutations
total_permuts <- 100

# CNV datasets analysed
all_CNV_sets <- c("Phase3_DELS",
                  "Handsaker_DELS",
                  "Zarrei_DELS",
                  "Abyzov_RANGES",
                  "Sud15_DELS")

# Filtering CNV datasets
for(CNV_set in all_CNV_sets) {
  CNV_ranges <- get(CNV_set)
  # Pass all ranges to UCSC seqname style
  seqlevelsStyle(CNV_ranges)      <- "UCSC"
  # I remove CNVs overlapping with masked regions:
  CNV_ranges <- CNV_ranges[!overlapsAny(CNV_ranges, human.mask.autosomal)]
  # Only CNVs in autosomes
  CNV_ranges <- filterChromosomes(CNV_ranges, chr.type = "autosomal")
  assign(CNV_set, CNV_ranges)
}


## GLOBAL - random intervals
### Randomizations done using random.intervals()

### ENRICHMENT CODING NUMBER OF CNVs
# Folder for results
folder <- "global_random_intervals"
for(CNV_set in all_CNV_sets) {
  ## RANDOM STARTS and CHROMOSOMES
  chromosomes <- read.table(paste0("global_random_intervals/output/global_random_intervals_random_chrs_", CNV_set, "_", total_permuts, ".txt"),
                            stringsAsFactors = F, header = F)
  print("chromosomes loaded")
  
  starts <- read.table(paste0("global_random_intervals/output/global_random_intervals_random_starts_", CNV_set, "_", total_permuts, ".txt"),
                       stringsAsFactors = F, header = F)
  print("starts loaded")
  
  # CNV_ranges
  load(paste0("global_random_intervals/output/global_random_intervals_original_", CNV_set, "_GRanges.RData"))
  
  print("CNV_ranges loaded")
  
  obs <- sum(overlapsAny(CNV_ranges, all_coding_RANGES, type = "any" ))
  
  # Random values
  ran <- c()
  
  while(ncol(chromosomes) > 0) {
    random_cnvs <- GRanges(seqnames = chromosomes[, 1],
                           ranges   = IRanges(start = starts[, 1],
                                              width = width(CNV_ranges)))
    
    ran <- c(ran, sum(sum(overlapsAny(random_cnvs, all_coding_RANGES, type = "any" ))))
    
    
    chromosomes <- chromosomes[, -1]
    starts      <- starts[, -1]
    
    if(class(starts) == "integer") {
      chromosomes <- matrix(chromosomes, ncol = 1)
      starts      <- matrix(starts, ncol = 1)
    }
    if(ncol(chromosomes)%%100 == 0) print(ncol(chromosomes))
    # print(ncol(chromosomes))
  }
  
  results <- calculate_pval_randomizations(obs, ran)
  
  save(results, file = paste0("enrichment_coding/", folder, "/general_enrichment_any", CNV_set, "_", total_permuts, "permuts.RData"))
  
  print(CNV_set)
  
}

print("End of global enrichment coding losses")


## LOCAL random intervals

# Folder for results
folder <- "local_random_intervals"
for(CNV_set in all_CNV_sets) {
  
  files <- list.files("local_random_intervals/output/")
  files <- files[grep(CNV_set, files)]
  files <- files[grep(paste0("_", total_permuts, "_permuts"), files)]
  ## RANDOM STARTS (for each fragment)
  all_starts <- c()
  for(segment_num in 1:278) {
    
    if(sum(grepl(paste0("_segment", segment_num, "_"), files)) == 1) {
      all_starts <-rbind(all_starts,
                         read.table(paste0("local_random_intervals/output/local_random_intervals_random_starts_", CNV_set, 
                                           "_subset_segment", segment_num, "_", 
                                           total_permuts, "_permuts.txt"),
                                    stringsAsFactors = F, header = F))
    }
    print(segment_num)
  }
  
  # CNV_ranges
  load(paste0("local_random_intervals/output/local_random_intervals_original_", CNV_set, "_GRanges.RData"))
  CNV_ranges <- whole_granges
  rm(whole_granges)
  
  obs <- sum(overlapsAny(CNV_ranges, all_coding_RANGES, type = "any" ))
  
  # Random values
  ran <- c()
  while(ncol(all_starts) > 0) {
    random_cnvs <- GRanges(seqnames = seqnames(CNV_ranges),
                           ranges   = IRanges(start = all_starts[, 1],
                                              width = width(CNV_ranges)))
    
    
    ran <- c(ran, sum(sum(overlapsAny(random_cnvs, all_coding_RANGES, type = "any" ))))
    
    all_starts      <- all_starts[, -1]
    
    if(class(all_starts) == "integer") {
      all_starts      <- matrix(all_starts, ncol = 1)
    }
    if(ncol(all_starts)%%500 == 0) print(ncol(all_starts))
    
  }
  
  results <- calculate_pval_randomizations(obs,ran)
  
  save(results, file = paste0("enrichment_coding/", folder, "/general_enrichment_any", CNV_set, "_", total_permuts, "permuts.RData"))
  
  print(CNV_set)
  
}

print("End of local enrichment coding losses")




