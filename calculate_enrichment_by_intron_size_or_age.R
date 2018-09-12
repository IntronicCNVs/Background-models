# THIS SCRIPT IS THE NEW VERSION OF "INFO PER INTRON CORRECTED.R"
library( "regioneR" )
library( "BSgenome.Hsapiens.UCSC.hg19.masked" )
library( "seqbias" )
library( "gtools")
library( "xlsx")

# #### Data loading and preparation
setwd( "/home/mrigau/_projects/A_FINAL_VERSION_CNVs/Revision_plos_genetics/extension_randomizations/" )
objects <- load("../data_for_randomizations_revision2.RData")
seqlevelsStyle(intronic_RANGES) <- "UCSC"
seqlevelsStyle(protein_coding_RANGES)      <- "UCSC"
source("../calculate_pvalue_randomizations.R")

### Number of permutations

total_permuts <- 10000

# Datasets
all_CNV_sets <- c("Phase3_DELS",
                  "Handsaker_DELS",
                  "Zarrei_DELS",
                  "Abyzov_RANGES",
                  "Sud15_DELS")

### CNV datasets preparation
# Filtering CNV datasets
for(CNV_set in all_CNV_sets) {
  CNV_ranges <- get(CNV_set)
  # Pass all ranges to UCSC seqname style
  seqlevelsStyle(CNV_ranges)      <- "UCSC"
  # Only CNVs in autosomes
  CNV_ranges <- filterChromosomes(CNV_ranges, chr.type = "autosomal")
  assign(CNV_set, CNV_ranges)
}


# ```
#
# ***********
#   # NEW RANDOMIZATIONS
#   ```{r message = FALSE, warning= F, echo = T, eval = T, fig.width = 10}
## Retrieving genome and mask
# Get genome and mask
human.genome    <- getGenomeAndMask(genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$genome
human.autosomal <- filterChromosomes(human.genome, organism="hg", chr.type="autosomal")

# Filter genome and mask, only autosomes
human.mask <- getGenomeAndMask(genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$mask
human.mask.autosomal <- filterChromosomes(human.mask, organism="hg", chr.type="autosomal")

# Mask with protein-coding genes (for randomizations excluding overlaps with genes)
protein_coding_RANGES <- filterChromosomes(protein_coding_RANGES, organism="hg", chr.type="autosomal")
masked.genes.and.gaps <- reduce(append(human.mask.autosomal, protein_coding_RANGES))

# Filtering CNV datasets(I don't want any variant overlapping low mappability regions, so I remove further  analysis)
for(CNV_set in all_CNV_sets) {
  CNV_ranges <- get(CNV_set)
  # I remove CNVs overlapping with masked regions:
  length_before <- length(CNV_ranges)

  CNV_ranges <- CNV_ranges[!overlapsAny(CNV_ranges, human.mask.autosomal)]
  print(c(CNV_set, length_before - length(CNV_ranges)))
  assign(CNV_set, CNV_ranges)
}

ages <- levels(intronic_RANGES$most_recent_age)
intronic_RANGES$size_quant20 <- quantcut(width(intronic_RANGES), 20)
protein_coding_RANGES$Age <- Ages[names(protein_coding_RANGES), "GeneAge"]
protein_coding_RANGES <- protein_coding_RANGES[!is.na(protein_coding_RANGES$Age)]

#######################################################################

## Function to calculate overlaps    
by_age_and_size <- function(granges, overlap_type) {
  num_cnvs_by_age     <- c()
  num_introns_by_age  <- c()
  num_genes_by_age    <- c()
  num_introns_by_size <- c()
  num_cnvs_by_size    <- c()
  num_genes_by_size   <- c()
  
  # Observed number of deletions per intronic range
  for(age in ages) {
    
    fo <- findOverlaps(granges, intronic_RANGES[intronic_RANGES$most_recent_age == age],
                       type = overlap_type)
    
    num_cnvs_by_age    <- c(num_cnvs_by_age, length(fo))
    num_introns_by_age <- c(num_introns_by_age, length(unique(subjectHits(fo))))
    num_genes_by_age   <- c(num_genes_by_age, sum(overlapsAny(protein_coding_RANGES[protein_coding_RANGES$Age == age],
                                                              granges[queryHits(fo)])))
  }
  
  for(intron_size in levels(intronic_RANGES$size_quant20)) {
    fo <- findOverlaps(granges, intronic_RANGES[intronic_RANGES$size_quant20 == intron_size],
                       type = overlap_type)
    
    
    num_cnvs_by_size <- c(num_cnvs_by_size, length(fo))
    num_introns_by_size <- c(num_introns_by_size, length(unique(subjectHits(fo))))
    num_genes_by_size <- c(num_genes_by_size,  sum(overlapsAny(protein_coding_RANGES[protein_coding_RANGES$Age == age],
                                                               granges[queryHits(fo)])))
  }
  list(num_cnvs_by_age     = num_cnvs_by_age,
       num_introns_by_age  = num_introns_by_age,
       num_genes_by_age    = num_genes_by_age,
       num_introns_by_size = num_introns_by_size,
       num_cnvs_by_size    = num_cnvs_by_size,
       num_genes_by_size   = num_genes_by_size)
} 
   

#### RANDOM.INTERVALS()
### Randomizations done using random.intervals()
## GLOBAL - random intervals

# Folder for results
folder <- "global_random_intervals"
for(CNV_set in all_CNV_sets[4]) {
  ## RANDOM STARTS
  chromosomes <- read.table(paste0("global_randomization/output/randomintervals/global_random_intervals_random_chrs_", CNV_set, "_", total_permuts, ".txt"),
                            stringsAsFactors = F, header = F)
  print("chromosomes loaded")
  
  starts <- read.table(paste0("global_randomization/output/randomintervals/global_random_intervals_random_starts_", CNV_set, "_", total_permuts, ".txt"),
                       stringsAsFactors = F, header = F)
  print("starts loaded")
  
  # CNV_ranges
  load(paste0("global_randomization/output/randomintervals/global_random_intervals_original_", CNV_set, "_GRanges.RData"))
  
  print("CNV_ranges loaded")
  
  obs_within <-  by_age_and_size(CNV_ranges, overlap_type  = "within")
  obs_any    <-  by_age_and_size(CNV_ranges, overlap_type = "any")
  
  
  # Random values
  ran_within <- c()
  ran_any    <- c()
  while(ncol(chromosomes) > 0) {
    random_cnvs <- GRanges(seqnames = chromosomes[, 1],
                           ranges   = IRanges(start = starts[, 1],
                                              width = width(CNV_ranges)))
    
    x_within <-  by_age_and_size(random_cnvs, overlap_type  = "within")
    x_any    <-  by_age_and_size(random_cnvs, overlap_type = "any")
    
    for(i in 1:length(x_within)) {
      ran_within[[names(x_within)[i]]] <- rbind(ran_within[[names(x_within)[i]]],
                                                x_within[[names(x_within)[i]]])
      ran_any[[names(x_any)[i]]]       <- rbind(ran_any[[names(x_any)[i]]],
                                                x_any[[names(x_any)[i]]])
    }
    
    chromosomes <- chromosomes[, -1]
    starts      <- starts[, -1]
    
    if(class(starts) == "integer") {
      chromosomes <- matrix(chromosomes, ncol = 1)
      starts      <- matrix(starts, ncol = 1)
    }
    if(ncol(chromosomes)%%1000 == 0) print(ncol(chromosomes))
    # print(ncol(chromosomes))
  }
  
  for(i in 1:3) {
    colnames(ran_within[[i]]) <- ages
    colnames(ran_any[[i]])    <- ages
  }
  for(i in 4:6) {
    colnames(ran_within[[i]]) <- levels(intronic_RANGES$size_quant20)
    colnames(ran_any[[i]])    <- levels(intronic_RANGES$size_quant20)
    
  }
  
  for(overlap_type in c("within", "any")) {
    obs_num_cnvs <- c(get(paste0("obs_", overlap_type))$num_cnvs_by_age,     get(paste0("obs_", overlap_type))$num_cnvs_by_size  )
    obs_num_introns <- c(get(paste0("obs_", overlap_type))$num_introns_by_age,     get(paste0("obs_", overlap_type))$num_introns_by_size)
    obs_num_genes <- c(get(paste0("obs_", overlap_type))$num_genes_by_age,     get(paste0("obs_", overlap_type))$num_genes_by_size)
    
    names(obs_num_cnvs)    <- c(ages, levels(intronic_RANGES$size_quant20))
    names(obs_num_introns) <- c(ages, levels(intronic_RANGES$size_quant20))
    names(obs_num_genes)   <- c(ages, levels(intronic_RANGES$size_quant20))
    
    all_ran_num_cnvs    <- cbind(get(paste0("ran_", overlap_type))$num_cnvs_by_age,     get(paste0("ran_", overlap_type))$num_cnvs_by_size)
    all_ran_num_introns <- cbind(get(paste0("ran_", overlap_type))$num_introns_by_age,  get(paste0("ran_", overlap_type))$num_introns_by_size)
    all_ran_num_genes   <- cbind(get(paste0("ran_", overlap_type))$num_genes_by_age,    get(paste0("ran_", overlap_type))$num_genes_by_size)
    
    results <- list()
    for(x in names(obs_num_cnvs)) {
      results[[x]]$num_cnvs    <- calculate_pval_randomizations(obs_num_cnvs[x],
                                                                unname(all_ran_num_cnvs[, x]))
      results[[x]]$num_introns <- calculate_pval_randomizations(obs_num_introns[x],
                                                                unname(all_ran_num_introns[, x]))
      results[[x]]$num_genes   <- calculate_pval_randomizations(obs_num_genes[x],
                                                                unname(all_ran_num_genes[, x]))
    }
    
    results_ages <- results[ages]
    save(results_ages, file = paste0("enrichment_by_age/", folder, "/", overlap_type, "_", CNV_set, "_", total_permuts, "permuts.RData"))
    results_sizes <- results[levels(intronic_RANGES$size_quant20)]
    save(results_sizes, file = paste0("enrichment_by_size/", folder, "/",  overlap_type, "_", CNV_set, "_", total_permuts, "permuts.RData"))
  }
  
print(CNV_set)
  
}

print("End of global enrichment by age and size")

## LOCAL random intervals

# Folder for results
folder <- "local_random_intervals"
for(CNV_set in all_CNV_sets[4]) {
  
  files <- list.files("local_randomization/output/randomintervals/")
  files <- files[grep(CNV_set, files)]
  files <- files[grep(paste0("_", total_permuts, "_permuts"), files)]
  ## RANDOM STARTS (for each fragment)
  all_starts <- c()
  for(segment_num in 1:278) {
    
    if(sum(grepl(paste0("_segment", segment_num, "_"), files)) == 1) {
      all_starts <-rbind(all_starts,
                         read.table(paste0("local_randomization/output/randomintervals/local_random_intervals_random_starts_", CNV_set, 
                                           "_subset_segment", segment_num, "_", 
                                           total_permuts, "_permuts.txt"),
                                    stringsAsFactors = F, header = F))
    }
    print(segment_num)
  }
  
  # CNV_ranges
  load(paste0("local_randomization/output/randomintervals/local_random_intervals_original_", CNV_set, "_GRanges.RData"))
  CNV_ranges <- whole_granges
  rm(whole_granges)
  
  obs_within <-  by_age_and_size(CNV_ranges, overlap_type  = "within")
  obs_any    <-  by_age_and_size(CNV_ranges, overlap_type = "any")
  
    ran_within <- c()
    ran_any    <- c()
    while(ncol(all_starts) > 0) {
      random_cnvs <- GRanges(seqnames = seqnames(CNV_ranges),
                             ranges   = IRanges(start = all_starts[, 1],
                                                width = width(CNV_ranges)))
      
      x_within <-  by_age_and_size(random_cnvs, overlap_type  = "within")
      x_any    <-  by_age_and_size(random_cnvs, overlap_type = "any")
      
      for(i in 1:length(x_within)) {
        ran_within[[names(x_within)[i]]] <- rbind(ran_within[[names(x_within)[i]]],
                                                  x_within[[names(x_within)[i]]])
        ran_any[[names(x_any)[i]]]       <- rbind(ran_any[[names(x_any)[i]]],
                                                  x_any[[names(x_any)[i]]])
      }
      
      
      all_starts      <- all_starts[, -1]
      
      if(class(all_starts) == "integer") {
        all_starts      <- matrix(all_starts, ncol = 1)
      }
      if(ncol(all_starts)%%500 == 0) print(ncol(all_starts))
      # print(ncol(chromosomes))
    }
    
    for(i in 1:3) {
      colnames(ran_within[[i]]) <- ages
      colnames(ran_any[[i]])    <- ages
    }
    for(i in 4:6) {
      colnames(ran_within[[i]]) <- levels(intronic_RANGES$size_quant20)
      colnames(ran_any[[i]])    <- levels(intronic_RANGES$size_quant20)
    }
    
    for(overlap_type in c("within", "any")) {
      obs_num_cnvs <- c(get(paste0("obs_", overlap_type))$num_cnvs_by_age,     get(paste0("obs_", overlap_type))$num_cnvs_by_size  )
      obs_num_introns <- c(get(paste0("obs_", overlap_type))$num_introns_by_age,     get(paste0("obs_", overlap_type))$num_introns_by_size)
      obs_num_genes <- c(get(paste0("obs_", overlap_type))$num_genes_by_age,     get(paste0("obs_", overlap_type))$num_genes_by_size)
      
      names(obs_num_cnvs)    <- c(ages, levels(intronic_RANGES$size_quant20))
      names(obs_num_introns) <- c(ages, levels(intronic_RANGES$size_quant20))
      names(obs_num_genes)   <- c(ages, levels(intronic_RANGES$size_quant20))
      
      all_ran_num_cnvs    <- cbind(get(paste0("ran_", overlap_type))$num_cnvs_by_age,     get(paste0("ran_", overlap_type))$num_cnvs_by_size)
      all_ran_num_introns <- cbind(get(paste0("ran_", overlap_type))$num_introns_by_age,  get(paste0("ran_", overlap_type))$num_introns_by_size)
      all_ran_num_genes   <- cbind(get(paste0("ran_", overlap_type))$num_genes_by_age,    get(paste0("ran_", overlap_type))$num_genes_by_size)
      
      results <- list()
      for(x in names(obs_num_cnvs)) {
        results[[x]]$num_cnvs    <- calculate_pval_randomizations(obs_num_cnvs[x],
                                                                  unname(all_ran_num_cnvs[, x]))
        results[[x]]$num_introns <- calculate_pval_randomizations(obs_num_introns[x],
                                                                  unname(all_ran_num_introns[, x]))
        results[[x]]$num_genes   <- calculate_pval_randomizations(obs_num_genes[x],
                                                                  unname(all_ran_num_genes[, x]))
      }
      
      results_ages <- results[ages]
      save(results_ages, file = paste0("enrichment_by_age/", folder, "/", overlap_type, "_", CNV_set, "_", total_permuts, "permuts.RData"))
      results_sizes <- results[levels(intronic_RANGES$size_quant20)]
      save(results_sizes, file = paste0("enrichment_by_size/", folder, "/",  overlap_type, "_", CNV_set, "_", total_permuts, "permuts.RData"))
    }
    
    print(CNV_set)
    
    
}


###### GLOBAL SPLITTING CNVS
## GLOBAL - random intervals

# Folder for results
folder <- "global_split_cnvs"
for(CNV_set in all_CNV_sets[c(1,3,4,5,2)]) {
  ## RANDOM STARTS
  chromosomes <- read.table(paste0("global_randomization/output/cutting_CNVs/global_random_intervals_random_chrs_", CNV_set, "_", total_permuts, ".txt"),
                            stringsAsFactors = F, header = F)#[, 1:1000]
  print("chromosomes loaded")
  
  starts <- read.table(paste0("global_randomization/output/cutting_CNVs/global_random_intervals_random_starts_", CNV_set, "_", total_permuts, ".txt"),
                       stringsAsFactors = F, header = F)#[, 1:1000]
  print("starts loaded")
  
  widths <- read.table(paste0("global_randomization/output/cutting_CNVs/global_random_intervals_random_widths_", CNV_set, "_", total_permuts, ".txt"),
                       stringsAsFactors = F, header = F)#[, 1:1000]
  print("widths loaded")
  
  # CNV_ranges
  load(paste0("global_randomization/output/cutting_CNVs/global_random_intervals_original_", CNV_set, "_GRanges.RData"))
  
  print("CNV_ranges loaded")
  
  obs_within <-  by_age_and_size(CNV_ranges, overlap_type  = "within")
  obs_any    <-  by_age_and_size(CNV_ranges, overlap_type = "any")
  
  
  # Random values
  ran_within <- c()
  ran_any    <- c()
  while(ncol(chromosomes) > 0) {
    random_cnvs <- GRanges(seqnames = chromosomes[, 1],
                           ranges   = IRanges(start = starts[, 1],
                                              width = widths[, 1]))
    
    x_within <-  by_age_and_size(random_cnvs, overlap_type  = "within")
    x_any    <-  by_age_and_size(random_cnvs, overlap_type = "any")
    
    for(i in 1:length(x_within)) {
      ran_within[[names(x_within)[i]]] <- rbind(ran_within[[names(x_within)[i]]],
                                                x_within[[names(x_within)[i]]])
      ran_any[[names(x_any)[i]]]       <- rbind(ran_any[[names(x_any)[i]]],
                                                x_any[[names(x_any)[i]]])
    }
    
    chromosomes <- chromosomes[, -1]
    starts      <- starts[, -1]
    widths      <- widths[, -1]
    
    if(class(starts) == "integer") {
      chromosomes <- matrix(chromosomes, ncol = 1)
      starts      <- matrix(starts, ncol = 1)
      widths      <- matrix(widths, ncol = 1)
      
    }
    if(ncol(chromosomes)%%1000 == 0) print(ncol(chromosomes))
    # print(ncol(chromosomes))
  }
  
  for(i in 1:3) {
    colnames(ran_within[[i]]) <- ages
    colnames(ran_any[[i]])    <- ages
  }
  for(i in 4:6) {
    colnames(ran_within[[i]]) <- levels(intronic_RANGES$size_quant20)
    colnames(ran_any[[i]])    <- levels(intronic_RANGES$size_quant20)
    
  }
  
  for(overlap_type in c("within", "any")) {
    obs_num_cnvs <- c(get(paste0("obs_", overlap_type))$num_cnvs_by_age,     get(paste0("obs_", overlap_type))$num_cnvs_by_size  )
    obs_num_introns <- c(get(paste0("obs_", overlap_type))$num_introns_by_age,     get(paste0("obs_", overlap_type))$num_introns_by_size)
    obs_num_genes <- c(get(paste0("obs_", overlap_type))$num_genes_by_age,     get(paste0("obs_", overlap_type))$num_genes_by_size)
    
    names(obs_num_cnvs)    <- c(ages, levels(intronic_RANGES$size_quant20))
    names(obs_num_introns) <- c(ages, levels(intronic_RANGES$size_quant20))
    names(obs_num_genes)   <- c(ages, levels(intronic_RANGES$size_quant20))
    
    all_ran_num_cnvs    <- cbind(get(paste0("ran_", overlap_type))$num_cnvs_by_age,     get(paste0("ran_", overlap_type))$num_cnvs_by_size)
    all_ran_num_introns <- cbind(get(paste0("ran_", overlap_type))$num_introns_by_age,  get(paste0("ran_", overlap_type))$num_introns_by_size)
    all_ran_num_genes   <- cbind(get(paste0("ran_", overlap_type))$num_genes_by_age,    get(paste0("ran_", overlap_type))$num_genes_by_size)
    
    results <- list()
    for(x in names(obs_num_cnvs)) {
      results[[x]]$num_cnvs    <- calculate_pval_randomizations(obs_num_cnvs[x],
                                                                unname(all_ran_num_cnvs[, x]))
      results[[x]]$num_introns <- calculate_pval_randomizations(obs_num_introns[x],
                                                                unname(all_ran_num_introns[, x]))
      results[[x]]$num_genes   <- calculate_pval_randomizations(obs_num_genes[x],
                                                                unname(all_ran_num_genes[, x]))
    }
    
    results_ages <- results[ages]
    save(results_ages, file = paste0("enrichment_by_age/", folder, "/", overlap_type, "_", CNV_set, "_", total_permuts, "permuts.RData"))
    results_sizes <- results[levels(intronic_RANGES$size_quant20)]
    save(results_sizes, file = paste0("enrichment_by_size/", folder, "/",  overlap_type, "_", CNV_set, "_", total_permuts, "permuts.RData"))
  }
  
  print(CNV_set)
  
}

print("End of global (splitting CNVs) enrichment by age and size")





#### CUSTOM functions for randomization
### Randomizations 
## LOCAL - custom

# Folder for results
folder <- "local_custom"
for(CNV_set in all_CNV_sets) {
  load(paste0("local_randomization/output/custom/local_custom_random_starts", 
              CNV_set, "_", total_permuts, "permuts.RData"))
  
  # CNV_ranges
  CNV_ranges <- random_starts[[1]]           
  starts     <- t(random_starts[[2]])
  
  # starts <- starts[, 1:100] #  This was just to have a smaller test
  # total_permuts <- 100   This was just to have a smaller test
  
  obs_within <- sum(overlapsAny(CNV_ranges, intronic_RANGES, type = "within")) 
  obs_any    <- sum(overlapsAny(CNV_ranges, intronic_RANGES, type = "any"))
  
  
  # Random values
  ran_within <- c()
  ran_any    <- c()
  while(ncol(starts) > 0) {
    random_cnvs <- GRanges(seqnames = seqnames(CNV_ranges),
                           ranges   = IRanges(start = starts[, 1],
                                              width = width(CNV_ranges)))
    
    
    ran_within <- c(ran_within, sum(overlapsAny(random_cnvs, intronic_RANGES, type = "within")))
    ran_any    <- c(ran_any   , sum(overlapsAny(random_cnvs, intronic_RANGES, type = "any")))
    
    starts      <- starts[, -1]
    
    if(class(starts) == "integer") {
      starts      <- matrix(starts, ncol = 1)
    }
    if(ncol(starts)%%100 == 0) print(ncol(starts))
  }
  
  
  
  # Calculating p.value as it would be calculated using permTest from regioneR
  results_within <- calculate_pval_randomizations(obs = obs_within,
                                                  ran = ran_within)
  results_any    <- calculate_pval_randomizations(obs = obs_any,
                                                  ran = ran_any)
  
  
  save(results_within, results_any, 
       file = paste0("enrichment_general_within_any/", folder, "/results_within_and_any_", 
                     CNV_set, "_",total_permuts, "permuts.RData"))
  
  # total_permuts <- 10000 
}





### PLOTS
col_and_name <- list(randomization_name =  c(local_random_intervals = "Local randomization",
                                             global_random_intervals = "Global randomization",
                                             global_split_cnvs = "Global (splitting CNV)",
                                             local_custom = "Local (custom)",
                                             global_custom = "Global (custom)"),
                     randomization_color = c(local_random_intervals = "royalblue",
                                             global_random_intervals = "gold1",
                                             global_split_cnvs = "orange3",
                                             local_custom = "firebrick2",
                                             global_custom = "seagreen3"))
overlap_type_names <- c(intronic_within = "Purely intronic deletions",
                        coding = "Coding-overlapping deletions")



### Figures (boxplot)
results_tables <- list()
for(age_or_size in c("enrichment_by_age", "enrichment_by_size")) {
  for(enrichment_in in c(#"num_genes", 
                         "num_cnvs", 
                         "num_introns")) {
    for(randomization_type in c(#"local_random_intervals",
                                "global_split_cnvs",
                                "global_random_intervals")) {
                            # "local_custom")) {
                            # global_custom)) {
  # if(randomization_type == "global_split_cnvs") {
  #   total_permuts <- 1000
  # } else {
  #   total_permuts <- 10000
  # }
  
  # pdf(paste0(age_or_size, "/", randomization_type, "/", randomization_type, "_enrichment_in_", enrichment_in, 
             # ".pdf"), width = 10, height = 5.5)
  
    par(mfrow = c(1,3), oma = c(1,0.5,0,0), mar = c(7,4,4,1))
    
    overlap_type <- "intronic_within" 
      for(CNV_set in all_CNV_sets[c(1,3,4)]) {
        
          # LOAD INTRONIC
          print(load(paste0(age_or_size, "/", randomization_type, "/within_", CNV_set, "_",  total_permuts, "permuts.RData")))
          
          if(age_or_size == "enrichment_by_age") results <- results_ages
          if(age_or_size == "enrichment_by_size") results <- results_sizes
          all_obs <- unlist(lapply(results, function(x) x[[enrichment_in]][["observed"]]))
          all_ran <- do.call("cbind", lapply(results, function(x) x[[enrichment_in]][["permuted"]]))
          alternatives <- unlist(lapply(results, function(x) x[[enrichment_in]][["alternative"]]))
          all_pvals <- unlist(lapply(results, function(x) x[[enrichment_in]][["pval"]]))
          ylim <- c(-3.3,3.3)
          
          percent_more_or_less <- (all_obs/apply(all_ran, 2, median)-1)*100
        
        
        
        color <- c(Abyzov_RANGES = "#FF9400",
                   Handsaker_DELS = "#93184E",
                   Zarrei_DELS = "#66A1D2",
                   Sud15_DELS = "#116A66",
                   Phase3_DELS = "#043C6B")[CNV_set]
        
        
        
        FCs <- list()
        for(age in names(all_obs)) {
          FCs[[age]] <-  log2(all_obs[age]/all_ran[,age])
        }
        
        bp <- boxplot(FCs, 
                      outline = F,  
                      main = paste(CNV_set, " - ", col_and_name$randomization_name[randomization_type],
                                   "\n", overlap_type_names[overlap_type]), 
                      ylim = ylim,
                      lwd = 0.7,
                      # border = border,
                      ylab= list(num_genes   = "log2(obs/exp number of genes)", 
                                 num_cnvs    = "log2(obs/exp number of deletions)",
                                 num_introns = "log2(obs/exp number of introns)")[[enrichment_in]],
                      col = color, names = NA,
                      cex.main = 1.35)
        abline(h = 0)
        
        if(age_or_size == "enrichment_by_age") {
          text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
               labels = names(results),
               srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
          mtext("Ancient --> Recent", side = 1, line = 5, xpd = T, cex = 0.8)
          mtext("---", side = 1, line = 5, xpd = T, cex = 0.8)
        }
        
        for(num_stars in c(1:3)) {
          if(num_stars == 1) {
            min_pval <- 0.005
            max_pval <- 0.05
          }
          if(num_stars == 2) {
            min_pval <- 0.0005
            max_pval <- 0.005
          }
          if(num_stars == 3) {
            min_pval <- 0
            max_pval <- 0.0005
          }
          # positive and significant FCs (add star)
          if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater") > 0) {
            trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater"
            yplus <- (bp$stats[5,])[trues] + 0.04
            text(x = (1:16)[trues], 
                 y = yplus,
                 rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                 xpd = T,
                 srt = 90,
                 adj = 0,
                 # pos = 3,
                 cex = 1.2, col = "red", xpd = T) # indianred1
          }
          
          
          # negative and significant FCs(add star)
          if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less") > 0) {
            trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less"
            yminus <- (bp$stats[1,])[trues] - 0.04
            text(x = (1:16)[trues], 
                 y = yminus,
                 rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                 xpd = T,
                 srt = 90,
                 adj = 1,
                 # pos = 3,
                 cex = 1.2, col = "red", xpd = T) # indianred1
          }
        }
        
      
        results_tables[[randomization_type]][[enrichment_in]][[age_or_size]][[CNV_set]] <- rbind(ratio = unlist(lapply(FCs, median)),
                                                                                                 percent_difference_from_expected = percent_more_or_less,
                                                                                                 p.val = all_pvals)
        
        }
    
   
    # dev.off()
    }
  }
}

for(randomization_type in c("local_random_intervals",
                            "global_random_intervals")) {  
  table_by_age_num_cnvs <- cbind(t(results_tables[[randomization_type]]$num_cnvs$enrichment_by_age$Phase3_DELS[c(1,3), ]),
                                 t(results_tables[[randomization_type]]$num_cnvs$enrichment_by_age$Zarrei_DELS[c(1,3), ]),
                                 t(results_tables[[randomization_type]]$num_cnvs$enrichment_by_age$Abyzov_RANGES[c(1,3), ]))
  
  table_by_age_num_cnvs <- as.data.frame(table_by_age_num_cnvs, stringsAsFactors = F)
  for(col_x in 1:6) {
    table_by_age_num_cnvs[, col_x] <- as.character(round(table_by_age_num_cnvs[, col_x], 4))
    table_by_age_num_cnvs[, col_x] <- gsub(table_by_age_num_cnvs[, col_x], pattern = ".", replacement = ",",fixed = T)
  }
  
  
  write.xlsx(table_by_age_num_cnvs,
              # sep = "\t",
              # quote = F,
              file = paste0("all_pvals_", randomization_type, "_num_cnvs.xls"))
}



## GROUPED AGES 

grouped_ages <- list(Ancient = c("FungiMetazoa", "Bilateria", "Chordata",
                                 "Euteleostomi",  "Sarcopterygii"),
                     Middle  = c("Tetrapoda", "Amniota", "Mammalia",
                                 "Theria","Eutheria" ),
                     Young   = c("Simiiformes", "Catarrhini", "Hominoidea",
                                 "Hominidae", "HomoPanGorilla", "HomoSapiens" ))

### Figures (boxplot)
results_tables_grouped <- list()
for(age_or_size in c("enrichment_by_age")) {
  for(enrichment_in in c("num_genes", "num_cnvs", "num_introns")) {
    for(randomization_type in c("local_random_intervals",
                                "global_random_intervals")) {
      # "local_custom")) {
      # global_custom)) {
      
      
      # pdf(paste0(age_or_size, "/", randomization_type, "/grouped_ages_", randomization_type, "_enrichment_in_", enrichment_in, 
      #            ".pdf"), width = 8, height = 4)
      
      
      
      
      # par(mfrow = c(1,2), oma = c(1,0.5,0,0), mar = c(7,4,4,2))
      par(mfrow = c(1,3), oma = c(1,0.5,0,0), mar = c(7,4,4,1))
      
      for(overlap_type in c("intronic_within")) {
        #, "coding")) {
        for(CNV_set in all_CNV_sets[c(1,3,4)]) {
          if(overlap_type == "intronic_within") {
            # LOAD INTRONIC
            load(paste0(age_or_size, "/", randomization_type, "/within_", CNV_set, "_",  total_permuts, "permuts.RData"))
            
            if(age_or_size == "enrichment_by_age") results <- results_ages
            if(age_or_size == "enrichment_by_size") results <- results_sizes
            all_obs <- unlist(lapply(results, function(x) x[[enrichment_in]][["observed"]]))
            all_ran <- do.call("cbind", lapply(results, function(x) x[[enrichment_in]][["permuted"]]))
            
            
            
            
            ## GROUP RESULTS
            grouped_obs <- c(ancient = sum(all_obs[grouped_ages$Ancient]),
                             middle  = sum(all_obs[grouped_ages$Middle]),
                             young   = sum(all_obs[grouped_ages$Young]))
            
            grouped_ran <- cbind(ancient = apply(all_ran[, grouped_ages$Ancient], 1, sum),
                                 middle  = apply(all_ran[, grouped_ages$Middle], 1, sum),
                                 young   = apply(all_ran[, grouped_ages$Young], 1, sum))
            
            
            percent_more_or_less <- (grouped_obs/apply(grouped_ran, 2, median)-1)*100
            
            
            alternatives <- c()
            all_pvals    <- c()
            for(col_x in 1:3) {
              alternatives <- c(alternatives, calculate_pval_randomizations(grouped_obs[col_x], grouped_ran[, col_x])$alternative)
              all_pvals    <- c(all_pvals, calculate_pval_randomizations(grouped_obs[col_x], grouped_ran[, col_x])$pval)
            }
            ylim <- c(-1.2,1.2)
            
            
          } 
          
          
          color <- c(Abyzov_RANGES = "#FF9400",
                     Handsaker_DELS = "#93184E",
                     Zarrei_DELS = "#66A1D2",
                     Sud15_DELS = "#116A66",
                     Phase3_DELS = "#043C6B")[CNV_set]
          
          
          
          FCs <- list()
          for(age in c("ancient", "middle", "young")) {
            FCs[[age]] <-  log2(grouped_obs[age]/grouped_ran[,age])
          }
          
           
          bp <- boxplot(FCs, 
                        outline = F,  
                        main = paste(CNV_set, " - ", col_and_name$randomization_name[randomization_type],
                                     "\n", overlap_type_names[overlap_type]), 
                        ylim = ylim,
                        lwd = 0.7,
                        # border = border,
                        ylab= list(num_genes   = "log2(obs/exp number of genes)", 
                                   num_cnvs    = "log2(obs/exp number of deletions)",
                                   num_introns = "log2(obs/exp number of introns)")[[enrichment_in]],
                        col = color, names = NA,
                        cex.main = 1.35)
          abline(h = 0)
          
         
            text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
                 labels = c("Ancient", "Middle", "Young"),
                 srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
         
          
          for(num_stars in c(1:3)) {
            if(num_stars == 1) {
              min_pval <- 0.005
              max_pval <- 0.05
            }
            if(num_stars == 2) {
              min_pval <- 0.0005
              max_pval <- 0.005
            }
            if(num_stars == 3) {
              min_pval <- 0
              max_pval <- 0.0005
            }
            # positive and significant FCs (add star)
            if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater") > 0) {
              trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater"
              yplus <- (bp$stats[5,])[trues] + 0.04
              text(x = (1:length(all_pvals))[trues], 
                   y = yplus,
                   rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                   xpd = T,
                   srt = 90,
                   adj = 0,
                   # pos = 3,
                   cex = 1.2, col = "red", xpd = T) # indianred1
            }
            
            
            # negative and significant FCs(add star)
            if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less") > 0) {
              trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less"
              yminus <- (bp$stats[1,])[trues] - 0.04
              text(x = (1:length(all_pvals))[trues], 
                   y = yminus,
                   rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                   xpd = T,
                   srt = 90,
                   adj = 1,
                   # pos = 3,
                   cex = 1.2, col = "red", xpd = T) # indianred1
            }
          }

            results_tables_grouped[[randomization_type]][[enrichment_in]][[age_or_size]][[CNV_set]] <- rbind(ratio = unlist(lapply(FCs, median)),
                                                                                                     percent_difference_from_expected = percent_more_or_less,
                                                                                                     p.val = all_pvals)
            
            
            
            
        }
        
      }
       
      dev.off()
    }
  }
}

# create excel files
for(randomization_type in c("local_random_intervals",
                            "global_random_intervals")) {  
  table_by_age_num_cnvs <- cbind(t(results_tables_grouped[[randomization_type]]$num_cnvs$enrichment_by_age$Phase3_DELS[c(1,3), ]),
                                 t(results_tables_grouped[[randomization_type]]$num_cnvs$enrichment_by_age$Zarrei_DELS[c(1,3), ]),
                                 t(results_tables_grouped[[randomization_type]]$num_cnvs$enrichment_by_age$Abyzov_RANGES[c(1,3), ]))
  
  table_by_age_num_cnvs <- as.data.frame(table_by_age_num_cnvs, stringsAsFactors = F)
  for(col_x in 1:6) {
    table_by_age_num_cnvs[, col_x] <- as.character(round(table_by_age_num_cnvs[, col_x], 4))
    table_by_age_num_cnvs[, col_x] <- gsub(table_by_age_num_cnvs[, col_x], pattern = ".", replacement = ",",fixed = T)
  }
  
  
  write.xlsx(table_by_age_num_cnvs,
             # sep = "\t",
             # quote = F,
             file = paste0("all_pvals_GROUPED_AGES_intronic_within", randomization_type, "_num_cnvs.xls"))
}






### PSEUDOGROUPED AGES
grouped_ages <- list(FungiMetazoa = "FungiMetazoa", 
                     Bilateria =     "Bilateria", 
                     Chordata =           "Chordata",
                     Euteleostomi =            "Euteleostomi",
                     "Sarcopterygii to Eutheria"  = c("Sarcopterygii", "Tetrapoda", "Amniota", "Mammalia",
                                 "Theria","Eutheria" ),
                     Primates   = c("Simiiformes", "Catarrhini", "Hominoidea",
                                 "Hominidae", "HomoPanGorilla", "HomoSapiens" ))

### Figures (boxplot)
for(age_or_size in c("enrichment_by_age")) {
  for(enrichment_in in c("num_genes", "num_cnvs", "num_introns")) {
    for(randomization_type in c("local_random_intervals",
                                "global_random_intervals")) {
      # "local_custom")) {
      # global_custom)) {
      
      
      pdf(paste0(age_or_size, "/", randomization_type, "/pseudogrouped_ages_", randomization_type, "_enrichment_in_", enrichment_in, 
                 ".pdf"), width = 8, height = 4)
      # ".pdf"), width = 8, height = 3.5)
      
      
      
      # par(mfrow = c(1,2), oma = c(1,0.5,0,0), mar = c(7,4,4,2))
      par(mfrow = c(1,3), oma = c(1,0.5,0,0), mar = c(7,4,4,1))
      
      for(overlap_type in c("intronic_within")) {
        #, "coding")) {
        for(CNV_set in all_CNV_sets[c(1,3,4)]) {
          if(overlap_type == "intronic_within") {
            # LOAD INTRONIC
            load(paste0(age_or_size, "/", randomization_type, "/within_", CNV_set, "_",  total_permuts, "permuts.RData"))
            
            if(age_or_size == "enrichment_by_age") results <- results_ages
            if(age_or_size == "enrichment_by_size") results <- results_sizes
            all_obs <- unlist(lapply(results, function(x) x[[enrichment_in]][["observed"]]))
            all_ran <- do.call("cbind", lapply(results, function(x) x[[enrichment_in]][["permuted"]]))
            
            ## GROUP RESULTS
            grouped_obs <- c(sum(all_obs[grouped_ages$FungiMetazoa]),
                             sum(all_obs[grouped_ages$Bilateria]),
                             sum(all_obs[grouped_ages$Chordata]),
                             sum(all_obs[grouped_ages$Euteleostomi]),
                             sum(all_obs[grouped_ages$`Sarcopterygii to Eutheria`]),
                             sum(all_obs[grouped_ages$Primates]))
            names(grouped_obs) <- names(grouped_ages)
            
          
            
            grouped_ran <- cbind(all_ran[,"FungiMetazoa"],
                                 all_ran[,"Bilateria"],
                                 all_ran[,"Chordata"],
                                 all_ran[,"Euteleostomi"],
                                 apply(all_ran[,grouped_ages$`Sarcopterygii to Eutheria`], 1, sum),
                                 apply(all_ran[,grouped_ages$Primates], 1, sum))
            colnames(grouped_ran) <- names(grouped_obs)
            
            alternatives <- c()
            all_pvals    <- c()
            for(col_x in 1:length(grouped_obs)) {
              alternatives <- c(alternatives, calculate_pval_randomizations(grouped_obs[col_x], grouped_ran[, col_x])$alternative)
              all_pvals    <- c(all_pvals, calculate_pval_randomizations(grouped_obs[col_x], grouped_ran[, col_x])$pval)
            }
            ylim <- c(-1.2,1.2)
            
            
          } 
          
          
          color <- c(Abyzov_RANGES = "#FF9400",
                     Handsaker_DELS = "#93184E",
                     Zarrei_DELS = "#66A1D2",
                     Sud15_DELS = "#116A66",
                     Phase3_DELS = "#043C6B")[CNV_set]
          
          
          
          FCs <- list()
          for(age in names(grouped_ages)) {
            FCs[[age]] <-  log2(grouped_obs[age]/grouped_ran[,age])
          }
          
          bp <- boxplot(FCs, 
                        outline = F,  
                        main = paste(CNV_set, " - ", col_and_name$randomization_name[randomization_type],
                                     "\n", overlap_type_names[overlap_type]), 
                        ylim = ylim,
                        lwd = 0.7,
                        # border = border,
                        ylab= list(num_genes   = "log2(obs/exp number of genes)", 
                                   num_cnvs    = "log2(obs/exp number of deletions)",
                                   num_introns = "log2(obs/exp number of introns)")[[enrichment_in]],
                        col = color, names = NA,
                        cex.main = 1.35)
          abline(h = 0)
          
          
          text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
               labels = names(grouped_ages),
               srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
          
          
          for(num_stars in c(1:3)) {
            if(num_stars == 1) {
              min_pval <- 0.005
              max_pval <- 0.05
            }
            if(num_stars == 2) {
              min_pval <- 0.0005
              max_pval <- 0.005
            }
            if(num_stars == 3) {
              min_pval <- 0
              max_pval <- 0.0005
            }
            # positive and significant FCs (add star)
            if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater") > 0) {
              trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater"
              yplus <- (bp$stats[5,])[trues] + 0.04
              text(x = (1:length(all_pvals))[trues], 
                   y = yplus,
                   rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                   xpd = T,
                   srt = 90,
                   adj = 0,
                   # pos = 3,
                   cex = 1.2, col = "red", xpd = T) # indianred1
            }
            
            
            # negative and significant FCs(add star)
            if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less") > 0) {
              trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less"
              yminus <- (bp$stats[1,])[trues] - 0.04
              text(x = (1:length(all_pvals))[trues], 
                   y = yminus,
                   rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                   xpd = T,
                   srt = 90,
                   adj = 1,
                   # pos = 3,
                   cex = 1.2, col = "red", xpd = T) # indianred1
            }
          }
          
        }
      }
      # results_tables[[randomization_type]][[overlap_type]] <- rbind(ratio = ratios,
      #                                                               percent_difference_from_expected = percent_more_or_less,
      #                                                               p.val = pvals)
      
      dev.off()
    }
  }
}


#### ONLY PRIMATES GROUPED
grouped_ages <- list(FungiMetazoa = "FungiMetazoa", 
                     Bilateria =    "Bilateria", 
                     Chordata =     "Chordata",
                     Euteleostomi = "Euteleostomi",
                     Sarcopterygii = "Sarcopterygii",
                     Tetrapoda = "Tetrapoda",
                     Amniota = "Amniota",
                     Mammalia = "Mammalia",
                     Theria = "Theria",
                     Eutheria = "Eutheria",
                     Primates   = c("Simiiformes", "Catarrhini", "Hominoidea",
                                    "Hominidae", "HomoPanGorilla", "HomoSapiens" ))

### Figures (boxplot)
for(age_or_size in c("enrichment_by_age")) {
  for(enrichment_in in c("num_genes", "num_cnvs", "num_introns")) {
    for(randomization_type in c("local_random_intervals",
                                "global_random_intervals")) {
      # "local_custom")) {
      # global_custom)) {
      
      
      pdf(paste0(age_or_size, "/", randomization_type, "/primates_grouped_ages_", randomization_type, "_enrichment_in_", enrichment_in, 
                 ".pdf"), width = 8, height = 4)
      # ".pdf"), width = 8, height = 3.5)
      
      
      
      # par(mfrow = c(1,2), oma = c(1,0.5,0,0), mar = c(7,4,4,2))
      par(mfrow = c(1,3), oma = c(0,0.5,2,0), mar = c(6,4,4,1))
      
      for(overlap_type in c("intronic_within")) {
        #, "coding")) {
        for(CNV_set in all_CNV_sets[c(1,3,4)]) {
          if(overlap_type == "intronic_within") {
            # LOAD INTRONIC
            load(paste0(age_or_size, "/", randomization_type, "/within_", CNV_set, "_",  total_permuts, "permuts.RData"))
            
            if(age_or_size == "enrichment_by_age") results <- results_ages
            if(age_or_size == "enrichment_by_size") results <- results_sizes
            all_obs <- unlist(lapply(results, function(x) x[[enrichment_in]][["observed"]]))
            all_ran <- do.call("cbind", lapply(results, function(x) x[[enrichment_in]][["permuted"]]))
            
            ## GROUP RESULTS
            grouped_obs <- c()
            for(gr in names(grouped_ages)) {
              grouped_obs <- c(grouped_obs, sum(all_obs[grouped_ages[[gr]]]))
            }
            names(grouped_obs) <- names(grouped_ages)
            
            grouped_ran <- c()
            for(gr in names(grouped_ages)) {
              if(length(grouped_ages[[gr]]) > 1) {
                grouped_ran <- cbind(grouped_ran, apply(all_ran[,grouped_ages[[gr]]], 1, sum))
              } else {
                grouped_ran <- cbind(grouped_ran, all_ran[, grouped_ages[[gr]]])
              }
            }
            colnames(grouped_ran) <- names(grouped_obs)
            
            alternatives <- c()
            all_pvals    <- c()
            for(col_x in 1:length(grouped_obs)) {
              alternatives <- c(alternatives, calculate_pval_randomizations(grouped_obs[col_x], grouped_ran[, col_x])$alternative)
              all_pvals    <- c(all_pvals, calculate_pval_randomizations(grouped_obs[col_x], grouped_ran[, col_x])$pval)
            }
            ylim <- c(-1.7,1.7)
            
            
          } 
          
          
          color <- c(Abyzov_RANGES = "#FF9400",
                     Handsaker_DELS = "#93184E",
                     Zarrei_DELS = "#66A1D2",
                     Sud15_DELS = "#116A66",
                     Phase3_DELS = "#043C6B")[CNV_set]
          
          
          
          FCs <- list()
          for(age in names(grouped_ages)) {
            FCs[[age]] <-  log2(grouped_obs[age]/grouped_ran[,age])
          }
          
          bp <- boxplot(FCs, 
                        outline = F,  
                        main = paste(CNV_set, " - ", col_and_name$randomization_name[randomization_type],
                                     "\n", overlap_type_names[overlap_type]), 
                        ylim = ylim,
                        lwd = 0.7,
                        # border = border,
                        ylab= list(num_genes   = "log2(obs/exp number of genes)", 
                                   num_cnvs    = "log2(obs/exp number of deletions)",
                                   num_introns = "log2(obs/exp number of introns)")[[enrichment_in]],
                        col = color, names = NA,
                        cex.main = 1.35)
          abline(h = 0)
          
          
          text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
               labels = names(grouped_ages),
               srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
          
          
          for(num_stars in c(1:3)) {
            if(num_stars == 1) {
              min_pval <- 0.005
              max_pval <- 0.05
            }
            if(num_stars == 2) {
              min_pval <- 0.0005
              max_pval <- 0.005
            }
            if(num_stars == 3) {
              min_pval <- 0
              max_pval <- 0.0005
            }
            # positive and significant FCs (add star)
            if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater") > 0) {
              trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater"
              yplus <- (bp$stats[5,])[trues] + 0.04
              text(x = (1:length(all_pvals))[trues], 
                   y = yplus,
                   rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                   xpd = T,
                   srt = 90,
                   adj = 0,
                   # pos = 3,
                   cex = 1.2, col = "red", xpd = T) # indianred1
            }
            
            
            # negative and significant FCs(add star)
            if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less") > 0) {
              trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less"
              yminus <- (bp$stats[1,])[trues] - 0.04
              text(x = (1:length(all_pvals))[trues], 
                   y = yminus,
                   rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                   xpd = T,
                   srt = 90,
                   adj = 1,
                   # pos = 3,
                   cex = 1.2, col = "red", xpd = T) # indianred1
            }
          }
          
        }
      }
      # results_tables[[randomization_type]][[overlap_type]] <- rbind(ratio = ratios,
      #                                                               percent_difference_from_expected = percent_more_or_less,
      #                                                               p.val = pvals)
      
      dev.off()
    }
  }
}




######### GROUPED SIZES

grouped_sizes <- matrix(levels(intronic_RANGES$size_quant20), ncol = 2, byrow = T)
names_grouped_sizes <- apply(grouped_sizes, 1, function(x) paste(unlist(strsplit(x, split = ","))[c(1,4)], collapse = ","))
xxx <- rep(names_grouped_sizes, each = 2)

grouped_sizes <- by(levels(intronic_RANGES$size_quant20), xxx, as.character)
grouped_sizes <- grouped_sizes[names_grouped_sizes]

results_tables_grouped_sizes <- list()
### Figures (boxplot)
for(age_or_size in c("enrichment_by_size")) {
  for(enrichment_in in c("num_genes", "num_cnvs", "num_introns")) {
    for(randomization_type in c("local_random_intervals",
                                "global_random_intervals")) {
    

      # pdf(paste0(age_or_size, "/", randomization_type, "/grouped_sizes_", randomization_type, "_enrichment_in_", enrichment_in,
      #            ".pdf"), width = 8, height = 4)
      # 
      # par(mfrow = c(1,3), oma = c(1,0.5,0,0), mar = c(7,4,4,1))
      
      overlap_type  <- "intronic_within"
        
        for(CNV_set in all_CNV_sets[c(1,3,4)]) {
          if(overlap_type == "intronic_within") {
            # LOAD INTRONIC
            load(paste0(age_or_size, "/", randomization_type, "/within_", CNV_set, "_",  total_permuts, "permuts.RData"))
            
            if(age_or_size == "enrichment_by_age") results <- results_ages
            if(age_or_size == "enrichment_by_size") results <- results_sizes
            all_obs <- unlist(lapply(results, function(x) x[[enrichment_in]][["observed"]]))
            all_ran <- do.call("cbind", lapply(results, function(x) x[[enrichment_in]][["permuted"]]))
            
            ## GROUP RESULTS
            ## GROUP RESULTS
            grouped_obs <- c()
            for(gr in names(grouped_sizes)) {
              grouped_obs <- c(grouped_obs, sum(all_obs[grouped_sizes[[gr]]]))
            }
            names(grouped_obs) <- names(grouped_sizes)
            
            grouped_ran <- c()
            for(gr in names(grouped_sizes)) {
              if(length(grouped_sizes[[gr]]) > 1) {
                grouped_ran <- cbind(grouped_ran, apply(all_ran[,grouped_sizes[[gr]]], 1, sum))
              } else {
                grouped_ran <- cbind(grouped_ran, all_ran[, grouped_sizes[[gr]]])
              }
            }
            colnames(grouped_ran) <- names(grouped_obs)
            
            percent_more_or_less <- (grouped_obs/apply(grouped_ran, 2, median)-1)*100
            
            alternatives <- c()
            all_pvals    <- c()
            for(col_x in 1:length(grouped_obs)) {
              alternatives <- c(alternatives, calculate_pval_randomizations(grouped_obs[col_x], grouped_ran[, col_x])$alternative)
              all_pvals    <- c(all_pvals, calculate_pval_randomizations(grouped_obs[col_x], grouped_ran[, col_x])$pval)
            }
            ylim <- c(-3,3)
            
            
          } 
          
          
          color <- c(Abyzov_RANGES = "#FF9400",
                     Handsaker_DELS = "#93184E",
                     Zarrei_DELS = "#66A1D2",
                     Sud15_DELS = "#116A66",
                     Phase3_DELS = "#043C6B")[CNV_set]
          
          
          
          FCs <- list()
          for(size in names(grouped_sizes)) {
            FCs[[size]] <-  log2(grouped_obs[size]/grouped_ran[,size])
          }
          
          bp <- boxplot(FCs, 
                        outline = F,  
                        main = paste(CNV_set, " - ", col_and_name$randomization_name[randomization_type],
                                     "\n", overlap_type_names[overlap_type]), 
                        ylim = ylim,
                        lwd = 0.7,
                        # border = border,
                        ylab= list(num_genes   = "log2(obs/exp number of genes)", 
                                   num_cnvs    = "log2(obs/exp number of deletions)",
                                   num_introns = "log2(obs/exp number of introns)")[[enrichment_in]],
                        col = color, names = NA,
                        cex.main = 1.35)
          abline(h = 0)
          
          
          text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
               labels = names(grouped_sizes),
               srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
          
          
          for(num_stars in c(1:3)) {
            if(num_stars == 1) {
              min_pval <- 0.005
              max_pval <- 0.05
            }
            if(num_stars == 2) {
              min_pval <- 0.0005
              max_pval <- 0.005
            }
            if(num_stars == 3) {
              min_pval <- 0
              max_pval <- 0.0005
            }
            # positive and significant FCs (add star)
            if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater") > 0) {
              trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "greater"
              yplus <- (bp$stats[5,])[trues] + 0.04
              text(x = (1:length(all_pvals))[trues], 
                   y = yplus,
                   rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                   xpd = T,
                   srt = 90,
                   adj = 0,
                   # pos = 3,
                   cex = 1.2, col = "red", xpd = T) # indianred1
            }
            
            
            # negative and significant FCs(add star)
            if(sum(all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less") > 0) {
              trues <- all_pvals >= min_pval & all_pvals < max_pval & alternatives == "less"
              yminus <- (bp$stats[1,])[trues] - 0.04
              text(x = (1:length(all_pvals))[trues], 
                   y = yminus,
                   rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                   xpd = T,
                   srt = 90,
                   adj = 1,
                   # pos = 3,
                   cex = 1.2, col = "red", xpd = T) # indianred1
            }
            
            results_tables_grouped_sizes[[randomization_type]][[enrichment_in]][[age_or_size]][[CNV_set]] <- rbind(ratio = unlist(lapply(FCs, median)),
                                                                                                                   percent_difference_from_expected = percent_more_or_less,
                                                                                                                   p.val = all_pvals)
            
            
          }
          
        }
     
        
      # dev.off()
    }
  }
}

# create excel file
for(randomization_type in c("local_random_intervals",
                            "global_random_intervals")) {  
  table_by_size_num_cnvs <- cbind(t(results_tables_grouped_sizes[[randomization_type]]$num_cnvs$enrichment_by_size$Phase3_DELS[c(1,3), ]),
                                 t(results_tables_grouped_sizes[[randomization_type]]$num_cnvs$enrichment_by_size$Zarrei_DELS[c(1,3), ]),
                                 t(results_tables_grouped_sizes[[randomization_type]]$num_cnvs$enrichment_by_size$Abyzov_RANGES[c(1,3), ]))
  
  table_by_size_num_cnvs <- as.data.frame(table_by_size_num_cnvs, stringsAsFactors = F)
  for(col_x in 1:6) {
    table_by_size_num_cnvs[, col_x] <- as.character(round(table_by_size_num_cnvs[, col_x], 4))
    table_by_size_num_cnvs[, col_x] <- gsub(table_by_size_num_cnvs[, col_x], pattern = ".", replacement = ",",fixed = T)
  }
  
  
  write.xlsx(table_by_size_num_cnvs,
             # sep = "\t",
             # quote = F,
             file = paste0("all_pvals_GROUPED_sizes_intronic_within", randomization_type, "_num_cnvs.xls"))
}


