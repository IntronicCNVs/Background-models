# Here we calculate the observed and expected number of deletions  per intron or gene  or the number of genes or introns with intronic deletions
# separating by gene/intron ages and sizes
# Overlap type "within": purely intronic introns
# Overlap type "any": any kind of overlap with introns

library( "regioneR" )
library( "BSgenome.Hsapiens.UCSC.hg19.masked" )
library( "seqbias" )
library( "gtools")
library( "xlsx")

# #### Data loading and preparation
# Set working directory to Background-models
setwd("Background-models/")
# Create necessary folders
if (!file.exists("enrichment_by_age_size")) dir.create(file.path(".", "enrichment_by_age_size"))
for(folder_name in c("enrichment_by_age", "enrichment_by_size")) {
  for(subfolder_name in c("global_random_intervals", "local_random_intervals")) {
    if (!file.exists(folder_name)) dir.create(file.path(".", folder_name))
    if (!file.exists(paste0(folder_name, "/", subfolder_name))) dir.create(file.path(folder_name, subfolder_name))
  }
}



# #### Data loading and preparation
load("deletions_introns_genes_ages.RData")
seqlevelsStyle(intronic_RANGES) <- "UCSC"
seqlevelsStyle(protein_coding_RANGES)      <- "UCSC"
source("calculate_pvalue_randomizations.R")


### Number of permutations
total_permuts <- 100

# Datasets
all_CNV_sets <- c("Phase3_DELS",
                  "Handsaker_DELS",
                  "Zarrei_DELS",
                  "Abyzov_RANGES",
                  "Sud15_DELS")

# Ages and size groups
ages <- levels(intronic_RANGES$most_recent_age)
intronic_RANGES$size_quant20 <- quantcut(width(intronic_RANGES), 20)
protein_coding_RANGES$Age <- Ages[names(protein_coding_RANGES), "GeneAge"]
protein_coding_RANGES <- protein_coding_RANGES[!is.na(protein_coding_RANGES$Age)]

#######################################################################

## Function to calculate overlaps    
by_age_and_size <- function(granges, overlap_type) {
  # by age
  num_cnvs_by_age     <- c() # Number of intronic CNVs in each age group
  num_introns_by_age  <- c() # Number of introns with intronic CNVs by age group
  num_genes_by_age    <- c() # Number of genes with intronic CNVs by age group
  # by intron size
  num_introns_by_size <- c() # Number of intronic CNVs in each size group
  num_cnvs_by_size    <- c() # Number of introns with intronic CNVs by size group
  num_genes_by_size   <- c() # Number of genes with intronic CNVs by size group
  
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
    # Putting together age and size numbers
    obs_num_cnvs    <- c(get(paste0("obs_", overlap_type))$num_cnvs_by_age,    get(paste0("obs_", overlap_type))$num_cnvs_by_size  )
    obs_num_introns <- c(get(paste0("obs_", overlap_type))$num_introns_by_age, get(paste0("obs_", overlap_type))$num_introns_by_size)
    obs_num_genes   <- c(get(paste0("obs_", overlap_type))$num_genes_by_age,   get(paste0("obs_", overlap_type))$num_genes_by_size)
    
    names(obs_num_cnvs)    <- c(ages, levels(intronic_RANGES$size_quant20))
    names(obs_num_introns) <- c(ages, levels(intronic_RANGES$size_quant20))
    names(obs_num_genes)   <- c(ages, levels(intronic_RANGES$size_quant20))
    
    all_ran_num_cnvs    <- cbind(get(paste0("ran_", overlap_type))$num_cnvs_by_age,     get(paste0("ran_", overlap_type))$num_cnvs_by_size)
    all_ran_num_introns <- cbind(get(paste0("ran_", overlap_type))$num_introns_by_age,  get(paste0("ran_", overlap_type))$num_introns_by_size)
    all_ran_num_genes   <- cbind(get(paste0("ran_", overlap_type))$num_genes_by_age,    get(paste0("ran_", overlap_type))$num_genes_by_size)
    
    # Calculate enrichment
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



## LOCAL background model

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


## PLOT RESULTS
### PLOTS
col_and_name <- list(randomization_name =  c(local_random_intervals = "Local randomization",
                                             global_random_intervals = "Global randomization"),
                     randomization_color = c(local_random_intervals = "royalblue",
                                             global_random_intervals = "gold1"))

overlap_type_names <- c(intronic_within = "Purely intronic deletions",
                        intronic_any    = "Intron-intersecting deletions") # Can be completely or partially overlapping  with introns
                        



### Figures (boxplot)
results_tables <- list()
for(age_or_size in c("enrichment_by_age", 
                     "enrichment_by_size")) {
  for(enrichment_in in c("num_genes",  
                         "num_cnvs", 
                         "num_introns")) {
    for(randomization_type in c("local_random_intervals",
                                "global_random_intervals")) {
      for(overlap_type in c("intronic_within" , "intronic_any")) {
        
        # pdf(paste0(age_or_size, "/", 
        #            randomization_type, "/", 
        #            randomization_type, "_enrichment_in_", enrichment_in,
        #            "_", overlap_type,
        #            ".pdf"), width = 10, height = 5.5)
        par(mfrow = c(1,3), oma = c(1,0.5,0,0), mar = c(8,4,4,1))
        
        for(CNV_set in all_CNV_sets[c(1,3,4)]) {
          # LOAD RESULTS
          print(load(paste0(age_or_size, "/", randomization_type, "/within_", CNV_set, "_",  total_permuts, "permuts.RData")))
          
          if(age_or_size == "enrichment_by_age")  results <- results_ages
          if(age_or_size == "enrichment_by_size") results <- results_sizes
          
          # all observed results
          all_obs <- unlist(lapply(results, function(x) x[[enrichment_in]][["observed"]]))
          # all random results
          all_ran <- do.call("cbind", lapply(results, function(x) x[[enrichment_in]][["permuted"]]))
          alternatives <- unlist(lapply(results, function(x) x[[enrichment_in]][["alternative"]]))
          all_pvals <- unlist(lapply(results, function(x) x[[enrichment_in]][["pval"]]))
          
          # % of deletions more or less than expected
          percent_more_or_less <- (all_obs/apply(all_ran, 2, median)-1)*100
          
          
          color <- c(Abyzov_RANGES = "#FF9400",
                     Handsaker_DELS = "#93184E",
                     Zarrei_DELS = "#66A1D2",
                     Sud15_DELS = "#116A66",
                     Phase3_DELS = "#043C6B")[CNV_set]
          
          # Calculate fold changes
          FCs <- list()
          for(age in names(all_obs)) {
            FCs[[age]] <-  log2(all_obs[age]/all_ran[,age])
          }
          
          bp <- boxplot(FCs, 
                        outline = F,  
                        main = paste(CNV_set, " - ", col_and_name$randomization_name[randomization_type],
                                     "\n", overlap_type_names[overlap_type]), 
                        ylim = c(-3.3,3.3),
                        lwd = 0.7,
                        # border = border,
                        ylab= list(num_genes   = "log2(obs/exp number of genes)", 
                                   num_cnvs    = "log2(obs/exp number of deletions)",
                                   num_introns = "log2(obs/exp number of introns)")[[enrichment_in]],
                        col = color, names = NA,
                        cex.main = 1.35)
          abline(h = 0)
          
          # x labels
          if(age_or_size == "enrichment_by_age") {
            text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
                 labels = names(results),
                 srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
            mtext("Ancient --> Recent", side = 1, line = 5.5, xpd = T, cex = 0.8)
          }
          if(age_or_size == "enrichment_by_size") {
            text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
                 labels = names(results),
                 srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
            mtext("Intron size (bp)", side = 1, line = 7, xpd = T, cex = 0.8)
           
          }
          
          # Significance (number of asterisks)
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
          
          
          results_tables[[randomization_type]][[enrichment_in]][[age_or_size]][[CNV_set]][[overlap_type]] <- rbind(ratio = unlist(lapply(FCs, median)),
                                                                                                                   percent_difference_from_expected = percent_more_or_less,
                                                                                                                   p.val = all_pvals)
      }
      
      
      # dev.off()
    }
  }
}
}

