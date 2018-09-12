library( "regioneR" )
library( "BSgenome.Hsapiens.UCSC.hg19.masked" )
library( "seqbias" )

# #### Data loading and preparation
# Set working directory to Background-models
setwd("Background-models/")
if (!file.exists("enrichment_general_within_any")) dir.create(file.path(".", "enrichment_general_within_any"))
if (!file.exists("enrichment_general_within_any/global_random_intervals")) dir.create(file.path("enrichment_general_within_any", "global_random_intervals"))
if (!file.exists("enrichment_general_within_any/local_random_intervals"))  dir.create(file.path("enrichment_general_within_any", "local_random_intervals"))

objects <- load("deletions_introns_genes_ages.RData")
seqlevelsStyle(intronic_RANGES) <- "UCSC"
seqlevelsStyle(protein_coding_RANGES)      <- "UCSC"
source("calculate_pvalue_randomizations.R")

### Number of permutations

total_permuts <- 10000

# Datasets
all_CNV_sets <- c("Phase3_DELS",
                  "Handsaker_DELS",
                  "Zarrei_DELS",
                  "Abyzov_RANGES",
                  "Sud15_DELS")


## GLOBAL RANDOMIZATION- random intervals

# Folder for results
folder <- "global_random_intervals"
for(CNV_set in all_CNV_sets) {
  ## RANDOM STARTS
  chromosomes <- read.table(paste0("global_random_intervals/output/global_random_intervals_random_chrs_", CNV_set, "_", total_permuts, ".txt"),
                            stringsAsFactors = F, header = F)
  print("chromosomes loaded")
  
  starts <- read.table(paste0("global_random_intervals/output/global_random_intervals_random_starts_", CNV_set, "_", total_permuts, ".txt"),
                            stringsAsFactors = F, header = F)
  print("starts loaded")
  
  # CNV_ranges
  load(paste0("global_random_intervals/output/global_random_intervals_original_", CNV_set, "_GRanges.RData"))
  
  print("CNV_ranges loaded")
  
  obs_within <- sum(overlapsAny(CNV_ranges, intronic_RANGES, type = "within")) 
  obs_any    <- sum(overlapsAny(CNV_ranges, intronic_RANGES, type = "any"))
  
  # Random values
  ran_within <- c()
  ran_any    <- c()
  while(ncol(chromosomes) > 0) {
    random_cnvs <- GRanges(seqnames = chromosomes[, 1],
                           ranges   = IRanges(start = starts[, 1],
                                              width = width(CNV_ranges)))
    
    
    ran_within <- c(ran_within, sum(overlapsAny(random_cnvs, intronic_RANGES, type = "within")))
    ran_any    <- c(ran_any   , sum(overlapsAny(random_cnvs, intronic_RANGES, type = "any")))
    
    chromosomes <- chromosomes[, -1]
    starts      <- starts[, -1]
   
    if(class(starts) == "integer") {
      chromosomes <- matrix(chromosomes, ncol = 1)
      starts      <- matrix(starts, ncol = 1)
    }
    if(ncol(chromosomes)%%100 == 0) print(ncol(chromosomes))
  }
  
  
  
  # Calculating p.value as it would be calculated using permTest from regioneR
  results_within <- calculate_pval_randomizations(obs = obs_within,
                                                  ran = ran_within)
  results_any    <- calculate_pval_randomizations(obs = obs_any,
                                                  ran = ran_any)
  
 
  save(results_within, results_any, 
       file = paste0("enrichment_general_within_any/", folder, "/results_within_and_any_", 
                     CNV_set, "_",
                     total_permuts, "permuts.RData"))

}




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
  
  obs_within <- sum(overlapsAny(CNV_ranges, intronic_RANGES, type = "within")) 
  obs_any    <- sum(overlapsAny(CNV_ranges, intronic_RANGES, type = "any"))
  
  
  # Random values
  ran_within <- c()
  ran_any    <- c()
  while(ncol(all_starts) > 0) {
    random_cnvs <- GRanges(seqnames = seqnames(CNV_ranges),
                           ranges   = IRanges(start = all_starts[, 1],
                                              width = width(CNV_ranges)))
    
    
    ran_within <- c(ran_within, sum(overlapsAny(random_cnvs, intronic_RANGES, type = "within")))
    ran_any    <- c(ran_any   , sum(overlapsAny(random_cnvs, intronic_RANGES, type = "any")))
    
    all_starts      <- all_starts[, -1]
    
    if(class(all_starts) == "integer") {
      all_starts      <- matrix(all_starts, ncol = 1)
    }
    if(ncol(all_starts)%%100 == 0) print(ncol(all_starts))
  }
  
  
  
  # Calculating p.value as it would be calculated using permTest from regioneR
  results_within <- calculate_pval_randomizations(obs = obs_within,
                                                  ran = ran_within)
  results_any    <- calculate_pval_randomizations(obs = obs_any,
                                                  ran = ran_any)
  
  
  save(results_within, results_any, 
       file = paste0("enrichment_general_within_any/", folder, "/results_within_and_any_", 
                     CNV_set, "_",
                     total_permuts, "permuts.RData"))
  
}




### PLOTS
col_and_name <- list(randomization_name =  c(local_random_intervals = "Local randomization",
                                             global_random_intervals = "Global randomization"),
                     randomization_color = c(local_random_intervals = "royalblue",
                                             global_random_intervals = "gold1"))


results_tables <- list()

for(randomization_type in c("local_random_intervals",
                            "global_random_intervals")) { 
  
  # pdf(file = paste0("enrichment_general_within_any/", randomization_type, "/", randomization_type,
  #            "_within_vs_any_general_",total_permuts, "_permuts.pdf"), width = 17, height = 9)
  titulo = "Purely intronic vs Intron intersecting deletions"
  par(mfrow = c(1,2))
  
  for(overlap_type in c("within", "any")) {
    ratios <- c()
    pvals  <- c()
    mads    <- c()
    alternatives <- c()
    percent_more_or_less <- c()
    for(CNV_set in all_CNV_sets[c(1,3,4,5,2)]) {
      load(paste0("enrichment_general_within_any/",
                  randomization_type, 
                  "/results_within_and_any_", CNV_set, "_", total_permuts, "permuts.RData"))

      
      results <- get(paste0("results_", overlap_type))  
     
      print(str(results))
      obs    <- results$observed  
      med    <- median(results$permuted)
      log2fc <- log2(obs/med)
      mad_x  <- mad(log2(obs/results$permuted))
      p.val  <- results$pval
      
        
      alternatives <- c(alternatives, results$alternative)
      ratios       <- c(ratios, log2fc)
      pvals        <- c(pvals, p.val)
      mads         <- c(mads, mad_x)
      percent_more_or_less <- c(percent_more_or_less, (obs/med-1)*100)
      
    }   
    
    
    names(ratios) <- all_CNV_sets[c(1,3,4,5,2)]
    
    names(ratios)[names(ratios) == "Phase3_DELS"] <- "Sudmant et al.\n(Nature)"
    names(ratios)[names(ratios) == "Handsaker_DELS"] <- "Handsaker et al.\n"
    names(ratios)[names(ratios) == "Zarrei_DELS"] <- "Zarrei et al.\n"
    names(ratios)[names(ratios) == "Abyzov_RANGES"] <- "Abyzov et al.\n"
    names(ratios)[names(ratios) == "Sud15_DELS"] <- "Sudmant et al.\n(Science)"
    
    
    
    
    par(mar = c(5,5.5,4,0.5))
    
    bp <- barplot(ratios, 
                  col = col_and_name$randomization_color[[randomization_type]],
                  
                  ylim = c(min(ratios - mads)- 0.4*dist(c(min(ratios + mads), max(ratios + mads))),
                           max(ratios + mads)+ 0.4*dist(c(min(ratios + mads), max(ratios + mads)))),
                  # 0.4),
                  main = paste0(col_and_name$randomization_name[[randomization_type]], "\n(",
                                list(within = "Purely intronic",
                                     any    = "Intron intersecting")[[overlap_type]], ")"),
                  ylab = paste("log2(obs/exp) deletions", overlap_type, "introns"))
    
    
    arrows(bp, (ratios - mads), 
           bp, (ratios + mads),  lwd = 1.25, angle = 90,
           code = 3, length = 0.025)
    
    
    # stars
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
      if(sum(pvals >= min_pval & pvals < max_pval & alternatives == "greater") > 0) {
        trues <- (pvals >= min_pval & pvals < max_pval & alternatives == "greater")
        yplus <- (ratios+ mads)[trues] + 0.02
        text(x = bp[trues], 
             y = yplus,
             rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
             xpd = T,
             srt = 90,
             adj = 0,
             cex = 2, col = "black", xpd = T) 
      }
      
      
      # negative and significant FCs(add star)
      if(sum(pvals >= min_pval & pvals < max_pval & alternatives == "less") > 0) {
        trues <- (pvals >= min_pval & pvals < max_pval & alternatives == "less")
        yminus <- (ratios - mads)[trues] - 0.02
        text(x = bp[trues], 
             y = yminus,
             rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
             xpd = T,
             srt = 90,
             adj = 1,
             cex = 2, col = "black", xpd = T) 
      }
    }
    results_tables[[randomization_type]][[overlap_type]] <- rbind(ratio = ratios,
                                                                  percent_difference_from_expected = percent_more_or_less,
                                                                  p.val = pvals)
  }
  
  # dev.off()
}

