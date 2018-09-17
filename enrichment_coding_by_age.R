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

protein_coding_RANGES$Age <- Ages[names(protein_coding_RANGES), "GeneAge"]
protein_coding_RANGES <- protein_coding_RANGES[!is.na(protein_coding_RANGES$Age)]

seqlevelsStyle(protein_coding_RANGES) <- "UCSC"
seqlevelsStyle(all_coding_RANGES)     <- "UCSC"

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


#### RANDOM.INTERVALS()
### Randomizations done using random.intervals()
## GLOBAL - random intervals

# Folder for results
folder <- "global_random_intervals"
for(CNV_set in all_CNV_sets) {
  ## RANDOM STARTS and chromosomes
  chromosomes <- read.table(paste0("global_random_intervals/output/global_random_intervals_random_chrs_", CNV_set, "_", total_permuts, ".txt"),
                            stringsAsFactors = F, header = F)
  print("chromosomes loaded")
  
  starts <- read.table(paste0("global_random_intervals/output/global_random_intervals_random_starts_", CNV_set, "_", total_permuts, ".txt"),
                       stringsAsFactors = F, header = F)
  print("starts loaded")
  
  # CNV_ranges
  load(paste0("global_random_intervals/output/global_random_intervals_original_", CNV_set, "_GRanges.RData"))
  
  print("CNV_ranges loaded")
  
  obs_genes <- unique(all_coding_RANGES$Gene[overlapsAny(all_coding_RANGES, 
                                                         CNV_ranges, type = "any" )])
  
  obs <- table(Ages[obs_genes, "GeneAge"])
  

  
  # Random values
  ran <- c()
  
  while(ncol(chromosomes) > 0) {
    random_cnvs <- GRanges(seqnames = chromosomes[, 1],
                           ranges   = IRanges(start = starts[, 1],
                                              width = width(CNV_ranges)))
    
    ran_genes <- unique(all_coding_RANGES$Gene[overlapsAny(all_coding_RANGES, 
                                                           random_cnvs, type = "any" )])
    
    ran <- rbind(ran,
                 table(Ages[ran_genes, "GeneAge"]))
    
    
    chromosomes <- chromosomes[, -1]
    starts      <- starts[, -1]
    
    if(class(starts) == "integer") {
      chromosomes <- matrix(chromosomes, ncol = 1)
      starts      <- matrix(starts, ncol = 1)
    }
    if(ncol(chromosomes)%%100 == 0) print(ncol(chromosomes))
    # print(ncol(chromosomes))
  }
  
  results <- list()
  
  for(x in 1:length(obs)) {
    results[[x]] <- calculate_pval_randomizations(obs[x],
                                                  ran[, x])
  }
  
  names(results) <- levels(Ages$GeneAge)
  
  save(results, file = paste0("enrichment_coding/", folder, "/by_age_losses_coding", CNV_set, "_", total_permuts, "permuts.RData"))
   
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
  
  obs_genes <- unique(all_coding_RANGES$Gene[overlapsAny(all_coding_RANGES, 
                                                         CNV_ranges, type = "any" )])
  
  obs <- table(Ages[obs_genes, "GeneAge"])
  
  
  
  # Random values
  ran <- c()
  while(ncol(all_starts) > 0) {
    random_cnvs <- GRanges(seqnames = seqnames(CNV_ranges),
                           ranges   = IRanges(start = all_starts[, 1],
                                              width = width(CNV_ranges)))
    
    ran_genes <- unique(all_coding_RANGES$Gene[overlapsAny(all_coding_RANGES, 
                                                           random_cnvs, type = "any" )])
    
    ran <- rbind(ran,
                 table(Ages[ran_genes, "GeneAge"]))
    
    
    all_starts      <- all_starts[, -1]
    
    if(class(all_starts) == "integer") {
      all_starts      <- matrix(all_starts, ncol = 1)
    }
    if(ncol(all_starts)%%500 == 0) print(ncol(all_starts))
    
  }
  
  results <- list()
  
  for(x in 1:length(obs)) {
    results[[x]] <- calculate_pval_randomizations(obs[x],
                                                  ran[, x])
  }
  
  names(results) <- levels(Ages$GeneAge)
  
  save(results, file = paste0("enrichment_coding/", folder, "/by_age_losses_coding", CNV_set, "_", total_permuts, "permuts.RData"))
  
  print(CNV_set)
  
  }

print("End of local enrichment coding losses")



### COMPARISON INTRONIC VS CODING (number of genes affected by each type)

### Figures (boxplot)

col_and_name <- list(randomization_name =  c(local_random_intervals = "Local randomization",
                                             global_random_intervals = "Global randomization"),
                     randomization_color = c(local_random_intervals = "royalblue",
                                             global_random_intervals = "gold1"))

overlap_type_names <- c(intronic_within = "Purely intronic deletions",
                        # intronic_any    = "Intron-intersecting deletions", # Can be completely or partially overlapping  with introns
                        coding = "Coding-overlapping deletions")


for(enrichment_in in c("num_genes")) {
    for(randomization_type in c("local_random_intervals",
                                "global_random_intervals")) {
         
      pdf(paste0("enrichment_coding/", randomization_type, "/", randomization_type, "_intronic_vs_coding_boxplot.pdf"), 
          width = 10, height = 4)
      par(mfrow = c(1,3), oma = c(1,0.5,0,0), mar = c(7,4,4,2))
      
      for(overlap_type in c("intronic_within", "coding")) {
        for(CNV_set in all_CNV_sets[c(1,3,4)]) {
          if(overlap_type == "intronic_within") {
            # LOAD INTRONIC
            print(load(paste0("enrichment_by_age/", randomization_type, "/within_", CNV_set, "_",  total_permuts, "permuts.RData")))
            results <- results_ages
            all_obs <- unlist(lapply(results, function(x) x[[enrichment_in]][["observed"]]))
            all_ran <- do.call("cbind", lapply(results, function(x) x[[enrichment_in]][["permuted"]]))
            alternatives <- unlist(lapply(results, function(x) x[[enrichment_in]][["alternative"]]))
            all_pvals <- unlist(lapply(results, function(x) x[[enrichment_in]][["pval"]]))
            ylim <- c(-3.3,3.3)
            
            
          } 
          if(overlap_type == "coding") {
            
            # LOAD CODING
            load(paste0("enrichment_coding/", randomization_type, "/by_age_losses_coding", CNV_set, "_", total_permuts, "permuts.RData"))
            
            all_obs <- unlist(lapply(results, function(x) x[["observed"]]))
            all_ran <- do.call("cbind", lapply(results, function(x) x[["permuted"]]))
            alternatives <- unlist(lapply(results, function(x) x[["alternative"]]))
            all_pvals <- unlist(lapply(results, function(x) x[["pval"]]))
            
            ylim <- list(Abyzov_RANGES = c(-3.5,3.5),
                         Handsaker_DELS = c(-3,3.5),
                         Zarrei_DELS = c(-4.4, 4.4),
                         Sud15_DELS = c(-3,4),
                         Phase3_DELS = c(-3,3))[[CNV_set]]
          }
          
          
          
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
                        ylab = "Log2( observed / expected )",
                        col = color, names = NA,
                        cex.main = 1.35)
          abline(h = 0)
          
         
            text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
                 labels = names(results),
                 srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
            mtext("Ancient --> Recent", side = 1, line = 5.2, xpd = T, cex = 0.8)
            mtext("---", side = 1, line = 5.2, xpd = T, cex = 0.8)
         
          
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
                   cex = 1.2, col = "red", xpd = T)
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
                   cex = 1.2, col = "red", xpd = T) 
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




## Figures (barplot % of genes)
for(enrichment_in in c("num_genes")) {
  for(randomization_type in c("local_random_intervals",
                              "global_random_intervals")) {
    
    pdf(paste0("enrichment_coding/", randomization_type, "/", randomization_type, "_intronic_vs_coding_barplot.pdf"),
        width = 10, height = 4)
    
    par(mfrow = c(1,3), oma = c(1,0.5,0,0), mar = c(7,4,4,1))
    
    for(overlap_type in c("intronic_within", "coding")) {
      for(CNV_set in all_CNV_sets[c(1,3,4)]) {
        if(overlap_type == "intronic_within") {
          # LOAD INTRONIC
          print(load(paste0("enrichment_by_age/", randomization_type, "/within_", CNV_set, "_",  total_permuts, "permuts.RData")))
          results <- results_ages
          all_obs <- unlist(lapply(results, function(x) x[[enrichment_in]][["observed"]]))
          all_ran <- do.call("cbind", lapply(results, function(x) x[[enrichment_in]][["permuted"]]))
          alternatives <- unlist(lapply(results, function(x) x[[enrichment_in]][["alternative"]]))
          all_pvals <- unlist(lapply(results, function(x) x[[enrichment_in]][["pval"]]))
          ylim <- c(-3.3,3.3)
          
          
        } 
        if(overlap_type == "coding") {
          
          # LOAD CODING
          load(paste0("enrichment_coding/", randomization_type, "/by_age_losses_coding", CNV_set, "_", total_permuts, "permuts.RData"))
          
          all_obs <- unlist(lapply(results, function(x) x[["observed"]]))
          all_ran <- do.call("cbind", lapply(results, function(x) x[["permuted"]]))
          alternatives <- unlist(lapply(results, function(x) x[["alternative"]]))
          all_pvals <- unlist(lapply(results, function(x) x[["pval"]]))
          
          ylim <- list(Abyzov_RANGES = c(-3.5,3.5),
                       Handsaker_DELS = c(-3,3.5),
                       Zarrei_DELS = c(-4.4, 4.4),
                       Sud15_DELS = c(-3,4),
                       Phase3_DELS = c(-3,3))[[CNV_set]]
        }
        
        
        
        color <- c(Abyzov_RANGES = "#FF9400",
                   Handsaker_DELS = "#93184E",
                   Zarrei_DELS = "#66A1D2",
                   Sud15_DELS = "#116A66",
                   Phase3_DELS = "#043C6B")[CNV_set]
        
        
        
 
        
        table_ages <- table(Ages[rownames(Ages) %in% names(protein_coding_RANGES), "GeneAge"])
       
        obs_val    <- all_obs/table_ages*100 
        median_ran <- apply(all_ran, 2, median)/table_ages*100
       
        all_ran    <- t(apply(all_ran, 1, function(x) x/table_ages*100))
        
        
        mads <- c()
        for(i in 1:length(obs_val)) {
          mads <- c(mads, mad(log2(obs_val[i]/all_ran[,i])))
        }
        
        
        ## FOLD CHANGE
        mads_no_NA <- mads
        mads_no_NA[is.na(mads_no_NA)] <- 0
        bp <- barplot(obs_val,
                      ylim = c(0, max(obs_val+mads_no_NA)*1.1),
                      beside = T, las = 1,
                      main = paste(CNV_set, " - ", randomization_type,
                                   "\n", overlap_type_names[overlap_type], "\n"), 
                      col = color,
                      border = NA,
                      names.arg = NA,
                      # legend = T, 
                      ylab = paste("% genes"),
                      las = 2,
                      args.legend = list(x = "bottom", bty = "n", horiz = F, cex = 1.3),
                      xpd = F)
        
        # lines(x = bp, median_ran, col = "white", lwd = 4)
        lines(x = bp, median_ran, col = "gray60", lwd = 2)
        
        text(x = bp, y =- max(obs_val+mads_no_NA)*0.05, 
             labels = names(obs_val),
             srt = 45,xpd = TRUE, adj = 1, cex = 1)
        
        
        arrows(bp, 
               obs_val,
               bp,
               obs_val + mads,
               xpd = T,
               lwd = 1.25, angle = 90,
               code = 2, length = 0.025)
       
        abline(h = 0)
        
        
        # text(x = 1:length(FCs), y = ylim[1]-0.1*(dist(ylim)), 
        #      labels = names(results),
        #      srt = 45,xpd = TRUE, adj = 1, cex = 0.85)
        # mtext("Ancient --> Recent", side = 1, line = 5.2, xpd = T, cex = 0.8)
        # mtext("---", side = 1, line = 5.2, xpd = T, cex = 0.8)
        # 
        
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
            yplus <- obs_val[trues] + 0.5
            text(x = bp[trues], 
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
            yminus <- obs_val[trues] +0.5
            text(x = bp[trues], 
                 y = yminus,
                 rep(paste0(rep("*", num_stars), collapse = ""),sum(trues)), 
                 xpd = T,
                 srt = 90,
                 adj = 0,
                 # pos = 3,
                 cex = 1.2, col = "black", xpd = T) # indianred1
          }
        }
        
      }
    }
   
    
    dev.off()
  }
}

