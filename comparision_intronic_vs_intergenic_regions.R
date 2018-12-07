## Comparison intronic - intergenic

#### Loading libraries
library( "regioneR" )
library( "BSgenome.Hsapiens.UCSC.hg19.masked" )
library( "beanplot")
library("gtools")

setwd("Background-models/")
if (!file.exists("comparison_intronic_intergenic_results")) dir.create(file.path(".", "comparison_intronic_intergenic_results"))

# Load and prepare data
load("data_for_comparison_intron_intergenic.RData")
seqlevelsStyle(intronic_RANGES) <- "UCSC"
seqlevelsStyle(protein_coding_RANGES)      <- "UCSC"
seqlevelsStyle(genes_without_UTR_RANGES) <- "UCSC"
seqlevelsStyle(Phase3_DELS) <- "UCSC"
source("calculate_pvalue_randomizations.R")

# Get genome and mask
human.genome    <- getGenomeAndMask(genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$genome
human.autosomal <- filterChromosomes(human.genome, organism="hg", chr.type="autosomal")

# Filter genome and mask, only autosomes
human.mask <- getGenomeAndMask(genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$mask
human.mask.autosomal <- filterChromosomes(human.mask, organism="hg", chr.type="autosomal")

# Mask with protein-coding genes (for randomizations excluding overlaps with genes)
protein_coding_RANGES <- filterChromosomes(protein_coding_RANGES, organism="hg", chr.type="autosomal")
masked.genes.and.gaps <- reduce(append(human.mask.autosomal, protein_coding_RANGES))

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


## Now with genes without UTR
genome_and_genes <- coverage(append(human.genome, reduce(genes_without_UTR_RANGES)))

coverages_RANGES <- GRanges()
for(i in paste0("chr", 1:22)) {
  coverage_Rle <- genome_and_genes[[i]]
  chr_RANGES <- GRanges(seqnames = as.character(i), ranges= ranges(coverage_Rle), coverage = coverage_Rle@values)
  coverages_RANGES <- append(coverages_RANGES, chr_RANGES)
}

intergenic_RANGES_no_UTR <- coverages_RANGES[coverages_RANGES$coverage == 1]

# Remove low mappability 
intergenic_RANGES_no_UTR <- intergenic_RANGES_no_UTR[!overlapsAny(intergenic_RANGES_no_UTR, human.mask.autosomal)]

# Make sure I didn't miss any exons

all_exons <- append(ppal_exons_RANGES, alt_exons_RANGES)
seqlevelsStyle(all_exons) <- "UCSC"
all_exons <- reduce(all_exons)

intergenic_RANGES_no_UTR <- intergenic_RANGES_no_UTR[!overlapsAny(intergenic_RANGES_no_UTR, all_exons)]


# Compare sizes
boxplot(list(intronic = width(intronic_RANGES),
             intergenic = width(intergenic_RANGES_no_UTR)), 
        ylab = "size (bp)",
        main = "Size of intronic and intergenic regions",
        outline = F)

table(cut(width(intronic_RANGES), 
          breaks = quantile(width(intronic_RANGES), seq(0, 1, by = 0.1))))



intergenic_regions_per_size_range <- table(cut(width(intergenic_RANGES_no_UTR), 
                                               breaks = quantile(width(intronic_RANGES), seq(0, 1, by = 0.1))))

par(mar = c(8,4,3,2))
bp <- barplot(intergenic_regions_per_size_range,
              ylab = "Number of intergenic regions",
              ylim = c(0, max(intergenic_regions_per_size_range)*1.1),
              las = 2,
              cex.names = 0.8,
              cex.main = 1,
              main = "Intergenic regions in each intronic region decile\n~19000 intronic regions per decile")

text(x = bp, y = intergenic_regions_per_size_range, labels = intergenic_regions_per_size_range,
     pos = 3, cex = 0.8)

## I take a sample of intronic regions, find intergenic regions  of same/similar length and check how many deletions within



sizes <- rbind(over2000bp = c(2000, max(width(intronic_RANGES))), 
               any_size  = c(50, max(width(intronic_RANGES))),
               sizes )
# save.image("comparison_intronic_intergenic/revision2/image_for_plots.RData")


# Number of permutations
total_permuts <- 100


Sys.time()

  intergenic <- intergenic_RANGES_no_UTR
  
  intron_intergenic <- c()
  deletions_sizes   <- c()
  for(perm in 1:total_permuts) {
    set.seed(123 + perm)
    # Select intronic ranges within size range
    sample_introns <- intronic_RANGES[width(intronic_RANGES) >= 50 &
                                        width(intronic_RANGES) <= max(width(intronic_RANGES))  ]
    # Sample of 700 (to start, but I will end up selecting 500 of them)
    sample_introns <- sample_introns[sample(1:length(sample_introns), 700, replace = F)]
    
    intergenic_sample <- GRanges()
    for(i in 1:length(sample_introns)) {
      min_intron_size <- width(sample_introns[i])
      max_intron_size <- width(sample_introns[i])
      # First try to find an integenic region of exactly the same size
      x <- intergenic[width(intergenic) %in% (min_intron_size:max_intron_size)]
      # If not possible, increase the size range by  20 (+10 and -10) and keep looking
      while(length(x) < 1) {
        min_intron_size <- min_intron_size-10
        max_intron_size <- max_intron_size+10
        x <- intergenic[width(intergenic) %in% (min_intron_size:max_intron_size)]
        }
      # Select one intergenic region of the size range
      intergenic_sample <- append(intergenic_sample, x[sample(1:length(x), 1)])
      }
    
    # Remove pairs with a duplicated intergenic region
    sample_introns     <- sample_introns[!duplicated(intergenic_sample)]
    intergenic_sample  <- intergenic_sample[!duplicated(intergenic_sample)]
    
    # take 500 regions 
    if(length(sample_introns) < 500) {
      intron_intergenic <- paste0("Not enough remaining intergenic regions (", length(sample_introns), ")")
      next
      } else {
        indice <- sample(1:length(sample_introns), 500, replace = F)
        sample_introns    <- sample_introns[indice]
        intergenic_sample <- intergenic_sample[indice]
       
        percent_lost_bp_intron      = sum(width(reduce(Phase3_DELS[overlapsAny(Phase3_DELS, sample_introns, type = "within")])))/sum(width(sample_introns))*100
        percent_lost_bp_intergenic  = sum(width(reduce(Phase3_DELS[overlapsAny(Phase3_DELS, intergenic_sample, type = "within")])))/sum(width(intergenic_sample))*100
          
        intron_intergenic <- rbind(intron_intergenic,
                                   c(intronic_within    = sum(overlapsAny(Phase3_DELS, sample_introns, type = "within")),
                                     intronic_any       = sum(overlapsAny(Phase3_DELS, sample_introns, type = "any")),
                                     intergenic_within  = sum(overlapsAny(Phase3_DELS, intergenic_sample, type = "within")),
                                     intergenic_any     = sum(overlapsAny(Phase3_DELS, intergenic_sample, type = "any")),
                                     total_intronic_bp  = sum(width(sample_introns)),
                                     total_intergenic_bp  = sum(width(intergenic_sample)),
                                     percent_lost_bp_intron = percent_lost_bp_intron,
                                     percent_lost_bp_intergenic = percent_lost_bp_intergenic,
                                     length             = length(sample_introns)))

        deletions_sizes <- rbind(deletions_sizes, 
                                 c(intr  = summary(width(Phase3_DELS)[overlapsAny(Phase3_DELS, sample_introns, type = "within")]),
                                   inter = summary(width(Phase3_DELS)[overlapsAny(Phase3_DELS, intergenic_sample, type = "within")])))
          
        }
        print(perm)
       }
    
    intron_intergenic <- list(size = c(50, max(width(Phase3_DELS))),
                              results = intron_intergenic)
    
    # Save data
    save(intron_intergenic, file = paste0("comparison_intronic_intergenic_results/INTRON_to_INTERGENIC_no_utr_",total_permuts, ".RData"))
    save(deletions_sizes,   file = paste0("comparison_intronic_intergenic_results/del_sizes_INTRON_to_INTERGENIC_no_utr_",total_permuts, ".RData"))
    
    
  
#### PLOTS
    
pdf("comparison_intronic_intergenic_results/Figure1_comparison_intronic_intergenic.pdf", 
    height = 7, width = 8)
par(mar = c(3,5,3,1), mfrow = c(2,3), oma = c(0,0,3,0))
layout(matrix(c(1,1,2,2,3,3,6,4,4,5,5,7), byrow = T, ncol = 6))

load(file = paste0("comparison_intronic_intergenic_results/INTRON_to_INTERGENIC_no_utr_", total_permuts, ".RData"))

lista <- intron_intergenic$results[, c("total_intronic_bp", "total_intergenic_bp", 
                                       "percent_lost_bp_intron", "percent_lost_bp_intergenic",
                                       "intronic_within", "intergenic_within")]
colnames(lista) <- c("Intronic (total bp)", "Intergenic (total bp)", 
                     "Intronic (% lost bp)", "Intergenic (% lost bp)",
                     "Intronic (num dels)", "Intergenic (num dels)")


for(i in c(5, 3, 1)) {
  j <- i+1
  if(i == 1) ylab = "Total content (Mb)" 
  if(i == 3) ylab = "% of region deleted"
  if(i == 5) ylab = "Number of deletions"
  
  bp <- boxplot(lista[, i:j], plot = F)
  
  boxplot(lista[, i:j], 
          las = 0, 
          outline = F,
          cex.axis = 1,
          names = c("Intronic", "Intergenic"),
          main = ylab,
          ylim = c(min(bp$stats), max(bp$stats)*1.1),
          ylab = ylab)
  
  tt <- t.test(lista[, i], lista[, j], paired = T)
  arrows(x0 = 1,
         x1 = 2,
         y0 = max(bp$stats)*1.05,
         y1 = max(bp$stats)*1.05,
         angle = 90, length = 0.05, code  = 3  )
  if(tt$p.value < 0.05) {
    pval <- formatC(tt$p.value, format = "e", digits = 2)
    if(pval == "0.00e+00") pval <- "<1.00e-100"
    text(x = 1.5, y = max(bp$stats)*1.1, labels =  pval,
         cex = 1)   
  } else {
    text(x = 1.5, y = max(bp$stats)*1.1, labels =  "n.s.",
         cex = 1)
  }
}





## sIZES OF DELETIONS
load( paste0("comparison_intronic_intergenic_results/del_sizes_INTRON_to_INTERGENIC_no_utr_",total_permuts, ".RData"))

deletions_sizes <- as.data.frame(deletions_sizes)

for(i in 3:4) {
  bp <- boxplot(deletions_sizes[, c(i, i+6)], plot = F)
  boxplot(deletions_sizes[, c(i, i+6)], 
          outline = F, ylim = c(min(bp$stats), max(bp$stats)*1.12),
          names = c("Intronic", "Intergenic"),
          ylab = "Deletion size",
          main = paste0("Deletion size (", c("min", "1st Qu", "median", "mean", "3rd Qu", "max")[i], ")"))
  tt <- t.test(deletions_sizes[, i], deletions_sizes[, i+6], paired = T)
  arrows(x0 = 1,
         x1 = 2,
         y0 = max(bp$stats)*1.05,
         y1 = max(bp$stats)*1.05,
         angle = 90, length = 0.05, code  = 3  )
  if(tt$p.value < 0.05) {
    text(x = 1.5, y = max(bp$stats)*1.1,  labels =  formatC(tt$p.value, format = "e", digits = 2),
         cex = 1)   
  } else {
    text(x = 1.5, y = max(bp$stats)*1.1, labels =  "n.s.",
         cex = 1)
  }
}

mtext(expression(bold("Comparison of deletion content in intronic and intergenic regions ")), 
      cex = 1.1, side = 3, outer = T, col = "gray50")

dev.off()

