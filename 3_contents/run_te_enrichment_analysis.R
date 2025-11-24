# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: script.R <tes.bed> <genome.bed> <elim.bed> <prefix>")
}

library(regioneR)
library(dplyr)

# Assign arguments to variables
te_bed_file <- args[1]  # TE bed
genome_bed <- args[2]  # genome bed (chrom start stop)
elim_bed <- args[3] # eliminated BED file
output_prefix <- args[4]  # output prefix 

# read TE bed
te_df <- read.table(te_bed_file, col.names=c("chr", "start", "stop", "dot", "score", "strand", "source", "family", "dot2", "col9"))

# convert genome bed (chromosome lengths) to granges
gen <- toGRanges(genome_bed) 

# convert eliminated regions to granges
elim_regions <- toGRanges(elim_bed)

# get list of TE family names
te_families <- unique(te_df$family)

# set seed
set.seed(123)

# initialize results data frame
results_list <- list()

# loop over TE families
for (family in te_families) {
  cat("Processing family:", family, "\n")
  
  # filter TEs for the current family
  te_family_df <- filter(te_df, family == !!family) %>% select(chr:stop)
  
  # convert to GRanges
  te_family <- toGRanges(te_family_df)
  
  # run permutation test
  te_enrichment <- permTest(A = elim_regions, 
                            B = te_family, 
                            genome = gen,
                            randomize.function = randomizeRegions,
                            evaluate.function = numOverlaps,
                            ntimes = 1000)
  
  # store results
  results_list[[family]] <- data.frame(
    family = family,
    z_score = te_enrichment$numOverlaps$zscore,
    p_value = te_enrichment$numOverlaps$pval,
    observed = te_enrichment$numOverlaps$observed,
    expected = mean(te_enrichment$numOverlaps$permuted)
  )
}

# add on a global analysis
te_family_df <- te_df %>% select(chr:stop)
te_family <- toGRanges(te_family_df)

te_enrichment <- permTest(A = elim_regions, 
                          B = te_family, 
                          genome = gen,
                          randomize.function = randomizeRegions,
                          evaluate.function = numOverlaps,
                          ntimes = 1000)

# store results
results_list[["All"]] <- data.frame(
  family = "All",
  z_score = te_enrichment$numOverlaps$zscore,
  p_value = te_enrichment$numOverlaps$pval,
  observed = te_enrichment$numOverlaps$observed,
  expected = mean(te_enrichment$numOverlaps$permuted))

# combine all results into a single data frame
results_df <- bind_rows(results_list)

# print and save results
output_filename <- paste0(output_prefix, "_te_family_enrichment_results.tsv")
write.table(results_df, output_filename, quote=FALSE, row.names = FALSE, sep='\t')
