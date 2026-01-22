library(tidyverse)

# get files
hic_cov_bed_file <- "nxCaeParv1_III_middle.hic.per-base.bed"
hic_delfies_bed_file <- "nxCaeParv1_III_middle.hic.delfies.bed"
genes_bed_file <- "nxCaeParv1_III_middle.genes.bed"

# load BED files into data frames
hic_cov_bed <- read.table(hic_cov_bed_file, col.names = c("chr", "start", "stop", "cov"))
hic_delfies_bed <- read.table(hic_delfies_bed_file, col.names = c("chr", "start", "stop", "name", "score", "strand"))
genes_bed <- read.table(genes_bed_file, col.names = c("chr", "start", "stop", "eleg_string"), fill = TRUE)
