library(tidyverse)

# get files
hifi_cov_bed_file <- "nxCaeMono1_IV_middle_1.hifi.per-base.bed"
hifi_delfies_bed_file <- "nxCaeMono1_IV_middle_1.hifi.delfies.bed"
genes_bed_file <- "nxCaeMono1_IV_middle_1.genes.bed"

# load BED files into data frames
hifi_cov_bed <- read.table(hifi_cov_bed_file, col.names = c("chr", "start", "stop", "cov"))
hifi_delfies_bed <- read.table(hifi_delfies_bed_file, col.names = c("chr", "start", "stop", "name", "score", "strand"))
genes_bed <- read.table(genes_bed_file, col.names = c("chr", "start", "stop", "eleg_string"), fill = TRUE)

