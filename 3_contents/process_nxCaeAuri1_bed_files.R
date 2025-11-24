library(tidyverse)

# get files
hifi_cov_bed_file <- "nxCaeAuri1_V_middle.hifi.per-base.bed"
hifi_delfies_bed_file <- "nxCaeAuri1_V_middle.hifi.delfies.bed"
hic_cov_bed_file <- "nxCaeAuri1_V_middle.hic.per-base.bed"
hic_delfies_bed_file <- "nxCaeAuri1_V_middle.hic.delfies.bed"
genes_bed_file <- "nxCaeAuri1_V_middle.genes.bed"

# load BED files into data frames
hifi_cov_bed <- read.table(hifi_cov_bed_file, col.names = c("chr", "start", "stop", "cov"))
hifi_delfies_bed <- read.table(hifi_delfies_bed_file, col.names = c("chr", "start", "stop", "name", "score", "strand"))
hic_cov_bed <- read.table(hic_cov_bed_file, col.names = c("chr", "start", "stop", "cov"))
hic_delfies_bed <- read.table(hic_delfies_bed_file, col.names = c("chr", "start", "stop", "name", "score", "strand"))
genes_bed <- read.table(genes_bed_file, col.names = c("chr", "start", "stop", "eleg_string"), fill = TRUE)

