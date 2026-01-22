#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)

#### II_right cov ####
# Assign arguments to variables
bed_file_1 <- "variable_sites/II_right.per-base.bed"  # First BED file (region.per-base.bed)
bed_file_2 <- "variable_sites/II_right.delfies.bed"  # Second BED file (region.delfies.bed)

# Load the BED files into data frames
cov_bed <- read.table(bed_file_1, col.names = c("chr", "start", "stop", "cov"))
delfies_bed <- read.table(bed_file_2, col.names = c("chr", "start", "stop", "name", "score", "strand"))

# Initialize a new score column in cov_bed with zeros
cov_bed$delfies <- 0

# Convert the start position of delfies_bed to a vector for matching
delfies_positions <- delfies_bed$start

# Match positions and update the new_score
for (i in 1:nrow(cov_bed)) {
  pos <- cov_bed$start[i]
  match_index <- which(delfies_positions == pos)
  if (length(match_index) > 0) {
    cov_bed$delfies[i] <- delfies_bed$score[match_index]
  }
}

# Plot
II_right_cov_p <- ggplot(data = cov_bed, aes(x = start / 1e6)) +
  geom_col(aes(y = cov, colour = "cov"), fill = "#c9c9c9") +
  geom_line(aes(y = delfies, colour = "delfies"), linewidth = 2) +
  theme_bw() +
  scale_colour_manual(
    values = c("cov" = "#c9c9c9", "delfies" = "black"),
    labels = c("All reads", "Tel-rep")
  ) +
  scale_x_continuous(expand = c(0, 0), limits=c(21.614, 21.618)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Position (Mb)") +
  ylab("Depth") +
  theme(legend.position = "None", 
        plot.margin = unit(c(0,0,0,0), "cm")) 

#### IV_right cov ####
# Assign arguments to variables
bed_file_1 <- "variable_sites/IV_right.per-base.bed"  # First BED file (region.per-base.bed)
bed_file_2 <- "variable_sites/IV_right.delfies.bed"  # Second BED file (region.delfies.bed)

# Load the BED files into data frames
cov_bed <- read.table(bed_file_1, col.names = c("chr", "start", "stop", "cov"))
delfies_bed <- read.table(bed_file_2, col.names = c("chr", "start", "stop", "name", "score", "strand"))

# Initialize a new score column in cov_bed with zeros
cov_bed$delfies <- 0

# Convert the start position of delfies_bed to a vector for matching
delfies_positions <- delfies_bed$start

# Match positions and update the new_score
for (i in 1:nrow(cov_bed)) {
  pos <- cov_bed$start[i]
  match_index <- which(delfies_positions == pos)
  if (length(match_index) > 0) {
    cov_bed$delfies[i] <- delfies_bed$score[match_index]
  }
}

# plot
IV_right_cov_p <- ggplot(data = cov_bed, aes(x = start / 1e6)) +
  geom_col(aes(y = cov, colour = "cov"), fill = "#c9c9c9") +
  geom_line(aes(y = delfies, colour = "delfies"), linewidth = 2) +
  theme_bw() +
  scale_colour_manual(
    values = c("cov" = "#c9c9c9", "delfies" = "black"),
    labels = c("All reads", "Tel-rep")
  ) +
  scale_x_continuous(expand = c(0, 0), limits=c(20.312, 20.315)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Position (Mb)") +
  ylab("Depth") +
  theme(legend.position = "None", 
        plot.margin = unit(c(0,0,0,0), "cm")) 

#### II_right dot ####
# read in bed file
dot_bed <- read.table("variable_sites/II:21604901-21624901.bed", header=T)

# mirror the values
dot_bed_mirrored <- dot_bed %>%
  mutate(
    temp_start = query_start,
    temp_end = query_end,
    query_start = reference_start,
    query_end = reference_end,
    reference_start = temp_start,
    reference_end = temp_end
  ) %>%
  select(-temp_start, -temp_end)

# combine into one dataframe
dot_bed_combined <- rbind(dot_bed, dot_bed_mirrored)

# calculate window size
window <- head(dot_bed_combined$query_end, 1)

II_right_dot <- ggplot(dot_bed_combined, aes(x=(query_start+21604901)/1e6, y=(reference_start+21604901)/1e6, fill=perID_by_events, height=window/1e6, width=window/1e6)) +
  geom_tile() +
  theme_bw() +
  xlab("Position (Mb)") + ylab("Position (Mb)") +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  geom_vline(xintercept=21.614, linetype=2) + geom_vline(xintercept=21.618, linetype=2) + 
  geom_hline(yintercept=21.614, linetype=2) + geom_hline(yintercept=21.618, linetype=2) + 
  scale_fill_gradient(low = "#b5d6e8", high = "#046cb3", name = "Percent\nidentity (%)", limits = c(90,100), oob=squish) +
  theme(legend.position = "none", 
        plot.margin = unit(c(0,0,0,0), "cm"))

#### IV_right dot ####
# read in bed file
dot_bed <- read.table("variable_sites/IV:20308549-20328549.bed", header=T)

# mirror the values
dot_bed_mirrored <- dot_bed %>%
  mutate(
    temp_start = query_start,
    temp_end = query_end,
    query_start = reference_start,
    query_end = reference_end,
    reference_start = temp_start,
    reference_end = temp_end
  ) %>%
  select(-temp_start, -temp_end)

# combine into one dataframe
dot_bed_combined <- rbind(dot_bed, dot_bed_mirrored)

# calculate window size
window <- head(dot_bed_combined$query_end, 1)

IV_right_dot <- ggplot(dot_bed_combined, aes(x=(query_start+20308549)/1e6, y=(reference_start+20308549)/1e6, fill=perID_by_events, height=window/1e6, width=window/1e6)) +
  geom_tile() +
  theme_bw() +
  xlab("Position (Mb)") + ylab("Position (Mb)") +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  geom_vline(xintercept=20.312, linetype=2) + geom_vline(xintercept=20.315, linetype=2) + 
  geom_hline(yintercept=20.312, linetype=2) + geom_hline(yintercept=20.315, linetype=2) + 
  scale_fill_gradient(low = "#b5d6e8", high = "#046cb3", name = "Percent\nidentity (%)", limits = c(90,100), oob=squish) +
  theme(legend.position = "none", 
        plot.margin = unit(c(0,0,0,0), "cm"))

#### PLOT #####
p <- (II_right_cov_p | IV_right_cov_p) /
  (II_right_dot | IV_right_dot) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18, face="bold"))

ggsave("variable_sites.pdf", p, height=8, width=9, units="in")