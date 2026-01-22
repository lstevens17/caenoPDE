#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
library(patchwork)

# Function to generate patchwork plot for a given TOLID
make_plot <- function(tolid) {
  # file paths
  cov_bed_file    <- paste0(tolid, "_orthosite.per-base.bed")
  delfies_bed_file <- paste0(tolid, "_orthosite.delfies.bed")
  genes_bed_file   <- paste0(tolid, "_genes.bed")
  
  # Load data
  cov_bed     <- read.table(cov_bed_file, col.names = c("chr", "start", "stop", "cov"))
  delfies_bed <- read.table(delfies_bed_file, col.names = c("chr", "start", "stop", "name", "score", "strand"))
  genes_bed   <- read.table(genes_bed_file, col.names = c("chr", "start", "stop", "name", "score", "strand"))
  
  # Add delfies score column
  cov_bed$delfies <- 0
  match_idx <- match(cov_bed$start, delfies_bed$start)
  cov_bed$delfies[!is.na(match_idx)] <- delfies_bed$score[match_idx[!is.na(match_idx)]]
  
  # bin to reduce size of plot
  bin_size <- 10
  
  cov_bed_binned <- cov_bed %>%
    mutate(bin = start %/% bin_size) %>%
    group_by(bin) %>%
    summarise(
      start = min(start),
      stop = max(stop),
      cov = mean(cov),
      delfies = max(delfies)
    )
  
  # Coverage + delfies plot
  p1 <- ggplot(data = cov_bed_binned, aes(x = start / 1e6)) + 
    geom_col(aes(y = cov, colour = "cov"), fill = "#c9c9c9") + 
    geom_line(aes(y = delfies, colour = "delfies"), linewidth = 2) +
    theme_bw() + 
    scale_colour_manual(
      values = c("cov" = "#c9c9c9", "delfies" = "black"),
      labels = c("All reads", "Tel-rep")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) + 
    xlab("Position (Mb)") + 
    ylab("Depth") + 
    theme(legend.position = "None",
          plot.margin = unit(c(0,0,0,0), "cm"), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  # Define xlim for gene track
  plot_start <- min(cov_bed$start) / 1e6
  plot_end   <- max(cov_bed$stop)  / 1e6
  
  # Gene track plot
  p2 <- ggplot(data=genes_bed) +
    geom_rect(aes(xmin=start/1e6, xmax=stop/1e6, ymin=2, ymax=5),
              fill="lightgrey", colour="black", linewidth=0.3) +
    geom_text(aes(x=(start+stop)/2/1e6, y=1, label=name)) + 
    geom_text(aes(x=(start+stop)/2/1e6, y=3.5, label=strand)) + 
    scale_x_continuous(expand = c(0, 0), limits=c(plot_start, plot_end)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,6)) +
    theme_bw() +
    xlab("Position (Mb)") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  # Return patchwork
  p1 / p2 + plot_layout(heights = c(4,1))
}

# List of TOLIDs
tolids <- c("nxCaeAuri1", "nxCaeMono1", "nxCaeParv1")

# Build plots and store in a named list
plots <- lapply(tolids, make_plot)
names(plots) <- tolids

# Optionally assign to objects in global env:
# for (tolid in tolids) assign(tolid, plots[[tolid]])

final_plot <- wrap_plots(plots, nrow = 3) 

ggsave("orthologous_site.pdf", final_plot, height=12, width=8, units="in")