library(ggplot2)
library(tidyverse)
library(patchwork)


# SET UP ####
# palette
cols <- c("I" = "#e40303", "II" = "#ff7900",
          "III" = "#ffda00", "IV" = "#04bf3b",
          "V" = "#004cff", "X" = "#a100c2")

# read in oxfrord table
oxford_table <- read.table("nxCaeMono1_vs_nxCaeAuri1.ce_oxford.scu.tsv", col.names=c("buscoID", "ce", "chr1", "position1", "chr2", "position2", "scu1", "scu2"))

# read in elimination sites
cm_elim_sites <- read.table("nxCaeMono1_elimination_sites.tsv", header=T)
ca_elim_sites <- read.table("nxCaeAuri1_elimination_sites.tsv", header=T)

# CAURI ####
## VIOLIN ##
# plot box plots
plot_boxplot <- function(chr) { 
  oxford_table %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot() +
    geom_boxplot(aes(x=scu1, y=position2 / 1e6), fill="#0092c1", alpha=0.5, outliers = FALSE) + 
    theme_bw() + 
    facet_grid(~chr2, scales="free_x") + 
    xlab("Somatic chromosome") + ylab("Position in C. monodelphis (Mb)") +
    theme(axis.title = element_blank())
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_boxplot)

# combine plots 
ca_violin_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect", nrow=1) 

## AVG POSITION ##
# plot avg position
plot_pos <- function(chr) {
  oxford_table %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot(aes(x = position1 / 1e6, y = position2 / 1e6)) +
    geom_smooth(span=0.1, method="loess") + 
    geom_vline(data = ca_elim_sites %>% filter(chromosome == chr), aes(xintercept = site / 1e6), linetype = 2) +
    #facet_grid(chr2 ~ chr1, scales = "free", space = "free") +
    scale_color_manual(values = cols, name = "Nigon\nelement") +
    xlab(expression("Position in " * italic("C. auriculariae") * " (Mb)")) +
    ylab(expression("Position in " * italic("C. monodelphis") * " (Mb)")) + 
    theme_bw() + 
    theme(legend.position="none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) 
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_pos)

# combine plots 
ca_pos_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect", nrow=1) 

## OXFORD ##
# plot oxford
plot_chromosome <- function(chr) {
  # Filter for chromosome once
  df <- oxford_table %>%
    filter(chr1 == chr & chr2 == chr)
  
  ggplot() +
    geom_point(data = df,
               aes(x = position1 / 1e6, y = position2 / 1e6, color = ce),
               alpha = 0.5) +
    geom_vline(data = ca_elim_sites %>% filter(chromosome == chr), aes(xintercept = site / 1e6), linetype = 2) +
    #facet_grid(chr2 ~ chr1, scales = "free", space = "free") +
    scale_color_manual(values = cols, name = "Nigon\nelement") +
    xlab(expression("Position in " * italic("C. auriculariae") * " (Mb)")) + 
    theme_bw() + 
    theme(legend.position="none",
          axis.title.y = element_blank()) 
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_chromosome)

# combine plots 
ca_ox_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect", nrow=1) 

# CMONO ####
## VIOLIN ##
# plot violin plots
plot_violin <- function(chr) { 
  oxford_table %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot() +
    geom_boxplot(aes(x=scu2, y=position1 / 1e6), fill="#0092c1", alpha=0.5, outliers = FALSE) + 
    theme_bw() + 
    facet_grid(~chr2, scales="free_x") + 
    xlab("Somatic chromosome") + ylab("Position in C. auriculariae (Mb)") +
    theme(axis.title = element_blank())
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_violin)

# combine plots 
cm_violin_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect", nrow=1) 

## AVG POSITION ##
# plot avg position
plot_pos <- function(chr) {
  oxford_table %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot(aes(x = position2 / 1e6, y = position1 / 1e6)) +
    geom_smooth(span=0.1, method="loess") + 
    geom_vline(data = cm_elim_sites %>% filter(chromosome == chr), aes(xintercept = site / 1e6), linetype = 2) +
    #facet_grid(chr1 ~ chr2, scales = "free", space = "free") +
    scale_color_manual(values = cols, name = "Nigon\nelement") +
    xlab(expression("Position in " * italic("C. monodelphis") * " (Mb)")) +
    ylab(expression("Position in " * italic("C. auriculariae") * " (Mb)")) + 
    theme_bw() + 
    theme(legend.position="none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) 
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_pos)

# combine plots 
cm_pos_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect", nrow=1) 

## OXFORD ##
# plot oxford
plot_chromosome <- function(chr) {
  # Filter for chromosome once
  df <- oxford_table %>%
    filter(chr1 == chr & chr2 == chr)
  
  ggplot() +
    geom_point(data = df,
             aes(x = position2 / 1e6, y = position1 / 1e6, color = ce),
             alpha = 0.5) +
    geom_vline(data = cm_elim_sites %>% filter(chromosome == chr), aes(xintercept = site / 1e6), linetype = 2) +
    #facet_grid(chr1 ~ chr2, scales = "free", space = "free") +
    scale_color_manual(values = cols, name = "Nigon\nelement") +
    xlab(expression("Position in " * italic("C. monodelphis") * " (Mb)")) +
    ylab(expression("Position in " * italic("C. auriculariae") * " (Mb)")) + 
    theme_bw() + 
    theme(legend.position="none",
          axis.title.y = element_blank()) 
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_chromosome)

# combine plots 
cm_ox_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect", nrow=1) 
# PLOT ####
p <- ca_violin_plot / ca_pos_plot / ca_ox_plot / cm_violin_plot / cm_pos_plot / cm_ox_plot

ggsave("box_pos_ox.pdf", plot = p, width=12, height=8, units="in")
ggsave("box_pos_ox.png", plot = p, width=12, height=8, units="in")
