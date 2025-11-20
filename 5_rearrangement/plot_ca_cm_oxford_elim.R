library(ggplot2)
library(tidyverse)
library(patchwork)


#### SET UP ####
# palette set up
cols <- c("I" = "#e40303", "II" = "#ff7900",
          "III" = "#ffda00", "IV" = "#04bf3b",
          "V" = "#004cff", "X" = "#a100c2")


# read in oxford table
oxford_table <- read.table("nxCaeAuri1_vs_nxCaeMono1.ce_oxford.tsv", col.names=c("buscoID", "ce", "chr1", "position1", "chr2", "position2"))
summary(oxford_table %>% filter(chr1 == chr2) %>% filter(ce == "unassigned"))
summary(oxford_table %>% filter(chr1 != chr2))
summary(oxford_table %>% filter(chr1 == chr2) %>% filter(ce != "unassigned"))
oxford_table <- oxford_table %>% filter(ce != "unassigned")

# read in elimination sites
cm_elim_sites <- read.table("nxCaeMono1_elimination_sites.tsv", header=T)
ca_elim_sites <- read.table("nxCaeAuri1_elimination_sites.tsv", header=T)

# reorder 
oxford_table$chr1 <- factor(oxford_table$chr1,
                     levels = c("I", "II", "III", "IV", "V", "X"))
oxford_table$chr2 <- factor(oxford_table$chr2,
                     levels = c("X", "V", "IV", "III", "II", "I"))

# plot oxford by chrom
plot_chromosome <- function(chr) {
  oxford_table %>%
    #filter(nigon != "unassigned") %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot() +
    geom_point(aes(x = position1 / 1e6, y = position2 / 1e6, color=ce), alpha=0.75) +
    geom_vline(data = ca_elim_sites %>% filter(chromosome == chr), aes(xintercept = site / 1e6), linetype = 2) +
    geom_hline(data = cm_elim_sites %>% filter(chromosome == chr), aes(yintercept = site / 1e6), linetype = 2) +
    facet_grid(chr2 ~ chr1, scales = "free", space = "free") +
    scale_colour_manual(values = cols) + 
    xlab(expression("Position in " * italic("C. auriculariae") * " (Mb)")) +
    ylab(expression("Position in " * italic("C. monodelphis") * " (Mb)")) + 
    theme_bw() + 
    theme(legend.position="none") 
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_chromosome)

# combine plots 
combined_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect") 

# save as pdf for further editing
ggsave("nxCaeMono1_vs_nxCaeAuri1.ce_oxford.elim_sites.pdf", plot = combined_plot, width=12, height=8, units="in")
ggsave("nxCaeMono1_vs_nxCaeAuri1.ce_oxford.elim_sites.png", plot = combined_plot, width=12, height=8, units="in")
