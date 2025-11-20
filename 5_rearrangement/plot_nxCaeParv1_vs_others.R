library(ggplot2)
library(tidyverse)
library(patchwork)


#### SET UP ####
# nigon set up
cols <- c("I" = "#e40303", "II" = "#ff7900",
          "III" = "#ffda00", "IV" = "#04bf3b",
          "V" = "#004cff", "X" = "#a100c2")

# read in oxford table
cp_vs_ca_oxford_table <- read.table("nxCaeParv1_vs_nxCaeAuri1.ce_oxford.tsv", col.names=c("buscoID", "ce", "chr1", "position1", "chr2", "position2"))
cp_vs_cm_oxford_table <- read.table("nxCaeParv1_vs_nxCaeMono1.ce_oxford.tsv", col.names=c("buscoID", "ce", "chr1", "position1", "chr2", "position2"))

# read in elimination sites
ca_elim_sites <- read.table("nxCaeAuri1_elimination_sites.tsv", header=T)
cm_elim_sites <- read.table("nxCaeMono1_elimination_sites.tsv", header=T)
cp_elim_sites <- read.table("nxCaeParv1_elimination_sites.tsv", header=T)

# reorder 
cp_vs_ca_oxford_table$chr1 <- factor(cp_vs_ca_oxford_table$chr1,
                            levels = c("I", "II", "III", "IV", "V", "X"))
cp_vs_cm_oxford_table$chr2 <- factor(cp_vs_cm_oxford_table$chr2,
                            levels = c("X", "V", "IV", "III", "II", "I"))

summary(cp_vs_ca_oxford_table %>%     
    filter(ce != "unassigned") %>%
    filter(chr1 == chr2))

#### PLOT BY CHROMOSOME CP VS CA ####
plot_chromosome <- function(chr) {
  cp_vs_ca_oxford_table %>%
    filter(ce != "unassigned") %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot() +
    geom_point(aes(x = position1 / 1e6, y = position2 / 1e6, color = ce), alpha=0.75) +
    geom_vline(data = cp_elim_sites %>% filter(chromosome == chr), aes(xintercept = site / 1e6), linetype = 2) +
    geom_hline(data = ca_elim_sites %>% filter(chromosome == chr), aes(yintercept = site / 1e6), linetype = 2) +
    facet_grid(chr2 ~ chr1, scales = "free", space = "free") +
    scale_color_manual(values = cols) +
    xlab(expression("Position in " * italic("C. parvicauda") * " (Mb)")) + 
    ylab(expression("Position in " * italic("C. auriculariae") * " (Mb)")) +
    theme_bw() + 
    theme(legend.position="none") 
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_chromosome)

# combine plots 
cp_vs_ca <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect") 

#### PLOT BY CHROMOSOME CP VS CM ####
plot_chromosome <- function(chr) {
  cp_vs_cm_oxford_table %>%
    filter(ce != "unassigned") %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot() +
    geom_point(aes(x = position1 / 1e6, y = position2 / 1e6, color = ce), alpha=0.75) +
    geom_vline(data = cp_elim_sites %>% filter(chromosome == chr), aes(xintercept = site / 1e6), linetype = 2) +
    geom_hline(data = cm_elim_sites %>% filter(chromosome == chr), aes(yintercept = site / 1e6), linetype = 2) +
    facet_grid(chr2 ~ chr1, scales = "free", space = "free") +
    scale_color_manual(values = cols) +
    xlab(expression("Position in " * italic("C. parvicauda") * " (Mb)")) + 
    ylab(expression("Position in " * italic("C. monodelphis") * " (Mb)")) +
    theme_bw() + 
    theme(legend.position="none") 
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_chromosome)

# combine plots 
cp_vs_cm <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect") 

# patchwork to panel plots
p <- cp_vs_ca / cp_vs_cm

# save
ggsave("nxCaeParv1_vs_others.nigon_oxford.elim_sites.png", plot = p, width=12, height=16, units="in")
ggsave("nxCaeParv1_vs_others.nigon_oxford.elim_sites.pdf", plot = p, width=12, height=16, units="in")
