library(ggplot2)
library(tidyverse)
library(patchwork)


# SET UP ####
# palette 
cols <- c("I" = "#e40303", "II" = "#ff7900",
          "III" = "#ffda00", "IV" = "#04bf3b",
          "V" = "#004cff", "X" = "#a100c2")

# read in oxfrord table
oxford_table <- read.table("sulstoni_vs_uteleia.ce_oxford.tsv", col.names=c("buscoID", "ce", "chr1", "position1", "chr2", "position2"))

# rename chromosomes
oxford_table <- oxford_table %>%
  mutate(chr1 = recode(chr1,
                              "OZ016475.1" = "I",
                              "OZ016476.1" = "II",
                              "OZ016477.1" = "III",
                              "OZ016478.1" = "IV",
                              "OZ016479.1" = "V",
                              "OZ016480.1" = "X"))

oxford_table <- oxford_table %>%
  mutate(chr2 = recode(chr2,
                              "OY752133.1" = "I",
                              "OY752134.1" = "II",
                              "OY752135.1" = "III",
                              "OY752136.1" = "IV",
                              "OY752137.1" = "V",
                              "OY752138.1" = "X"))

#### CSULS ####

## AVG POSITION ##
# plot avg position
plot_pos <- function(chr) {
  oxford_table %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot(aes(x = position1 / 1e6, y = position2 / 1e6)) +
    facet_grid(~chr2) + 
    geom_smooth(span=0.1, method="loess") +
    scale_color_manual(values = cols, name = "Nigon\nelement") +
    scale_x_continuous(breaks = scales::breaks_width(5)) +
    scale_y_continuous(breaks = scales::breaks_width(5)) +
    xlab(expression("Position in " * italic("C. sulstoni") * " (Mb)")) +
    ylab(expression("Position in " * italic("C. uteleia") * " (Mb)")) +
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
cs_pos_plot <- wrap_plots(plots) +
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
    scale_color_manual(values = cols, name = "Nigon\nelement") +
    scale_x_continuous(breaks = scales::breaks_width(5)) +
    scale_y_continuous(breaks = scales::breaks_width(5)) +
    xlab(expression("Position in " * italic("C. sulstoni") * " (Mb)")) +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y = element_blank())
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_chromosome)

# combine plots
cs_ox_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect", nrow=1)

#### CUTEL ####

## AVG POSITION ##
# plot avg position
plot_pos <- function(chr) {
  oxford_table %>%
    filter(chr1 == chr & chr2 == chr) %>%
    ggplot(aes(x = position2 / 1e6, y = position1 / 1e6)) +
    facet_grid(~chr2) + 
    geom_smooth(span=0.1, method="loess") +
    scale_color_manual(values = cols, name = "Nigon\nelement") +
    scale_x_continuous(breaks = scales::breaks_width(5)) +
    scale_y_continuous(breaks = scales::breaks_width(5)) +
    xlab(expression("Position in " * italic("C. uteleia") * " (Mb)")) +
    ylab(expression("Position in " * italic("C. sulstoni") * " (Mb)")) +
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
cu_pos_plot <- wrap_plots(plots) +
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
    scale_color_manual(values = cols, name = "Nigon\nelement") +
    scale_x_continuous(breaks = scales::breaks_width(5)) +
    scale_y_continuous(breaks = scales::breaks_width(5)) +
    xlab(expression("Position in " * italic("C. uteleia") * " (Mb)")) +
    ylab(expression("Position in " * italic("C. sulstoni") * " (Mb)")) +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y = element_blank())
}

# generate the individual plots
chromosomes <- c("I", "II", "III", "IV", "V", "X")
plots <- lapply(chromosomes, plot_chromosome)

# combine plots
cu_ox_plot <- wrap_plots(plots) +
  plot_layout(axis_titles = "collect", nrow=1)
# PLOT ####
p <- cs_pos_plot / cs_ox_plot / cu_pos_plot / cu_ox_plot

# save
ggsave("cs_cu_box_pos_ox.png", plot = p, width=12, height=8, units="in")
ggsave("cs_cu_box_pos_ox.pdf", plot = p, width=12, height=8, units="in")
