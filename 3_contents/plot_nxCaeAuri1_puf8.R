#!/usr/bin/env Rscript
library(ggplot2)
library(patchwork)
library(ggrepel)

source("process_nxCaeAuri1_bed_files.R")

# set limits 
plot_start <- 8462000/1e6
plot_end <- 8623000/1e6

# plot hifi cov
hifi_p <- ggplot() + 
  geom_rect(data = hifi_cov_bed, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0, ymax=cov), fill = "#c9c9c9") + 
  geom_rect(data = hifi_delfies_bed, aes(xmin=(start/1e6)-0.001, xmax=(stop/1e6)+0.001, ymin=0, ymax=score), fill = "black") + 
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits= c(0, 200)) + 
  xlab("Position (Mb)") + 
  ylab("Depth") + 
  coord_cartesian(xlim=c(plot_start, plot_end)) + 
  theme(legend.position = "None", 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0.2,0), "cm"))

# plot hic cov
hic_p <- ggplot() + 
  geom_col(data = hic_cov_bed, aes(x = start / 1e6, y = cov), color = "#c9c9c9") + 
  geom_rect(data = hic_delfies_bed, aes(xmin=(start/1e6)-0.001, xmax=(stop/1e6)+0.001, ymin=0, ymax=score), fill = "black") + 
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits= c(0, 180)) + 
  xlab("Position (Mb)") + 
  ylab("Depth") + 
  coord_cartesian(xlim=c(plot_start, plot_end)) +
  theme(legend.position = "None", 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))

# seed 
set.seed(121)

# plot genes 
gene_p <- ggplot(data=genes_bed) + 
  geom_rect(aes(xmin=start/1e6, xmax=stop/1e6, ymin=5, ymax=9),fill="lightgrey", colour="black", linewidth=0.3) + 
  #geom_text(aes(x = (start + stop) / (2 * 1e6), y = 2, label = eleg_string), size = 3) + 
  geom_text_repel(aes(x = (start + stop) / (2 * 1e6), y = 5, label = eleg_string), 
                  size = 3, 
                  nudge_y           = -3,
                  force=10, 
                  direction         = "x",
                  segment.size      = 0.3,
                  segment.curvature = 0) + 
  #geom_segment(data = genes_bed %>% filter(!is.na(eleg_string) & eleg_string != ""), aes(x = (start + stop) / (2 * 1e6), xend = (start + stop) / (2 * 1e6), y = 3, yend = 4.5), colour = "black", linewidth=0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10)) +
  theme_bw() +
  xlab("Position (Mb)") +
  coord_cartesian(xlim=c(plot_start, plot_end)) + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))

p <- hifi_p / hic_p / gene_p + plot_layout(heights = c(10, 10, 2))

ggsave("V_middle_hifi_hic_genes.png", p, height=180, width=150, units="mm")
ggsave("V_middle_hifi_hic_genes.pdf", p, height=180, width=150, units="mm")
