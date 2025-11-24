library(ggplot2)
library(png)
library(ggrepel)
library(tidyverse)
library(patchwork)
library(scales)
library(ggplotify)

#### a) Gene orthology summary ####
# Read data
table <- read.table("gene_orthology_summary.tsv", header = TRUE, sep="\t") 

# Custom facet labels (mapping old names to new)
species_labels <- c(
  "nxCaeAuri1" = "italic('C. auriculariae')",
  "nxCaeMono1" = "italic('C. monodelphis')",
  "nxCaeParv1" = "italic('C. parvicauda')"
)

# Plot with facets for species
go_p <- table %>% 
  ggplot(aes(x = copy, y = count, fill = status)) + 
  geom_col(position = "stack") + 
  facet_wrap(~species, labeller = as_labeller(species_labels, default = label_parsed)) +  
  theme_bw() + 
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#FDB863", "#E66101")) +
  theme(
    legend.position = c(0.97, 0.95), 
    legend.justification = c(1, 1), 
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.title = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"),
  ) +
  labs(x = element_blank(), 
       y = "Count",
       fill = "Category") 


#### b) C monodelphis puf-8 ####

source("process_nxCaeMono1_bed_files.R")

#set limits 
plot_start <- 10937000/1e6
plot_end <- 10967000/1e6

# plot hifi cov
cm_hifi_p <- ggplot() + 
  geom_rect(data = hifi_cov_bed, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0, ymax=cov), fill = "#c9c9c9") + 
  geom_rect(data = hifi_delfies_bed, aes(xmin=(start/1e6)-0.0002, xmax=(stop/1e6)+0.0002, ymin=0, ymax=score), fill = "black") + 
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits= c(0, 650)) + 
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
cm_gene_p <- ggplot(data=genes_bed) + 
  geom_rect(aes(xmin=start/1e6, xmax=stop/1e6, ymin=5, ymax=9),fill="lightgrey", colour="black", linewidth=0.3) + 
  #geom_text(aes(x = (start + stop) / (2 * 1e6), y = 2, label = eleg_string), size = 3) + 
  geom_text_repel(aes(x = (start + stop) / (2 * 1e6), y = 5, label = eleg_string), 
                  size = 3, 
                  nudge_y           = -3,
                  force=10, 
                  direction         = "x",
                  segment.size      = 0.3,
                  segment.curvature = 0) + 
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

#### c) C parvicauda III ####

source("process_nxCaeParv1_bed_files.R")

# set limits 
plot_start <- 9350000/1e6
plot_end <- 9630000/1e6

# plot hic cov
cp_hic_p <- ggplot() + 
  geom_col(data = hic_cov_bed, aes(x = start / 1e6, y = cov), color = "#c9c9c9") + 
  geom_rect(data = hic_delfies_bed, aes(xmin=(start/1e6)-0.001, xmax=(stop/1e6)+0.001, ymin=0, ymax=score), fill = "black") + 
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits= c(0, 190)) + 
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
cp_gene_p <- ggplot(data=genes_bed) + 
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


#### d-e) Dot plots ####
# read in bed file
IV_middle_bed <- read.table("nxCaeAuri1_IV_6924762-7184171.bed", header=T)

# mirror the values
IV_middle_bed_mirrored <- IV_middle_bed %>% 
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
IV_middle_combined <- rbind(IV_middle_bed, IV_middle_bed_mirrored) 

# calculate window size
window <- head(IV_middle_combined$query_end, 1)

IV_middle_combined_p <- ggplot(IV_middle_combined, aes(x=query_start/1e3, y=reference_start/1e3, fill=perID_by_events, height=window/1e3, width=window/1e3)) + 
  geom_tile() + 
  theme_bw() + 
  xlab("Position (kb)") + ylab("Position (kb)") + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
  scale_fill_continuous(name="Percent identity (%)") + 
  scale_fill_gradient(low = "#b5d6e8", high = "#046cb3", name = "Percent\nidentity (%)", limits = c(90,100), oob=squish) + 
  theme(legend.position = "none")

# read in bed file
IV_right_bed <- read.table("nxCaeAuri1_IV_18405029-18607180.bed", header=T)

# mirror the values
IV_right_bed_mirrored <- IV_right_bed %>% 
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
IV_right_combined <- rbind(IV_right_bed, IV_right_bed_mirrored) 

# calculate window size
window <- head(IV_right_combined$query_end, 1)

IV_right_combined_p <- ggplot(IV_right_combined, aes(x=query_start/1e3, y=reference_start/1e3, fill=perID_by_events, height=window/1e3, width=window/1e3)) + 
  geom_tile() + 
  theme_bw() + 
  xlab("Position (kb)") + ylab("Position (kb)") + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
  scale_fill_continuous(name="Percent identity (%)") + 
  scale_fill_gradient(low = "#b5d6e8", high = "#046cb3", name = "Percent\nidentity (%)", limits = c(90,100), oob=squish) + 
  theme(legend.position = "none")

# plot as raster and the reread in 
# Save dot plot panels as high-res raster
ggsave("IV_middle_combined_p.png", IV_middle_combined_p, dpi=600, width=100, height=100, units="mm")
ggsave("IV_right_combined_p.png", IV_right_combined_p, dpi=600, width=100, height=100, units="mm")

# Read back as raster grobs
IV_middle_combined_raster <- as.raster(readPNG("IV_middle_combined_p.png"))
IV_right_combined_raster  <- as.raster(readPNG("IV_right_combined_p.png"))

# Turn into ggplot objects
IV_middle_combined_pg <- as.ggplot(~grid::grid.raster(IV_middle_combined_raster))
IV_right_combined_pg  <- as.ggplot(~grid::grid.raster(IV_right_combined_raster))

#### Assemble plot ####
p <- (go_p | (cm_hifi_p / cm_gene_p + plot_layout(heights=c(6,1)))) / cp_hic_p / cp_gene_p / (IV_middle_combined_pg | IV_right_combined_pg) + plot_layout(nrow = 4, heights=c(3,3,0.5,4))

ggsave("figure_3.png", p, height=260, width=210, units="mm")
ggsave("figure_3.pdf", p, height=260, width=210, units="mm")

# for editing later
ggsave("nxCaeParv1_gene_plot.png", cp_gene_p, height=50, width=200, units="mm")

# get legend for dotplots
IV_right_combined_p <- ggplot(IV_right_combined, aes(x=query_start/1e3, y=reference_start/1e3, fill=perID_by_events, height=window/1e3, width=window/1e3)) + 
  geom_tile() + 
  theme_bw() + 
  xlab("Position (kb)") + ylab("Position (kb)") + 
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
  scale_fill_continuous(name="Percent identity (%)") + 
  scale_fill_gradient(low = "#b5d6e8", high = "#046cb3", name = "Percent\nidentity (%)", limits = c(90,100), oob=squish) 

ggsave("dotplot_legend.pdf", IV_right_combined_p, height=50, width=200, units="mm")
