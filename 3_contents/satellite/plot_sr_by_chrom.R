library(tidyverse)
library(ggplot2)
library(patchwork)

generate_trf_plot <- function(tolid, species_name) { 
  # Read the samtools faidx file
  faidx <- read_tsv(paste0(tolid, ".1.primary.fa.fai"), 
                    col_names = c("chr", "length", "offset", "linebases", "linewidth")) %>%
    filter(chr %in% c("I", "II", "III", "IV", "V", "X")) %>%
    select(chr, length) 
  
  all_clusters <- read.table(paste0(tolid, ".1.primary.fa_50kb.trf_cov.bed"), 
                             col.names=c("chr", "start", "stop", "count", "span", "size", "cov", "cluster")) %>%
    filter(chr %in% c("I", "II", "III", "IV", "V", "X"))
  
  p <- ggplot() +
    geom_rect(data = faidx,
              aes(xmin = 0, xmax = length, ymin = 0, ymax = 1),
              color = "black", fill = "#ebebeb", alpha = 1) +
    geom_col(data=all_clusters, aes(x=start, y=cov, fill=cluster)) + 
    facet_grid(chr ~ .) + 
    xlab("Position (Mb)") + ylab("Proportion") + 
    scale_x_continuous(labels=function(x)x/1e6) + 
    scale_y_continuous(limits=c(0,1), breaks = c(0.0, 0.5, 1.0)) + 
    scale_fill_discrete(name = "Tandem repeat ID") +
    coord_cartesian(ylim=c(0,1)) + 
    theme_bw() + 
    theme(
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.key.size = unit(0.8, "lines")
    ) + 
    ggtitle(bquote(italic(.(species_name))))  
  
  # Read and parse elimination sites
  elimination_sites_list <- read_tsv(paste0(tolid, "_elimination_sites.tsv"),
                                     col_names = c("chr", "site")) %>%
    split(.$chr) %>% lapply(function(x) x$site)
  
  for (chr in names(elimination_sites_list)) {
    p <- p + geom_vline(
      data = tibble(chr = chr, xintercept = elimination_sites_list[[chr]]),
      aes(xintercept = xintercept),
      color = "black", linetype = "21", linewidth = 0.75
    ) +
      scale_linetype_manual(name = "Annotations", values = c("Break site" = "21"))
  }
  
  return(p)
}

# Define mapping of ToLID to species name
species_map <- c(
  "nxCaeAuri1" = "C. auriculariae",
  "nxCaeMono1" = "C. monodelphis",
  "nxCaeParv1" = "C. parvicauda"
)

# Generate plots using the full species name
plot_list <- setNames(
  lapply(names(species_map), function(tolid) {
    generate_trf_plot(tolid, species_map[[tolid]])
  }),
  names(species_map)
)

# Combine plots
p <- plot_list[["nxCaeAuri1"]] / plot_list[["nxCaeMono1"]] / plot_list[["nxCaeParv1"]] +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 16, face="bold"))

# Save the combined plot
ggsave("sr_by_chrom.png", p, height=300, width=250, units="mm")
ggsave("sr_by_chrom.pdf", p, height=300, width=250, units="mm")
