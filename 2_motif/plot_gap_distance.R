library(ggplot2)
library(tidyverse)
library(patchwork)

motif_df <- read.table("motif_details.tsv", header=T, sep='\t') 

motif_df <- motif_df %>%
  select(
    species = Species,
    distance = Distance.from.break.site..bp.,
    gap = Gap.length..bp.
  )

gap_df <- motif_df %>%
  filter(species != "C. parvicauda") %>% 
  mutate(gap = as.numeric(gap))

gap_plot <- ggplot(data=gap_df, aes(x=gap)) + 
  geom_histogram(binwidth=1) + 
  facet_grid(~species) +  
  theme_bw() + 
  xlab("Gap length (bp)") + 
  ylab("Count") + 
  coord_cartesian(xlim=c(0,70)) + 
  theme(strip.text = element_text(face = "italic"))  + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

dist_plot <- ggplot(data=motif_df, aes(x=distance)) + 
  geom_vline(xintercept=0, linetype=2) + 
  geom_histogram(binwidth=1) + 
  scale_x_reverse() +
  facet_grid(~species) +  
  theme_bw() + 
  xlab("Distance from telomere addition site (bp)") + 
  ylab("Count") + 
  coord_cartesian(xlim=c(35,-10)) + 
  theme(strip.text = element_text(face = "italic"))  + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

p <- gap_plot / dist_plot + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size=20))

ggsave("gap_length_distance.png", p, height=6, width=8, units="in")
ggsave("gap_length_distance.pdf", p, height=6, width=8, units="in")

