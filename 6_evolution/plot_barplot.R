library(ggplot2)
library(tidyverse)
library(patchwork)

df <- read.table("0.25_matrix.tsv", col.names=c("species", "bs_count"))

phylo_order <- c(
  "Caenorhabditis_briggsae",
  "Caenorhabditis_nigoni",
  "Caenorhabditis_sinica",
  "Caenorhabditis_zanzibari",
  "Caenorhabditis_latens",
  "Caenorhabditis_remanei",
  "Caenorhabditis_tropicalis",
  "Caenorhabditis_wallacei",
  "Caenorhabditis_doughertyi",
  "Caenorhabditis_brenneri",
  "Caenorhabditis_elegans",
  "Caenorhabditis_inopinata",
  "Caenorhabditis_kamaaina",
  "Caenorhabditis_afra",
  "Caenorhabditis_sulstoni",
  "Caenorhabditis_sp36",
  "Caenorhabditis_macrosperma",
  "Caenorhabditis_japonica",
  "Caenorhabditis_imperialis",
  "Caenorhabditis_uteleia",
  "Caenorhabditis_sp8", 
  "Caenorhabditis_quiockensis",
  "Caenorhabditis_angaria",
  "Caenorhabditis_virilis",
  "Caenorhabditis_drosophilae",
  "Caenorhabditis_bovis",
  "Caenorhabditis_plicata",
  "Caenorhabditis_parvicauda",
  "Caenorhabditis_monodelphis",
  "Caenorhabditis_auriculariae",
  "Diploscapter_coronatus"
)

df$species <- factor(df$species, levels = phylo_order)
df_ordered <- df[order(df$species), ]

p <- ggplot() + 
  geom_col(data=df_ordered, aes(x=rev(species),y=bs_count)) + 
  theme_bw() + 
  coord_flip() + 
  ylab("Break site count (n)") + 
  scale_y_continuous(expand=c(0,0)) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid = element_blank())
 
ggsave("barplots.pdf", p, width=6, height=8, units="in")