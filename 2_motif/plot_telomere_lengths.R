library(ggplot2)
library(tidyverse)
library(patchwork)

nxCaeAuri1 <- read.table("nxCaeAuri1_telomere_lengths.tsv", col.names=c("count", "category")) %>% mutate(species = "C. auriculariae")
nxCaeMono1 <- read.table("nxCaeMono1_telomere_lengths.tsv", col.names=c("count", "category")) %>% mutate(species = "C. monodelphis")
nxCaeParv1 <- read.table("nxCaeParv1_telomere_lengths.tsv", col.names=c("count", "category")) %>% mutate(species = "C. parvicauda") %>% filter(category == "somatic")

df <- rbind(nxCaeAuri1, nxCaeMono1, nxCaeParv1)

p1 <- ggplot() + 
  geom_density(data=df, aes(y=count*6, fill=category)) + 
  facet_grid(~species) +
  theme_bw() + 
  ylab("Telomere length (bp)") + 
  scale_y_log10() + 
  coord_flip() + 
  scale_fill_manual(values=c("#ef4136", "#1c75bc")) + 
  theme(legend.position = c(0.92,0.82), legend.background = element_rect(colour = 'black'))

nxCaeParv1_all_counts <- read.table("nxCaeParv1_telomere_repeat_counts.tsv", col.names=c("count", "readname"))

nxCaeParv1_all_counts <- nxCaeParv1_all_counts %>% select(count) %>% mutate(label = "all reads") %>% mutate(species="C. parvicauda")
nxCaeParv1 <- nxCaeParv1 %>% select(count) %>% mutate(label = "somatic reads") %>% mutate(species="C. parvicauda")

df <- rbind(nxCaeParv1_all_counts, nxCaeParv1)

p2 <- ggplot() +
  geom_histogram(data=df, aes(x=count*6, fill=label)) +
  facet_grid(~species) + 
  theme_bw() + 
  xlab("Telomere length (bp)") + 
  scale_fill_manual(values=c("grey", "#1c75bc")) + 
  theme(legend.position = c(0.92,0.82), legend.background = element_rect(colour = 'black'))

p <- p1 / p2 

ggsave("telomere_lengths.pdf", p, height=8, width=10, units="in")