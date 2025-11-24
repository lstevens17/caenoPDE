library(ggplot2)
library(tidyverse)

# read in tables
ca <- read.table("nxCaeAuri1_te_family_enrichment_results.tsv", header=T) %>%
  mutate(species = "nxCaeAuri1")
cm <- read.table("nxCaeMono1_te_family_enrichment_results.tsv", header=T) %>%
  mutate(species = "nxCaeMono1")
cp <- read.table("nxCaeParv1_te_family_enrichment_results.tsv", header=T) %>%
  mutate(species = "nxCaeParv1")

# combine dataframes
df <- rbind(ca, cm, cp) %>% mutate(classification = str_extract(family, "^[^/]+")) 

# facet labels
species_labels <- c(
  "nxCaeAuri1" = "italic('C. auriculariae')",
  "nxCaeMono1" = "italic('C. monodelphis')",
  "nxCaeParv1" = "italic('C. parvicauda')"
)

# keep only those families that show enrichment and ignore all summary
df <- df %>%
  group_by(family) %>%
  filter(any(abs(z_score) >= 0.5)) %>%
  ungroup() %>%
  filter(family != "Satellite") %>% 
  filter(family != "Low_complexity") %>%
  filter(family != "Simple_repeat") %>%
  filter(family != "Unknown")

# create a new column for significance level
df <- df %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

# adjust y position for significance labels
df <- df %>%
  mutate(y_position = if_else(z_score >= 0, z_score + 2, z_score - 2))


# plot
p <- ggplot(data=df, aes(x=reorder(family, z_score), y=z_score, fill=classification)) + 
  geom_col() + 
  facet_wrap(~species, labeller = as_labeller(species_labels, default = label_parsed)) +  
  geom_hline(yintercept=0, linetype=2) +
  coord_flip() + 
  labs(x = "Repeat superfamily", y="Z-score") + 
  theme_bw() + 
  ylim(-19, 19) + 
  theme(legend.position = "none",
        strip.text = element_text(size="14")) + 
  geom_text(aes(label=significance, y=y_position), size=5, vjust=0.8, fontface="bold")

# save
ggsave("te_zscores.png", p, height=250, width=225, units="mm")
