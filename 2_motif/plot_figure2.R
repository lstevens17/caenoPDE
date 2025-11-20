library(ggseqlogo)
library(universalmotif)
library(ggplot2)
library(patchwork)
library(tidyverse)

#### MOTIF PLOTS ####
# Ca motif
meme <- read_meme("nxCaeAuri1_meme1_50bp_meme.txt")

ca_m1 <- ggplot() +
  geom_logo(meme@motif) +
  theme_bw()  +
  scale_x_continuous(expand=c(0,0), breaks = seq(1, 14, by = 1), limits=c(-0.5, 14.5)) +
  theme(panel.grid = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"),
        text = element_text(size = 14)) 

meme <- read_meme("nxCaeAuri1_meme2_40bp_meme.txt")

ca_m2 <- ggplot() +
  geom_logo(meme@motif) +
  theme_bw()  +
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 16, by = 1), labels = function(x) x + 14) +
  theme(panel.grid = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"),
        text = element_text(size = 14),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

# Cm motif
meme <- read_meme("nxCaeMono1_meme1_80bp_meme.txt")

cm_m1 <- ggplot() +
  geom_logo(meme@motif) +
  theme_bw()  +
  scale_x_continuous(expand=c(0,0), breaks = seq(1, 13, by = 1), limits=c(-2.5, 13.5)) +
  theme(panel.grid = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"),
        text = element_text(size = 14)) 

meme <- read_meme("nxCaeMono1_meme2_40bp_meme.txt")

cm_m2 <- ggplot() +
  geom_logo(meme@motif) +
  theme_bw()  +
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 16, by = 1), labels = function(x) x + 13) +
  theme(panel.grid = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"),
        text = element_text(size = 14),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

# Cp motif
meme <- read_meme("nxCaeParv1_110bp_meme.txt")

cp_motif <- ggplot() +
  geom_logo(meme@motif) +
  theme_bw() +
  scale_x_continuous(expand=c(0,0), breaks = seq(0, 24, by = 1)) + 
  xlab("Position (bp)") +
  theme(panel.grid = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"),
        text = element_text(size = 14)) 

#### COV PLOT ####
# Load the BED files into data frames
cov_bed <- read.table("V_right.per-base.bed", col.names = c("chr", "start", "stop", "cov"))
delfies_bed <- read.table("V_right.delfies.bed", col.names = c("chr", "start", "stop", "name", "score", "strand"))


# filter
cov_bed <- cov_bed %>% filter(start > 20.9440*1e6 & start < 20.9442*1e6)
delfies_bed <- delfies_bed %>% filter(start > 20.9440*1e6 & start < 20.9442*1e6)

p1 <- ggplot() +
  geom_rect(data = cov_bed, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0, ymax=cov), fill = "#c9c9c9") +
  geom_rect(data = delfies_bed, aes(xmin=(start/1e6)-0.000001, xmax=(stop/1e6)+0.000001, ymin=0, ymax=score), linewidth = 2, fill="black") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,600)) +
  ylab("Depth") +
  theme(legend.position = "None",
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 14),
        plot.margin = unit(c(0.2,0,0,0), "cm"))

p2 <- ggplot(data = cov_bed, aes(x = start / 1e6)) +
  geom_rect(aes(xmin=20944023/1e6, xmax=20944099/1e6, ymin=2, ymax=8), fill="#1F78B4") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10)) +
  theme_bw() +
  xlab("Position (Mb)") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(0,0,0.2,0), "cm"),
        text = element_text(size = 14))

covplot <- p1 / p2 + plot_layout(heights=c(8,1))


#### SPECIFICITY PLOT ####
nxCaeAuri1 <- read.table("nxCaeAuri1_motif_pairs.filtered.dedup.high_scoring.labelled.no_telomere.bed", col.names=c("chrom", "start", "stop", "orientation", "score", "strand", "gap_len",
                                                                                                        "m2_start", "m2_end", "m2_name", "m2_score", "m2_strand",
                                                                                                        "m1_start", "m1_stop", "m1_name", "m1_score", "m1_strand", "location")) %>% mutate(species = "nxCaeAuri1")

nxCaeMono1 <- read.table("nxCaeMono1_motif_pairs.filtered.dedup.high_scoring.labelled.no_telomere.bed", col.names=c("chrom", "start", "stop", "orientation", "score", "strand", "gap_len",
                                                                                                        "m2_start", "m2_end", "m2_name", "m2_score", "m2_strand",
                                                                                                        "m1_start", "m1_stop", "m1_name", "m1_score", "m1_strand", "location")) %>% mutate(species = "nxCaeMono1")

nxCaeParv1 <- read.table("nxCaeParv1_fimo_out.high_scoring.labelled.bed", col.names=c("chrom", "start", "stop", "dot", "score", "strand", "location")) %>% mutate(species = "nxCaeParv1")

species_labels <- c(
  "nxCaeAuri1" = "italic('C. auriculariae')",
  "nxCaeMono1" = "italic('C. monodelphis')",
  "nxCaeParv1" = "italic('C. parvicauda')"
)

df <- bind_rows(nxCaeAuri1, nxCaeMono1, nxCaeParv1)

df <- df %>%
  mutate(location = factor(location,
                           levels = c("real", "eliminated", "retained"),
                           labels = c("Break site", "Eliminated", "Retained")))


location_counts <- df %>%
  group_by(location, species) %>%
  summarise(n = n(), .groups = "drop")

spec <- df %>%
  ggplot(aes(x=location, y=score)) +
  geom_point(position="jitter", pch=21, aes(fill=location)) +
  geom_boxplot(aes(fill=location), alpha=0.5, outlier.shape = NA) +
  theme_bw() +
  labs(x = "Location", y = "Score") +
  facet_grid(~species, labeller = as_labeller(species_labels, default = label_parsed)) +
  geom_text(data = location_counts, aes(x = location, y = 143, label = n), size = 5) +
  theme(legend.position = "none", text = element_text(size = 14))

# compile plot and save
p <- (((ca_m1 | plot_spacer() | ca_m2) + plot_layout(widths=c(15,1.5,11))) /
  ((cm_m1 | plot_spacer() | cm_m2) + plot_layout(widths=c(16,0.5,11)))) /
  ((cp_motif | plot_spacer()) +  plot_layout(widths=c(24,0))) / 
  covplot / 
  spec + plot_layout(heights=c(1,1,1,5,4)) 

ggsave("Figure_2.pdf", p, height=14, width=11, units="in")