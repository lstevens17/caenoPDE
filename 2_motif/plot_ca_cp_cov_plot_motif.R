library(ggplot2)
library(patchwork)
library(tidyverse)

#### CA PLOT ####
# Load the BED files into data frames
cov_bed <- read.table("I_right.per-base.bed", col.names = c("chr", "start", "stop", "cov"))
delfies_bed <- read.table("I_right.delfies.bed", col.names = c("chr", "start", "stop", "name", "score", "strand"))

left_lim = 18031700
right_lim = 18031950

# filter
cov_bed <- cov_bed %>% filter(start > left_lim & start < right_lim)
delfies_bed <- delfies_bed %>% filter(start > left_lim & start < right_lim)

ca_p1 <- ggplot() +
  geom_rect(data = cov_bed, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0, ymax=cov), fill = "#c9c9c9") +
  geom_rect(data = delfies_bed, aes(xmin=(start/1e6)-0.000001, xmax=(stop/1e6)+0.000001, ymin=0, ymax=score), linewidth = 2, fill="black") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,130)) +
  ylab("Depth") +
  theme(legend.position = "None",
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 14),
        plot.margin = unit(c(0.2,0,0,0), "cm"))

ca_p2 <- ggplot(data = cov_bed, aes(x = start / 1e6)) +
  geom_rect(aes(xmin=18031803/1e6, xmax=18031840/1e6, ymin=2, ymax=8), fill="#1F78B4") +
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

left_lim - right_lim

#### CP PLOT ####
# Load the BED files into data frames
cov_bed <- read.table("II_right.per-base.bed", col.names = c("chr", "start", "stop", "cov"))
delfies_bed <- read.table("II_right.delfies.bed", col.names = c("chr", "start", "stop", "name", "score", "strand"))


# filter
left_lim = 17202753
right_lim = 17203003

# filter
cov_bed <- cov_bed %>% filter(start > left_lim & start < right_lim)
delfies_bed <- delfies_bed %>% filter(start > left_lim & start < right_lim)

cp_p1 <- ggplot() +
  geom_rect(data = cov_bed, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0, ymax=cov), fill = "#c9c9c9") +
  geom_rect(data = delfies_bed, aes(xmin=(start/1e6)-0.000001, xmax=(stop/1e6)+0.000001, ymin=0, ymax=score), linewidth = 2, fill="black") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,65)) +
  ylab("Depth") +
  theme(legend.position = "None",
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 14),
        plot.margin = unit(c(0.2,0,0,0), "cm"))

cp_p2 <- ggplot(data = cov_bed, aes(x = start / 1e6)) +
  geom_rect(aes(xmin=17202842/1e6, xmax=17202866/1e6, ymin=2, ymax=8), fill="#1F78B4") +
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

left_lim - right_lim

covplot <- ca_p1 / ca_p2 / cp_p1 / cp_p2 + plot_layout(heights=c(8,1,8,1))

ggsave("ca_cp_cov_plot_motif.pdf", covplot, height=10, width=9, units="in")