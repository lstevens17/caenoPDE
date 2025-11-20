library(Biostrings)
library(ggplot2)
library(ggseqlogo)
library(patchwork)

# read alignments
ca_aln <- readDNAStringSet("gg_dinucleotide/nxCaeAuri1_40bp.fa")
cm_aln <- readDNAStringSet("gg_dinucleotide/nxCaeMono1_40bp.fa")
cp_aln <- readDNAStringSet("gg_dinucleotide/nxCaeParv1_40bp.fa")

# get seqs 
ca_seqs <- as.character(ca_aln)
cm_seqs <- as.character(cm_aln)
cp_seqs <- as.character(cp_aln)

# ggseqlogo 
ca_p <- ggplot() +
  geom_logo(ca_seqs, seq_type = "dna") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(1, 40, by = 5),
                     labels = seq(1, 40, by = 5) - 21) +
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) + 
  ggtitle(expression(italic("C. auriculariae"))) + 
  geom_vline(xintercept=20.5, linetype=2)

cm_p <- ggplot() +
  geom_logo(cm_seqs, seq_type = "dna") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(1, 40, by = 5),
                     labels = seq(1, 40, by = 5) - 21) +
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) + 
  ggtitle(expression(italic("C. monondelphis"))) + 
  geom_vline(xintercept=20.5, linetype=2)

cp_p <- ggplot() +
  geom_logo(cp_seqs, seq_type = "dna") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(1, 40, by = 5),
                     labels = seq(1, 40, by = 5) - 21) +
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) + 
  ggtitle(expression(italic("C. parvicauda"))) + 
  geom_vline(xintercept=20.5, linetype=2)
  
# plot 
p <- ca_p / cm_p / cp_p +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face="bold"))

# save
ggsave("gg_dinucleotide.png", p, height=150, width=150, units="mm")
ggsave("gg_dinucleotide.pdf", p, height=150, width=150, units="mm")