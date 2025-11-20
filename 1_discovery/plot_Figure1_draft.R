library(ggplot2)
library(ggtree)
library(treeio)
library(patchwork)
library(grid)
library(png)
library(tidyverse)

#### TREE ####
# read in rooted tree
tree <- read.tree("astral_tree_supermatrix.rooted.noDcoro_noCbovi.treefile.txt")

# remove underscores and italicise
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- paste0("italic('", tree$tip.label, "')")

# plot
tree <- ggtree(tree) + 
  geom_tiplab(parse=TRUE, size = 4) + 
  geom_treescale(x = 0.725, y = 14, width = 0.1, fontsize = 4) + 
  xlim(0,2.2) + 
  theme(plot.margin = unit(c(0,1,0,1), "cm"),)

#### HIC PLOTS ####

# Read the PNG image
ca_png <- readPNG("nxCaeAuri1_yahs.hic_screenshot.png")
cm_png <- readPNG("nxCaeMono1_yahs.hic_screenshot.png")
cp_png <- readPNG("nxCaeParv1_yahs.hic_screenshot.png")

# Convert to raster object for ggplot2
ca_g <- rasterGrob(ca_png, interpolate = TRUE)
cm_g <- rasterGrob(cm_png, interpolate = TRUE)
cp_g <- rasterGrob(cp_png, interpolate = TRUE)

# Plot in ggplot2
ca_hic <- ggplot() +
  annotation_custom(ca_g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()  
cm_hic <- ggplot() +
  annotation_custom(cm_g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() 
cp_hic <- ggplot() +
  annotation_custom(cp_g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() 


#### BUSCO PLOTS
# Function to generate plots for different species
generate_nigon_paint_plot <- function(tolid) {
  
  # Read the location file
  location_table <- read_tsv(paste0(tolid, "_location.tsv")) %>%
    group_by(query_chr) %>%
    filter(n() > 10)
  
  # Read the samtools faidx file
  faidx <- read_tsv(paste0(tolid, ".1.primary.fa.fai"),
                    col_names = c("query_chr", "length", "offset", "linebases", "linewidth")) %>%
    filter(query_chr %in% c("I", "II", "III", "IV", "V", "X")) %>%
    select(query_chr, length) %>%
    mutate(length = as.numeric(length) / 1e6)
  
  # Set palette
  cols <- c("I" = "#E40303", "II" = "#FF8C00",
            "III" = "#FFED00", "IV" = "#04bf3b",
            "V" = "#004CFF", "X" = "#a100c2")
  
  # only plot x-axis if tolid is nxCaeParv1
  x_axis_title <- if (tolid == "nxCaeParv1") {
    element_text()
  } else {
    element_blank()
  }
  
  # Plot
  p <- ggplot(location_table) +
    geom_rect(aes(xmin = (position / 1e6) - 0.05,
                  xmax = (position / 1e6) + 0.05,
                  ymax = 0, ymin = 10,
                  fill = assigned_chr), alpha = 0.4) +
    geom_rect(data = faidx,
              aes(xmin = 0, xmax = length, ymin = 0, ymax = 10),
              color = "black", fill = "white", alpha = 0.0) +
    facet_grid(query_chr ~ ., switch = "y") +  
    xlab("Position (Mb)") +
    #scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_continuous(breaks = NULL, expand = c(0,0)) +
    scale_fill_manual(values = cols, name = "Nigon") +
    theme_bw() +  #
    theme(
      panel.border = element_blank(),           # Remove the boxes around facets
      strip.background = element_blank(),       # Remove the background for facet strips
      panel.grid = element_blank(),
      strip.text.y.left = element_text(size=14, angle=0),
      text = element_text(size = 14),
      legend.position = "none", 
      axis.title.y = element_blank(), 
      axis.title.x = x_axis_title,
    )
  
  # Read and parse elimination sites into a list
  elimination_sites_list <- read_tsv(paste0(tolid, "_elimination_sites.tsv"),
                                     col_names = c("chr", "site")) %>%
    mutate(site = site / 1e6) %>%
    split(.$chr) %>%           # Split by chromosome
    lapply(function(x) x$site) # Extract the site positions for each chromosome
  
  # Add vertical dotted lines for elimination sites
  for (chr in names(elimination_sites_list)) {
    p <- p + geom_vline(
      data = tibble(query_chr = chr, xintercept = elimination_sites_list[[chr]]),
      aes(xintercept = xintercept),
      color = "black",
      linewidth = 1,
      linetype = "21"
    ) + geom_point(
      data = tibble(query_chr = chr, x = elimination_sites_list[[chr]]),
      aes(x=x, y=9),
      shape=25, 
      color="black",
      fill="grey",
      stroke=1,
      size=3
    )
  }
  
  # Save the plot
  assign(paste0(tolid, "_nigon"), p, envir = .GlobalEnv)
}

# List of tolids to process
tolids <- c("nxCaeAuri1", "nxCaeMono1", "nxCaeParv1")

# Loop through each species and generate the plots
for (tolid in tolids) {
  generate_nigon_paint_plot(tolid)
}

#### COVERAGE PLOTS ####
### region 1 ###
# Load the BED files into data frames
cov_bed <- read.table("nxCaeMono1_II_middle_2.per-base.bed", col.names = c("chr", "start", "stop", "cov"))
delfies_bed <- read.table("nxCaeMono1_II_middle_2.delfies.bed", col.names = c("chr", "start", "stop", "name", "score", "strand"))

# Plot
II_middle2_cov <- ggplot() + 
  geom_rect(data = cov_bed, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0, ymax=cov), fill = "#c9c9c9") + 
  geom_rect(data = delfies_bed, aes(xmin=(start/1e6)-0.0005/1.5, xmax=(stop/1e6)+0.0005/1.5, ymin=0, ymax=score), linewidth = 2, fill="black") +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,500)) + 
  xlab("Position (Mb)") + 
  ylab("Depth") + 
  theme(legend.position = "None")

### region 2 ###
# Load the BED files into data frames
cov_bed <- read.table("nxCaeMono1_III_right.per-base.bed", col.names = c("chr", "start", "stop", "cov"))
delfies_bed <- read.table("nxCaeMono1_III_right.delfies.bed", col.names = c("chr", "start", "stop", "name", "score", "strand"))


# Plot
III_right <- ggplot() + 
  geom_rect(data = cov_bed, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0, ymax=cov), fill = "#c9c9c9") + 
  geom_rect(data = delfies_bed, aes(xmin=(start/1e6)-0.0009/1.5, xmax=(stop/1e6)+0.0009/1.5, ymin=0, ymax=score), linewidth = 2, fill="black") +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,500)) + 
  xlab("Position (Mb)") + 
  ylab("Depth") + 
  theme(legend.position = "None")

#### PLOT ####
p <- (tree | (ca_hic / cm_hic / cp_hic) | (nxCaeAuri1_nigon / nxCaeMono1_nigon / nxCaeParv1_nigon)) / (III_right + II_middle2_cov) + 
  plot_layout(heights=c(4,1))

ggsave("Figure_1.png", p, height=14, width=11, units="in")
ggsave("Figure_1.pdf", p, height=14, width=11, units="in")
