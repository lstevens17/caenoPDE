library(ggplot2)
library(patchwork)
library(grid)
library(png)


# Read the PNG image
ca_png <- readPNG("nxCaeAfra1_yahs.hic_screenshot.png")
cd_png <- readPNG("nxCaeDros1_yahs.hic_screenshot.png")
cp_png <- readPNG("nxCaePlic1_yahs.hic_screenshot.png")
cw_png <- readPNG("nxCaeWall1_yahs.hic_screenshot.png")

# Convert to raster object for ggplot2
ca_g <- rasterGrob(ca_png, interpolate = TRUE)
cd_g <- rasterGrob(cd_png, interpolate = TRUE)
cp_g <- rasterGrob(cp_png, interpolate = TRUE)
cw_g <- rasterGrob(cw_png, interpolate = TRUE)

# Plot in ggplot2
ca_hic <- ggplot() +
  annotation_custom(ca_g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() + 
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
cd_hic <- ggplot() +
  annotation_custom(cd_g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() + 
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
cp_hic <- ggplot() +
  annotation_custom(cp_g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() + 
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
cw_hic <- ggplot() +
  annotation_custom(cw_g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() + 
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

# plot 
p <- (ca_hic | cd_hic) / (cp_hic | cw_hic) 

# save
ggsave("other_hic_plots.png", p, height=14, width=14, units="in")
ggsave("other_hic_plots.pdf", p, height=14, width=14, units="in")