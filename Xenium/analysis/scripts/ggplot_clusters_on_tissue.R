library(ggplot2)
library(viridis)  # For a nice color palette
library(Polychrome)

### create polychrome color palette ----
p30 = createPalette(30, c("#FF0000", "#00FF00", "#0000FF"), range = c(20, 80))
p30 <- sortByHue(p30)
p30 <- as.vector(t(matrix(p30, ncol=3)))
swatch(p30)
names(p30) = NULL

# Create plot data
plot_data <- data.frame(
  x = xenium.obj@meta.data$x,
  y = xenium.obj@meta.data$y,
  # cluster = xenium.obj@meta.data$SCT_snn_res.1.4
  # cluster = xenium.obj@meta.data$SCT_snn_res.0.8
  cluster = xenium.obj@meta.data$banksy_cluster
)

# Create the plot
ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.75, alpha = 0.7) +
  # scale_color_viridis_d(option = "turbo") +  # Using viridis color palette
  scale_color_manual(values=p30) +
  theme_minimal() +
  labs(
    x = "Spatial X",
    y = "Spatial Y",
    color = "Cluster",
    title = "Spatial Distribution of Clusters"
  ) +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1,  # Make plot square
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.key.size = unit(1.5, "lines")  # Makes the legend dots bigger
  ) +
  coord_fixed() +  # Maintain aspect ratio
  guides(color = guide_legend(override.aes = list(size = 3)))  # This specifically makes the legend dots bigger

# Save the plot if needed
ggsave("spatial_distribution.png", width = 12, height = 10, dpi = 300)