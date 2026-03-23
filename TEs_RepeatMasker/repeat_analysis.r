library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(svglite)

df <- read.csv('repeats_updated.csv')

gene_locus <- list(
  'Actin' = list(start = 2143197, end = 2145226),
  'α-Tubulin' = list(start = 523332, end = 525367),
  'β-Tubulin' = list(start = 1245400, end = 1247425),
  'Elongation factor 1' = list(start = 1412453, end = 1414989),
  'PST130_P495001' = list(start = 2674372, end = 2674809)
)

all_plot_data <- list()
gene_counter <- 0

for (gene in unique(df$Gene)) {
  gene_counter <- gene_counter + 1
  gene_group <- df[df$Gene == gene, ]
  
  # Filter out repeats that are beyond 20kb from CDS boundaries
  gene_group <- gene_group[gene_group$proximity_to_CDS <= 20000, ]
  
  if (nrow(gene_group) == 0) next
  
  # Get actual CDS length for this gene
  if (gene %in% names(gene_locus)) {
    actual_cds_length <- gene_locus[[gene]]$end - gene_locus[[gene]]$start
  } else {
    actual_cds_length <- 2000  # Default fallback
  }
  
  
  gene_group$relative_position <- ifelse(
    gene_group$stream == "Upstream",
    -gene_group$proximity_to_CDS,
    gene_group$proximity_to_CDS
  )
  
  # Create start and end positions for each repeat
  # Use original repeat lengths but ensure minimum visibility
  original_lengths <- gene_group$End - gene_group$Start
  scaled_lengths <- pmax(original_lengths, 100)  # Minimum 100bp for visibility
  
  gene_group$scaled_start <- gene_group$relative_position - scaled_lengths/2
  gene_group$scaled_end <- gene_group$relative_position + scaled_lengths/2
  
  region_start <- -20000
  region_end <- 20000
  
  cds_start <- -100
  cds_end <- 100
  
  gene_group$gene_name <- gene
  gene_group$gene_order <- gene_counter
  gene_group$cds_start <- cds_start
  gene_group$cds_end <- cds_end
  gene_group$region_start <- region_start
  gene_group$region_end <- region_end
  gene_group$actual_cds_length <- actual_cds_length
  
  unique_repeats_gene <- unique(gene_group$Repeats)
  n_repeats_gene <- length(unique_repeats_gene)
  
  if (n_repeats_gene <= 12) {
    colors_gene <- brewer.pal(max(3, min(n_repeats_gene, 12)), "Paired")
  } else {
    colors_gene <- c(brewer.pal(12, "Paired"), 
                     brewer.pal(min(n_repeats_gene - 12, 12), "Set3"))
  }
  
  gene_color_mapping <- setNames(colors_gene[1:n_repeats_gene], unique_repeats_gene)
  gene_group$repeat_color <- gene_color_mapping[gene_group$Repeats]
  
  all_plot_data[[gene]] <- gene_group
}

plot_data <- do.call(rbind, all_plot_data)

gene_order <- unique(plot_data$gene_name)[order(unique(plot_data$gene_order))]
plot_data$gene_name <- factor(plot_data$gene_name, levels = rev(gene_order))

text_data <- plot_data
text_data$x_center <- (text_data$scaled_start + text_data$scaled_end) / 2
text_data$y_center <- 1.05

p <- ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = region_start, xmax = region_end, 
                ymin = 0.75, ymax = 1), fill="white", color = "black", linewidth = 0.8) +
  
  # CDS region
  geom_rect(data = plot_data,
            aes(xmin = cds_start, xmax = cds_end, 
                ymin = 0.75, ymax = 1),
            fill = "lightgrey", color = "darkgrey", linewidth = 0.5) +
  
  # Repeat regions
  geom_rect(data = plot_data,
            aes(xmin = scaled_start, xmax = scaled_end, 
                ymin = 0.75, ymax = 1, fill = repeat_color)) +
  
  geom_text_repel(data = text_data,
                  aes(x = x_center, y = y_center, 
                      label = new_repeats, color = repeat_color),
                  size = 3,
                  fontface = "bold",
                  direction = "y",
                  nudge_y = 0.3,
                  box.padding = 0.2,
                  point.padding = 0.1,
                  segment.size = 0.5,
                  segment.alpha = 0.7,
                  max.overlaps = Inf,
                  min.segment.length = 0.1,
                  force = 3,
                  force_pull = 0.2) +
  
  scale_fill_identity() +
  scale_color_identity() +
  
  facet_wrap(~ gene_name, ncol = 1, scales = "free_x", 
             strip.position = "left") +

  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.ticks.x = element_line(color = "black"),
    
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),

    plot.margin = margin(5, 5, 5, 5),
    panel.spacing = unit(1, "lines"),
  ) +
  
  labs(x = "Distance from CDS Start (bp)") +
  
  coord_cartesian(ylim = c(0.75, 1.5)) +
  
  scale_x_continuous(
    breaks = c(-20000, -10000, -2000, 0, 2000, 10000, 20000),
    labels = c("-20kb", "-10kb", "-2kb", "CDS", "+2kb", "+10kb", "+20kb"),
    limits = c(-21000, 21000)
  )

n_genes <- length(unique(plot_data$gene_name))
plot_height <- 3 * n_genes + 2

ggsave("all_genes_repeat_diagrams.svg", plot = p, 
       width = 10, height = plot_height, 
       units = "in", dpi = 300,
       bg = "white")

create_gene_legend <- function(gene_name, plot_data) {
  gene_data <- plot_data[plot_data$gene_name == gene_name, ]
  unique_repeats <- unique(gene_data$Repeats)
  unique_colors <- unique(gene_data$repeat_color)
  
  legend_df <- data.frame(
    repeat_type = unique_repeats,
    color = unique_colors,
    x = 1,
    y = 1:length(unique_repeats),
    stringsAsFactors = FALSE
  )
  
  ggplot(legend_df, aes(x = x, y = y)) +
    geom_tile(aes(fill = color), width = 0.8, height = 0.8, 
              color = "white", linewidth = 1) +
    geom_text(aes(label = repeat_type), hjust = 0, x = 1.5, 
              size = 4, fontface = "bold") +
    scale_fill_identity() +
    theme_void() +
    xlim(0.5, max(nchar(legend_df$repeat_type)) * 0.2 + 2) +
    ylim(0.5, max(legend_df$y) + 0.5) +
    labs(title = paste(gene_name, "Repeat Types")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
}

for (gene in unique(plot_data$gene_name)) {
  legend_plot <- create_gene_legend(gene, plot_data)
  filename <- paste0(gsub("[^A-Za-z0-9α-β]", "_", gene), "_legend.svg")
  
  ggsave(filename, plot = legend_plot, 
         width = 8, height = length(unique(plot_data[plot_data$gene_name == gene, "Repeats"])) * 0.5 + 1, 
         units = "in", dpi = 300, bg = "white")
}
