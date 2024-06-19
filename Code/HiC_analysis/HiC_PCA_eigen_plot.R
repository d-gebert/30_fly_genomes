# Load required libraries
library(ggplot2)
library(reshape2)
library(gridExtra)

# Load the Hi-C contact matrix and ensure values are treated as numeric
matrix_data <- read.csv("/Volumes/Seagate5TB/Daniel/32spec_project/Analyses/Genome_structure/HiC_remap/genome_matrix_balanced_log10.csv", header=TRUE, row.names=1, check.names=FALSE)

# Load the PCA values
pca_data <- read.csv("/Volumes/Seagate5TB/Daniel/32spec_project/Analyses/Genome_structure/HiC_remap/pca_values.csv", header=TRUE, row.names=1)

# Extract bin names
bin_names <- colnames(matrix_data)

# Identify chromosome transitions
chromosomes <- sapply(strsplit(bin_names, ":"), `[`, 1)
chromosome_transitions <- which(diff(as.numeric(as.factor(chromosomes))) != 0) + 0.5

# Melt the matrix for ggplot
melted_matrix <- melt(as.matrix(matrix_data))

# Reverse the y-axis by changing the order of Var2
melted_matrix$Var2 <- factor(melted_matrix$Var2, levels = rev(levels(melted_matrix$Var2)))

# Adjust the chromosome transitions for the reversed y-axis
chromosome_transitions_adjusted <- nrow(matrix_data) - chromosome_transitions + 0.5

# Create the heatmap
heatmap <- ggplot(data = melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.ticks = element_blank(),
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_blank()
  )

# Add lines to separate chromosomes
for (pos in chromosome_transitions) {
  heatmap <- heatmap + geom_hline(yintercept = nrow(matrix_data) - pos + 0.5, color = "black", size = 0.25) +
    geom_vline(xintercept = pos, color = "black", size = 0.25)
}

# Add border around the heatmap
heatmap <- heatmap + geom_segment(aes(x = 0.25, y = 0.25, xend = 0.25, yend = nrow(matrix_data) + 0.25), size = 0.25) +
  geom_segment(aes(x = 0.25, y = 0.25, xend = ncol(matrix_data) + 0.25, yend = 0.25), size = 0.25) +
  geom_segment(aes(x = ncol(matrix_data) + 0.25, y = 0.25, xend = ncol(matrix_data) + 0.25, yend = nrow(matrix_data) + 0.25), size = 0.25) +
  geom_segment(aes(x = 0.25, y = nrow(matrix_data) + 0.25, xend = ncol(matrix_data) + 0.25, yend = nrow(matrix_data) + 0.25), size = 0.25)

# Create a bar plot for the PCA values
pca_data$bin <- factor(rownames(pca_data), levels = rownames(matrix_data))
pca_data$PC1_color <- ifelse(pca_data$PC1 > 0, "red", "blue")

pca_barplot <- ggplot(pca_data, aes(x = bin, y = PC1, fill = PC1_color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  theme_minimal() +
  theme(
    axis.text = element_blank(), 
    axis.title = element_blank(), 
    axis.ticks = element_blank(),
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_blank()
  )

# Combine the heatmap and the PCA bar plot
combined_plot <- grid.arrange(heatmap, pca_barplot, ncol = 1, heights = c(3, 1))

# Save the combined plot to a file
ggsave("/Volumes/Seagate5TB/Daniel/32spec_project/Analyses/Genome_structure/HiC_remap/genome_heatmap_with_pca.pdf", plot = combined_plot, width = 12, height = 16)
