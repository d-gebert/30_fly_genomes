# Load required packages
library(ggplot2)
library(dplyr)

# Read the data from the CSV file
breakpoints_data <- read.csv("/Volumes/Seagate5TB/Daniel/32spec_project/Analyses/Genome_synteny/gene_order_breakpoints_Dmel.csv")

# Filter out scaffolds, assuming scaffolds have names different from A, B, C, D, E, F
filtered_data <- breakpoints_data %>%
  filter(grepl("^[A-F]$", Chromosome))  # This regex matches only chromosomes labeled A to F

# Define custom colors for each chromosome
chromosome_colors <- c("A" = "#C00000", "B" = "#0AA84E", "C" = "#003C9B", "D" = "#F7B300", "E" = "#702FA0", "F" = "#F07500")

# Create the plot
p <- ggplot(filtered_data, aes(x = Bin, y = Count, fill = Chromosome)) +
    geom_bar(stat = "identity", position = "dodge", width = 1) +
    scale_fill_manual(values = chromosome_colors, name = "Chromosome") +  # Set custom colors
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      plot.background = element_blank()    # Remove background grid
    ) +
    labs(title = "Breakpoints across Chromosomes",
         x = "Bin (genes per bin)",
         y = "Number of breakpoints") +
    theme(plot.title = element_text(hjust = 0.5),  # Centering the plot title
          axis.text.x = element_text(angle = 90, hjust = 1))  # Vertical x labels for clarity +

# Output the plot to a PDF
ggsave("Chromosome_gene_order_breakpoints.pdf", plot = p, width = 12, height = 4)

# Alternatively, display the plot
print(p)