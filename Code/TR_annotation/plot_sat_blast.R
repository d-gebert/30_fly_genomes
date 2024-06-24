# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)

# Read the BLAST output CSV file
blast_results <- read.csv("Dbip.cen.A.arrays.csv", stringsAsFactors = FALSE)

# Rename the query names by dropping all characters starting with a dash
blast_results$query.acc.ver <- str_remove(blast_results$query.acc.ver, "-.*")

# Filter the data to include only alignments with identity >= 98%
filtered_results <- blast_results %>% filter(identity >= 0)

# Create a plot for each subject sequence
unique_subjects <- unique(filtered_results$subject.acc.ver)

for (subject in unique_subjects) {
  # Filter the data for the current subject sequence
  subject_data <- filtered_results %>% filter(subject.acc.ver == subject)
  
  # Create the plot with rectangles for each alignment
  p <- ggplot(subject_data) +
    geom_rect(aes(xmin = s.start, xmax = s.end, ymin = 0, ymax = 1, fill = query.acc.ver), color = NA) +
    labs(title = paste("Alignments for Subject:", subject), 
         x = "Subject Sequence Position", 
         y = NULL, 
         fill = "Query Sequence") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.title.y = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),  # Set panel background to white
      plot.background = element_rect(fill = "white", color = NA)    # Set plot background to white
    ) +  
    facet_wrap(~ query.acc.ver, scales = "free_y", ncol = 1)  # Separate panels for each query sequence in a single column
  
  # Display the plot
  print(p)
  
  # Save the plot as a file with a white background
  ggsave(filename = paste("plot_", subject, ".pdf", sep=""), plot = p, width = 10, height = 3, bg = "white")
}