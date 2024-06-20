library(ggplot2)

# Extract data from files
file1 <- paste0("/Volumes/Seagate5TB/Daniel/32spec_project/Analyses/Genome_synteny/muller_block_sizes.txt")
data1 <- read.table(text = gsub("#", "\t", readLines(file1)), header=T, comment.char = "")

# Extract vectors
A <- ifelse(is.na(data1$A), NA, data1$A / 1000)
B <- ifelse(is.na(data1$B), NA, data1$B / 1000)
C <- ifelse(is.na(data1$C), NA, data1$C / 1000)
D <- ifelse(is.na(data1$D), NA, data1$D / 1000)
E <- ifelse(is.na(data1$E), NA, data1$E / 1000)

median(na.omit(data1$A))
median(na.omit(data1$B))
median(na.omit(data1$C))
median(na.omit(data1$D))
median(na.omit(data1$E))

# Boxplot of block sizes per muller element
ggplot() + 
  geom_boxplot(aes(x="A", y=A), fill="#C00000", outlier.alpha = 0.1, outlier.size = 0.75) +
  geom_boxplot(aes(x="B", y=B), fill="#0AA84E", outlier.alpha = 0.1, outlier.size = 0.75) +
  geom_boxplot(aes(x="C", y=C), fill="#003C9B", outlier.alpha = 0.1, outlier.size = 0.75) +
  geom_boxplot(aes(x="D", y=D), fill="#F7B300", outlier.alpha = 0.1, outlier.size = 0.75) +
  geom_boxplot(aes(x="E", y=E), fill="#702FA0", outlier.alpha = 0.1, outlier.size = 0.75) +
  labs(x ="Muller element", y = "Block size [kb]") + 
  scale_y_log10() +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")


# Extract data from files
file2 <- paste0("/Volumes/Seagate5TB/Daniel/32spec_project/Analyses/Genome_synteny/muller_breaks_per_dist.txt")
data2 <- read.table(text = gsub("#", "\t", readLines(file2)), header=T, comment.char = "")

# Remove muller F
my_data <- subset(data2, mull != "F")

# Custom function to calculate mean and IQR
mean_with_iqr <- function(x) {
  data.frame(y = mean(x), ymin = quantile(x, 0.25), ymax = quantile(x, 0.75))
}

# Line plot of mean number of breaks per evolutionary distance
# Plot using custom IQR function (instead of std.dev: mean_cl_normal)
ggplot(my_data, aes(x = dist, y = brks, group = mull, colour = mull, fill = mull)) + 
  stat_summary(fun.data = mean_with_iqr, geom = "ribbon", alpha = 0.2, colour = NA) +
  stat_summary(fun.y = mean, geom = "line", size = 1) +
  labs(x = "Evolutionary distance", y = "Breaks per Mb") +
  scale_color_manual(values = c("A" = "#C00000", "B" = "#0AA84E", "C" = "#003C9B", "D" = "#F7B300", "E" = "#702FA0")) +
  scale_fill_manual(values = c("A" = "#C00000", "B" = "#0AA84E", "C" = "#003C9B", "D" = "#F7B300", "E" = "#702FA0")) +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")


# Extract data from files
file3 <- paste0("/Volumes/Seagate5TB/Daniel/32spec_project/Analyses/Genome_synteny/muller_block_sizes_per_dist.txt")
data3 <- read.table(text = gsub("#", "\t", readLines(file3)), header=T, comment.char = "")

# Remove muller F
my_data <- subset(data3, mull != "F")
my_data$brks <- my_data$brks / 1000

# Custom function to calculate mean and IQR
mean_with_iqr <- function(x) {
  data.frame(y = mean(x), ymin = quantile(x, 0.25), ymax = quantile(x, 0.75))
}

# Line plot of mean number of breaks per evolutionary distance
# Plot using custom IQR function (instead of std.dev: mean_cl_normal)
ggplot(my_data, aes(x = dist, y = brks, group = mull, colour = mull, fill = mull)) + 
  stat_summary(fun.data = mean_with_iqr, geom = "ribbon", alpha = 0.2, colour = NA) +
  stat_summary(fun.y = mean, geom = "line", size = 1) +
  labs(x = "Evolutionary distance", y = "Block size [kb]") +
  scale_color_manual(values = c("A" = "#C00000", "B" = "#0AA84E", "C" = "#003C9B", "D" = "#F7B300", "E" = "#702FA0")) +
  scale_fill_manual(values = c("A" = "#C00000", "B" = "#0AA84E", "C" = "#003C9B", "D" = "#F7B300", "E" = "#702FA0")) +
  scale_y_log10() +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")


ggplot(my_data, aes(x = dist, y = brks, group = mull, colour = mull, fill = mull)) + 
  stat_summary(fun.data = mean_with_iqr, geom = "ribbon", alpha = 0.2, colour = NA) +
  stat_summary(fun.y = mean, geom = "line", size = 1) +
  labs(x = "Evolutionary distance", y = "Block size [kb]") +
  scale_color_manual(values = c("A" = "#C00000", "B" = "#0AA84E", "C" = "#003C9B", "D" = "#F7B300", "E" = "#702FA0")) +
  scale_fill_manual(values = c("A" = "#C00000", "B" = "#0AA84E", "C" = "#003C9B", "D" = "#F7B300", "E" = "#702FA0")) +
  scale_y_log10() +
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(1, 1),  # Position at the top right of the plot area
        legend.justification = c(1, 1),  # Anchor the legend at the top-right corner
        legend.margin = margin(t = 25, r = 25, b = 10, l = 10)) # Adjust margin to fit inside the plot
