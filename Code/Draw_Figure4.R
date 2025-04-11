# Load required libraries
library(ggplot2)
library(cowplot)
library(dplyr)

# Read the data
data <- read.csv("combined_summary_with_sd_corrected.csv", stringsAsFactors = FALSE)

# Function to create a barplot with a box matching grid lines
create_barplot <- function(data, y_var, se_var, phase_val, ld_val) {
  filtered_data <- data %>%
    filter(Phase == phase_val, LDThres == ld_val)
  
  # Construct title with new format
  title <- paste("Phase", phase_val, ", LDThres", ld_val)
  
  p <- ggplot(filtered_data, aes(x = factor(chr), y = get(y_var), fill = NetType)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = get(y_var) - get(se_var), 
                      ymax = get(y_var) + get(se_var)),
                  position = position_dodge(0.9), width = 0.25) +
    labs(x = "Chromosome", 
         y = ifelse(y_var == "AUROC", "AUROC", "AUPR"),
         title = title) +
    scale_fill_manual(values = c("M" = "#d95f02", "H" = "#1b9e77")) +
    theme_minimal() +
    theme(legend.position = "right",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 14, hjust = 0.5),  # Center the title
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          panel.background = element_rect(fill = "white", colour = NA),  # No border here
          panel.border = element_rect(colour = "grey80", size = 0.5, fill = NA),  # Box matching grid
          plot.background = element_rect(fill = "white", colour = "white"))
  
  return(p)
}

# Create individual plots
# Row 1: AUROC, Phase 1
plot_A <- create_barplot(data, "AUROC", "AUROC.se", 1, "0.2")
plot_B <- create_barplot(data, "AUROC", "AUROC.se", 1, "0.8")

# Row 2: AUPR, Phase 1
plot_C <- create_barplot(data, "AUPR", "AUPR.se", 1, "0.2")
plot_D <- create_barplot(data, "AUPR", "AUPR.se", 1, "0.8")

# Row 3: AUROC, Phase 3
plot_E <- create_barplot(data, "AUROC", "AUROC.se", 3, "0.2")
plot_F <- create_barplot(data, "AUROC", "AUROC.se", 3, "0.8")

# Row 4: AUPR, Phase 3
plot_G <- create_barplot(data, "AUPR", "AUPR.se", 3, "0.2")
plot_H <- create_barplot(data, "AUPR", "AUPR.se", 3, "0.8")

# Adjust subfigures
plot_A <- plot_A + theme(legend.position = "none", plot.margin = unit(c(1, 3, 1, 1), "lines"))
plot_B <- plot_B + labs(y = "")

plot_C <- plot_C + theme(legend.position = "none", plot.margin = unit(c(1, 3, 1, 1), "lines"))
plot_D <- plot_D + labs(y = "")

plot_E <- plot_E + theme(legend.position = "none", plot.margin = unit(c(1, 3, 1, 1), "lines"))
plot_F <- plot_F + labs(y = "")

plot_G <- plot_G + theme(legend.position = "none", plot.margin = unit(c(1, 3, 1, 1), "lines"))
plot_H <- plot_H + labs(y = "")

# Combine plots using plot_grid
combined_plot <- plot_grid(plot_A, plot_B,
                           plot_C, plot_D,
                           plot_E, plot_F,
                           plot_G, plot_H,
                           labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                           ncol = 2, nrow = 4,
                           align = "v", axis = "l") +
  theme(plot.background = element_rect(fill = "white", colour = "white"))

# Save the plot
ggsave("Figure4.pdf", combined_plot, width = 14, height = 16, bg = "white")

# Display the plot
print(combined_plot)