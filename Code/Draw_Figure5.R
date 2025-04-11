# Load required libraries
library(ggplot2)
library(cowplot)
library(dplyr)

# Read the data
data <- read.csv("combined_summary_with_sd_corrected.csv", stringsAsFactors = FALSE)

# Filter data for NetType = "H"
data_H <- data %>% filter(NetType == "H")

# Function to create a barplot
create_barplot <- function(data, y_var, se_var, phase_val) {
  filtered_data <- data %>%
    filter(Phase == phase_val)
  
  fontsize=14
  # Construct title with new format
  title <- paste("Phase", phase_val)
  filtered_data$LDThres = as.factor(filtered_data$LDThres)
  p <- ggplot(filtered_data, aes(x = factor(chr), y = get(y_var), fill = LDThres)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = get(y_var) - get(se_var), 
                      ymax = get(y_var) + get(se_var)),
                  position = position_dodge(0.9), width = 0.25) +
    labs(x = "Chromosome", 
         y = ifelse(y_var == "AUROC", "AUROC", "AUPR"),
         title = title) +
    theme_minimal() +
    theme(legend.position = "right",
          axis.title.x = element_text(size = fontsize),
          axis.title.y = element_text(size = fontsize),
          axis.text.x = element_text(size = fontsize),
          axis.text.y = element_text(size = fontsize),
          plot.title = element_text(size = fontsize, hjust = 0.5),  # Center the title
          legend.text = element_text(size = fontsize),
          legend.title = element_text(size = fontsize),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(colour = "grey80", size = 0.5, fill = NA),  # Box matching grid
          plot.background = element_rect(fill = "white", colour = "white"))
  
  return(p)
}

# Create individual plots for NetType = "H"
# Row 1: AUROC
plot_A <- create_barplot(data_H, "AUROC", "AUROC.se", 1)  # Phase 1, AUROC
plot_B <- create_barplot(data_H, "AUROC", "AUROC.se", 3)  # Phase 3, AUROC

# Row 2: AUPR
plot_C <- create_barplot(data_H, "AUPR", "AUPR.se", 1)    # Phase 1, AUPR
plot_D <- create_barplot(data_H, "AUPR", "AUPR.se", 3)    # Phase 3, AUPR

# Adjust subfigures (similar to previous layout)
plot_A <- plot_A + theme(legend.position = "none", plot.margin = unit(c(1, 4, 1, 1), "lines"))
plot_B <- plot_B + labs(y = "")

plot_C <- plot_C + theme(legend.position = "none", plot.margin = unit(c(1, 4, 1, 1), "lines"))
plot_D <- plot_D + labs(y = "")

# Combine plots using plot_grid
combined_plot <- plot_grid(plot_A, plot_B,
                           plot_C, plot_D,
                           labels = c("A", "B", "C", "D"),
                           ncol = 2, nrow = 2,
                           align = "v", axis = "l") +
  theme(plot.background = element_rect(fill = "white", colour = "white"))

# Save the plot
ggsave("Figure5.pdf", combined_plot, width = 14, height = 8, bg = "white")

# Display the plot
print(combined_plot)