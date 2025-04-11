# Load required libraries
library(ggplot2)
library(cowplot)
library(dplyr)

# Read the data
data <- read.csv("combined_summary_with_sd_corrected.csv", stringsAsFactors = FALSE)

# Preprocess LDThres to factor
data <- data %>%
  mutate(LDThres = as.factor(LDThres))

# Filter data for NetType = "H" and LDThres = "0.8"
data_H_08 <- data %>% 
  filter(NetType == "H", LDThres == "0.8")

# Function to create a barplot
create_barplot <- function(data, y_var, se_var) {
  # No title since we'll use a single title for the figure
  p <- ggplot(data, aes(x = factor(chr), y = get(y_var), fill = factor(Phase))) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = get(y_var) - get(se_var), 
                      ymax = get(y_var) + get(se_var)),
                  position = position_dodge(0.9), width = 0.25) +
    labs(x = "Chromosome", 
         y = ifelse(y_var == "AUROC", "AUROC", "AUPR")) +
    scale_fill_manual(values = c("1" = "#1b0077", "3" = "#d90002"), 
                      name = "Phase") +  # Specify colors and legend title
    theme_minimal() +
    theme(legend.position = "right",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(colour = "grey80", size = 0.5, fill = NA),  # Box matching grid
          plot.background = element_rect(fill = "white", colour = "white"))
  
  return(p)
}

# Create individual plots for NetType = "H" and LDThres = "0.8"
# Column 1: AUROC
plot_A <- create_barplot(data_H_08, "AUROC", "AUROC.se")  # AUROC

# Column 2: AUPR
plot_B <- create_barplot(data_H_08, "AUPR", "AUPR.se")    # AUPR

# Combine plots using plot_grid with a single title
combined_plot <- plot_grid(plot_A, plot_B,
                           labels = c("A", "B"),
                           ncol = 2, nrow = 1,
                           align = "v", axis = "l") +
  theme(plot.background = element_rect(fill = "white", colour = "white"))

# Save the plot
ggsave("Figure6.pdf", combined_plot, width = 14, height = 4, bg = "white")

# Display the plot
print(combined_plot)