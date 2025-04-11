#Step 1: Load Required Libraries and Data
# Install packages if not already installed
# install.packages(c("ggplot2", "dplyr", "tidyr", "ggrepel", "pheatmap"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(pheatmap)

# Set working directory to where your files are stored (adjust path as needed)
setwd("~/Manuscripts/70MeSHSNPAssoc/Results/Prediction")

topk=30
# List all chromosome files
file_list <- paste0("H_byDisease_Phase3_Chr", 1:22, "_r208.ld_predict_top",topk,"_evid.txt")

# Function to read and add chromosome column
read_chrom_file <- function(file) {
  chrom <- gsub("chr|.txt", "", file)  # Extract chromosome number
  data <- read.delim(file, header = FALSE, sep = "\t", 
                     col.names = c("Disease_ID", "Disease_Name", "Description", 
                                   "EFO_ID", "SNP_ID", "P_value", "Evidence")) %>%
    mutate(Chromosome = as.numeric(chrom))
  return(data)
}

# Load and combine all files
all_data <- lapply(file_list, read_chrom_file) %>%
  bind_rows()  # Combine into one data frame

# Convert p-values to numeric and calculate -log10(p-value)
all_data <- all_data %>%
  mutate(P_value = as.numeric(P_value),
         LogP = -log10(P_value))

# Remove duplicates if desired
all_data_unique <- all_data %>% distinct(Disease_Name, SNP_ID, .keep_all = TRUE)


rsids <- unique(all_data_unique$SNP_ID)
MESHids <- unique(all_data_unique$Disease_ID)
print(length(MESHids))

# Load library
library(biomaRt)
# Connect to Ensembl BioMart (GRCh38/hg38)
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "https://www.ensembl.org")
# # Query Ensembl for SNP positions

## For SMALL list of rsids
# snp_positions_original <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
#                        filters = "snp_filter",
#                        values = rsids,
#                        mart = ensembl)

## For LARGE list of rsids
chunk_size <- 1000
rsid_chunks <- split(rsids, ceiling(seq_along(rsids) / chunk_size))
snp_positions_list <- lapply(rsid_chunks, function(chunk) {
  getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
        filters = "snp_filter", values = chunk, mart = ensembl)
})
snp_positions_original <- bind_rows(snp_positions_list)

saveRDS(snp_positions_original,paste0("snp_positions_original",topk,".rdata"))

###################
topk=30
snp_positions_original = readRDS(paste0("snp_positions_original",topk,".rdata"))

library(dplyr)

# Assuming snp_positions is from biomaRt (Option 1)
snp_positions <- snp_positions_original %>%
  rename(SNP_ID = refsnp_id, Chromosome_hg38 = chr_name, Position_hg38 = chrom_start) %>%
  select(SNP_ID, Chromosome_hg38, Position_hg38) %>%
  filter(Chromosome_hg38 %in% c(1:22)) %>%
  mutate(Chromosome_hg38 = as.numeric(Chromosome_hg38))

# Merge with original data, keeping only hg38 columns
all_data_with_positions <- all_data_unique %>%
  select(-Chromosome) %>%  # Drop original Chromosome column to avoid conflict
  left_join(snp_positions, by = "SNP_ID")

# Rename for simplicity
all_data_with_positions <- all_data_with_positions %>%
  rename(Chromosome = Chromosome_hg38, Position = Position_hg38)

#######################
#Step 2: Visualization 1 - Genome-Wide Manhattan Plot
library(ggplot2)
library(dplyr)
library(ggrepel)

# Filter out rows where Chromosome is NA
all_data_with_positions <- all_data_with_positions %>%
  filter(!is.na(Chromosome)) %>%  # Remove NA chromosomes
  mutate(Chromosome = as.numeric(Chromosome),
         Position = as.numeric(Position))

# Check for remaining missing positions
missing_positions <- all_data_with_positions %>% filter(is.na(Position) | is.na(Chromosome))
if (nrow(missing_positions) > 0) {
  message("Warning: ", nrow(missing_positions), " SNPs lack position data.")
} else {
  message("All SNPs have valid chromosome and position data.")
}

# Manhattan plot with chromosomes in one row, numbers at bottom
manhattan_plot <- ggplot(all_data_with_positions, aes(x = Position, y = LogP, color = as.factor(Chromosome))) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
  facet_wrap(~Chromosome, scales = "free_x", nrow = 1, switch = "x") +  # Move labels to bottom
  labs(title = "",
       x = "",  # No x-axis label
       y = "-log10(P-value)", 
       color = "Chromosome") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),  # Increase y-axis text (if visible)
        axis.title.y = element_text(size = 14),  # Increase y-axis label
        plot.title = element_text(size = 16, hjust = 0.5),  # Increase title size
        strip.text = element_text(size = 12),  # Increase chromosome labels
        strip.background = element_blank()) +
  geom_text_repel(data = subset(all_data_with_positions, P_value < 5e-8),
                  aes(label = SNP_ID), size = 3, max.overlaps = 20)  # Increase SNP label size

# Display plot
print(manhattan_plot)


# Save with a wider width to accommodate 22 facets
ggsave("../../Figure/manhattan_plot.png", width = 18, height = 6, units = "in")


#Step 3: Visualization 2 - Bar Plot of SNP Counts per Disease Across Chromosomes
# Count SNPs per disease
disease_counts <- all_data_unique %>%
  group_by(Disease_Name) %>%
  summarise(SNP_Count = n()) %>%
  arrange(desc(SNP_Count))

# Bar plot
bar_plot = ggplot(disease_counts, aes(x = reorder(Disease_Name, SNP_Count), y = SNP_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "",
       x = "", y = "Number of SNPs") +
  theme_minimal()+
  theme(plot.title = element_text(size = 16, hjust = 0.5),  # Increase title size
        axis.title.x = element_text(size = 14),  # Increase x-axis label
        axis.title.y = element_text(size = 14),  # Increase y-axis label
        axis.text.x = element_text(size = 12),  # Increase x-axis text
        axis.text.y = element_text(size = 12))  # Increase y-axis text (disease names)

write.csv(disease_counts,paste0("disease_counts_top",topk,".csv"))
# Display plot
print(bar_plot)

ggsave("../../Figure/bar_plot.png", width = 6, height = 18)


# library('cowplot')
# plot_grid(bar_plot, manhattan_plot, labels=c("A", "B"), ncol = 2, nrow = 1)
# 
# ggsave("../../Figure/Summary_topKEvidence_Manhattan.pdf", width = 10, height = 20)

library(magick)

# Read the Manhattan plot and rotate it
manhattan_img <- image_read("../../Figure/manhattan_plot.png")
manhattan_rotated_img <- image_rotate(manhattan_img, 90)

# Save the rotated image
image_write(manhattan_rotated_img, "../../Figure/manhattan_plot_rotated.png")

library(cowplot)

# Load the rotated Manhattan plot as a drawable object
manhattan_rotated <- ggdraw() + draw_image("../../Figure/manhattan_plot_rotated.png")

# Combine the plots
combined_plot <- plot_grid(bar_plot, manhattan_rotated, labels = c("A", "B"), ncol = 2, nrow = 1)

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave(paste0("../../Figure/Summary_topKEvidence_Manhattan_top",topk,".pdf"), plot = combined_plot, width = 12, height = 18)


