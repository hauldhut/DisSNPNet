# Read the input files
summary_data <- read.csv("/Users/hauldhut/Manuscripts/70MeSHSNPAssoc/Figure/combined_summary_data.csv", stringsAsFactors = FALSE)
mapping_results <- read.csv("/Users/hauldhut/Manuscripts/70MeSHSNPAssoc/Figure/Map_MeSH2SNP_2_M-H.csv", stringsAsFactors = FALSE)

# Merge the data frames based on Method, chr, ld, and phase
combined_data <- merge(summary_data, 
                       mapping_results[, c("Method", "chr", "ld", "phase", "NuDisease3")],
                       by = c("Method", "chr", "ld", "phase"),
                       all.x = TRUE)

# Calculate standard deviations using corrected formula
combined_data$auROCbyD.sd <- combined_data$auROCbyD.se * sqrt(combined_data$NuDisease3)
combined_data$auPRbyD.sd <- combined_data$auPRbyD.se * sqrt(combined_data$NuDisease3)

# Write to new file
write.csv(combined_data, "/Users/hauldhut/Manuscripts/70MeSHSNPAssoc/Figure/combined_summary_with_sd.csv", row.names = FALSE)

# Print summary
print(head(combined_data))