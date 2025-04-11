# Read the two result files
results_H <- read.csv("/Users/hauldhut/Manuscripts/70MeSHSNPAssoc/Figure/Map_MeSH2SNP_2_H.csv", stringsAsFactors = FALSE)
results_M <- read.csv("/Users/hauldhut/Manuscripts/70MeSHSNPAssoc/Figure/Map_MeSH2SNP_2_M.csv", stringsAsFactors = FALSE)

# Add Method column to each
results_H$Method <- "H"
results_M$Method <- "M"

# Combine the two data frames
combined_results <- rbind(results_H, results_M)

# Reorder columns to put Method first
combined_results <- combined_results[, c("Method", "chr", "ld", "phase", 
                                         "NuDisease", "NuSNP", "NuDisease3")]

# Write to new file
write.csv(combined_results, "/Users/hauldhut/Manuscripts/70MeSHSNPAssoc/Figure/Map_MeSH2SNP_2_M-H.csv", row.names = FALSE)

# Print summary
print(head(combined_results))
print(tail(combined_results))