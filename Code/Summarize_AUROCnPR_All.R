
setwd("~/Manuscripts/70MeSHSNPAssoc/Figure")

# Define the parameter vectors
Phases <- c(1, 3)
LDs <- c("r208", "r2def")
Methods <- c("M", "H")

# Create an empty list to store data frames
all_data <- list()

# Loop through all combinations
for (method in Methods) {
  for (phase in Phases) {
    for (ld in LDs) {
      # Construct filename
      filename <- paste0("Summary_AUROCnPR_", method, "_Phase", phase, "_", ld, ".ld.csv")
      
      # Check if file exists
      if (file.exists(filename)) {
        # Read the csv file
        df <- read.csv(filename)
        
        # Add the three columns
        df$Method <- method
        df$phase <- phase
        df$ld <- ld
        
        # Add to the list
        all_data[[length(all_data) + 1]] <- df
      } else {
        warning(paste("File not found:", filename))
      }
    }
  }
}

# Combine all data frames into one
combined_data <- do.call(rbind, all_data)

# Optionally, reorder columns to have Method, phase, ld first
combined_data <- combined_data[, c("Method", "phase", "ld", "chr", 
                                   "auROCbyD", "auROCbyD.se", 
                                   "auPRbyD", "auPRbyD.se")]

# Write the combined data to a new CSV file (optional)
write.csv(combined_data, "combined_summary_data.csv", row.names = FALSE)

# View the first few rows
head(combined_data)