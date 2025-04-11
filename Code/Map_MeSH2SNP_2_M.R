# Define parameters
Phases <- c(1, 3)
LDs <- c("r208", "r2def")
chrs <- 1:22

# Initialize results data frame
results <- data.frame(chr = integer(),
                      ld = character(),
                      phase = integer(),
                      NuDisease = integer(),
                      NuSNP = integer(),
                      NuDisease3 = integer(),
                      stringsAsFactors = FALSE)

# Loop through all combinations
for (chr in chrs) {
  # Read bipartite network
  bipartite_file <- paste0("/Users/hauldhut/Data/GWAS/CAUSALdb/Chr_", chr, "_Assoc.txt_BinaryInteraction.csv")
  if (!file.exists(bipartite_file)) next
  
  bipartite <- read.csv(bipartite_file, stringsAsFactors = FALSE)
  
  # Get unique SNPs and diseases from bipartite network
  bipartite_snps <- unique(bipartite$SNP)
  bipartite_diseases <- unique(bipartite$disease)
  
  # NuDisease is simply the number of unique diseases in bipartite network
  nudisease <- length(bipartite_diseases)
  
  for (phase in Phases) {
    for (ld in LDs) {
      # Read SNP LD network
      snp_file <- paste0("/Users/hauldhut/Data/1KGP/1kg_phase", phase, "_chr", chr, "_", ld, ".ld_Net.txt")
      if (!file.exists(snp_file)) next
      
      snp_net <- read.table(snp_file,
                            sep = "\t",
                            col.names = c("rsID1", "LDweight", "rsID2"),
                            stringsAsFactors = FALSE)
      
      # Get unique SNPs from LD network
      unique_snps <- unique(c(snp_net$rsID1, snp_net$rsID2))
      
      # Calculate mapped SNPs
      mapped_snps <- intersect(bipartite_snps, unique_snps)
      
      # Calculate diseases with at least 3 SNPs
      # Create a data frame of SNP-disease mappings with mapped SNPs
      mapped_mappings <- bipartite[bipartite$SNP %in% mapped_snps, ]
      
      # Count SNPs per disease
      disease_snp_count <- table(mapped_mappings$disease)
      nudisease3 <- sum(disease_snp_count >= 3)
      
      # Add results
      results <- rbind(results, data.frame(
        chr = chr,
        ld = ld,
        phase = phase,
        NuDisease = nudisease,
        NuSNP = length(mapped_snps),
        NuDisease3 = nudisease3,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Write results to file
write.csv(results, "bipartite_snp_mapping_results.csv", row.names = FALSE)

# Print summary
print(head(results))