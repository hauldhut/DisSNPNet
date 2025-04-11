# Load required libraries
library(dplyr)

perf = "AUPR"

# Read the data
data <- read.csv("combined_summary_with_sd_corrected.csv", stringsAsFactors = FALSE)

# Filter data for NetType "M" and "H" (Figure4)
compare = "NetType"
figure = "Figure4"
data1 <- data %>% filter(NetType == "M")
data2 <- data %>% filter(NetType == "H")

# data1 <- data %>% filter(LDThres == "0.2")
# data2 <- data %>% filter(LDThres == "0.8")

# data1 <- data %>% filter(Phase == 1)
# data2 <- data %>% filter(Phase == 3)

# # Select "H" (Figure5)
# compare = "LDThres"
# figure = "Figure5"
# data_H <- data %>% filter(NetType == "H")
# # Filter data for LDThres "0.2" and "0.8"
# data1 <- data_H %>% filter(LDThres == "0.2")
# data2 <- data_H %>% filter(LDThres == "0.8")

# # Select NetType = "H" and LDThres = "0.8" (Figure6)
# compare = "Phase"
# figure = "Figure6"
# data_H_08 <- data %>% filter(NetType == "H", LDThres == "0.8")
# # Filter data for LDThres "0.2" and "0.8"
# data1 <- data_H_08 %>% filter(Phase == 1)
# data2 <- data_H_08 %>% filter(Phase == 3)

# Perform Welch's t-test (unequal variances)
t_test_result <- t.test(data1[[perf]], data2[[perf]],
                        alternative = "less",  # Two-tailed test
                        var.equal = FALSE,          # Unequal variances
                        conf.level = 0.95)          # 95% confidence interval

# Print the results
cat("Welch's Two-Sample T-Test for AUROC between ",compare,"\n")
print(t_test_result)

# Extract and display key results
p_value <- t_test_result$p.value
t_stat <- t_test_result$statistic
df <- t_test_result$parameter
mean1 <- mean(data1[[perf]], na.rm = TRUE)
mean2 <- mean(data2[[perf]], na.rm = TRUE)
ci <- t_test_result$conf.int

cat("Welch's Two-Sample T-Test for AUROC between ",compare,"\n")
cat("\nSummary of Results:\n")
cat("Mean AUROC for data1:", round(mean1, 4), "\n")
cat("Mean AUROC for data2:", round(mean2, 4), "\n")
cat("T-statistic:", round(t_stat, 4), "\n")
cat("Degrees of Freedom:", round(df, 4), "\n")
cat("P-value:", format.pval(p_value, digits = 4), "\n")
cat("95% Confidence Interval for Difference in Means:", round(ci[1], 4), "to", round(ci[2], 4), "\n")

# # Optional: Save results to a file
# sink("t_test_AUROC_NetType_M_vs_H.txt")
# print("Welch's Two-Sample T-Test for AUROC between NetType M and H")
# print(t_test_result)
# cat("\nSummary of Results:\n")
# cat("Mean ",perf," for data1:", round(mean1, 4), "\n")
# cat("Mean ",perf," for data2:", round(mean2, 4), "\n")
# cat("T-statistic:", round(t_stat, 4), "\n")
# cat("Degrees of Freedom:", round(df, 4), "\n")
# cat("P-value:", format.pval(p_value, digits = 4), "\n")
# cat("95% Confidence Interval for Difference in Means:", round(ci[1], 4), "to", round(ci[2], 4), "\n")
# sink()

print(perf)
increase = (mean2-mean1)*100/mean1
cat(mean1, mean2, increase,"\n")
data1_ordered <- data1[order(data1[[perf]]),]
data1_first_last <- rbind(data1_ordered[1, ], data1_ordered[nrow(data1_ordered), ])
print(data1_first_last)

data2_ordered <- data2[order(data2[[perf]]),]
data2_first_last <- rbind(data2_ordered[1, ], data2_ordered[nrow(data2_ordered), ])
print(data2_first_last)

#####
if (figure=="Figure4"){
  data1_sel = data1[order(data1$NetType, data1$Phase, data1$LDThres, data1$chr),][,c("NetType", "LDThres", "Phase","chr","AUROC", "AUPR")]
  data2_sel = data2[order(data2$NetType, data2$Phase, data2$LDThres, data2$chr),][,c("NetType", "LDThres", "Phase","chr","AUROC", "AUPR")]
}
if (figure=="Figure5"){
  data1_sel = data1[order(data1$NetType, data1$Phase, data1$LDThres, data1$chr),][,c("NetType","LDThres", "Phase","chr","AUROC", "AUPR")]
  data2_sel = data2[order(data2$NetType, data2$Phase, data2$LDThres, data2$chr),][,c("NetType","LDThres", "Phase","chr","AUROC", "AUPR")]
}

if (figure=="Figure6"){
  data1_sel = data1[order(data1$NetType, data1$Phase, data1$LDThres, data1$chr),][,c("NetType","LDThres", "Phase","chr","AUROC", "AUPR")]
  data2_sel = data2[order(data2$NetType, data2$Phase, data2$LDThres, data2$chr),][,c("NetType","LDThres", "Phase","chr","AUROC", "AUPR")]
}


data1_sel$AUROC = round(data1_sel$AUROC,3)
data1_sel$AUPR = round(data1_sel$AUPR,3)

data2_sel$AUROC = round(data2_sel$AUROC,3)
data2_sel$AUPR = round(data2_sel$AUPR,3)

data_cbind = cbind(data1_sel, data2_sel)
write.csv(data_cbind,paste0(figure,"_data.csv"),row.names = FALSE)
data_cbind