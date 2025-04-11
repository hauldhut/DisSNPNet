# This folder contains all source code used in the manuscript

## Performance assessment
- **InvestigateData_MeSH2SNP.R**: Investigate distribution of Known Disease-SNP Associations Across Populations and Chromosomes from CAUSALdb --> Figure 2
- **InvestigateData_NetProps.R**: Investigate structural properties of SNP LD Networks constructed with 1KGP Phase 1 and 3 Data --> Figure 3
  
## Performance Evaluation
- **SNP_M_KFold_ROC.R**: Run KFold cross-validation for SNP LD Networks for Phases <- c(1, 3), LDs <- c("r208", "r2def")
  - Will generate files: "Summary_AUROCnPR_", method, "_Phase", phase, "_", ld, ".ld.csv"
- **SNP_H_KFold_ROC.R**: Run KFold cross-validation for Heterogeneous Networks for Phases <- c(1, 3), LDs <- c("r208", "r2def")

## Performance Assessment
- **Summarize performance assessment**: Follow steps below
    1. Summarize AUROC and AUPR: Summarize_AUROCnPR_All.R --> combined_summary_data.csv
    2. Estimate number of diseases and SNPs mapped to SNP LD Networks and Heterogeneous Networks: Map_MeSH2SNP_2_M.R and Map_MeSH2SNP_2_H.R (output: Map_MeSH2SNP_2_M.csv and Map_MeSH2SNP_2_H.csv, respectively), then use Map_MeSH2SNP_Combine_M-H.R to combine these .csv files --> Map_MeSH2SNP_2_M-H.csv 
    3. Combine 2 files: combined_summary_data.csv and Map_MeSH2SNP_2_M-H.csv by Summarize_AUROCnPR_All_w_sd.R --> combined_summary_with_sd.csv
    4. combined_summary_with_sd.csv --> combined_summary_with_sd_corrected.csv (replaced r2def and r208 by 0.2 and 0.8 in column ld, respectively)
- **Performance comparison**: Input: combined_summary_with_sd_corrected.csv
  - **Draw_Figure4.R**: Compare Performance of DisSNPNet on SNP LD and Heterogeneous Networks Across Chromosomes --> Figure 4
  - **Draw_Figure5.R**: Compare Performance of DisSNPNet on Heterogeneous Networks Across LD Thresholds --> Figure 5
  - **Draw_Figure6.R**: Compare Performance of DisSNPNet on Heterogeneous Networks Across Datasets with LD threshold r² ≥ 0.8 --> Figure 6
  - **Compare_Mean_by_tTest.R**: Perform Statistical Tests
 
## Prediction of Novel Disease-SNP Associations
  - **SNP_MH_Predict_Final.R**: Predict, select top k and collect evidence
  - **Summarize_topKEvidence_Final.R**: Visualize Number of Supporting Evidence for Top-Ranked SNPs Predicted by DisSNPNet on Heterogeneous Networks --> Figure 7
  - **Visualize_Evidence_Final.R**: Visualize Distribution of Top 30 Predicted SNPs Across Diseases and Chromosomes  --> Figure 8
  - **Enrich_Pathways_Final.R**: KEGG Pathway Enrichment for Genes Associated with Predicted SNPs --> Figure 9

