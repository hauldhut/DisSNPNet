# DisSNPNet: Predicting Disease-Associated SNPs Using Linkage Disequilibrium, Disease Similarity, and 1KGP Datasets with Evidence-Based Validation

![Workflow of DisSNPNet for Predicting Disease-Associated SNPs](https://github.com/hauldhut/DisSNPNet/blob/main/Figure1.png)

## Repo structure
- **Data**: Contains all data 
- **Code**: Contains all source code to reproduce all the results
- **Results**: To store simulation results
  - **Prediction**: To store top predicted SNPs and their evidence
- **Figure**: To store figures for performance comparison and evidence collection

## How to run
- Install R packages
  - *RandomWalkRestartMH, ontologyIndex, MeSHSim, phenoscanner, doParallel, ROCR, ggplot2, Metrics, hash*
- Download the repo
- Follow instructions in the folder **Code** to run
  
- *Note*: For large and large networks (i.e., heterogeneous networks for large-size chromosomes and with LD threshold >= 0.2), it is recommended to run on a multi-core and at least 16 GB RAM computer

