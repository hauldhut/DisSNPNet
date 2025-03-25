![image](https://github.com/user-attachments/assets/318488f0-fe68-44dc-8881-1c278098adaa)![image](https://github.com/user-attachments/assets/92785278-ea84-46d9-b87a-e37f658eecfe)# DisSNPNet
Predicting Disease-Associated SNPs Using Linkage Disequilibrium, Disease Similarity, and 1KGP Datasets with Evidence-Based Validation

![Workflow of DisSNPNet for Predicting Disease-Associated SNPs](https://github.com/hauldhut/DisSNPNet/blob/main/Figure1.png)

## Repo structure
- **Data**: Contains all data 
- **Code**: Contains all source code to reproduce all the results
- **Results**: To store simulation results
  - **Evidence**: To store collected evidence

## How to run
- Install R packages
  - *RandomWalkRestartMH, ontologyIndex, MeSHSim, phenoscanner, doParallel, ROCR, ggplot2, Metrics, hash*
- Download the repo
- Follow instructions in the folder **Code** to run
  
- *Note*: For large and large networks (i.e., heterogeneous networks for large-size chromosomes and with LD threshold >= 0.2), it is recommended to run on a multi-core and at least 16 GB RAM computer

