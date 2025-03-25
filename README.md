# DisSNPNet
Predicting Disease-Associated SNPs Using Linkage Disequilibrium, Disease Similarity, and 1KGP Datasets with Evidence-Based Validation

![Network construction](https://github.com/hauldhut/MHDR/blob/main/Figure1.png)

## Repo structure
- **Data**: Contains all data 
- **Code**: Contains all source code to reproduce all the results
- **Results**: To store simulation results
  - **Evidence**: To store collected evidence

## How to run
- Install R packages
  - *RandomWalkRestartMH, igraph, foreach, doParallel, ROCR, ggplot2, Metrics, hash*
- Download the repo
- Follow instructions in the folder **Code** to run
  
- *Note*: For large and complex networks (i.e., multiplex disease networks, and multiplex-hetergeneous networks of drugs and diseases), it is recommended to run on a multi-core and at least 16 GB RAM computer

