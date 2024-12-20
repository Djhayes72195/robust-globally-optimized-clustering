# A Hybrid Approach to Globally Optimal Clustering

### Overview

This repository contains a novel clustering algorithm implemented in R. Its purpose is to find globally optimal solutions while minimizing computational costs. The algorithm balances speed and optimality by combining elements of the Expectation-Maximization (EM) algorithm for local optimization and the Gibbs sampling algorithm for global search.

For a detailed description of the algorithm, including theoretical insights and comparisons with classical clustering methods like K-means, soft clustering, and methods for selecting K, please refer to "MastersReport.pdf" included in this repository.

### Usage

To run the algorithm, you will need an installation of R. By default, the implementation is configured to use the wine dataset from the UCI Machine Learning Repository. You can easily adapt the configuration to use a different dataset, provided it meets the following conditions:

- All features must be numerical.
- There should be no class labels.

The following steps would be required to adapt the code to a new dataset:

- Modify this line to point to the .data file you wish to use.
    - `df <- read.table(file='wine.data', sep=",", header=FALSE)`
- Modify or eliminate this line to drop any class labels if they are present.
    - `df$V1 = NULL`
- Set an appropriate value for K. See "MastersReport.pdf" for guidance on selecting a suitable value for K.



### Contact Information

If you any any questions, reccomendations, or opportunities for collaboration, I would love to hear from you. Please feel free to reach out via email at Djhayes72195@gmail.com.


