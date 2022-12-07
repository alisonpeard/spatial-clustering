# Clustering Methods for Spatial Networks
Code for dissertation _Uncovering Unbiased Meso-scale Structures in Spatial Networks_ as part of the MSc in Mathematical Modelling and Scientific Computing 2020, Universtity of Oxford. This repository contains Python packages for performing spatially-corrected clustering on networks and methods for plotting the results spatially.<br>

Clustering methods include modularity optimisation by spectral methods, and spatial backbone extraction as in [2] followed by standard modularity optimisaiton and core-periphery detection. Additionally asymmetric core-periphery detection code for directed graphs is included [3]. Spatial benchmark graphs are generated using a number of different spatial null models [1].

## Files
1. `dissertation_noappendix.pdf`: Dissertation for MSc Mathematical Modelling and Scientific Computing (2020-2021)
2. `spatial-clustering/spatial_graphs/spatial_graphs.p`: Module for spatial graph classes which are children of NetworkX's Graph and DiGraph classes. Includes class methods for instantiating spatial graph instances from various spatial null models and node partitions.
3. `spatial-clustering/notebooks/Classic Modularity Visualisations.ipynb`: notebook for generating figures used in dissertation



## References
1. Spatial Correlations in Attribute Communities, Cerina et al. (2011), doi:10.1371/ journal.pone.0037507
2. [Community Detection in Customer Store Networks](https://www.maths.ox.ac.uk/study-here/postgraduate-study/industrially-focused-mathematical-modelling-epsrc-cdt/infomm-resear-46), Leal Cervantes (2022)
3. Core–periphery structure in directed networks, Elliot et al. (2020), http://dx.doi.org/10.1098/rspa.2019.0783
