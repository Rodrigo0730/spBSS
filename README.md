# Spatial Independent Component Analysis for Heteroskedastic random fields
This project contains the code to reproduce the results from the paper "Spatial independent component analysis for heteroskedastic random fields" by Morales-Martínez, R., Nordhausen, K. and Ruiz-Gazen, A.

The goal of the simulations was to show how the novel methods introduced in the paper, spFOBI and spJADE, outperform the current models for independent component analysis in a spatial setting.

### Functionality
All the code was written in R and requires the packages found in requirements.txt file. Note that some packages are only needed if the MPI frameworks are done in Slurm. The structure is as follows:
* Scripts: contains the main scripts for the generation and analysis of the data.
* R: contains the helper R files used in the scripts.

### Important
The code uses a modified version of the R library spGARCH which can be installed by using the provided tarball. It can be downloaded under releases i the right sidebar.

### Authors
Morales-Martínez, R., Nordhausen, K. and Ruiz-Gazen, A.
### License
GNU GPLv3
### References
Otto, P., Schmid, W. & Garthoff, R. (2018). Generalised spatial and spatiotemporal autoregressive conditional heteroscedasticity. Spatial Statistics, 26, 125–145. 

Muehlmann C, Sipila M, Cappello C, De Iaco S, Nordhausen K, Taskinen S, Virta J (2025). SpatialBSS: Blind Source Separation for Multivariate Spatial Data.
