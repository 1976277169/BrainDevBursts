# BrainDevBursts
Code for paper "Temporal ordering of input modulates connectivity formation in a developmental neuronal network model of the cortex" by C Hartley, SF Farmer and L Berthouze. 
BiorXiv reference to follow. 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1297410.svg)](https://doi.org/10.5281/zenodo.1297410)


### List of files: 

#### Main files: 

* **runSim.m**: Main simulation file. Runs one simulation for a given set of parameters

* **cluster_sweep.m**: Performs a systematic sweep of key parameters (Hurst exponents, thresholds, LTP/LTD parameters). 

* **calculateNormalisedNetworkStatistics.m**: Calculate normalised network metrics. Requires prior execution of runSim. 


#### Ancillary files: 

* **createRandomNetwork.m**: Generates a random (directed, no self-loop) network with a given number of links.

* **FDN.m**: Generates fractional differencing noise with parameter d

* **generateFDN.m**: A wrapping function for DFN specifying Hurst exponent. 

* **generateIBI_from_FDN.m**: Generates sequences of IBIs from FDN with given Hurst. 

* **generateControlIBIseqs.m**: Generates sequences of IBIs of various Hurst using noise values from a reference IBIseq with given Hurst. 

* **[progressbar](https://uk.mathworks.com/matlabcentral/fileexchange/6922-progressbar)**: Displays a progress bar. Code redistributed. See license.txt in directory. 


### Requirements: 

This code requires access to functions from the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/). Code was developed and tested with version 2017_01_15_BCT.  
