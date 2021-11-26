# Structural Behavioral Economics with Python and Julia

### Authors

- Massimiliano Pozzi (Bocconi University, pozzi.massimiliano@studbocconi.it)
- Salvatore Nunnari (Bocconi University, https://snunnari.github.io, salvatore.nunnari@unibocconi.it)

### Introduction

Here you will find Python and Julia notebooks to replicate the structural estimates of a set of papers in the field of Behavioral Economics. Following DellaVigna (*Handbook of Behavioral Economics*, 2018), we define structural as *"the estimation of a model on data that recovers estimates (and confidence intervals) for some key behavioral parameters"*. [Python](https://www.python.org) and [Julia](https://julialang.org) are two programming languages that offer some perks compared to others often used by economists, since they are:

1. Free
2. Open source 
3. Fast!

We hope this work can help and ease the transition for future researchers' moving from Stata, R and Matlab to Python and Julia.

### Organization

This repository contains one folder for each paper whose structural estimates our notebooks are intended to replicate. The naming convention we use for each of these folders is authors-year where "year" is the year of publication. Each folder contains five sub-folders named: paper, original_replication_files, code, input and output. The sub-folder "paper" contains .pdf files for the paper and its appendices. The sub-folder "original_replication_files" contains the replication datasets and codes downloaded from the journal's website. The sub-folder "code" contains the notebooks we created to replicate the structural estimates in the paper with Python and Julia; the Python notebook has extension .ipynb (and can be opened with [Jupyter](https://jupyter.org)) while the Julia notebook has extension .jl (and can be opened with [Neptune](https://github.com/compleathorseplayer/Neptune.jl) and [Pluto](https://github.com/fonsp/Pluto.jl)). The sub-folder "input" contains all the data imported in the notebooks to generate the estimates. The sub-folder "output" contains .csv files with the estimation results obtained running the notebooks.

### Papers

As of right now this repository contains only the notebooks to replicate the aggregate estimates of Bruhin, Fehr, and Schunk (*Journal of the European Economic Association*, 2019) and Augenblick, and Rabin (*Review of Economic Studies*, 2019). However, more notebooks and papers will be uploaded soon. We have already written and tested the code to replicate the structural estimates in the following papers: 

- Ned Augenblick and Matthew Rabin, 2019, "An Experiment on Time Preference and Misprediction in Unpleasant Tasks", *Review of Economic Studies*, 86(3): 941&ndash;975
	- Estimation Methods: Maximum Likelihood
	- Modeling of Heterogeneity: Implementation Errors
- Adrian Bruhin, Ernst Fehr, and Daniel Schunk, 2019, "The Many Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences", *Journal of the European Economic Association*, 17(4): 1025&ndash;1069
	- Estimation Methods: Maximum Likelihood and SEAM/EM-Algorithm for Mixture Model
	- Modeling of Heterogeneity: Random Utility, Finite Mixture
- Stefano DellaVigna and Devin Pope, 2018, "What Motivates Effort? Evidence and Expert Forecasts", *Review of Economic Studies*, 85(2): 1029&ndash;1069
	- Estimation Methods: GMM (with Bootstrapped Standard Errors) and Non-Linear-Least-Squares
	- Modeling of Heterogeneity: Random Utility for Non-Linear-Least-Squares
- Stefano DellaVigna, John A. List, Ulrike Malmendier, and Gautam Rao, 2017, "Voting To Tell Others", *Review of Economic Studies*, 84(1): 143&ndash;181
	- Estimation Methods: Simulated Method of Moments
	- Modeling of Heterogeneity: Random Parameters			

We chose these papers to replicate because (a) they were recently published in general interest journals and, as such, they constitute the state-of-the-art in structural behavioral economics, (b) the replication files are available on the journal's website, and (c) they allow us to explore a wide array of different estimation methods and modeling of heterogeneity (the two building blocks of a structural model, as discussed in Sections 4.1 and 4.2 of DellaVigna 2018 and detailed in the list above).
 
 ### References on Structural Models in Economics
 
 - DellaVigna, Stefano, 2018, ["Structural Behavioral Economics"](http://snunnari.github.io/dellavigna), *Handbook of Behavioral Economics Vol. 1*, North-Holland, 613&ndash;723
