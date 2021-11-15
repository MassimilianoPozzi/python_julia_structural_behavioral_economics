# Structural Behavioral Economics with Python and Julia

### Introduction

Here you will find Python and Julia notebooks to replicate the structural estimates of a set of papers in the field of Behavioral Economics. Following DellaVigna (*Handbook of Behavioral Economics*, 2018), we define structural as *"the estimation of a model on data that recovers estimates (and confidence intervals) for some key behavioral parameters"*. [Python](https://www.python.org) and [Julia](https://julialang.org) are two programming languages that offer some perks compared to others often used by economists, since they are:

1. Free
2. Open source 
3. Fast!

We hope this work can help and ease the transition for future researchers' moving from Stata, R and Matlab to Python and Julia.

### Organization

This repository contains one folder for each paper whose structural estimates our notebooks are intended to replicate. The naming convention we use for each of these folders is authors-journal-year where "journal" is the journal where the paper was published and "year" is the year of publication. Each folder contains three sub-folders named: code, input and output. Code contains the notebooks to replicate the main estimates of the  papers. The Python notebook has extension .IPYNB while the Julia notebook has extension .jl. Input contains all the data used in the notebooks to replicate the estimates. Output contains the results obtained running the notebooks in .csv format.

### Papers

As of right now this repository contains only the notebooks to replicate the aggregate estimates of Bruhin, Fehr, and Schunk (*Journal of the European Economic Association*, 2019). However, more papers will be uploaded soon. We have already written and tested the code for the following papers: 

<br>

- Augenblick, Ned and Matthew Rabin, 2019, "An Experiment on Time Preference and Misprediction in Unpleasant Tasks", *Review of Economic Studies*, 86(3): 941&ndash;975; Estimation Methods: Maximum Likelihood (with Bootstrapped Standard Errors)
- Bruhin, Adrian, Ernst Fehr, and Daniel Schunk, 2019, "The mMny Faces of Human Sociality: Uncovering the Distribution and Stability of Social Preferences", *Journal of the European Economic Association*, 17(4): 1025&ndash;1069; Estimation Methods: Maximum Likelihood and EM-Algorithm for Mixture Model
- DellaVigna, Stefano and Devin Pope, 2018, "What Motivates Effort? Evidence and Expert Forecasts", *Review of Economic Studies*, 2018, 85(2): 1029&ndash;1069; Estimation Methods: GMM (with Bootstrapped Standard Errors) and Non-Linear-Least-Squares
- DellaVigna, Stefano, John A. List, Ulrike Malmendier, and Gautam Rao, 2017, "Voting To Tell Others", *Review of Economic Studies*, 84(1): 143&ndash;181; Estimation Methods: Simulated Method of Moments			

We chose these papers to replicate because (a) they were recently published in general interest journals and, as such, they constitute the state-of-the-art in structural behavioral economics, (b) the replication files are available on the journal's website, and (c) they allow us to explore a wide array of different estimation methods and modeling of heterogeneity (the two building blocks of a structural model, as discussed in Sections 4.1 and 4.2 of DellaVigna 2018).
 
