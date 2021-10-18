# Replicated-Spectral-Variance-Estimator
 This repository contains simulations attesting to the better performance of globally-centered autocovariance and limiting variance estimators while using parallel Markov chains. The improvements are mostly studied graphically in form of ACF plots and ESS disgnostics. The following examples from Agarwal and Vats (2020) have been thorougly explored for a comparitative study:
 
 #### Gaussian mixture model
 #### Vector autoregressive model (VAR-1)
 #### Boomerang distribution
 #### Sensor network model
 #### Poisson change-point model
 #### Network crawling
 
 All the figures and tables in the paper are produced serially using the script "one_run.R" and saved in the folder "AllOut". This script uses 
pre-obtained R objects to save simulation time. Each example folder contains the code for producing these objects, figures, and tables. Read the readme.txt for detailed explanation. 
 
