####################################################################################
###################  Code for producing plots and tables in  #######################
##################  Globally-centered autocovariance in MCMC   #####################
########################### Agarwal and Vats (2020)  ###############################
####################################################################################

All the figures and tables in the paper are produced serially using the script "one_run.R" and saved in the folder "AllOut". This script uses 
pre-obtained R objects to save simulation time. Each example folder contains the code for producing these objects, figures, and tables. 
These R objects have been copied to the "AllOut" folder. Following is the information regarding Figures (1-10), Tables (1-3), and their associated R objects.

######### Figure-1 uses:

- "AllOut/gaussian-two_chains.Rdata" produced by "Examples/GaussianMix/output.R"

######### Figure-2 uses:

- "AllOut/var-five_chains.Rdata" produced by "Examples/VAR-1/acf_plots.R"

######### Figure-3 uses:

- "AllOut/var-truth.Rdata" produced by "Examples/VAR-1/output.R"
- "AllOut/var-conv_data_min500_max50000.Rdata" produced by "Examples/VAR-1/simulations.R"

######## Table-1 uses:

- "AllOut/var-out_check.pts_freq1000.Rdata" produced by "Examples/VAR-1/simulations.R"

####### Figure-4 uses:

- "AllOut/boom-two_chains_sp_1_3_8.Rdata" produced by "Examples/Boomerang/output.R"
- "AllOut/boom-five_chains_1_3_8.Rdata" produced by "Examples/Boomerang/acf_plots.R"
- "AllOut/boom-two_chains_sp_1_10_7.Rdata" produced by "Examples/Boomerang/output.R"
- "AllOut/boom-five_chains_1_10_7.Rdata" produced by "Examples/Boomerang/acf_plots.R"

###### Table-2 uses:
- "AllOut/boom1-out_check.pts_m2_freq1000.Rdata" produced by "Examples/Boomerang/simulations.R"
- "AllOut/boom1-out_check.pts_m5_freq1000.Rdata" produced by "Examples/Boomerang/simulations.R"

###### Table-3 uses:
- "AllOut/boom2-out_check.pts_m2_freq1000.Rdata" produced by "Examples/Boomerang/simulations.R"
- "AllOut/boom2-out_check.pts_m5_freq1000.Rdata" produced by "Examples/Boomerang/simulations.R"

###### Figure-5 uses:
- "AllOut/boom1-conv_data_m5_min500_max1e+05.Rdata" produced by "Examples/Boomerang/simulations.R"
- "AllOut/boom2-conv_data_m5_min500_max50000.Rdata" produced by "Examples/Boomerang/simulations.R"

##### Figure-6 uses:
- "AllOut/sensor-five_chains.Rdata" produced by "Examples/SensorNetwork/acf_plots.R"

##### Figure-7 uses:
- "AllOut/sensor-conv_data_m5_min500_max2e+05.Rdata" produced by "Examples/SensorNetwork/simulations.R"

##### Supp: Figure-1 uses:
 - "AllOut/poisson-two_chains.Rdata" produced by "Examples/PoissonChange/acf_plots.R"

##### Supp: Figure-2 uses:
 - "AllOut/poisson-two_chains.Rdata" produced by "Examples/PoissonChange/acf_plots.R"

##### Supp: Figure-3 uses:
 - "AllOut/magnolia-two_chains.Rdata" produced by "Examples/MagnoliaNetwork/acf_plots.R"

##### Supp: Figure-4 uses:
 - "AllOut/var-neg_five_chains.Rdata" produced by "Examples/negVAR-1/acf_plots.R"

##### Supp: Figure-5 uses:
 - "AllOut/var-set1_truth.Rdata", "AllOut/conv_data_min500_max50000_set1.Rdata", and "AllOut/conv_data_min500_max50000_set2.Rdata" produced by "Examples/highdim_VAR-1/simulations.R"

##### Supp: Figure-6 uses:
 - "AllOut/var-100_chains.Rdata" produced by "Examples/morechainsVAR/acf_plots.R"