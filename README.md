## R0-methods-comparison

Code and data to reproduce the analysis in "A comparative analysis of statistical methods to estimate the reproduction number in emerging epidemics with implications for the current COVID-19 pandemic", https://doi.org/10.1093/cid/ciaa1599.

***Code***
- 'report_simulations.Rmd' contains the code to run the simulation experiments.
- 'real_data.R' containst the code to run the analysis on COVID data.
- 'methods.R' contains functions for a variety of methods that will return R0 estimates and associated 95% confidence intervals.
- 'produce_simulations.R' makes use of two other scripts in order to generate simulated epidemic curves, save them in the data folder and visualise them.
   - 'simulation.R' contains the code to produce stochastic SEIR simulation.
   - 'BP_simulation.R' produces simulations according to a Branching Process.
 

***Datasets***
- The COVID data used is available at: https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
