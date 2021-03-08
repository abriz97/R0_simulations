# Want to run similar analyses having missed the first generations
# Also, want to do make sure we use daily incidence as to 

library(dplyr)
library(bbmle)
library(EpiEstim)
library(R0)
library(ggplot2)
library(ggridges)
library(deSolve)
library(truncnorm)

# specify repository path
repo_path <- '~/Documents/mini_project_2/R0-methods-comparison/'
setwd(repo_path)

# source functions for estimating R0 & data simulation
source('methods.R')
source('fit_R0_seq.R')
source('data_simulation.R')
source('assess_performance.R')


#===== 1. Fit to simulated data with varying levels of random noise =====#

# simulate 250 datasets, assuming average infectious and incubation periods of 6 and 14 days
simulations <- simulate_data(N=20, gamma=1/6, sigma=1/14)

# Can add some options in simulate data, for example the number of missed generations of infections,
# and then something that allows for increasing reporting rate.
#simulations <- simulate_data_missed_generations(N=250, gamma=1/6, sigma=1/14)

# define the generation time distribution (parameters are mean and standard deviation)
GTd <- R0::generation.time("gamma", c(20/7, 7.4/7))
# GTd <- R0::generation.time("gamma", c(5.2, 1.72))

# list to store results from simulations
Sim_results <- list()

get_methods

# fit each method at sequential time points of epidemic growth phase
# note: running this locally in a loop will take a long time.
# where possible, we recommend using parallel computing for this.

# for each of the 3 noise levels
for(i in 1:3){
  
  # list for storage
  results_i <- list()
  
  # for each of the 250 simulations
  for(j in 1:ncol(simulations$sims)){
    
    # data for fitting - simulation j with noise level i
    sim_ij <- data.frame(week=1:nrow(simulations[[i]]), cases=simulations[[i]][ ,j],
                         country=paste("sim", j, sep='_'))
    
    # maybe remove the initial weeks here
    sim_ij <- sim_ij %>% filter(week > 6)
    
    
    # fit & store results
    results_i[[j]] <- fit_R0_seq(data=sim_ij, mean_GT=20/7, sd_GT=7.4/7, GTd=GTd, GT_week=3,
                                 methods = c("EG_Lin","EG_P", "EG_MLE","EpiEstim", "WP" , "WT", "BR"))
  }
  
  # store results from each noise level
  Sim_results[[i]] <- results_i
}

i = 2
# calculate performance metrics
metrics <- performance_metrics(trueR0=simulations$pars, estimates=Sim_results[[i]], max_weeks=15, min_peak=15)

# choose metrics to plot
get_metrics()
want <- c('Bias','Coverage','Uncertainty','RMSE')

# metrics summary plot & bias plot
SummPlot <- summary_plot(metrics$metrics_summ, include=want)
BPlot <- bias_plot(metrics$metrics_indiv)
