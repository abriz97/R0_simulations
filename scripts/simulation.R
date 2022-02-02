#' is.wholenumber
#'
#' For parameter checking
is.wholenumber <- function(x) {
  x == round(x)
}

check_args <- function(
  arnaught, t_E, t_I, Npv, S_init, E_init, I_init, n_t, n_steps_per_t
) {
  # Check all input parameters
  stopifnot(is.wholenumber(Npv) && length(Npv) == 1 && Npv >= 1)
  stopifnot(is.wholenumber(n_t) && length(n_t) == 1 && n_t >= 1)
  stopifnot(
    is.wholenumber(n_steps_per_t) && length(n_steps_per_t) == 1 &&
      n_steps_per_t >= 1
  )
  stopifnot(
    is.numeric(arnaught) && arnaught > 0 &&
      (length(arnaught) == 1 || length(arnaught) == n_t + 1)
  )
  stopifnot(
    is.numeric(t_E) && length(t_E) == 1 && t_E >= 0
  )
  stopifnot(
    is.numeric(t_I) && length(t_I) == 1 && t_I > 0
  )
}

construct_beta <- function(arnaught, t_I, n_t) {
  beta_t_all <- arnaught / t_I
  if(length(arnaught) == 1) {
    function(t) beta_t_all
  } else {
    approxfun(0:n_t, beta_t_all)
  }
}

simulate_seir <- function(
  arnaught, t_E, t_I,
  Npv, S_init, E_init, I_init,
  n_t, n_steps_per_t = 1,
  method = 'stochastic'
) {
  func <- if(method == 'stochastic') simulate_seir_stochastic else simulate_seir_ode
  func(
    arnaught, t_E, t_I,
    Npv, S_init, E_init, I_init,
    n_t, n_steps_per_t
  )
}



# For more info on simulating a discrete-time approximation of a S(E)IR model, follow the link:
# https://cran.r-project.org/web/packages/odin/vignettes/discrete.html


#' Simulate a discrete-time approximation of a continuous-time, discrete-state stochastic S(E)IR model.
#' 
#' Ed Baskerville
#' 15 April 2020
#' 
#' No age structure.
#' 
#' @param n_t Number of units of time to simulate.
#' @param n_steps_per_t Number of discrete simulation steps to take per
#' unit of simulation time.
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0. Average number of new infections produced by an infection in a susceptible population. A scalar or a vector of length `n_t + 1`, which specifies R0 at the start of each timestep. R0 is linearly interpolated between timesteps.
#' @param t_E Mean latent period. If set to 0, the model reduces to an SIR.
#' @param t_I Mean duration of infectiousness.
# # for debugging
# arnaught
# Npv = Npv
# E_init = 0          # Initial in E
# I_init = 1e-5 * Npv     # Initial in I
# t_E = 3               # time from exposed to infected
# t_I = 5               # time from infected to recovery
# n_t = 1000            # total simulation time
# max_R0 = 2.0          # Initial, max R0
# min_R0 = 0.9         # Minimum with intervention
# end_max_time = 50     # time of intervention
# start_min_time = 70  #time when R0 hits lowest value
# n_steps_per_t = 10
# CONTINUOUS = TRUE       # If TRUE, R decreases according to arctan after intervention. If FALSE, R


simulate_seir_stochastic <- function(
  N,
  gamma,
  sigma,
  arnaught,
  alpha = 1,
  n_t = 105, 
  n_steps_per_t = 10,
  reporting_fraction = 1,
  min_peak = 0,
  N0_par = 10^6,
  I0_par = 10
) {
  
  if (N %% length(arnaught) != 0){return("N must be divisible by the lenght of arnaught")}
  if (alpha > 1 | alpha < 0){return("alpha must be a real number between 0 and 1")}
  
  # Model is given by the schema:
  # S   --(beta)->   E   --(sigma)->    I   --(gamma)->   R
  
  # Precompute a few things
  delta_t <- 1 / n_steps_per_t
  
  # Draws a binomial based on a rate
  draw <- function(n, rate) {
    p <- 1 - exp(-rate * delta_t)
    rbinom(1, n, p)
  }
  
  # Function to compute beta at a particular time
  beta <- construct_beta(arnaught, 1/gamma, n_t)
  
  
  # Step forward from t to t + delta_t
  step <- function(t, S_prev, E_prev, I_prev) {
    
    # SEIR model
    dS <- draw(S_prev, beta(t) * I_prev^(alpha) / Npv)
    dIR <- draw(I_prev, gamma)
    dEI <- draw(E_prev, sigma)
    list(
      S = S_prev - dS,
      E = E_prev + dS - dEI,
      I = I_prev + dEI - dIR,
      dS = dS,
      dEI = dEI,
      dIR = dIR
    )
  }
  
  ## NOW WE CAN START WITH THE SIMULATION

  # dataframes to store simulations and parameters data.
  sim_data <- data.frame(matrix(ncol=N, nrow = n_t - 1))
  colnames(sim_data) <- paste("sim", seq(1:N), sep="")
  pars_data <- data.frame(sim=1, R0=arnaught, Np=NA, I0=NA)
  
  index = 1
  n_tries = 0
  
  while(index != N + 1){
    
    # initial conditions
    Npv = N0_par
    I_init = I0_par
    S_init = Npv - I_init
    
    
    #count tries
    n_tries = n_tries + 1
    
    # Define parameters
    pars <- c(R0=arnaught, gamma=gamma, sigma=sigma)
    
    # Set up state vectors over time
    S <- numeric(n_t + 1)
    E <- numeric(n_t + 1)
    I <- numeric(n_t + 1)
    
    S[1] <- S_init
    E[1] <- 0
    I[1] <- I_init
    
    # Track transitions over time
    dS <- rep(NA, n_t)
    dEI <- rep(NA, n_t)
    dIR <- rep(NA, n_t)
    
    # Simulate
    for(tt in 1:(n_t-1)) {
      S_prev <- S[tt]
      E_prev <- E[tt]
      I_prev <- I[tt]
      
      # Do we really need to set these to 0?
      dS[tt+1] <- 0
      dEI[tt+1] <- 0
      dIR[tt+1] <- 0
      
      for(i in 1:n_steps_per_t) {
        state_next <- step(tt + delta_t * (i - 1), S_prev, E_prev, I_prev)
        S_prev <- state_next$S
        E_prev <- state_next$E
        I_prev <- state_next$I
        dS[tt+1] <- dS[tt+1] + state_next$dS
        dEI[tt+1] <- dEI[tt+1] + state_next$dEI
        dIR[tt+1] <- dIR[tt+1] + state_next$dIR
      }
      
      S[tt+1] <- S_prev
      E[tt+1] <- E_prev
      I[tt+1] <- I_prev
    }
    
    # remove first NA
    dEI <- dEI[-1]
    # extract the incidence of infections and multiply by the reporting fraction
    incidence <- rbinom(n = length(dEI), size = dEI , prob = reporting_fraction)
    
    if( max(incidence) != 0 && min_peak <= which.max(incidence) ){
      
      # fill in output dataframe with incidence and true parameters
      sim_data[1:length(incidence) ,index] <- incidence
      pars_data[index, ] <- c(index, arnaught, Npv, I_init)
      
      index = index + 1 
    }
    
  } 
  # return a list of simulated data, with varying noise levels, & parameter values used
  return(list(sims=sim_data, I=I, pars=pars_data, n_tries = n_tries))

}













#' Simulate a deterministic ODE approximation of a continuous-time, discrete-state stochastic S(E)IR model.
#' 
#' Ed Baskerville
#' 15 April 2020
#' 
#' No age structure.
#' 
#' @param arnaught Basic reproduction number (ratio), squiggly-R-0. Average number of new infections produced by an infection in a susceptible population. A scalar or a vector of length `n_t + 1`, which specifies R0 at the start of each timestep. R0 is linearly interpolated between timesteps.
#' @param t_E Mean latent period. If set to 0, the model reduces to an SIR.
#' @param t_I Mean duration of infectiousness.
#' @param n_t Number of units of time to simulate.
simulate_seir_ode <- function(
  arnaught, t_E, t_I,
  Npv, S_init, E_init, I_init,
  n_t,
  n_steps_per_t = 1 # Ignored; included so the function signature matches stochastic version
) {
  library(deSolve)
  
  check_args(
    arnaught, t_E, t_I, Npv, S_init, E_init, I_init, n_t, n_steps_per_t
  )
  
  beta <- construct_beta(arnaught, t_I, n_t)
  d_dt <- function(t, y, params) {
    dS <- y['S'] * beta(t) * y['I'] / Npv
    dIR <- y['I'] / t_I
    
    if(t_E > 0) {
      # SEIR model
      dEI <- y['E'] / t_E
      list(c(
        S = -dS,
        E = dS - dEI,
        I = dEI - dIR,
        R = dIR,
        cum_dS = dS,
        cum_dEI = dEI
      ), NULL)
    }
    else {
      # SIR model
      list(c(
        S = -dS,
        E = 0,
        I = dS - dIR,
        R = dIR,
        cum_dS = dS,
        cum_dEI = dS
      ), NULL)
    }
  }
  
  y_init <- c(
    S = S_init,
    E = if(t_E > 0) E_init else 0,
    I = if(t_E > 0) I_init else E_init + I_init,
    R = 0,
    cum_dS = 0,
    cum_dEI = 0
  )
  #automatic ode solver is lsoda, an "automatic stiff/non-stiff solver"
  as.data.frame(ode(y_init, 0:n_t, d_dt, NULL)) %>%
    mutate(dS = cum_dS - lag(cum_dS, 1)) %>%
    mutate(dEI = cum_dEI - lag(cum_dEI, 1)) %>%
    mutate(dIR = R - lag(R, 1))
}



##########################################################
#=== Simulate according to a Poisson Branching Process
##########################################################

# Each infected individual produces a random number of secondary infections simultaneously.
# The generation time is distributed according to a discrete distribution
# The distribution of offsprings is given by a Poisson(R_n)

# == Input

# R0
# Generation distribution


# == Settings

# 3 columns matrix: time - parents - incidence.
# - time denotes day
# - parents denotes the number of primary infections which will cause secondary infections on a given day
# - incidence denotes the number of new infections produced by parents on a given day

# == How?
# Parents ---(Poisson)--> Incidence ---(Generation time)--> Parents ...

# For day in 1: end:

# count the parents, and produce the new infections for the same day via Poisson realizations
# incidence[t] <- sum(rpoisson(n_parents[t], R0))

# take the new cases, and decide when they will cause new infections:
# n_parents <- n_parents + rgeneration(incidence[t], t)

rgeneration <- function(n, t, par, mean = 6.5, sd = 4.61)
{
  # reparametrise in terms of shape and scale
  sh = mean^2/sd^2
  sc = sd^2/mean
  tmp <- table(t + round(rgamma(n, shape = sh, scale = sc)))
  
  for(i in 1:length(tmp)){
    index <- as.numeric(names(tmp[i]))
    count <- as.numeric(tmp)[i]
    par[index] <- par[index] + count
  }
  
  # return summary of determined days of infections
  return(par)
}


simulation_BP_individual <- function(n_days, R0, seeds, reporting_fraction, mean, sd)
{
  
  # prepare vectors to store results
  n_parents <- rep(0, n_days)
  incidence <- rep(NA, n_days)
  
  #decide day of infections for the seeds
  n_parents <- rgeneration(n = seeds, t = 0, par = n_parents)
  
  for (day in 1:n_days){
    incidence[day] <- sum(rpois(n_parents[day], R0))
    n_parents <- rgeneration(n = incidence[day], t = day, par = n_parents, mean=mean, sd=sd)
  }
  
  n_parents <- n_parents[1:n_days]
  obs_inc <- rbinom(n = n_days, size = incidence, prob = reporting_fraction)
  return(data.frame(day = 1:n_days, n_parents = n_parents, obs_incidence = obs_inc))
}


simulation_BP <- function(N, n_days = 50 , R0 = 3, seeds = 5, reporting_fraction = 1, mean = 6.5, sd = 4.61)
{
  
  # Prepare list to store all the simulations
  output <- data.frame()
  
  for (i in 1:N){
    output <- rbind(output, simulation_BP_individual(n_days = n_days, R0 = R0, seeds = seeds, reporting_fraction = reporting_fraction, mean=mean, sd=sd)$obs_incidencei)
  }
  
  output <- t(output)
  colnames(output) <- paste0("sim", 1:N)
  rownames(output) <- c()
  
  return(output)
  
}






