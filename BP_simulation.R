#=== Simulate according to a Poisson Branching Process?

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


rgeneration <- function(n, t, mean = 6.5, sd = .62, par){
  
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


simulation_BP_individual <- function(n_days, R0, seeds, Reporting_fraction = 1){
  
  # prepare vectors to store results
  n_parents <- rep(0, n_days)
  incidence <- rep(NA, n_days)
  
  #decide day of infections for the seeds
  n_parents <- rgeneration(n = seeds, t = 0, par = n_parents)
  
  for (day in 1:n_days){
    incidence[day] <- sum(rpois(n_parents[day], R0))
    n_parents <- rgeneration(n = incidence[day], t = day, par = n_parents)
  }
  
  n_parents <- n_parents[1:n_days]
  obs_inc <- rbinom(n = n_days, size = incidence, prob = Reporting_fraction)
  return(data.frame(day = 1:n_days, n_parents = n_parents, obs_incidence = obs_inc))
}


simulation_BP <- function(N, n_days = 50 , R0 = 3, seeds = 5, Reporting_fraction = 1){
  
  # Prepare list to store all the simulations
  output <- data.frame()
  
  for (i in 1:N){
    output <- rbind(output, simulation_BP_individual(n_days = n_days, R0 = R0, seeds = seeds, Reporting_fraction = Reporting_fraction)$obs_incidence)
  }
  
  output <- t(output)
  colnames(output) <- paste0("sim", 1:N)
  rownames(output) <- c()
  
  return(output)
  
}
