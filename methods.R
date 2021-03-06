#=== Define the moment generating function of the GT distribution, g(a)
Mz_ga <- function(t, r, mean_GT, sd_GT){
  dgamma(t, shape=(mean_GT^2/sd_GT^2), scale=(sd_GT^2/mean_GT))*exp(-t*r)
}

#=== Helper function to avoid infinities in integrations
# doesn't seem to work as I want, let s forget about it for the moment.
integrate_ifcan <- function(method, r){
  tryCatch(
    {1/integrate(Mz_ga, r=r, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    return(TRUE)},
    error=function(e){
    print(paste0(method,": Unable to estimate R0 with current data (tryCatch)"))
    return(FALSE)
    })
}


#=== 1. Linear exponential growth method
EG_Lin <- function(data, mean_GT, sd_GT){
  
  # log transform cases, add small constant to time-series if zeros are present
  if(any(data$cases==0)){
    data$casesC <- data$cases+0.5
    data$logcases <- log(data$casesC)
  }else{
    data$logcases <- log(data$cases)
  }
  
  # define growth rate, r
  lm_mod <- lm(logcases~ day, data=data)
  r <- coef(lm_mod)[2]
  rciL <- confint(lm_mod, level=0.95)[2,1]
  rciU <- confint(lm_mod, level=0.95)[2,2]

  #if(integrate_ifcan("EG_Lin", rciL)){ 
  if (rciL > -.37){
    
    # calculate R by integrating over the GT distribution, g(a)
    R <- 1/integrate(Mz_ga, r=r, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    
    # calculate 95% CI for R
    RciL <- 1/integrate(Mz_ga, r=rciL, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    RciU <- 1/integrate(Mz_ga, r=rciU, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    
    results <- c(R, RciL, RciU)
    
  }else{
    print("EG_P: Unable to estimate R0 with current data")
    results <- c(NA, NA, NA)
  }
  return(results)
}



#=== 2. Poisson exponential growth method
EG_P <- function(data, mean_GT, sd_GT){
  
  # define growth rate, r
  P_mod <- glm(cases~ day, data=data, family="poisson")
  r_p <- coef(P_mod)[2]
  rciL <- confint(P_mod, level=0.95)[2,1]
  rciU <- confint(P_mod, level=0.95)[2,2]
  
  if(rciL > -.37){
    
    # calculate R by integrating over the GT distribution, g(a)
    R_p <- 1/integrate(Mz_ga, r=r_p, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    
    # calculate 95% CI for R
    RciL <- 1/integrate(Mz_ga, r=rciL, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    RciU <- 1/integrate(Mz_ga, r=rciU, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    
    results <- c(R_p, RciL, RciU)
  
  }else{
    print("EG_P: Unable to estimate R0 with current data")
    results <- c(NA, NA, NA)
  }
  return(results)
}


#=== 3. Maximum likelihood exponential growth method

EG_MLE <- function(data, mean_GT, sd_GT){
  
  # define log-likelihood (negative required by mle2)
  ll <- function(I0, r){
    -sum(data$cases*(r*data$day+log(I0))-I0*exp(r*data$day)-lfactorial(data$cases))
  }
        
  tryCatch({
    est <- mle2(ll, start=list(I0=1, r=0.2))
    r_mle <- coef(est)[2]
    ci_mle <- c(confint(est, quietly=TRUE)[2,1], confint(est, quietly=TRUE)[2,2])
  }, error=function(e){print(paste("EG_MLE: Unable to estimate R0 with current data"))})
  
  # check the validity of results (i.e. not NA & are above min_r)
  if(is.na(ci_mle[1])+is.na(ci_mle[2])==0 & ci_mle[1]>-0.38){
    
    # calculate R by integrating over the GT distribution, g(a)
    R <- 1/(integrate(Mz_ga, r=r_mle, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT))$value
    R_CIL <- 1/(integrate(Mz_ga, r=ci_mle[1], lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT))$value
    R_CIU <- 1/(integrate(Mz_ga, r=ci_mle[2], lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT))$value
    
    results <- c(R, R_CIL, R_CIU)
    
  }else{
    results <- c(NA, NA, NA)
  }
  return(results)
}


#=== 4. EpiEstim method

epiestim <- function(data, mean_GT, sd_GT){
  
  tryCatch({
    R_EpiEstim <- EpiEstim::estimate_R(data$cases, method="parametric_si", 
                                       config=make_config(list(mean_si=mean_GT, std_si=sd_GT,
                                                               t_start=2, t_end=as.numeric(nrow(data)))))
  }, error=function(e){print(paste("EpiEstim: Unable to estimate R0 with current data"))})
  
  if(exists("R_EpiEstim")==TRUE){
    results <- c(R_EpiEstim$R$`Mean(R)`, R_EpiEstim$R$`Quantile.0.025(R)`, R_EpiEstim$R$`Quantile.0.975(R)`)
  } else {
    results <- c(NA, NA, NA)
  }
  return(results)
}


#=== 5. White and Pagano method

WP <- function(data, GTd){
  
  tryCatch({
    R0_WP <- R0::est.R0.ML(data$cases, GT=GTd, begin=1, end=as.numeric(nrow(data)))
  }, error=function(e){print(paste("White & Pagano: Unable to estimate R0 with current data"))})
  
  if(exists("R0_WP")==TRUE){
    results <- c(R0_WP$R, R0_WP$conf.int)
  }else{
    results <- c(NA, NA, NA)
  }
  return(results)
}

#=== 6. Wallinga and Teunis method (time-varying)

WT <- function(data, GTd){
  
  tryCatch({
    
    # estimate time-dependant R0
    R0_WT <- R0::est.R0.TD(data$cases, GT=GTd, correct=TRUE, begin=1, end=as.numeric(nrow(data)))
    
    # smooth the time-dependant estimates across period of interest
    R0_WTav <- R0::smooth.Rt(R0_WT, time.period=(nrow(data)-1))
    
  }, error=function(e){print(paste("Wallinga & Teunis: Unable to estimate R0 with current data"))})
  
  if(exists("R0_WTav")==TRUE){
    results <- c(R0_WTav$R, R0_WTav$conf.int)
  }else{
    warnings("Wallinga & Teunis: R0 estimate not inside confidence interval")
    results <- c(NA, NA, NA)
  }
  return(results)
}

#== 7. Bettencourt & Ribeiro method

BR <- function(data, GTd){
  
  tryCatch({
    
    # estimate time-dependant R0
    R0_BR <- R0::est.R0.SB(data$cases, GT=GTd, begin=1, end=as.numeric(nrow(data)))

    # smooth the time-dependant estimates across period of interest
    R0_BRav <- R0::smooth.Rt(R0_BR, time.period=(nrow(data)-1))
    
  }, error=function(e){print(paste("Bettencourt & Ribeiro: Unable to estimate R0 with current data"))})
  
  if(exists("R0_BRav")==TRUE){
    results <- c(R0_BRav$R, R0_BRav$conf.int)
  }else{
    results <- c(NA, NA, NA)
  }
  return(results)
}


#=== Helper function to get methods names 
get_methods <- function() c("EG_Lin", "EG_P", "EG_MLE", "EpiEstim", "WP", "WT", 'BR')


