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

integrate_ifcan2 <- function(r, mean_GT, sd_GT){
  
  shape=(mean_GT^2/sd_GT^2); scale=(sd_GT^2/mean_GT)
  
  if( r > - (sd_GT^2/mean_GT)^(-1))
  {
    return((1 + scale*r)^shape)
  }else{print("integrate_ifcan2: Unable to estimate R0 with current r")}
}


#=== Linear exponential growth method
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
  
  # as we are always using a Gamma distributed generation interval, let s avoid to do numerical integration
  # and let s just use closed formed identity for MGF
  if (rciL > - (sd_GT^2/mean_GT)^-1 ){
    
    R <- integrate_ifcan2(r, mean_GT, sd_GT)
    RciL <- integrate_ifcan2(rciL, mean_GT, sd_GT)
    RciU <- integrate_ifcan2(rciU, mean_GT, sd_GT)
    
    # calculate R by integrating over the GT distribution, g(a)
    # R <- 1/integrate(Mz_ga, r=r, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    
    # calculate 95% CI for R
    # RciL <- 1/integrate(Mz_ga, r=rciL, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    # RciU <- 1/integrate(Mz_ga, r=rciU, lower=0, upper=Inf, subdivisions=1e8, mean_GT=mean_GT, sd_GT=sd_GT)$value
    
    results <- c(R, RciL, RciU)
    
  }else{
    print("EG_P: Unable to estimate R0 with current data")
    results <- c(NA, NA, NA)
  }
  return(results)
}



#=== Standard EpiEstim method

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

#get estimated R0

#=== My EpiEstim
epiestim_adjust <- function(data, mean_GT, sd_GT, window_length = 14, instantaneous = TRUE, adjust = TRUE){
  
  # How does the data need to be formatted?
  if(!all(colnames(data) %in% c('cases', 'day'))){stop('data may only contain columns "cases" and "day"')} 
  
  # save first day of data for sanity checks later
  min_day <- min(data$day)
  
  # log transform cases, add small constant to time-series if zeros are present
  if(any(data$cases==0)){
    data$casesC <- data$cases+0.5
    data$logcases <- log(data$casesC)
  }else{
    data$logcases <- log(data$cases)
  }
  
  # Fit line through logged cases
  exp_fit <- lm(data=data, logcases ~ day)
  
  
  # Find estimated start ( ie cases = 1  => logcases = 0 or 0.5 depending on correction)
  # done by finding intersection of fitted line with corresponding horizontal line
  if ("casesC" %in% names(data)){
    start <- round((.5 -exp_fit$coefficients[1])/exp_fit$coefficients[2])
  }else{
    start <- -round(exp_fit$coefficients[1]/exp_fit$coefficients[2])
  }
  start <- as.numeric(start)
  
  if(start < -10000){
    warnings("Estimated start too far into past: R0 estimation would be inappropriate")
    results <- c(NA, NA, NA)
    return(list(results = results, estimated_start = NA))
  }
  
  #== Megan's extension
  
  # define t_start for EpiEstim purposes
  t_start <- 1
  
  if (start < min_day){
    
    # introduce extrapolated data and add them to data
    tmp <- data.frame(day = start:min_day)
    tmp$logcases <- predict.lm(exp_fit, newdata = data.frame(day = start:min_day ))
    
    if ("casesC" %in% names(data)){
      tmp$casesC <- exp(tmp$logcases)
      tmp$cases <- tmp$casesC - 0.5
    }else{
      tmp$cases <- exp(tmp$logcases)
    }
    
    data <- rbind(tmp, data)
    
    # and inform EpiEstim we have introduced data prior to what the days we want to estimate for.
    if (instantaneous){t_start <- dim(tmp)[1] + t_start}
  }
  
  I_t <- data$cases
  
  #== adjust's correction :get the first day, extend it 4 generations times to be sure, and dthen apply eE
  if (adjust){
    
    fake_start <- round(min(start, min_day) - 4*mean_GT)
    tmp <- data.frame(day = fake_start:min(start, min_day-1))
    # Previous cases are x_0*gamma^{t} (t being negative)
    tmp$logcases <- predict.lm(exp_fit, newdata = tmp)
    
    if ("casesC" %in% names(data)){
      tmp$casesC <- exp(tmp$logcases)
      tmp$cases <- tmp$casesC - 0.5
    }else{
      tmp$cases <- exp(tmp$logcases)
    }
    
    tmp <- filter(tmp, cases > 0)
    
    #get rea
    t_start <- dim(tmp)[1] + t_start
    I_t <- c(tmp$cases, I_t)
  }
  
  
  
  #== apply epiEstim

  # t_start can t be 1 for some reason, so define 
  t_end <- t_start + window_length - 1
  if (t_start == 1){t_start = 2}
  

  
  # Apply method
  tryCatch({
    R_EpiEstim <- EpiEstim::estimate_R(I_t, method="parametric_si", 
                                       config=make_config(list(mean_si=mean_GT, std_si=sd_GT,
                                                               t_start=t_start, t_end=t_end)))
  }, error=function(e){print(paste("EpiEstim_adjust: Unable to estimate R0 with current data"))})
  
  if(exists("R_EpiEstim")==TRUE){
    results <- c(R_EpiEstim$R$`Mean(R)`, R_EpiEstim$R$`Quantile.0.025(R)`, R_EpiEstim$R$`Quantile.0.975(R)`)
  } else {
    results <- c(NA, NA, NA)
  }
  return(list(results = results, estimated_start = start))
}



