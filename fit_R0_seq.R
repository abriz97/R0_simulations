#=== Function to estimate time-constant R0 with each method at sequential ===#
#=== time points in the epidemic growth phase and output results as csv's ===#

#FOT TESTING:
# data=sim_ij
# mean_GT = 6.5
# sd_GT = 0.62
# GTd=GTd
# GT_week= 1
# methods = get_methods


fit_R0_seq <- function(data, mean_GT, sd_GT, GTd, GT_days,
                       methods = c("EG_Lin", "EG_P", "EG_MLE", "EpiEstim", "WP", "WT", 'BR')){
  
  # list to store results
  store <- list()
  
  # define peak - max time point of highest reported cases
  peak_idx <- which.max(data$cases)[length(which.max(data$cases))] 
  peak <- data$day[peak_idx]
  
  # What are we doing exactly here? Can this be moved to the simulating part?
  if(peak > min(data$day) + GT_days* 2){
    
    # define sections: from 2 generation times (approximated in terms of weeks) on, up to peak
    # maybe change: from first observation.
    sections <- seq(from= min(data$day) + GT_days* 2, to=peak, by = GT_days)
    
    # Loop to fit each method to each section with increasing number of time points
    for(t in 1:length(sections)){
      
      # extract section of the epidemic curve for fitting
      # Maybe insert something here, substituting the 1 with some other parameter.
      s_t <- filter(data, day <= sections[t])
      
      # fit each method
      R_est <- c()
      if("EG_Lin" %in% methods){
        eglin <- EG_Lin(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT);
        R_est <- rbind(R_est, eglin)
      }
      if("EG_P" %in% methods){
        egp <- EG_P(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT);
        R_est <- rbind(R_est, egp)
        }
      if("EG_MLE" %in% methods){
        egmle <- EG_MLE(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT);
        R_est <- rbind(R_est, egmle)
        }
      if("EpiEstim" %in% methods){
        epiest <- epiestim(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT);
        R_est <- rbind(R_est, epiest)
        }
      if("WP" %in% methods){
        wp <- WP(data=s_t, GTd=GTd);
        R_est <- rbind(R_est, wp)
        }
      if("WT" %in% methods){
        wt <- as.numeric(WT(data=s_t, GTd=GTd));
        if (!is.na(wt) & wt[1] <= wt[3] && wt[1] >= wt[2]){
          R_est <- rbind(R_est, wt)
        }else{
          R_est <- rbind(R_est, c(NA, NA, NA))
          }
        }
      if("BR" %in% methods){
        br <- BR(data=s_t, GTd=GTd);
        R_est <- rbind(R_est, br)}
      
      # store results
      store[[t]] <- cbind(paste(data$country[1]), sections[t], t + 1,methods, R_est, peak)
      colnames(store[[t]]) <- c("Country", "Ndays","#GTofData", "method", "R0", "CI_L", "CI_U", "peak")
    }
    
    # bind results by country and output as csv
    results_all <- do.call("rbind", store)
    return(results_all)
  }
}


