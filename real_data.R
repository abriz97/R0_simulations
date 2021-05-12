require(data.table)
require(EpiEstim)
require(lubridate)
require(ggplot2)
require(ggpubr)

# set up directories' paths
repo_path <- '~/Documents/mini_project_2/R0_simulations'
setwd(repo_path)
data_path <- '~/Documents/mini_project_2/COVID-19/csse_covid_19_data/csse_covid_19_time_series/'
plots_path <- '~/Documents/mini_project_2/plots'



global_data <- read.csv(file.path(data_path, 'time_series_covid19_confirmed_global.csv'))
global_data <- as.data.table(global_data)
dim(global_data) # 275 468

colnames(global_data)
# we have 2 columns concerning Country and Regions
# 2 with Lat and Long which i can forget about
# The rest are dates, written as X'month'.'day'.'year'
global_data[, c("Lat", "Long") := NULL,]

# == extracts legend

library(cowplot)
library(ggpubr)



# == EXCTRACT COUNTRY-SPECIFIC DATA FROM GLOBAL DASTASET ==

extract_provinces <- function(country){
  global_data[Country.Region == country]$Province.State
}

extract_by_country <- function(country, province = ''){
  
  # Extract data dropping Province/State & Country/region columns
  if (province != ''){
    data <- global_data[ Country.Region ==  country & Province.State == province, -c(1,2)] 
  }else{
    data <- global_data[ Country.Region ==  country & Province.State == '', -c(1,2)]
  }
  
  if (dim(data)[1] != 1){
    print('Wrong dimensions: specify province')
    return(data)
  }
  
  # Get data in the right format
  data <- melt(data,
               variable.name = "date",
               value.name = 'cases')
  
  # fix dates 
  data$date <- gsub(pattern = '\\.', replacement = '/',data$date)
  data$date <- gsub(pattern = 'X', replacement = '',data$date)
  mdy(data$date)
  
  # and fix incidence 
  f <- function(x){
    y = c(0,x[1:(length(x)-1)])
    return(x - y)
  }
  data$cases <- f(data$cases)
  
  return(data)
}


# == MAKE FUNCTION THAT COMPARES WITH AND WITHOUT CORRECTION ==

EpEstimComparison <- function(data, WINDOW = 7, eg_window = 7, sd = SD, mean = MEAN){
  
  last_date <- tail(data$date,1)
  data$day <- seq_along(data$cases)
  data$date <- NULL
  
  # apply normal EpiEstim
  t_end <- (1 + WINDOW):(dim(data)[1])
  t_start <- t_end - WINDOW + 1
  
  epestim <- EpiEstim::estimate_R(data$cases, method = "parametric_si",
                                  config = make_config(
                                    list(
                                      mean_si = MEAN,
                                      std_si = sd,
                                      t_start = t_start,
                                      t_end = t_end
                                    )))
  
  # extract S
  S <- max(which(epestim$si_distr > 0))
  
  #and then get only what we want:
  epestim <- epestim$R[,c("Mean(R)", "Quantile.0.025(R)", "Quantile.0.975(R)")]
  
  # apply correction
  tmp <- data[1:eg_window,]
  
  if(any(tmp$cases==0)){
    tmp$casesC <- tmp$cases+0.5
    tmp$logcases <- log(tmp$casesC)
  }else{
    tmp$logcases <- log(tmp$cases)
  }
  
  # EG_growth
  exp_fit <- lm(logcases ~ day, data = tmp)
  #1/MGF_integral(exp_fit$coefficients[2], mean_GT = mean_si, sd_GT = sd_si)
  #1/MGF_integral(confint(exp_fit, level=0.95)[2, c(1,2)], mean_GT = mean_si, sd_GT = sd_si)
  
  # introduce extrapolated data
  day0 <- data$day[1] - 1
  tmp <- data.frame(day =  (-S):0 + day0,
                    logcases = predict.lm(exp_fit, newdata = data.frame(day = (-S):0  + day0 )))
  
  if (any(data[1:eg_window,]$cases==0)){
    tmp$casesC <- exp(tmp$logcases)
    tmp$cases <- tmp$casesC - 0.5
    tmp$cases[tmp$cases < 0] <- 0
  }else{
    tmp$cases <- exp(tmp$logcases)
  }
  
  data <- rbind(tmp[,c('cases', 'day')], data)
  
  # adjust windows: use the fact that tail data is the same:
  t_end <- rev(dim(data)[1] - seq_along(t_end) + 1)
  t_start <- t_start <- t_end - WINDOW + 1
  
  epestim_adj <- EpiEstim::estimate_R(c(data$cases), method = "parametric_si",
                                      config = make_config(
                                        list(
                                          mean_si = MEAN,
                                          std_si = sd,
                                          t_start = t_start,
                                          t_end = t_end
                                        )))$R[,c("Mean(R)", "Quantile.0.025(R)", "Quantile.0.975(R)")]
  
  epestim$date <- mdy(last_date) - rev(seq_along(t_end)-1)
  epestim_adj$date <- mdy(last_date) - rev(seq_along(t_end)-1)
  
  epestim$adjustment = FALSE
  epestim_adj$adjustment = TRUE
  
  tmp <- rbind(epestim, epestim_adj)
  colnames(tmp)[c(1,2,3)] <- c('M', 'l', 'u')
  
  return(tmp)
  
}


# Parametrisations:
# Challen et al.  'generation interval with a mean 4.8 (95% CrI 4.3–5.41) and SD1.7 (95% CrI 1.0–2.6)'

SOURCE = 'Challen et al';MEAN = 4.8; SD = 1.7

SOURCE = 'Ganyani et al';MEAN = 5.20; SD = 1.72
SOURCE = 'Ganyani et al';MEAN = 3.95; SD = 1.51

SOURCE = 'Megan par';MEAN = 6.5; SD = 2.6


SOURCE = 'CV1';MEAN = 6; SD = 3
SOURCE = 'CV2';MEAN = 4; SD = 2


title <- paste0(SOURCE, ', mean = ', MEAN, ', sd = ', SD)




# ITALY

data <- extract_by_country('Italy')
data <- data[31:50,]
pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

italy_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)

gg <- ggplot(italy_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') + xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "bottom")

my_legend <- get_legend(gg)

gg <- gg + theme(legend.position = "none")

italy <- ggarrange(pp, gg,
                   widths = c(1, 2))
italy <- annotate_figure(italy, top = 'Italy')

# UK 

data <- extract_by_country('United Kingdom')[33:60]

pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() + 
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

UK_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)


gg <- ggplot(UK_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') + xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none")

UK <- ggarrange(pp, gg,
                widths = c(1, 2))
UK <- annotate_figure(UK, top = 'United Kingdom')

# FRANCE

data <- extract_by_country('France')
data <- data[35:60]
pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

france_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)


gg <- ggplot(france_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') + xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none")

France <- ggarrange(pp, gg,
                    widths = c(1, 2))
France <- annotate_figure(France, top = 'France')

# US

data <- extract_by_country('US')
data <- data[39:69]
pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')


US_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)


gg <- ggplot(US_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') +  xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none")

US <- ggarrange(pp, gg,
                widths = c(1, 2))
US <- annotate_figure(US, top = 'United States')




# Germany
data <- extract_by_country('Germany')
data <- data[35:70]
pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

Ger_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)

gg <- ggplot(Ger_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') +  xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none")

Germany <- ggarrange(pp, gg,
                widths = c(1, 2))
Germany <- annotate_figure(Germany, top = 'Germany')

# Sweden
data <- extract_by_country('Sweden')
data <- data[42:80]
pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

Swe_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)

gg <- ggplot(Swe_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') +  xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none")

Sweden <- ggarrange(pp, gg,
                     widths = c(1, 2))
Sweden <- annotate_figure(Sweden, top = 'Sweden')


# Spain 

data <- extract_by_country('Spain')
data <- data[35:65]
pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

Spain_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)

gg <- ggplot(Spain_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') +  xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none")

Spain <- ggarrange(pp, gg,
                    widths = c(1, 2))
Spain <- annotate_figure(Spain, top = 'Spain')           


# Netherlands

data <- extract_by_country('Netherlands')
data <- data[39:69]
pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

Neth_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)

gg <- ggplot(Neth_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') +  xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none")

Netherlands <- ggarrange(pp, gg,
                   widths = c(1, 2))
Netherlands <- annotate_figure(Netherlands, top = 'Netherlands')  


# Denmark
data <- extract_by_country('Denmark')
data <- data[47:77]
pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

Denmark_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)

gg <- ggplot(Denmark_res, aes(x = date, y = M, color = adjustment, fill = adjustment)) +
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  ylab('Rt') +  xlab('') +
  theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position = "none")

Denmark <- ggarrange(pp, gg,
                         widths = c(1, 2))
Denmark <- annotate_figure(Denmark, top = 'Denmark') 



## EVERYTHING:

everything <- ggarrange(italy, France, UK, US, Germany, Sweden, Spain, Netherlands, Denmark,  common.legend = T) 
everything <- ggarrange(everything, as_ggplot(my_legend), heights = c(97,3), nrow = 2) 
annotate_figure(everything, top = title)


final <- ggarrange(italy, UK, Sweden, US, common.legend = T) 
final <- ggarrange(final, as_ggplot(my_legend), heights = c(95,5), nrow = 2)



# Change in reporting correction.
require(rootSolve)

mean_estimates <- italy_res$M
mean_estimates <- seq(from = 1.1, to = 8, by = .1)

alpha <- (.10 - .05)/30
rho_t <- .5 + (7 + seq_along(mean_estimates)) * alpha

A1 <- mean_estimates^(-1)
A2 <- alpha * MEAN / rho_t 
power <- -(SD/MEAN)^2

fun <- function(x, i){y <- A1*x + A2*x^power - 1; return(y[i])} 

new_estimates <- rep(0, length(mean_estimates))
for(i in seq_along(mean_estimates)){
  f <- function(x){fun(x,i)} 
  new_estimates[i] <- as.numeric(uniroot(f, lower = 0.1, upper = 10))
}

