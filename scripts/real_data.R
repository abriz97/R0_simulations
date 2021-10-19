require(data.table)
require(EpiEstim)
require(lubridate)
require(ggplot2)
# == extracts legend
require(ggpubr)
library(cowplot)


# set up directories' paths
repo_path <- '~/Documents/mini_project_2/R0_simulations'
setwd(repo_path)
data_path <- file.path(repo_path, 'data')
plots_path <- '~/Documents/mini_project_2/plots'



######################################################
#############    COVID-19 Analysis       #############
######################################################

global_data <- read.csv(file.path(data_path, 'time_series_covid19_confirmed_global.csv'))
global_data <- as.data.table(global_data)
dim(global_data) # 275 468

colnames(global_data)
# we have 2 columns concerning Country and Regions
# 2 with Lat and Long which i can forget about
# The rest are dates, written as X'month'.'day'.'year'
global_data[, c("Lat", "Long") := NULL,]


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
                                      std_si = SD,
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
                                          std_si = SD,
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

# SOURCE = 'Challen et al';MEAN = 4.8; SD = 1.7

SOURCE = 'Ganyani et al';MEAN = 5.70; SD = 1.72 # (5.7, 1.72)
# SOURCE = 'Ganyani et al';MEAN = 3.95; SD = 1.51
# SOURCE = 'Megan par';MEAN = 6.5; SD = 2.6

# SOURCE = 'CV0';MEAN = 4; SD = 2
# SOURCE = 'CV1';MEAN = 6; SD = 3
# SOURCE = 'CV2';MEAN = 8; SD = 4


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

italy_res$loc <- 'Italy'
tmp <- which(italy_res$date == min(italy_res$date))
tmp1 <- which(italy_res$adjustment)
R0_final <- italy_res[intersect(tmp, tmp1),]

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


UK_res$loc <- 'UK'
tmp <- which(UK_res$date == min(UK_res$date))
tmp1 <- which(UK_res$adjustment)
R0_final <- rbind(R0_final, UK_res[intersect(tmp, tmp1),])

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

france_res$loc <- 'France'
tmp <- which(france_res$date == min(france_res$date))
tmp1 <- which(france_res$adjustment)
R0_final <- rbind(R0_final, france_res[intersect(tmp, tmp1),])

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

US_res$loc <- 'USA'
tmp <- which(US_res$date == min(US_res$date))
tmp1 <- which(US_res$adjustment)
R0_final <- rbind(R0_final, US_res[intersect(tmp, tmp1),])


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

Ger_res$loc <- 'Germany'
tmp <- which(Ger_res$date == min(Ger_res$date))
tmp1 <- which(Ger_res$adjustment)
R0_final <- rbind(R0_final, Ger_res[intersect(tmp, tmp1),])

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

Swe_res$loc <- 'Sweden'
tmp <- which(Swe_res$date == min(Swe_res$date))
tmp1 <- which(Swe_res$adjustment)
R0_final <- rbind(R0_final, Swe_res[intersect(tmp, tmp1),])


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

Spain_res$loc <- 'Spain'
tmp <- which(Spain_res$date == min(Spain_res$date))
tmp1 <- which(Spain_res$adjustment)
R0_final <- rbind(R0_final, Spain_res[intersect(tmp, tmp1),])


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

Neth_res$loc <- 'Netherlands'
tmp <- which(Neth_res$date == min(Neth_res$date))
tmp1 <- which(Neth_res$adjustment)
R0_final <- rbind(R0_final, Neth_res[intersect(tmp, tmp1),])


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

Denmark_res$loc <- 'Denmark'
tmp <- which(Denmark_res$date == min(Denmark_res$date))
tmp1 <- which(Denmark_res$adjustment)
R0_final <- rbind(R0_final, Denmark_res[intersect(tmp, tmp1),])


## EVERYTHING:

everything <- ggarrange(italy, France, UK, US, Germany, Sweden, Spain, Netherlands, Denmark,  common.legend = T) 
everything <- ggarrange(everything, as_ggplot(my_legend), heights = c(97,3), nrow = 2) 
annotate_figure(everything, top = title)


final <- ggarrange(italy, UK, Sweden, US, common.legend = T) 
final <- ggarrange(final, as_ggplot(my_legend), heights = c(95,5), nrow = 2)

ggsave(file.path(plots_path, paste0('realdata_M', MEAN,'_SD', SD, '.png')), final, w=10,h=5)
ggsave(file.path(plots_path, paste0('realdata_M', MEAN,'_SD', SD, '.png')), final, w=10,h=5)

R0_final

# ITALY

# data <- extract_by_country('Italy')
# data <- data[3:50,]
# pp <- ggplot(data, aes(x = mdy(date), y = log(cases))) + geom_point() +
#   geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) +
#   geom_vline(mapping = aes(xintercept = mdy('2/21/20')), color= 'green') + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')
# 
# italy_res <- EpEstimComparison(data, eg_window = 7, WINDOW = 7)
# 
# gg <- ggplot(dplyr::filter(italy_res, adjustment == FALSE), aes(x = date, y = M)) +
#   geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
#   ylab('Rt') + xlab('') +
#   theme_bw() +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = 'none')
# 
# my_legend <- get_legend(gg)
# 
# gg <- gg + theme(legend.position = "none")
# 
# italy <- ggarrange(pp, gg,
#                    widths = c(1, 1))
# italy <- annotate_figure(italy, top = 'Italy')




######################################################
###############      Other stuff       ###############
######################################################


### MEASLES ###

require(data.table)
data("Measles1861")

data <- data.frame(date = seq_along(Measles1861$incidence),
                   cases = Measles1861$incidence)

MEAN <- sum(seq_along(Measles1861$si_distr)*Measles1861$si_distr)
SD <- sqrt(sum(seq_along(Measles1861$si_distr)^2*Measles1861$si_distr) - MEAN^2)

res <- EpEstimComparison(data, WINDOW = 7, eg_window = 10)

L <- length(Measles1861$incidence)

res$date <- L + 1 - rev(seq_along(res$date))
res$date[1:(L-7)] <- L + 1 - rev(seq_along(res$date[1:(L-7)]))

pp <-  ggplot(data, aes(x = date, y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')

gg <- ggplot(dplyr::filter(res, !is.na(M)), aes(x = date, y = M, color = adjustment, fill = adjustment)) + 
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  labs(y = 'R_t', x ='') +
  theme_bw() + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))

p <- ggarrange(pp, gg, widths = c(1, 2))
p <- annotate_figure(p, top = 'Measles 1861')
ggsave(filename = file.path(plots_path,'measles1861.png'), plot = p, w=8, h=6)

### SARS ####

data("SARS2003")
data <- data.frame(date = seq_along(SARS2003$incidence),
                   cases = SARS2003$incidence)

data <- dplyr::filter(data, date > 17 & date < 40)
data$date <- seq_along(data$cases)

MEAN <- sum(seq_along(SARS2003$si_distr)*SARS2003$si_distr)
SD <- sqrt(sum(seq_along(SARS2003$si_distr)^2*SARS2003$si_distr) - MEAN^2)

res <- EpEstimComparison(data, WINDOW = 7, eg_window = 7)

L <- dim(data)[1]
res$date <- L + 1 - rev(seq_along(res$date))
res$date[1:(L-7)] <- L + 1 - rev(seq_along(res$date[1:(L-7)]))


pp <- ggplot(data, aes(x = date, y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')


gg <- ggplot(dplyr::filter(res, !is.na(M)), aes(x = date, y = M, color = adjustment, fill = adjustment)) + 
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  labs(y = 'R_t', x ='') +
  theme_bw() + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))

p <- ggarrange(pp, gg, widths = c(1, 2))
p <- annotate_figure(p, top = 'SARS 2003')
ggsave(filename = file.path(plots_path,'sars2003.png'), plot = p, w=8, h=6)



### SMALLPOX ####

data("Smallpox1972")
data <- data.frame(date = seq_along(Smallpox1972$incidence),
                   cases = Smallpox1972$incidence)

data <- dplyr::filter(data, date > 31)
data$date <- seq_along(data$cases)

MEAN <- sum(seq_along(Smallpox1972$si_distr)*Smallpox1972$si_distr)
SD <- sqrt(sum(seq_along(Smallpox1972$si_distr)^2*Smallpox1972$si_distr) - MEAN^2)

res <- EpEstimComparison(data, WINDOW = 7, eg_window = 7)

L <- dim(data)[1]
res$date <- L + 1 - rev(seq_along(res$date))
res$date[1:(L-7)] <- L + 1 - rev(seq_along(res$date[1:(L-7)]))


pp <- ggplot(data, aes(x = date, y = log(cases))) + geom_point() +
  geom_smooth(method='lm', formula= y~x, data = data[1:7,], se=F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('')


gg <- ggplot(dplyr::filter(res, !is.na(M)), aes(x = date, y = M, color = adjustment, fill = adjustment)) + 
  geom_line() + geom_ribbon(aes(ymin = l, ymax = u), alpha = .3) + 
  labs(y = 'R_t', x ='') +
  theme_bw() + 
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))

p <- ggarrange(pp, gg, widths = c(1, 2))
p <- annotate_figure(p, top = 'Smallpox 1972')
ggsave(filename = file.path(plots_path,'smallpox1972.png'), plot = p, w=8, h=6)
