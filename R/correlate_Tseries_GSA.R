# The purpose of this script is to calculate correlation coefficients
# between pairs of time series (e.g. range size vs. species count).
# The initial data are output from the 'bootstrap_assemblages_GSA.R' script,
# and are available as 'DS5_assemblage_data_subsampled.csv', which is read in here.

# for stationarity and cross-correlation tests on time series
library(fUnitRoots) 
library(TSA)

# packages that speed up computation
# WARNING: this will give R the power to compute on all cores.
# If you'd rather run things in a loop (1 core):
# change 'foreach X %dopar% Y' to 'for(x in X){ Y(x) }'
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

dat <- read.csv('Data/DS5_assemblage_data_subsampled.csv', stringsAsFactors=FALSE)
bins <- unique(dat$bin)

# remove subsample iterations where all species were singletons
# this removes only 1 subsample (of 500) from the published data
dat <- dat[!is.na(dat$mst_med),]

# determine the number of successful bootstrap replications
n <- sapply(bins, function(t) length(which(dat$bin==t)))
nruns <- min(n)

resp_vars <- c('oc_wo.singl_mean','mst_med','species_count')
pred_vars <- c('species_count','mst_sample') 

cors <- matrix(nrow = length(resp_vars), ncol = length(pred_vars),
                      dimnames = list(resp_vars, pred_vars))

# set up implicit clusters (run time 3.5x faster than in series)
ncores <- detectCores() - 1
registerDoParallel(ncores) 

# Non-parametric correlation for a single time series at a time
# This function will run in parallel across bootstrap replicates,
# within a loop across time bins
get_cor <- function(){ 
  cols <- c(pred, r_var)
  series_list <- lapply(bins, function(t) dat[which(dat$bin==t)[iter], cols])
  series <- matrix(unlist(series_list), ncol=length(cols), 
                   dimnames=list(bins, cols), byrow=TRUE)
  
  # fit AR1 model to the predictor, non-stationary series
  AR <- arima(series[,pred], order=c(1, 0, 0)) 
  series[,pred] <- as.numeric(AR$residuals) 
  
  # return only the point estimate of the coefficient
  # error bars come from bootstrapping, so no need to save significance values here
  cor(series[,c(r_var,pred)], method='kend')[1,2] 
} 

for (k in 1:length(pred_vars)){
  pred <- pred_vars[k]

  for (i in 1:length(resp_vars)){ 
    r_var <- resp_vars[i]
    
    # calculate across all iterations (in parallel, otherwise this is slowest step)
    cors_i <- foreach(iter=1:nruns, .combine=c) %dopar% get_cor() 
    
    # 95% CI for 2-sided tests
    cors_avg <- mean(cors_i, na.rm = TRUE)
    cors_CI <- quantile(cors_i, probs=c(0.025, 0.975), na.rm=TRUE)
    
    # round and combine into 1 string
    rounded <- sapply(c(cors_avg, cors_CI), round, digits=2) 
    cors[i,k] <- paste(rounded[1], ' (', paste(rounded[2:3], collapse=', '), ')', sep='')
  } # end loop through range metrics
} # loop through each predictor variable

stopImplicitCluster()

# return mean and 95% confidence interval values
cors


################################################
# Test for stationarity in representative series
# (representative series = median over iterations)

# get representative sp count time series
count <- sapply(bins, function(t){ 
  rows <- which(dat$bin==t)
  median(dat[rows, 'species_count'], na.rm=TRUE)
} )

# test for stationarity in original time series
adfTest(count) 

# modify series and test for stationarity again
count_AR <- arima(count, order=c(1, 0, 0)) 
count_e <- as.numeric(count_AR$residuals)
adfTest(count_e) 

# return estimate of 1st order autocorrelation strength
count_AR$coef

# get representative sample aggregation time series
samp_agg <- sapply(bins, function(t){ 
  rows <- which(dat$bin==t)
  median(dat[rows, 'mst_sample'], na.rm=TRUE)
} )
adfTest(samp_agg)
samp_agg_AR <- arima(samp_agg, order=c(1, 0, 0))
adfTest(samp_agg_AR$residuals)
samp_agg_AR$coef

# calculate cross-correlations at a range of time lags

# get representative mst med time series
span <- sapply(bins, function(t){ 
  rows <- which(dat$bin==t)
  median(dat[rows, 'mst_med'], na.rm=TRUE)
} )

# get representative occupancy time series
occ <- sapply(bins, function(t){ 
  rows <- which(dat$bin==t)
  median(dat[rows, 'oc_wo.singl_mean'], na.rm=TRUE)
} )

# cross-correlation plot inspection
prewhiten(count, span, x.model=count_AR) 
prewhiten(count, occ, x.model=count_AR) 
