# Decide whether to run comparison on all species (TRUE)
# or only on those mid-existence (FALSE; main text analysis)
all_spp <- FALSE

# for weighted Wilcoxon test
library(survey) 

##################
#### data prep ###
##################

# read in helper functions
source('R/function_distribution_summary_GSA.R')

# read in regionally subsampled range data
# NB that preferred environment is calculated along only 2 axes (grain size omitted)
dat <- read.csv('Data/DS6_species_data_subsampled.csv',
                stringsAsFactors = FALSE)
dat <- dat[!is.na(dat$n_oc),]

# get all the data associated with each species (e.g. higher taxa)
tax_level <- 'species'
source('R/read_data_for_bootstrap_subsampling.R')

###################################
#### calculate survivor change ####
###################################

# First and last appearance dates (FADs and LADs)
fad_row <- lapply(unique(dat$sp), function(sp){
  sp_rows <- which(pbd$unique_name==sp)
  fad <- pbd$stage[min(sp_rows)]
  which(dat$sp==sp & dat$bin==fad)
})
lad_row <- lapply(unique(dat$sp), function(sp){
  sp_rows <- which(pbd$unique_name==sp)
  lad <- pbd$stage[max(sp_rows)]
  last_delta <- stages$name[which(stages$name==lad)-1] # bin before lad
  which(dat$sp==sp & dat$bin==last_delta)
})
doubles <- intersect(unlist(lad_row), unlist(fad_row)) # only observed over 1 boundary

# OPTIONAL:
# keep species that cross only 1 boundary; leave out FAD and LAD crossers otherwise
ends <- setdiff(union(unlist(fad_row), unlist(lad_row)), doubles) 
if (!all_spp) {
  dat <- dat[-ends,]
}

var_cols <- c('n_oc','mst','n_suitable','species_count')
v <- var_cols[1] # pick the variable to summarise, occupancy or MST
if (v=='mst'){
  keep <- !is.na(dat$mst)
  dat <- dat[keep,]
}

# find boundaries across which more than 1 species survives
any_surv <- function(bin, bins, dat, bin_col, name_col){
  next_bin <- bins[which(bins==bin)+1]
  dat_t1 <- dat[dat[,bin_col]==bin,]
  dat_t2 <- dat[dat[,bin_col]==next_bin,]
  survivors <- intersect(dat_t1[,name_col], dat_t2[,name_col])
  if (length(survivors) > 1) { 
    return(bin) 
  } else { return(vector()) }
}
surv_bounds_l <- sapply(stages$name, any_surv, bins=stages$name, dat=dat, bin_col='bin', name_col='sp')
surv_bounds <- unlist(surv_bounds_l)
# note that Pleistocene is excluded because no modern data are considered

calculate_delta <- function(bin, bins, dat, bin_col, name_col, var_cols){
  next_bin <- bins[which(bins==bin)+1]
  dat_t1 <- dat[dat[,bin_col]==bin,]
  dat_t2 <- dat[dat[,bin_col]==next_bin,]
  survivors <- intersect(dat_t1[,name_col], dat_t2[,name_col])
  
  delta <- lapply(survivors, function(sp) {
    t1row <- which(dat_t1[,name_col]==sp)
    t2row <- which(dat_t2[,name_col]==sp)
    var_delta <- dat_t2[t2row, var_cols] - dat_t1[t1row, var_cols]
    
    output <- data.frame(var_delta) # do this step separately so numeric class is preserved
    output[,c('species','bin1')] <- c(sp, bin)
    return(output)
  })
  
  delta_df <- data.frame(matrix(unlist(delta), ncol= length(var_cols)+2, byrow=TRUE))
  return(delta_df)
}

delta_l <- lapply(surv_bounds, calculate_delta, bins=stages$name, dat=dat, 
                  bin_col='bin', name_col='sp', var_cols=var_cols)
delta_dat <- do.call('rbind', delta_l)

colnames(delta_dat) <- c(var_cols,'sp','bin')
to_numeric <- function(x){ as.numeric(as.character(x))  }
delta_dat[,1:length(var_cols)] <- 
  apply(delta_dat[,1:length(var_cols)], 2, to_numeric)


###################################
### summarise survivor response ###
###################################

# first determine whether each bin sees more or fewer species

# get summary change change of all species in a time bin
sumry_delta_l <- lapply(unique(delta_dat$bin), function(b){
  summarize_distrib(delta_dat[delta_dat$bin==b, 'species_count'], return_quarts=TRUE) 
})
tmp_spp_count <- data.frame(do.call('rbind', sumry_delta_l))
sumry_spp_count <- cbind(tmp_spp_count, unique(delta_dat$bin))
colnames(sumry_spp_count) <- c('Q1','Q2','Q3','avg','sd','bin') 

fall_bins <- which(sumry_spp_count$Q2 <= 0) 
# Rhaetian shows ca. no change

# summarise range metric and compare between groups

# get summary change over all species in a time bin
sumry_delta_l <- lapply(unique(delta_dat$bin), function(b){
  summarize_distrib(delta_dat[delta_dat$bin==b, v], return_quarts=TRUE) 
})
tmp <- data.frame(do.call('rbind', sumry_delta_l))
sumry_delta <- cbind(tmp, unique(delta_dat$bin))
colnames(sumry_delta) <- c('Q1','Q2','Q3','avg','sd','bin') 

# get n-values
sumry_delta$n <- sapply(unique(delta_dat$bin), 
                        function(b) length(which(delta_dat$bin==b)))

sumry_delta$event <- 'Rise'
sumry_delta$event[fall_bins] <- 'Fall'

sumry_delta$event <- as.factor(sumry_delta$event)
table(sumry_delta$event) # inspect sample size

###########################
# The table of range size shifts is presented as Table S4. Feel free to export.

tabl <- sumry_delta[,c('bin','event','avg','sd','n')]
colnames(tabl) <- c('Early bin','Species count','Range shift','Shift SD','Survivor n')

###########################
# plot survivor range size shift

wtd <- sapply(levels(sumry_delta$event), function(type){
				x <- sumry_delta[sumry_delta$event==type,]
				weighted.mean(x$avg, x$n/sum(x$n))
})
line_seg <- data.frame(horiz=1:length(wtd), vert=wtd)

sumry_delta <- sumry_delta[order(sumry_delta$avg, decreasing=TRUE),]

# Figure 3
boxplot(avg ~ event, data=sumry_delta)
segments(x0=c(0.6,1.6), x1=c(1.4,2.4), y0=wtd, col='blue', lwd=2)


###########################
# nonparametric test

# weight by n spp
design <- svydesign(ids = ~0, data = sumry_delta, weights = ~n) 
svyranktest(avg~event, design, test = 'wilcoxon') 

# print the weighted mean of each group
wMean <- function(vals, w){
  w <- w/sum(w)
  wVals <- vals*w
  sum(wVals)
}

r <- sumry_delta$event=='Rise'
riseMean <- wMean(sumry_delta$avg[r], sumry_delta$n[r])
fallMean <- wMean(sumry_delta$avg[!r], sumry_delta$n[!r])
sapply(c(riseMean,fallMean), round, 2)
