# The purpose of this script is to build linear mixed-effects models
# on species that cross stage boundaries (i.e. survivors).
# The initial data are output from the 'bootstrap_across_species_GSA.R' script,
# and are available as 'DS6_species_data_subsampled.csv', which is read in here.

library(lme4)
library(arm)

# function to return slected distribution moments
source('R/function_distribution_summary_GSA.R')

# read in regionally subsampled range data
dat <- read.csv('Data/DS6_species_data_subsampled.csv', stringsAsFactors=FALSE)
dat <- dat[!is.na(dat$n_oc),]

# read in all data associated with a species
source('R/read_data_for_bootstrap_subsampling.R')

var_cols <- c('n_oc', 'mst','n_suitable','species_count') 
n_col <- length(var_cols)*2

#############################################################
# calculate variables

#### survivor change in occupancy and community species count ###

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
surv_bounds_l <- sapply(stages$name, any_surv, bins=stages$name, 
                        dat=dat, bin_col='bin', name_col='sp')
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
    standing_t2 <- (dat_t2[t2row, var_cols])
    var_delta <- standing_t2 - dat_t1[t1row, var_cols]
    
    # combine data separately so numeric class is preserved
    output <- data.frame(cbind(var_delta, standing_t2)) 
    output[,c('species','bin')] <- c(sp, bin)
    return(output)
  })
  
  delta_df <- data.frame(matrix(unlist(delta), ncol = n_col+2, byrow=TRUE))
  return(delta_df)
}

# calculate delta for species in each time bin
delta_l <- lapply(surv_bounds, calculate_delta, bins=stages$name, dat=dat, 
                  bin_col='bin', name_col='sp', var_cols=var_cols)
delta_dat <- do.call('rbind', delta_l)

# structure data appropriately (otherwise all variables are factors)
colnames(delta_dat) <- c(var_cols, paste0('t2',var_cols), 'sp','bin') 
delta_dat[,1:n_col] <- apply(delta_dat[,1:n_col], 2, function(x) as.numeric(as.character(x)) )
delta_dat$sp <- as.character(delta_dat$sp)

#### add other variables associated with each species ####
new_vars <- c('phylum','order','family','pref_lith','pref_bath','time_mid')
get_vars <- function(sp, dat){
  row_pos <- which(pbd$unique_name==sp)[1]
  pbd[row_pos,new_vars]
}
extra_dat <- lapply(delta_dat$sp, get_vars, dat=delta_dat)
dat2add <- do.call('rbind', extra_dat)
delta_dat <- cbind(delta_dat, dat2add)

#### determine where each record is within a species' duration ###

# First and last appearance dates (FADs and LADs)
fad_row <- lapply(unique(delta_dat$sp), function(sp, dat){
  sp_rows <- which(pbd$unique_name==sp)
  fad <- pbd$stage[min(sp_rows)]
  which(dat$sp==sp & dat$bin==fad)
}, dat=delta_dat)
lad_row <- lapply(unique(delta_dat$sp), function(sp, dat){
  sp_rows <- which(pbd$unique_name==sp)
  lad <- pbd$stage[max(sp_rows)]
  last_delta <- stages$name[which(stages$name==lad)-1] # bin before lad
  which(dat$sp==sp & dat$bin==last_delta)
}, dat=delta_dat)

# only observed over 1 boundary
doubles <- intersect(unlist(lad_row), unlist(fad_row)) 

delta_dat$dur_pt <- 'm'
delta_dat$dur_pt[unlist(fad_row)] <- 'f'
delta_dat$dur_pt[unlist(lad_row)] <- 'l'

# doubleton species could be a separate class, but they behave like species mid existence
delta_dat$dur_pt[doubles] <- 'm' 

# convert and order factor for easier interpretation
delta_dat$dur_pt <- relevel(as.factor(delta_dat$dur_pt), 'm')

resp <- 'n_oc' # 'mst'

# remove singletons if MST is response
if (resp=='mst'){
  keep <- !is.na(delta_dat$mst)
  delta_dat <- delta_dat[keep,]
}

##############################################
# build models

# baseline (simplest) model
old_fmla <- as.formula(paste(resp, '~ 1 + (1 | bin)')) 

# write out possible fixed effects to add
fixd <- c('species_count','dur_pt','phylum','pref_lith','pref_bath','time_mid') 
# note - we include stage as a random effect to account for the unavoidable part 
# of study design: species living at the same time tend to be similar to one another

# streamline model-building code
legit_lmer <- function(basic, fancy, dataset){
  lme_new <- lmer(update(basic, fancy), data=dataset)
}

# custom function to save model fit information 
mod_sumry <- function(mod){
  y <- paste(formula(mod)[2])
  preds <- paste(formula(mod)[-(1:2)], collapse = ' ')
  mod_aic <- extractAIC(mod) 
  c(y, preds, round(mod_aic,2))
}

# count the number of ways to pick terms (i.e. number of models fit)
combos <- function(n, r){
  factorial(n) / (factorial(r) * factorial(n-r))
} 
nmods <- sum(sapply(0:length(fixd), combos, n=length(fixd))) 

# empty matrix to store model AICs and formulae
mod_fits <- matrix(NA, nrow=nmods, ncol=4) 

# empty list to store fitted models
mods <- vector(mode='list', length=nmods) 

# build models with combinations of 1 or 5 fixed effects
for (i in 1:length(fixd)){
  term <- fixd[i]
  
  # add each single fixed effect by itself
  new_fmla <- as.formula(paste('~ . +', term))
  mods[[i]] <- legit_lmer(basic=old_fmla, fancy=new_fmla, dataset=delta_dat)
  mod_fits[i,] <- mod_sumry(mods[[i]])
  i <- i + length(fixd) # i is the row index in the model summary matrix
  
  # combination of (n-1) fixed effects
  new_fmla <- paste('~ . +', paste(fixd[fixd!=term], collapse = '+') ) 
  mods[[i]] <- legit_lmer(basic=old_fmla, fancy=new_fmla, dataset=delta_dat)
  mod_fits[i,] <- mod_sumry(mods[[i]])
}

# combinations of 2 or 4 fixed effects
pred_pairs <- combn(fixd, 2)
n_pair <- ncol(pred_pairs)
for (i in 1:n_pair){
  new_fmla <- as.formula(paste('~ . +', paste(pred_pairs[,i], collapse = '+')))
  # j is the row index in the model summary matrix
  j <- i + length(fixd)*2 
  mods[[j]] <- legit_lmer(basic=old_fmla, fancy=new_fmla, dataset=delta_dat)
  mod_fits[j,] <- mod_sumry(mods[[j]])
  j <- j + n_pair
  
  # (n-2) fixed effects
  n2less <- fixd[! fixd %in% pred_pairs[,i]]
  new_fmla <- paste('~ . +', paste(n2less, collapse = '+') ) 
  mods[[j]] <- legit_lmer(basic=old_fmla, fancy=new_fmla, dataset=delta_dat)
  mod_fits[j,] <- mod_sumry(mods[[j]])
}

# combinations of 3 fixed effects
triads <- combn(fixd, 3)
for (i in 1:ncol(triads)){
  new_fmla <- as.formula(paste('~ . +', paste(triads[,i], collapse = '+')))
  # j is the row index in the model summary matrix
  j <- i + length(fixd)*2 + n_pair*2
  mods[[j]] <- legit_lmer(basic=old_fmla, fancy=new_fmla, dataset=delta_dat)
  mod_fits[j,] <- mod_sumry(mods[[j]])
}

# simplest model: no fixed effects
mods[[nmods-1]] <- legit_lmer(basic=old_fmla, fancy=old_fmla, dataset=delta_dat)
mod_fits[nmods-1,] <- mod_sumry(mods[[nmods-1]])

# most complex model: all possible effects
new_fmla <- paste('~ . +', paste(fixd, collapse = '+') ) 
mods[[nmods]] <- legit_lmer(basic=old_fmla, fancy=new_fmla, dataset=delta_dat)
mod_fits[nmods,] <- mod_sumry(mods[[nmods]])

##############################################
# calculate model weights and inspect coefficient details

mod_fits <- as.data.frame(mod_fits, stringsAsFactors = FALSE)
colnames(mod_fits) <- c('response','predictors','n_params','AIC')
mod_fits$AIC <- as.numeric(mod_fits$AIC)

mod_fits <- mod_fits[order(mod_fits$AIC, decreasing = FALSE),] 
mod_fits$deltaAIC <- mod_fits$AIC - mod_fits$AIC[1] 

# calculate AIC weights
mod_fits$weight <- sapply(-0.5*mod_fits$deltaAIC, exp)
mod_fits$weight <- mod_fits$weight / sum(mod_fits$weight)

# tables S5 is a formatted versions of the mod_fits object

# 5 best models
mod_indices <- as.numeric(rownames(mod_fits[1:5,]))

# cumulative weight
sum(mod_fits$weight[1:5]) 

# save confidence intervals for coefficients on a log-odds scale
get_odds <- function(ranking){
  refit <- update(mods[[ranking]], data=delta_dat, verbose=FALSE)
  cc <- confint(refit, parm='beta_', method='profile') 
  ltab <- cbind(est=fixef(mods[[ranking]]), cc) 
  tab_round <- round(ltab, 3)
  mod_odds <- cbind(Predictor=rownames(cc), tab_round)
}
odds_l <- lapply(mod_indices, get_odds)
odds_df <- data.frame( do.call('rbind', odds_l) )

# add a column to specify model rank
cur_cols <- colnames(odds_df)
odds_df$Rank <- ''
at_intrcpt <- grep('Intercept', odds_df$Predictor)
odds_df$Rank[at_intrcpt] <- 1:5

# tables S5 is a formatted versions of the mod_fits object
mod_fits

# table S5 is a formatted version of the odds_df object
odds_df
