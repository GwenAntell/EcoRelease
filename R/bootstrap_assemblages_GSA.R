# The purpose of this script is to bootstrap 'assemblage-based' samples,
# starting from the output of the bin_occurrences_from_PBDB_records.R script
# This script contains the options to:
# * specify the subsampling spatial constraints (radius and area)
# * weight cell sampling by proximity
# * remove species that occur once in an assemblage
# * remove the single most common species in an assemblage
# * do classical rarefaction or SQS or neither

# set subsampling parameters
grid_res <- 1 # raster resolution in degrees
ncells <- 12 # number of grid cells in each iteration subsample
r_km <- 1500 # radius of maximum sampling distance from centroid in km
nreps <- 500 # number of subsamples taken for each time bin
coverage <- NULL # specify a value from 0-1; 0.7 is the limit for the supplied dataset
size <- NULL # can set this from 0-50 for the supplied dataset; the limit is in Turonian
equal_probs <- TRUE # set to FALSE to weight cell sampling probability by proximity
remove_singles <- FALSE # set to TRUE to remove species that occur once in an assemblage
remove_dominant <- FALSE # set to TRUE to remove single most common species in each assemblage
eck.proj <- "+proj=eck6 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"

# packages crucial to analysis
library(vegan)
library(iNEXT)
library(sp)
library(raster)
library(rgeos)
# save names of all of the above to put packages on all cores later
pkgs <- c('vegan','iNEXT','sp','raster','rgeos') 

# packages that speed up computation
# WARNING: this will give R the power to compute on all cores.
# If you'd rather run things in a loop (1 core), 
# change 'foreach X %dopar% Y' to 'for(x in X){ Y(x) }' below
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

# read in occurrence and time bin data:
source('R/read_data_for_bootstrap_subsampling.R')

# If using classical rarefaction, omit Tremadocian (Lower Ordovician) - 
# not enough fossils to sample
if (! is.null(size)){ 
  stages <- stages[-1,] 
}

# read in relevant functions
source('R/function_distribution_summary_GSA.R')
source('R/function_find_seed_cells_GSA.R')
source('R/function_subsample_assemblage_approach_GSA.R')

# calculate weights for sampling
# (the same values are good for any circle of same radius and resolution)
if (equal_probs==FALSE){
  # define gridded world
  x.raster <- raster(ncol=360/grid_res, nrow=180/grid_res, xmn=-180, xmx=180, ymn=-90, ymx=90)
  proj4string(x.raster) <- CRS("+proj=longlat")
  raster_proj <- projectRaster(x.raster, crs = eck.proj)
  raster_proj <- setValues(raster_proj, c(1:ncell(raster_proj))) 
  center <- matrix(c(0, 0), 1, 2)
  pt <- SpatialPoints(center, proj4string = CRS("+proj=longlat"))
  pt_eck6 <- spTransform(pt, eck.proj) 
  seed_cell <- extract(raster_proj, pt_eck6) 
  seed_pt <- xyFromCell(raster_proj, seed_cell, spatial = TRUE) 
  
  # units set to m
  seed_buf <- buffer(seed_pt, width=r_km*1000) 
  
  cr <- crop(raster_proj, seed_buf)
  circle <- mask(cr, seed_buf)
  inCircle <- which(is.na(values(circle))==FALSE)
  centroids <- xyFromCell(circle, inCircle, spatial = TRUE)
  gcdists <- spDistsN1(centroids, seed_pt) 
  
  # deal with special case of the center cell
  gcdists[which.min(gcdists)] <- NA
  
  # squared inverse weight because inverse alone is too weak an effect
  weightFun <- function(x){ 
    gcdists[x] ^(-2) 
  }
  probs <- vapply(1:length(gcdists), FUN.VALUE = 1.1, FUN = weightFun)
  cumProb <- sum(probs, na.rm = TRUE) 
  probs <- probs/cumProb
} else { probs <- NULL }

### RUN ON ACTUAL DATA ###

# find seed cells in PBDB dataset
seed_list <- lapply(unique(pbd$stage), find_seed_cells, dat=pbd, bin_col='stage',
                    coord_cols=c('centroid_long', 'centroid_lat'), cell_col='cell_number', 
                    prj=CRS(eck.proj), r=r_km, n=ncells)
names(seed_list) <- unique(pbd$stage)

# set up implicit clusters (run time 3.5x faster than in series)
ncores <- detectCores() - 1
registerDoParallel(ncores) 

# subsample over all time bins
# this step takes 5-10 min without rarefaction, on an 8 core, 3.6 GHz processor
res_l <- foreach(bin=stages$name, .packages=pkgs) %dopar% 
  subsample_bin(bin=bin, dat=pbd, bin_col='stage', seeds=seed_list, nreps=nreps, 
                cell_col='cell_number', coord_cols=c('centroid_long','centroid_lat'), 
                name_col='unique_name', r=r_km, ncells=ncells, prj=CRS(eck.proj), 
                weight=probs, coverage=coverage, size=size, 
                remove_singles=remove_singles, remove_dominant=remove_dominant) 
stopImplicitCluster()
# NB: the subsampling procedure contains randomized steps. Some iterations of a step
# may fail if data are too sparse for a specified rarefaction of species count. 
# The function times out after a certain # of attempts.
# If by bad luck so many attempts fail that the max # of attempts is reached,
# this line can be run again.
# The threshold for timeout can be changed in line 166 of the subsampling function.

# convert list of output from each stage to datafame
n_col <- ncol(res_l[[1]])
res_df <- data.frame(matrix(nrow = nreps*nrow(stages), ncol=n_col))
for (i in 1:length(res_l)){ 
  stgRows <- ((i-1)*nreps+1) : ((i-1)*nreps+nreps)
  res_df[stgRows,] <- res_l[[i]] 
}

range_names <- lapply(c('oc_w.singl','oc_wo.singl','mst'), 
                      function(x) paste(x, c('min','med','max','mean'), sep='_'))
raref_names <- 'species_count'
suffix <- c(sapply(0:2, function(x) paste0('Hill', x)), 'SC')					  
if (!is.null(coverage)){   
  raref_names <- c(raref_names, paste0('SQS_', suffix)) 
}  
if (!is.null(size)){   
  raref_names <- c(raref_names, paste0('CR_', suffix)) 
}

colnames(res_df) <- 
  c(unlist(range_names), raref_names,'pool_n','mst_sample','seed_cell','bin')

# save output (descriptive file name recommended so there's a record of data attributes)
file_name <- 'Data/DS5_assemblage_data_subsampled_v2.csv'
write.csv(res_df, file_name, row.names=FALSE)
