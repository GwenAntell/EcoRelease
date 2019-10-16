# The purpose of this script is to bootstrap 'species-based' samples,
# starting from the output of the bin_occurrences_from_PBDB_records.R script

# set subsampling parameters 
grid_res <- 1 # raster resolution in degrees
n_cells <- 12 # number of grid cells in each iteration subsample
r_km <- 1500 # radius of maximum sampling distance from centroid in km
n_reps <- 500 # number of subsamples taken for each time bin
eck.proj <- "+proj=eck6 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"

# packages crucial to analysis
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(vegan)
# save package names as object to put on all cores later
pkgs <- c('sp','rgdal','raster','rgeos','vegan') 

# packages that speed up computation
# WARNING: this will give R the power to compute on all cores.
# If you'd rather run things in a loop (1 core), 
# change 'foreach X %dopar% Y' to 'for(x in X){ Y(x) }' below
library(foreach)
library(iterators)
library(doParallel)

# read in function that performs subsampling
source('R/function_subsample_species_approach_GSA.R')

# If this script is run immediately after bin_occurrences_from_PBDB_records.R,
# the following line is unnecessary (all relevant objects will be loaded).
# But if jumping in to the analysis from a clean workspace, 
# read in occurrence and time bin data:
source('R/read_data_for_bootstrap_subsampling.R')

# define gridded world
x.raster <- raster(ncol=360/grid_res, nrow=180/grid_res, 
                   xmn=-180, xmx=180, ymn=-90, ymx=90)
proj4string(x.raster) <- CRS("+proj=longlat")
raster_proj <- projectRaster(x.raster, crs = eck.proj)
raster_proj <- setValues(raster_proj, c(1:ncell(raster_proj))) 

# subset pbd occurrence data to species that have known environmental preferences
occ_env <- pbd[, c('pref_lith','pref_bath')]
pref_unk <- apply(occ_env, 1, function(x) any(x == 'uk') ) 
pbd <- pbd[!pref_unk,]

# row indices must perfectly sequential for later steps
row.names(pbd) <- 1:nrow(pbd) 

# load environmental raster stacks
# NB: these should be unzipped and the tif placed in the /Data folder
r_names <- paste0('Data/', c('lithology','bathymetry'), '_from_all_PBDB_marine.tif') 
env_r <- lapply(r_names, brick)

# remove the Cambrian
env_r <- lapply(env_r, dropLayer, 1) 

env_r <- lapply(env_r, setNames, stages$name)

# all cells that have environmental data have occurrences in them
# but some cells have incomplete environmental data, so:
# cross-reference across variables
occs_wo_dat <- vector(mode='list', length=nrow(stages))
for (lyr in 1:nrow(stages)){
  lyr_vals <- sapply(env_r, function(r) getValues(r[[lyr]]) )
  
  cells_w_dat <- apply(lyr_vals, 1, function(x){ ! any(is.na(x)==TRUE) } )
  
  # not all PBDB occurrence records fall in cells with environmental data
  # need to remove such records
  pbd_bin <- pbd[pbd$stage==stages$name[lyr], ]
  occs2toss <- which(! pbd_bin$cell_number %in% which(cells_w_dat==TRUE) )
  occs_wo_dat[[lyr]] <- row.names(pbd_bin[occs2toss,])
  
  lyr_vals[cells_w_dat==FALSE,] <- NA
  env_r <- lapply(1:length(env_r), 
                  function(r){ 
                    this_r <- env_r[[r]]
                    r <- setValues(this_r, lyr_vals[,r], layer=lyr) 
                  }
  )
}
# all non-NA cells now have data for all enviro variables

# take out occurrence records in cells where environmental conditions are unknown
# these cannot be drawn for subsamples
pbd <- pbd[-as.numeric(unlist(occs_wo_dat)),]

# subsample over all time bins
# WARNING: this take 1 hour for 500 iterations (8 core, 3.6 GHz processor)
ncores <- detectCores()
registerDoParallel(ncores) 
res_l <- foreach(bin=stages$name, .packages=pkgs) %dopar% 
  subsample_bin(bin=bin, occ_dat=pbd, env_dat=env_r, bin_col='stage', 
                cell_col='cell_number', coord_cols=c('centroid_long','centroid_lat'), 
                name_col='unique_name', r_km=r_km, n_cells=n_cells, n_reps=n_reps,
                env_pref=c('pref_lith','pref_bath'), prj=CRS(eck.proj) 
  ) 
stopImplicitCluster()

# convert list of output from each stage to datafame
res_df <- data.frame()
for (i in 1:length(res_l)){ 
  res_df <- rbind(res_df, res_l[[i]], stringsAsFactors=FALSE) 
}
colnames(res_df) <- c('sp','bin','n_oc','mst','n_suitable','species_count')

# save output
file_name <- 'Data/DS6_species_data_subsampled_v2.csv'
write.csv(res_df, file_name, row.names = FALSE)
