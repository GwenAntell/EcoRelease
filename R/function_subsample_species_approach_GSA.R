
# FUNCTION OF 1 TIME BIN

subsample_bin <- function(occ_dat, env_dat, bin, bin_col, cell_col, coord_cols, name_col, 
                          r_km, prj, env_pref, n_cells, n_reps){
  
  # FUNCTION OF 1 SEED CELL ITERATION
  # this must be nested so that subsample_bin is the parent environment
  iteration <- function(seed_cell){
    
    # put a circle around the focal cell (i.e. seed)
    seed_pt <- occ_centroids[ which(occ_centroids@data[,cell_col]==seed_cell)[1], ]
    seed_buf <- buffer(seed_pt, width=r_km*1000) # radius in meters
    
    occs_in_region <- unique(occ_centroids[seed_buf,]@data[,cell_col]) 
    suitable_in_region <- intersect(occs_in_region, suitable)
	
	# sample the seed cell by default, to guarantee at least 1 occurrence of the focal species
	suitable_in_region <- suitable_in_region[-which(suitable_in_region==seed_cell)] 
    n_suitable <- length(suitable_in_region) + 1
    if (n_suitable < n_cells){ stop('nowhere to live') } 
	
	# random, unweighted subsampling
	# include seed cell in subsample so >=1 occurrence is guaranteed
    sampled_cells <- c(seed_cell, sample(sample(suitable_in_region), n_cells-1, replace=FALSE) )
    sampled_sp_rows <- sp_rows[ which(occ[sp_rows, cell_col] %in% sampled_cells) ]
	
	sampled_spp <- occ[ occ[,cell_col] %in% sampled_cells, name_col]
	n_spp <- length(unique(sampled_spp))
	
    n_oc <- length(occ[sampled_sp_rows, cell_col]) 
    
    # minimum spanning tree length, in km, from great circle distance matrix
    if (n_oc < 2) { mst <- NA } else {
      pts_sp <- SpatialPoints(occ[sampled_sp_rows, coord_cols], proj4string = prj)
      gcdists <- spDists(pts_sp, longlat = FALSE) 
      mst_sp <- spantree(gcdists)
      mst_length_sp <- sum(mst_sp$dist)/1000 
      mst <- sqrt(mst_length_sp) 
    }
    
    return( c('n_oc'=n_oc, 'mst'=mst, 'n_suitable'=n_suitable, 'n_spp'=n_spp) )
  } 
  
  occ <- occ_dat[occ_dat[,bin_col]==bin, c(coord_cols, cell_col, bin_col, name_col, env_pref)]
  env <- lapply(env_dat, function(r) r[[bin]] )
  
  occ_centroids <- SpatialPointsDataFrame(occ[, coord_cols], occ, proj4string = prj)
  spp <- unique(occ[, name_col])
  
  cells_w_dat <- intersect(unique(occ[,cell_col]), which(is.na(getValues(env[[1]]))==FALSE))
  n_spp <- length(spp)
  n_output_col <- 4
  output_spp <- matrix(nrow=n_spp, ncol=n_output_col)
  
  for (i in 1:n_spp){
    sp_rows <- which(occ[,name_col]==spp[i])
    
    #### find all cells in the world that are suitable for focal species ###
    
    # start with assumption that all occupied cells are suitable
    suitable <- cells_w_dat
    
    for (j in 1:length(env_pref)){
      pref <- occ[sp_rows[1], env_pref[j]]
      env_j <- env[[j]]
      rat <- env_j@data@attributes[[1]][-1,] # raster value ID linked to each pref
      
      # get raster IDs that represent enviro suitable for the species
      # raster cells with enviro conditions of 'both' are fine for any species
      if (pref=='both'){ pref_ID <- rat$ID } else { pref_ID <- rat$ID[ rat$category %in% c(pref,'both')] }
      
      env_values <- getValues(env_j)[cells_w_dat]
      suitable <- intersect(suitable, cells_w_dat[env_values %in% pref_ID])
    }
    
    if (length(suitable) < n_cells){ output_spp[i,] <- rep(NA, n_output_col); next }
    
    #### iterate across all regions that contain the species ###
    
    occs_in_suitable <- sp_rows[ which(occ[sp_rows, cell_col] %in% suitable) ]
    seeds_sp <- occ[occs_in_suitable, cell_col]
    if (length(seeds_sp)==0) { output_spp[i,] <- rep(NA, n_output_col); next } 
    
    success <- 0 # keep count of how many subsamples achieve sufficient habitat coverage
    output_sp <- matrix(nrow=n_reps, ncol=n_output_col)

    while(success < n_reps){
      
      # select seed cell at random; NB: 'sample' fcn makes mistake if only 1 item to pick
      if (length(seeds_sp)>1){ seed_cell <- sample(sample(seeds_sp), 1) } else { seed_cell <- seeds_sp }
      
      # try subsampling; return a matrix with as many rows as species or an empty vector
      sampled_dat <- tryCatch({ iteration(seed_cell) }, error = function(err) return(vector()) )
      
      # case where subsample attempt was successful
      if (length(sampled_dat) > 0){ 
        success <- success + 1
        output_sp[success,] <- sampled_dat
        
        # insufficient habitat in region around given seed
      } else { seeds_sp <- seeds_sp[-which(seeds_sp==seed_cell)] } 
      
      if (length(seeds_sp)==0) { break }
    } # end while loop to achieve n samples	  	  
    
    #### summarise across iterations ###
    
    # output matrix of while loop could be incompletely filled if habitat is sparse
    output_sp <- output_sp[!is.na(output_sp[,1]),]
    
    if (nrow(output_sp)==0) { output_spp[i,] <- rep(NA, n_output_col); next } 
    if (nrow(output_sp)==1) { output_spp[i,] <- sp_i; next } 
    
    # save one iteration where the prop occ is median among iterations
    
    # force median to pick a realised value (otherwise would return mean):
    # check if even or odd number of rows; take off a row if needed so there is an odd number
    if (nrow(output_sp) > 1) {  
      test_val <- nrow(output_sp)/2 
      if ( test_val==round(test_val) ){ output_sp <- output_sp[-which.min(output_sp[,1]),] }
    }
    
    med_row <- which(output_sp[,1]==median(output_sp[,1])) 
    if (length(med_row)==1){ representative <- med_row } else { representative <- sample(med_row, 1) }
    output_spp[i,] <- output_sp[representative,] 
  } # end loop through species
  
  output_spp <- cbind(spp, rep(bin, n_spp), output_spp)
  return(output_spp) # species does not occur in any samples if this is NA
}
