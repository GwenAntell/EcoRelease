
# FUNCTION OF 1 SUBSAMPLE (1 try on 1 starting cell within 1 bin):

subsample <- function(dat, seed_cell, occ_centroids, cell_col, coord_cols, name_col, 
                      r, ncells, prj, weight=NULL, coverage=NULL, size=NULL, 
                      remove_singles=FALSE, remove_dominant=FALSE){
  
  # take a subsample of cells
  seed_pt <- occ_centroids[which(dat[,cell_col]==seed_cell)[1],]
  seed_buf <- buffer(seed_pt, width=r*1000) # radius in meters
  
  if (is.null(weight)){  
    pool <- unique(occ_centroids[seed_buf,]@data[,cell_col]) # cells within radius of seed cell
    sampled_cells <- sample(sample(pool), ncells, replace=FALSE) # random, unweighted subsampling
  } else {
    cr <- crop(raster_proj, seed_buf)
    circle <- mask(cr, seed_buf) # constrain sampling area
    pool <- sort(unique(dat$cell_number[dat$cell_number %in% values(circle)])) # list of occupied cells
    pool <- pool[!pool==seed_cell] # seed cell will be included in sample anyway
    x <- which(values(circle) %in% pool==TRUE) # get raster positions of occupied cells
    
    # check circle isn't on world edge
    if (length(values(circle)[!is.na(values(circle))]) == length(weight)) { 
      values(circle)[!is.na(values(circle))] <- weight
      sampled_cells <- c(seed_cell, sample(pool, ncells-1, prob = values(circle)[x]))
    } else { stop('fell off edge of the Earth') }
  }
  
  pres_in_sample <- dat[dat[,cell_col] %in% sampled_cells,]
  
  if (! is.null(size)){  if (nrow(pres_in_sample) < size) { stop('no good') } }
  if (remove_singles==TRUE){
    name_freq <- table(pres_in_sample[,name_col])
    singles <- names(name_freq[name_freq==1])
    pres_in_sample <- pres_in_sample[! pres_in_sample[,name_col] %in% singles,]
  }
  if (remove_dominant==TRUE){
    name_freq <- table(pres_in_sample[,name_col])
    dominatrix <- names(which.max(name_freq))
    pres_in_sample <- pres_in_sample[pres_in_sample[,name_col] != dominatrix,]
  }
  spp <- unique(pres_in_sample[,name_col])
  nspp <- length(spp)
  
# RANGE SIZE METRICS
  
  # occupancy for each species
  oc_spp <- as.numeric(table(pres_in_sample[,name_col]))
  
  # minimum spanning tree length for each species, in km, based on great circle distances
  mst_spp <- sapply(spp, function(x){
    rows <- which(pres_in_sample[,name_col]==x)
    if (length(rows) < 2) { return(NA) } else {
      pts_sp <- SpatialPoints(pres_in_sample[rows, coord_cols], proj4string = prj)
      gcdists <- spDists(pts_sp, longlat = FALSE) 
      mst_sp <- spantree(gcdists)
      mst_length_sp <- sum(mst_sp$dist)/1000 
      return(sqrt(mst_length_sp)) 
    }
  } # end function to calculate species site dispersion
  )
  
  # summarise range metrics as min, med, max, avg
  oc_w.singles <- t(summarize_distrib(oc_spp, include_singles=TRUE, return_quarts=FALSE))
  oc_wo.singles <- t(summarize_distrib(oc_spp, include_singles=FALSE, return_quarts=FALSE))
  mst <- t(summarize_distrib(mst_spp, na.rm=TRUE, return_quarts=FALSE))
 
  # combine into output
  metrics <- c('w.singl','wo.singl','mst')
  nm_paste <- function(prefix, vect){
    paste(prefix, vect, sep='_')
  }
  range_names <- lapply(metrics, nm_paste, vect = dimnames(oc_w.singles)[[2]])
  range_dat <- cbind(oc_w.singles, oc_wo.singles, mst) 
  dimnames(range_dat)[[2]] <- unlist(range_names) 
  
# DIVERSITY METRICS
  
  richness_dat <- matrix(nspp)
  # note: there will be other species not counted if singletons or dominant species were removed
  
  if (! is.null(coverage)){
   # suppress warning about extrapolation since any coverage < required will be discarded later
  SQS_full <- suppressWarnings(estimateD(list(c(ncells, oc_spp)), datatype='incidence_freq',
                                         base='coverage', level=coverage, conf=NULL))
										 
  # q0=richness, q1=exp of Shannon's entropy index, q2=inverse of Simpson's concentration index
  SQS <- as.matrix(SQS_full[c('q = 0', 'q = 1', 'q = 2', 'SC')])
  if (SQS[,'SC'] < coverage) { stop('go get more rock') }
  richness_dat <- cbind(richness_dat, SQS)
  }
  
  if (! is.null(size)){
  CR_full <- suppressWarnings(estimateD(list(c(ncells, oc_spp)), datatype='incidence_freq',
                                        base='size', level=size, conf=NULL))
  CR <- as.matrix(CR_full[c('q = 0', 'q = 1', 'q = 2', 'SC')])
  richness_dat <- cbind(richness_dat, CR)
  }
 
# SAMPLE METRICS
  
  pool_n <- length(pool)
  
  # minimum spanning tree length of sample cells, in km, based on great circle distances
  # subset to unique cells for faster processing 
  rows2keep <- sapply(sampled_cells, function(cell) which(dat[,cell_col]==cell)[1] )
  dispersion_dat <- dat[rows2keep,] 
  pts_sample <- SpatialPoints(dispersion_dat[, coord_cols], proj4string = prj)
  gcdists_sample <- spDists(pts_sample, longlat = FALSE)
  spant_sample <- spantree(gcdists_sample)
  mst_sample <- sqrt(sum(spant_sample$dist)/1000) 
  
  output_sample <- cbind(range_dat, richness_dat, pool_n, mst_sample, seed_cell)
  return(output_sample)
} # end subsample function

# FUNCTION OF ONE TIME BIN

subsample_bin <- function(dat, bin, bin_col, seeds, nreps, cell_col, coord_cols, name_col, 
                          r, ncells, prj, weight=NULL, coverage=NULL, size=NULL, 
                          remove_singles=FALSE, remove_dominant=FALSE){
  
  dat <- dat[dat[,bin_col]==bin, c(coord_cols, cell_col, bin_col, name_col)]
  seeds <- seeds[[bin]]
  occ_centroids <- SpatialPointsDataFrame(dat[, coord_cols], dat, proj4string = prj)
  
  # subsample cells within radius of each center/seed cell as many times as nreps
  
  # keep count of how many subsamples achieve sufficient coverage
  success <- 0 
  # keep track of subsampling attempted so can break out of while loop
  attempts <- 0 
  
  n_col <- 16
  if(! is.null(coverage)){ n_col <- n_col + 4 }
  if(! is.null(size)){ n_col <- n_col + 4 }
  output <- matrix(nrow=nreps, ncol=n_col) 
  
  while(success < nreps){
    attempts <- attempts + 1
	
	# select seed cell at random
    if (length(seeds)>1){
      seed_cell <- sample(sample(seeds), 1) 
	# NB: 'sample' fcn makes mistake if only 1 item to pick
    } else { seed_cell <- seeds } 
    
    # try subsampling; return a matrix with as many rows as species or an empty vector
    sampled_dat <- tryCatch({
      subsample(dat=dat, occ_centroids=occ_centroids, seed_cell=seed_cell,
                cell_col=cell_col, coord_cols=coord_cols, name_col=name_col, 
                ncells=ncells, r=r, prj=prj, weight=weight, coverage=coverage, size=size, 
                remove_singles=remove_singles, remove_dominant=remove_dominant)
      }, error = function(err){ 
	    return(vector()) 
	  }
    )
    
	# case where subsample met coverage threshold
    if (length(sampled_dat) > 0){ 
      success <- success + 1
      output[success,] <- sampled_dat
    }
    
	# safety catch
    if (attempts > nreps*200) { stop(paste(bin,'maximum attempts reached')) } 
  } # end while loop to achieve n samples
  
  output <- data.frame(output)
  colnames(output) <- dimnames(sampled_dat)[[2]]
  output$bin <- bin
  return(output)
}
