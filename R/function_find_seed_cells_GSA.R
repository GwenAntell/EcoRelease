
# function to try all possible starting cells in bin (i.e. all occupied cells)
# save the position of any cells that contain a given pool size within the given buffer

find_seed_cells <- function(dat, bin, bin_col, coord_cols, cell_col, prj, r, n){
  sl.oc <- dat[dat[,bin_col]==bin, c(coord_cols, cell_col)]
  possible_seeds <- unique(sl.oc[, cell_col])
  
  # subset sl.oc to unique cells for faster (?) processing
  rows2keep <- sapply(possible_seeds, function(cell) which(sl.oc[,cell_col]==cell)[1] )
  sl.oc <- sl.oc[rows2keep,]
  
  # convert coordinates of cell centroids to spatial object
  occ_centroids <- SpatialPointsDataFrame(sl.oc[, coord_cols], sl.oc, proj4string=prj)
  
  # test whether each occupied cell could be used for subsampling
  test_cell <- function(seed_cell){
    seed_pt <- occ_centroids[which(sl.oc[,cell_col]==seed_cell),]
    seed_buf <- buffer(seed_pt, width=r*1000) # radius in meters
    pool <- occ_centroids[seed_buf,]@data[,cell_col] # cells within radius of seed cell
    ifelse(length(pool)>=n, return(seed_cell), return(NA))
  }
  
  # return vector of cells that lie within radius of at least n cells
  viable <- sapply(possible_seeds, test_cell)
  viable <- viable[!is.na(viable)]
  return(viable) 
}

