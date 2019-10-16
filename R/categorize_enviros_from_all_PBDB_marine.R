# This script takes ALL non-terrestrial PBDB occurrence data and uses them to characterise grid cells by environment.
# All taxa, at all identification levels, are included - but fossils must be 'regular' not e.g. trace.
# The outputs of this script are 2 raster files (tif format), which are called in 'bin_occurrences_from_PBDB_records.R'
# The authors already provide these raster files, so it is unneccessary to run this script when replicating analysis.
# This script is provided anyway for full transparency about how the tif files were generated.

library(sp)
library(raster)
library(rgeos)

library(foreach)
library(iterators)
library(doParallel)

# for ease of downloading time scale
library(velociraptr) 

pkgs <- c('sp','raster','rgeos') # save to put on all cores

# This is a large download. Consider increasing timeout settings or downloading from the GUI instead of API.
#pbdUrl <- paste0(
#  'http://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Brachiopoda,',
#  '%20Bivalvia&interval=Cambrian,Holocene&time_rule=overlap&envtype=!terr',
#  '&show=full,classext,ident,img,etbasis,strat,lith,env,timebins,timecompare,resgroup,crmod')
#all_occs <- read.csv(pbdUrl, header=TRUE, skip=19, stringsAsFactors=FALSE)

# Alternatively, use the same play data as in bin_occurrences_from_PBDB_records.R:
pbdUrl <- paste0('https://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Pectinida',
                 '&taxon_reso=genus&idmod=!af,cf,sl,if,eg,qm,qu&interval=Cambrian,Holocene&time_rule=overlap',
                 '&envtype=!terr&show=attr,class,ident,img,abund,coll,coords,loc,paleoloc,strat,stratext,lith,lithext,',
                 'geo,timebins,resgroup,ref,refattr,ent')
all_occs <- read.csv(pbdUrl, header=TRUE, skip=21, stringsAsFactors=FALSE)

# the time-binning and basic cleaning follows that of bin_occurrences_from_PBDB_records.R

# some records are missing coordinates
all_occs <- subset(all_occs, !is.na(all_occs$paleolat) & !is.na(all_occs$paleolng)) 
# coordinates too imprecise:
all_occs <- all_occs[-which(all_occs$geogscale=='basin' & all_occs$latlng_precision=='minutes'),] 
# reduce object size
all_occs <- all_occs[, c('max_ma','min_ma','paleolng','paleolat','lithology1','environment')] 

################### assign stage affiliation ##################

# The ICS chronostratigraphic chart is available at www.stratigraphy.org
# The same dates can be downloaded through the Macrostrat API, directly or via {velociraptr}
stages <- downloadTime('international ages')
stages$name <- row.names(stages)
stages <- stages[order(stages$b_age, decreasing=TRUE), ]
stages$b_round <- stages$t_round <- 0

stages2omit <- c('Stage 2','Stage 3','Stage 4','Wuliuan',
                 'Drumian','Guzhangian','Paibian','Jiangshanian',
                 'Stage 10',
                 'Floian','Darriwilian',
                 'Katian', # otherwise no seed cells for Sandbian
                 'Aeronian', # otherwise no Rhuddanian or Aeronian seed cells
                 'Homerian','Ludfordian',
                 'Pragian', 
                 'Eifelian',
                 'Bashkirian',
                 'Kasimovian', 
                 'Sakmarian','Kungurian', # no Artinskian species records or Sakmarian/Artinskian grain size
                 'Olenekian', # otherwise only 1 abundance datum for Olenekian
                 'Sinemurian', # Hettangian is too poorly sampled
                 'Bajocian', # otherwise no Aalenian seed cells
                 'Hauterivian','Barremian', # no seed cells for Haut., Barremian or Valanginian alone
                 'Santonian', # otherwise nothing survives from Coniacian
                 'Thanetian',
                 'Bartonian', # otherwise no environmental data for Bartonian
                 'Aquitanian', # otherwise no seeds here or in Chattian
                 'Serravallian', # otherwise no seed cells for Langhian
                 'Messinian', # otherwise no seed cells for Messinian
                 'Calabrian','Middle Pleistocene','Late Pleistocene', # otherwise weird extinction rates
				         'Northgrippian','Meghalayan') # lump all Holocene records so they're easy to remove later 
stages_trunc <- stages[!(stages$name %in% stages2omit),] #remove lumped stages

# round boundary ages, with more digits for more recent boundaries

round_age <- function(vect, digits, round_up=TRUE){
  vect <- (10^digits)*vect
  ifelse(round_up==TRUE, vect <- ceiling(vect), vect <- floor(vect))
  return( vect / (10^digits) )
}

u10 <- which(stages_trunc$b_age < 10)
u150 <- which(stages_trunc$b_age < 150 & stages_trunc$b_age > 10)
old <- which(stages_trunc$b_age > 150)

for (group in 1:3){
  bins <- list(u10, u150, old)[[group]]
  digits <- c(2, 1, 0)[group]
  
  # round down younger boundary (terminus) and up older boundary (beginning) per stage
  for (i in bins){
    b <- stages_trunc$b_age[i]
    t <- stages_trunc$b_age[i+1]
    stages_trunc$b_round[i] <- round_age(b, digits=digits, round_up=TRUE)
    stages_trunc$t_round[i] <- round_age(t, digits=digits, round_up=FALSE)
  }
}

binned <- data.frame()

# loosely based on velociraptr constrainAges fcn
for (i in 1:(nrow(stages_trunc)-1)) { # Assign stages_trunc by myr
  stage <- stages_trunc$name[i]
  next_stage <- stages_trunc$name[i+1] # name of more recent stage
  # include names of all bins that are lumped with the focal time bin
  included <- stages$name[which(stages$name==stage):(which(stages$name==next_stage)-1)]
  
  # define the buffered time window for the focal bin
  b <- stages_trunc$b_round[i]
  t <- stages_trunc$t_round[i]
  
  earlypos <- union(which(all_occs$max_ma<=b & all_occs$max_ma>=t), which(all_occs$early_interval %in% included)) 
  latepos <- union(which(all_occs$min_ma<=b & all_occs$min_ma>=t), which(all_occs$late_interval %in% included))
  keepers <- intersect(earlypos, latepos)
  
  if (length(keepers)>0){
    cleanbin <- cbind(all_occs[keepers,], stage)
    binned <- rbind(binned, cleanbin)
  }
}

################### find grid occurrences ##################

# define gridded world
grid_res <- 1 
x.raster <- raster(ncol=360/grid_res, nrow=180/grid_res, xmn=-180, xmx=180, ymn=-90, ymx=90)
proj4string(x.raster) <- CRS("+proj=longlat")
eck.proj <- "+proj=eck6 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
raster_proj <- projectRaster(x.raster, crs = eck.proj)
raster_proj <- setValues(raster_proj, c(1:ncell(raster_proj))) 

# get grid cells from occurrence points
pts <- SpatialPointsDataFrame(binned[c('paleolng','paleolat')], data = binned, 
                              coords.nrs = which(colnames(binned) %in% c('paleolng','paleolat')==T), 
                              proj4string = CRS("+proj=longlat"))
pts_eck6 <- spTransform(pts, eck.proj)
binned$cell_number <- extract(raster_proj, pts_eck6) 
binned[c('centroid_long', 'centroid_lat')] <- xyFromCell(raster_proj, binned$cell_number)

############## derive environment from each record's lithology #############

# do not make inference for mixed lithologies
# otherwise, mixed lithology will be overwritten as siliciclastic
binned$lithology1[grep('mixed', binned$lithology1)] <- NA 

# The following section is modified from Nurnberg and Aberhan (2013)

# Define variables for search 
carb <- c("\"carbonate\"", "\"limestone\"", "\"reef rocks\"", "bafflestone", "bindstone", "dolomite", 
          "framestone", "grainstone", "lime mudstone", "packstone", "rudstone", "floatstone",
          "wackestone")
clast <- c("\"shale\"", "\"siliciclastic\"", "\"volcaniclastic\"", "claystone", "conglomerate",  
           "mudstone", "phyllite", "quartzite", "sandstone", "siltstone", "slate", "schist")
shallow <- c("coastal indet.", "delta front", "delta plain",
             "deltaic indet.", "estuary/bay", "foreshore", "interdistributary bay",
             "lagoonal", "lagoonal/restricted shallow subtidal",
             "marginal marine indet.", "open shallow subtidal", "fluvial-deltaic indet.",
             "paralic indet.", "peritidal", "prodelta", "sand shoal",
             "shallow subtidal indet.", "shoreface", "transition zone/lower shoreface",
             "intrashelf/intraplatform reef", "reef, buildup or bioherm",
             "perireef or subreef", "platform/shelf-margin reef") # last 2 lines include reefs
deep <- c("basinal (carbonate)", "basinal (siliceous)", "basinal (siliciclastic)",
          "deep-water indet.", "deep subtidal indet.", "deep subtidal ramp",
          "deep subtidal shelf", "offshore", "offshore indet.",
          "offshore shelf", "slope", "submarine fan", "offshore ramp", 
          "basin reef", "slope/ramp reef") # last line includes deep water reefs

# categorise each record for each enviro type 
l <- b <- rep(NA, nrow(binned))

# Lithology
l[binned$lithology1 %in% carb] <- "c"  
l[binned$lithology1 %in% clast] <- "s"

# Bathymetry
b[binned$environment %in% shallow] <- "p" # p = proximal = shallow
b[binned$environment %in% deep] <- "d"

now_names <- c('lithnow','bathnow')
binned[,now_names] <- cbind(l,b)

############## rasterize environmental conditions through time #############

find_dominant <- function(dat, bin, env){
  slc <- dat[which(dat$stage==bin),]
  slc <- slc[!is.na(slc[,env]),] 
  cell_numbers <- unique(slc$cell_number)
  dominant_env <- sapply(cell_numbers, function(x){
    env_tabl <- table(slc[slc$cell_number==x, env]) 
    if (length(env_tabl)==1) { 
	  dominant <- names(env_tabl[1]) 
    } else {
      env_tabl <- env_tabl[1:2]/sum(env_tabl[1:2])
      if (any(env_tabl[1:2] > 0.8)) {
        if (env_tabl[1] > 0.8){ 
		  dominant <- names(env_tabl[1]) 
		} else { 
		  dominant <- names(env_tabl[2]) 
		} 
      } else { 
	    dominant <- 'both' 
	  }
    }
    dominant
  })
  cbind(cell_numbers, dominant_env)
}

# iterate over all environmental axes

ncores <- detectCores() - 1
registerDoParallel(ncores) # set up implicit clusters 

for (env in now_names){
raster_env <- foreach(bin=stages_trunc$name, .packages=pkgs) %dopar% 
  find_dominant(dat=binned, bin=bin, env=env)

# stack rasters with conditions for each time bin
stack_env <- stack()

# If using full download:
#for (lay in 1:length(raster_env[-1])){ # leave off Holocene

# If using Pectinid data, restrict to time periods when Pectinids existed:
pect_bins <- 11:64
for (lay in pect_bins){

    # template raster to fill
	r <- raster(raster_proj) 
	
	# numbers of the cells with data
	layer_cells <- as.numeric(raster_env[[lay]][,1]) 
	layer_env <- as.factor(raster_env[[lay]][,2]) 
	
	# add dominant env values to corresponding cells
	r[layer_cells] <- layer_env 
	
	# troubleshoot
	if (length(layer_cells)==0) { print(stages_trunc$name[lay]) } 
	
	# format raster attribute table (RAT)
	r <- ratify(r) 
	rat <- levels(r)[[1]]
	rat[,env] <- levels(layer_env)
	levels(r) <- rat
	
	stack_env <- addLayer(stack_env, r)
}
names(stack_env) <- 
  # stages_trunc$name[-nrow(stages_trunc)]
  stages_trunc$name[pect_bins]
raster_name <- paste('Data/',env,'_from_all_PBDB_marine_Pectinida.tif', sep='')
writeRaster(stack_env, filename=raster_name, format='GTiff', bylayer=FALSE, 
            overwrite=TRUE, options="COMPRESS=NONE")

} # end loop through all environmental axes
stopImplicitCluster()
