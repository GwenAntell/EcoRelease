library(stringr)
library(divDyn)
library(velociraptr)
library(vegan)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
# save package names as object to put on all cores later
pkgs <- c('sp','raster','rgeos') 

# packages that speed up computation
# WARNING: this will give R the power to compute on all cores.
# If you'd rather run things in a loop (1 core), change 'foreach' to 'for(){}' below
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

# This script links to a dataset that is 10% of the size of the full one used for the paper,
# because users may experience memory timeouts attempting to download the full file.
# To download the data run in the paper, set 'base_name=Brachiopoda,Bivalvia' instead of '=Pectinida'.
# A suggestion for changing memory settings, if necessary:
  # getOption("timeout")
  # options(timeout=1200) # set to 20 Minutes
pbdUrl <- paste0('https://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Pectinida',
                 '&taxon_reso=genus&idmod=!af,cf,sl,if,eg,qm,qu&interval=Cambrian,Holocene&time_rule=overlap',
                 '&envtype=!terr&show=attr,class,ident,img,abund,coll,coords,loc,paleoloc,strat,stratext,lith,lithext,',
                 'geo,timebins,resgroup,ref,refattr,ent')
pbd <- read.csv(pbdUrl, header=TRUE, skip=21, stringsAsFactors=FALSE)
pbd <- subset(pbd, !is.na(pbd$paleolat) & !is.na(pbd$paleolng)) 

#################################
#### assign stage affiliation ###
#################################

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
  'Homerian','Ludfordian','Pragian','Eifelian','Bashkirian','Kasimovian', 
  'Sakmarian','Kungurian', # no Artinskian species records
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

# assign stages to records
# loosely based on {velociraptr} constrainAges function
for (i in 1:(nrow(stages_trunc)-1)) { 
  stage <- stages_trunc$name[i]
  next_stage <- stages_trunc$name[i+1] 
  # include names of all bins that are lumped with the focal time bin
  included <- stages$name[which(stages$name==stage):(which(stages$name==next_stage)-1)]
  
  # define the buffered time window for the focal bin
  b <- stages_trunc$b_round[i]
  t <- stages_trunc$t_round[i]
  
  earlypos <- union(which(pbd$max_ma<=b & pbd$max_ma>=t), which(pbd$early_interval %in% included)) 
  latepos <- union(which(pbd$min_ma<=b & pbd$min_ma>=t), which(pbd$late_interval %in% included))
  keepers <- intersect(earlypos, latepos)
  
  if (length(keepers)>0){
    cleanbin <- cbind(pbd[keepers,], stage)
    binned <- rbind(binned, cleanbin)
  }
}

# Delete records assigned to more than 1 time bin
# This can happen for records with uncertain time intervals spanning boundaries
list_occ_no <- lapply(bins, 
  function(t) unique(binned$occurrence_no[binned$stage==t])
  )
dupes <- unlist(list_occ_no)[duplicated(unlist(list_occ_no))]
dupe_rows <- which(binned$occurrence_no %in% dupes)

# In the full dataset, there is duplication for 3 intervals. Modify this section for other data.
if (length(dupe_rows) > 0){
# A check of original references indicates that these should be Toarcian:
not_Pliens <- intersect(which(binned$max_ma==183 & binned$min_ma==182), dupe_rows)
binned <- binned[-intersect(not_Pliens, which(binned$stage=='Pliensbachian')),]

# A check of original references indicates that these should be Induan:
not_Chang <- intersect(which(binned$max_ma==252.17 & binned$min_ma==251.2), dupe_rows)
binned <- binned[-intersect(not_Chang, which(binned$stage=='Changhsingian')),]

# When using the Pectinida dataset example in this script, these lines are unnecessary:
  # not_Fortun <- intersect(which(binned$max_ma==485.4 & binned$min_ma==485), dupe_rows)
  # binned <- binned[-intersect(not_Fortun, which(binned$stage=='Fortunian')),]
}

# set cutoff for acceptable spatial precision
binned <- binned[-which(binned$geogscale=='basin' & binned$latlng_precision=='degrees'),] 

bins <- as.character(unique(binned$stage))

####################################
#### taxonomic name vetting  #######
####################################

# NB: the authors used a preliminary version of the {divDyn} function 'cleansp'.
# Now that the package has an official release, the published function is used instead.
# The arguments 'misspells' and 'stems' were added to the function by modifying GSA code below.
# Thus, these arguments are set to FALSE in this script because the steps are implemented already.

binned <- binned[binned$identified_rank=="species",]
binned$short_name <- cleansp(binned$identified_name, misspells=FALSE, stems=FALSE)

# remove records with 'sp.', 'spp.', or single letters
binned <- binned[!is.na(binned$short_name),] 
binned <- binned[! binned$difference %in% c('nomen dubium', 'nomen vanum', 'nomen nudum'),] 

# remove specimens without higher taxonomy entered - need that info for mixed effects models
no_higher_dat <- c(which(binned$class==''), which(binned$order==''), which(binned$family==''))
if (length(no_higher_dat) > 0){
  binned <- binned[-unique(no_higher_dat),]  
}
binned$accepted_name <- cleansp(binned$accepted_name, misspells=FALSE, stems=FALSE)

# if accepted name isn't sp level, use identified name instead
binned$unique_name <- paste(binned$phylum, binned$accepted_name)
unac <- which(! binned$accepted_rank %in% c('species','subspecies')) 
binned$unique_name[unac] <- paste(binned$phylum[unac], binned$short_name[unac])

#### catch -us/-um/-a/-ae mismatches (remove nominative/genitive suffix) ###
# NB: these lines can now be ommitted by setting stems=TRUE in cleansp function above
get_root <- function(s){
  last2 <- substr(s, nchar(s)-1, nchar(s))
  last1 <- substr(s, nchar(s), nchar(s))
  root <- s
  if (last2 %in% c('us', 'um', 'ae', 'is', 'es', 'ii')){ # check for 2-letter suffix
    root <- substr(s, 1, nchar(s)-2)
  } else { # check for 1-letter suffix
    if (last1 %in% c('a', 'i', 'e')){ root <- substr(s, 1, nchar(s)-1) }
  }
  return(root)
}

root_names <- sapply(binned$unique_name, get_root)
all_unique <- unique(binned$unique_name)
root_unique <- unique(root_names) 
# root names shorter than all names, so suffix synonyms exist

# split out Phylum/genus combination from species epithet

split_space <- strsplit(binned$unique_name, split='_')
gen <- sapply(split_space, function(n) n[1] )

# find discrepancies in number of species within genus before and after removing roots
for (g in names(table(gen))){
  full_names <- all_unique[grep(g, all_unique)]
  trunc_names <- root_unique[grep(g, root_unique)]
  if (length(full_names) == length(trunc_names)) { next } # no synonyms found
  
  for (a in trunc_names){
    syn_rows <- which(root_names==a) # rows of duplicate/synonymous names
    syn_names <- unique(binned$unique_name[syn_rows])
    if (length(syn_names) >1){ 
      
      # replace synonym with first alphabetically 
      replacement <- syn_names[1]
      be_replaced <- syn_names[-1]
      binned$unique_name[which(binned$unique_name %in% be_replaced)] <- replacement
    }
  } # loop through species epithets 
  
} # loop through genus epithets

#### replace common dipthongs with alternative or easily misspelled versions ###
# NB: these lines can now be ommitted by setting misspells=TRUE in cleansp function above

misspell <- function(term){
  term <- gsub('ue','u', term, fixed=TRUE)
  term <- gsub('ae','e', term, fixed=TRUE)
  term <- gsub('ll','l', term, fixed=TRUE)
  term <- gsub('ss','s', term, fixed=TRUE)
  term <- gsub('y', 'i', term, fixed=TRUE)
  term
}
altered <- unique(misspell(binned$unique_name))

# number of pairs assumed synonymous
length(unique(binned$unique_name))-length(altered) 

for (g in names(table(gen))){
  originals <- all_unique[grep(g, all_unique)]
  revisions <- altered[grep(misspell(g), altered)]
  
  if (length(originals) == length(revisions)) { next } # no synonyms found
  
  tab_syns <- table(misspell(originals)) # table value = 1 if spelling is ok
  issue_names <- names(tab_syns[tab_syns>1])
  
  for (s in issue_names){
    syns <- originals[ misspell(originals)==s ] # entries in list that are synonymous
    
    # replace synonym with first alphabetically 
    replacement <- syns[1]
    be_replaced <- syns[-1]
    binned$unique_name[which(binned$unique_name %in% be_replaced)] <- replacement
  } # loop through species epithets 
} # loop through genus epithets

# Notes:
# There are lines ommitted here that change individual taxon names
# (synonymize or correct misspellings) and delete individual records (where invalid).
# There are 3 reasons for the omission:
# 1. This script is general for a PBDB download of any taxon,
# and aims to be concise while retaining all necessary lines
# 2. Authors will implement corrections in the PBDB directly, 
# so post hoc corrections will become redundant
# 3. Some of the sources used to vet individual species' records are under embargo

#################################################
#### convert points to raster grid cells ########
#################################################

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
binned$centroid_lat <- binned$centroid_long <- NA
binned[c('centroid_long', 'centroid_lat')] <- xyFromCell(raster_proj, binned$cell_number)

#################################################
#### determine environment for each grid cell ###
#################################################
# This section is not necessary unless the cross-species analysis will be run,
# and will account for environmental preferences.
# If running assemblage analysis (time series analysis) only, skip to next section.

# Do not make inference for mixed lithologies
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
             # include reefs ('basin' and 'slope/ramp' are deep water reefs):
             "intrashelf/intraplatform reef", "reef, buildup or bioherm",
             "perireef or subreef", "platform/shelf-margin reef") 
deep <- c("basinal (carbonate)", "basinal (siliceous)", "basinal (siliciclastic)",
          "deep-water indet.", "deep subtidal indet.", "deep subtidal ramp",
          "deep subtidal shelf", "offshore", "offshore indet.",
          "offshore shelf", "slope", "submarine fan", "offshore ramp", 
          "basin reef", "slope/ramp reef") 

# categorise each record for each enviro type 
l <- b <- rep(NA, nrow(binned))

# Lithology
l[binned$lithology1 %in% carb] <- "c"  
l[binned$lithology1 %in% clast] <- "s"

# Bathymetry (p = proximal = shallow)
b[binned$environment %in% shallow] <- "p" 
b[binned$environment %in% deep] <- "d"

env_axes <- c('lithnow','bathnow')
binned <- cbind(binned, data.frame('lithnow'=l,'bathnow'=b, stringsAsFactors=TRUE) )

#### species' environmental preferences ###

# Infer environmental conditions for records that lack them.
# Try to use data from all PBDB occurrences in the same cell.
# If both environmental conditions are found in the cell, make no inference:
# the species could have had any of the possible preferences.

# Load environmental raster stacks.
# The code that generated these rasters is available: 
# categorize_enviros_from_all_PBDB_marine.R
r_names <- paste0('Data/',c('lithology','bathymetry'),'_from_all_PBDB_marine.tif')
env_r <- lapply(r_names, brick)
env_r <- lapply(env_r, setNames, stages_trunc$name[-1])

fill_dat <- function(bin, occ_dat, env_dat, bin_col, cell_col, env_axes){
	occ_bin <- occ_dat[ occ_dat[,bin_col]==bin, ]
	env_bin <- data.frame(sapply(env_dat, function(r) getValues(r[[bin]]) ))
	colnames(env_bin) <- env_axes
	
	new_dat <- sapply(1:nrow(occ_bin), function(i){
		record <- occ_bin[i,]
		if (any(is.na(record[env_axes])==TRUE)){
		  env2infer <- env_axes[which(is.na(record[env_axes])==TRUE)]
		  cell <- as.numeric(record[cell_col])
		  record[env2infer] <- env_bin[cell,env2infer]
		  output <- record[env_axes]
		
		# all conditions aleady known; nothing to replace
		} else { output <- record[env_axes] } 
		
		output <- as.matrix(output)
	}) # apply across all rows
	t(new_dat)
} # end function

ncores <- detectCores() - 1
registerDoParallel(ncores) 
inferred_env <- foreach(bin=bins, .packages=pkgs) %dopar%	
				fill_dat(bin=bin, occ_dat=binned, env_dat=env_r, 
				bin_col='stage', cell_col='cell_number', env_axes=env_axes)
stopImplicitCluster()

dat2add <- data.frame(do.call('rbind', inferred_env), stringsAsFactors=FALSE)
colnames(dat2add) <- env_axes

for (e in 1:length(env_axes)){
	rat <- levels(env_r[[e]][[1]])[[1]]
	rat$category <- as.character(rat$category)

	numeric2fix <- intersect( which(! dat2add[,e] %in% letters), which(is.na(dat2add[,e])==FALSE) )
	dat2add[numeric2fix,e] <- vapply(dat2add[numeric2fix,e], FUN.VALUE='', FUN=switch, 
		'0'=NA,'1'=rat$category[rat$ID==1],'2'=rat$category[rat$ID==2],'3'=rat$category[rat$ID==3])
		
	# make no inference if both enviro conditions in cell
	dat2add[which(dat2add[,e]=='both'), e] <- NA 
	
	# (necessary in order to tabulate frequency below)
	dat2add[,e] <- as.factor(dat2add[,e]) 
}

binned[-which(binned$stage=='Fortunian'), env_axes] <- dat2add

# determine sp environmental preference (if any) from all occurrences
pref <- function(sp, name_vect, env_vect){
  elements <- which(name_vect==sp)
  env_tabl <- table(env_vect[elements])
  if (sum(env_tabl[1:2])==0) { pref <- 'uk' } else {
    env_tabl <- env_tabl[1:2]/sum(env_tabl[1:2])
    if (any(env_tabl[1:2] > 0.8)) {
      if (env_tabl[1] > 0.8){ 
        pref <- names(env_tabl[1]) 
      } else { 
        pref <- names(env_tabl[2]) 
      } 
    } else { pref <- 'both' }
  }
  list(pref, elements)
}

# list each species' preference along each enviro axis
registerDoParallel(ncores) 
lith_pref_l <- foreach(sp = unique(binned$unique_name)) %dopar% 
	pref(sp=sp, name_vect=binned$unique_name, env_vect=binned$lithnow)
bath_pref_l <- foreach(sp = unique(binned$unique_name)) %dopar% 
	pref(sp=sp, name_vect=binned$unique_name, env_vect=binned$bathnow)
stopImplicitCluster()

# add each environmental preference of a focal species to all its occurrence records
# (rows are not in order so this is best in a loop)
binned[,c('pref_lith','pref_bath')] <- NA
for (x in 1:length(lith_pref_l)){ 
  prefs <- sapply(list(lith_pref_l, bath_pref_l),
    function(pref_list) pref_list[[x]][[1]] )
  
  # rows where species occurs
  x_rows <- lith_pref_l[[x]][[2]] 
  
  binned[x_rows, c('pref_lith','pref_bath')] <- 
    matrix(rep(prefs, length(x_rows)), ncol = 2, byrow = TRUE)
}


##############################
#### shorten dataframe #######
##############################

# Remove occurences that are singular in space-time
nm_freq <- table(binned$unique_name)
singles <- names( which(nm_freq == 1) )
binned <- binned[! binned$unique_name %in% singles,]

# Many nested functions:
# subset by stage, then species, then find unique species-cells combinations

# omit duplicate cell-species combinations
st_sp_summary <- function(dat, sp){
  sp_st <- dat[dat$unique_name==sp,] 
  temp_list <- lapply(unique(sp_st$cell_number), 
    function(x) sp_st[which(sp_st$cell_number==x)[1],]
    )
  do.call('rbind', temp_list)
}

# apply st_sp_summary over species in a stage
st_summary <- function(dat, bin){
  slc <- dat[dat$stage==bin,]
  temp_list <- lapply(unique(slc$unique_name), st_sp_summary, dat=slc)
  do.call('rbind', temp_list)
}

registerDoParallel(ncores)
stage_dat_list <- foreach(bin=bins) %dopar% st_summary(bin=bin, dat=binned)
stopImplicitCluster()
fin <- do.call('rbind', stage_dat_list)

# reduce object size so save memory space during later analysis
cols2keep <- c("phylum","class","order","family","genus","unique_name","stage","cell_number",
               "centroid_long","centroid_lat","pref_lith","pref_bath")
fin <- fin[, cols2keep]

# export dataframe
df_name <- paste('Data/DS1_unique_taxonLocationTime_occurrences_Pectinida.csv', sep='')
write.csv(fin, file = df_name, row.names = FALSE)

# If running subsequent scripts directly from this one, such that all objects are still loaded,
# rename the lumped stage object for later use.
stages <- stages_trunc
