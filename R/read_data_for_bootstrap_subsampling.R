pbd <- read.csv('Data/DS1_unique_taxonLocationTime_occurrences.csv', stringsAsFactors=FALSE)

# remove Cambrian occurrences - sampling too sparse
pbd <- pbd[pbd$stage!='Fortunian',] 

# prepare stage data - lump certain stages together
library(velociraptr)
stages <- downloadTime('international ages')
stages$name <- row.names(stages)
stages <- stages[order(stages$b_age, decreasing=TRUE), ]
# The ICS chronostratigraphic chart is available at www.stratigraphy.org
# The same dates can be downloaded through the Macrostrat API, 
# directly or via {velociraptr} as written

stages2omit <- 
  c('Stage 2','Stage 3','Stage 4','Wuliuan',
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
  'Sakmarian','Kungurian', # no Artinskian species records
  'Olenekian', # otherwise only 1 abundance datum for Olenekian
  'Sinemurian', # Hettangian is too poorly sampled
  'Bajocian', # otherwise no Aalenian seed cells
  'Hauterivian','Barremian', # no seed cells for Haut., Barremian or Valanginian alone
  'Santonian', # otherwise nothing is sampled across the Coniacian
  'Thanetian',
  'Bartonian', # otherwise no environmental data for Bartonian
  'Aquitanian', # otherwise no seeds here or in Chattian
  'Serravallian', # otherwise no seed cells for Langhian
  'Messinian', # otherwise no seed cells for Messinian
  'Calabrian','Middle Pleistocene','Late Pleistocene', 
  'Northgrippian','Meghalayan') # lump all Holocene records so they're easy to remove  
stages <- stages[!(stages$name %in% stages2omit),] 

# remove Cambrian and Holocene
stages <- stages[-c(1, nrow(stages)),] 

rownames(stages) <- as.character(1:nrow(stages))

# ensure chronological sorting of occurrences, oldest to youngest
getAge <- function(x){ 
  stgBool <- stages$name==pbd$stage[x]
  stages$Midpoint[stgBool]
}
pbd$time_mid <- sapply(1:nrow(pbd), getAge)
pbd <- pbd[order(pbd$time_mid, decreasing = TRUE), ] 

rownames(pbd) <- as.character(1:nrow(pbd))
