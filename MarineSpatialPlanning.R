####Example SDM####
#install.packages("maxnet")
#install.packages("terra")
#install.packages("predicts")

#load needed packages
library(maxnet) #for making models
library(terra) #for working w raster files
library(predicts) #for analyzing models

coral<- read.csv("coral.csv") #access coral locations
coral <- coral[,2:3] #subset the x and y coordinates, omit ID column

#crs(coral) <- "epsg:32629"

#next, list the input files as batch process
env_vars <- rast(list.files("/Users/kevinbuck/Desktop/AdvancedGIS/Lab5/CoralMaxent/Input",
                            pattern= ".asc" , #only extracting .asc files
                            full.names = T))
#plot(env_vars[[4]]) #plot one of the environmental input rasters
set.seed(123) #be able to randomly sample background same way 
backgound <- spatSample(env_vars, #spatially sampling the input rasters
                        #next line takes 1000 random points, eliminates NA values
                        #and returns the result as a raster (points)
                        size= 1000, "random", na.rm=TRUE,as.points=T
) #note the output is a raster
#points(coral) #plots the points of known coral presence
presvals <- extract(env_vars,coral) #take input conditions at presence locations
presvals <- presvals[,-1] #taking off the id column
backvals <- values(backgound) #save data of the random bkgd points

#making a df with presence absence associating input conditions w presence and
#absence from the presvals and random backgd vals data
pb <- c(rep(1,nrow(presvals)), rep(0,nrow(backvals)))
sdmdata <- data.frame(cbind(pb,rbind(presvals,backvals)))

#head(sdmdata).    #look at all the data
#tail(sdmdata)
#summary(sdmdata)
#pairs(sdmdata[,2:10], cex=0.1)

sdmdata <- na.omit(sdmdata) #get rid of NA data
#make a maxent model with the data
#syntax is maxnet([species occurrence], [inputs] )
mod <- maxnet(sdmdata[,1], sdmdata[,-1])

#look at plot of the model response curves
plot(mod)

#create a prediction from the model
prediction <- predict(env_vars, #input raster stack
                      mod, #model to use
                      clamp=F,#extrapolate, dont just interpolate
                      type="cloglog", #prediction type is essentially a  probability surface
                      na.rm=T)
plot(prediction) #look at predition
#save prediciton
writeRaster(prediction,"Maxnet_prediction_example.tif")

#checking how model fits given vector of presence, and vector of abs vals
model_evaluation <- pa_evaluate(sdmdata[sdmdata==1,],
                                sdmdata[sdmdata==0,],
                                mod)
#looking at receiver operator curve of model
plot(model_evaluation, "ROC")



####New SDMs####

#bring in environmental vars again
env_vars <- rast(list.files("/Users/kevinbuck/Desktop/AdvancedGIS/Lab5/CoralMaxent/Input",
                            pattern= ".asc" , #only extracting .asc files
                            full.names = T))
#fixing previous spelling mistakes
# background <- backgound
#saving background points for future use
# saveRDS(background,"background.rds")

background <- readRDS(file="background.rds")
backvals <- values(background)

#installing necessary packages
#install.packages("SDMtune")
library(SDMtune)

#Read in Vulnerable Marine Ecosystem Occurrence Data
VME_data <- read.csv("VMEs_database.csv")

#convert the data into spatial points to allow spatial operations
VME_data <- vect(VME_data, geom=c("MiddleLongitude","MiddleLatitude"),
                 crs="epsg:4362" , #correct coordinate ref system
                 keepgeom=T)
#assign the right CRS (UTM 29N) to environmental variables
crs(env_vars) <- "epsg:32629"

#crop occurrence data to extent of environmental data
VME_data <- crop(VME_data,project(ext(env_vars),
                                  #project env_vars to latlong
                                  from = "epsg:32629", 
                                  to = "epsg:4326"))
#convert back to simple dataframe
VME_data <- values(VME_data)
#peek at the df structure to make sure
str(VME_data)

#summarize how many of each VME indicator we have
summary_vme <- table(VME_data$VME_Indicator,useNA = "ifany")
summary_vme
#get rid of the NA classes (function takes out rows that have empty column spaces)
VME_data <- VME_data[complete.cases(VME_data[,'VME_Indicator']),]
#Standardize capitalization of "C" in Coral; grep like function basically
VME_data$VME_Indicator <- gsub('Black Coral','Black coral', VME_data$VME_Indicator)
#check to see it worked again
summary_vme <- table(VME_data$VME_Indicator,useNA = "ifany")
summary_vme

#specify min cutoff for number of records thats acceptable
cutoff <- 100
#create a character vect of names of the VME Indicator species
categories_over_100 <- names(summary_vme[summary_vme > cutoff])
categories_over_100
#running a for loop to repeat maxent modeling for all VME Indicator species over 100 records

# make file structure to save outputs
Occ_points <- list()
MX_mods <- list()
raster_predictions <- list()
PA_predictions <- list()
Evaluations <- list()

#here comes the loop
for (category in categories_over_100){
  #create species object for the category we are iterating over
  #this is making a temp table of what species is, its long and lat
  species <- na.omit(VME_data[VME_data$VME_Indicator==category,
                              c("VME_Indicator","MiddleLongitude","MiddleLatitude")])
  #convert to vect
  species <- vect(species, geom=c("MiddleLongitude", 
                                  "MiddleLatitude"), 
                  crs = "epsg:4326" , keepgeom=F)
  #project to UTM 29N CRS
  species <- project(species, "epsg:32629")
  #remove NAs and duplicate values in a given raster cell
  species <- thinData(as.data.frame(species, geom="XY"),env = env_vars)
  #save occurrence points used to train model
  Occ_points[[category]] <- species
  #Extract environmental values at presence locations
  presvals <- extract(env_vars,species[,-1])
  #get rid of useless id column
  presvals <- presvals[,-1]
  #Create Presence-Background column, as many 1s as presence, 0s as bkgd points
  pb <- c(rep(1,nrow(presvals)),rep(0,nrow(backvals)))
  #Create a combined SDM dataset, get rid of NAs
  sdmdata <- data.frame(cbind(pb,rbind(presvals,backvals)))
  sdmdata <- na.omit(sdmdata)
  #train Maxent model and save it
  mod <- maxnet(sdmdata[,1], sdmdata[,-1])
  MX_mods[[category]] <- mod
  #make a prediction raster, save it
  prediction <- predict(
    env_vars, #input
    mod, #model to use
    clamp=F,#make a prediction, dont just interpolate
    type="cloglog", #create a probability surface
    na.rm=T #get rid of NAs
    )
  raster_predictions[[category]] <- prediction
  #create an evaluation object, save it
  evaluation <- pa_evaluate(sdmdata[sdmdata==1,],
                            sdmdata[sdmdata==0,],
                            mod, type="cloglog")
  Evaluations[[category]] <- Evaluations
  #create a binary presence-absence raster based on determiend
  # max specificity and sensitivity threshold
  presence_absence_raster <- prediction
  presence_absence_raster[presence_absence_raster>threshold(evaluation)$max_spec_sens] <- 1
  presence_absence_raster[presence_absence_raster!=1] <- 0
  #save that binary raster
  PA_predictions[[category]] <- presence_absence_raster
}

#look at outputs
category <- categories_over_100[8]
#close all plot windows
dev.off()
#set layout of plot window to have 2x2 panels
par(mfrow = c(2,2))
#plot the raster prediction
plot(raster_predictions[[category]])
#add the occurrence points to the plot
points(Occ_points[[category]][,2:3])
#plot the ROC curve
plot(Evaluations[[category]], "ROC")
plot(PA_predictions[[category]])
barplot(as.matrix(threshold(Evaluations[[category]])), beside =T,
         legent.text=T, args.legend = list(x ="topright"))

#save outputs and write to file
PA_pred_rast <- rast(PA_predictions)
writeRaster(PA_pred_rast, paste0("PA_",names(PA_pred_rast), ".tif"), overwrite=T)
Prob_rasters <- rast(raster_predictions)
writeRaster(Prob_rasters, paste0("Prob_",names(PA_pred_rast), ".tif"), overwrite=T)

####Marine SpatialPlanning#####
#install.packages("prioritizr")
#install.packages("Rsymphony")

#import the necessary packages
library(terra)
library(sf)
library(prioritizr)
library(Rsymphony)

#read in the environmental rasters weve been using
env_vars <- rast(list.files("/Users/kevinbuck/Desktop/AdvancedGIS/Lab5/CoralMaxent/Input",
                            pattern= ".asc" , #only extracting .asc files
                            full.names = T))
#bring in background data
background <- readRDS(file="background.rds")
backvals <- values(background)
#how big will the grid be (unit is m^2 bc its UTM 29N)
cellarea <- 2.5e8 
#make a a hexagonal grid based on extent of env variables
#code is making matrix of extent, making feutres from corner poitns of that extent
#and making Planning Unit grid from that
PU_Layer <- st_make_grid(
  st_multipoint(matrix(ext(env_vars),ncol=2)),
                         cellsize = sqrt(2*(cellarea)/sqrt(3)),
                         square=F)
#set correct CRS for grid
st_crs(PU_Layer) = 32629
#look at the grid
plot(PU_Layer)
#Clip the Panning Unit Layer to just marine area
#Using one of the marine data layers
#plot(env_vars) #looking for one with just marine area
Marine_area <- env_vars[[4]]
crs(Marine_area) <- "epsg:32629"
#set anywhere in the mask layer that is NOT NA to a value of 1
Marine_area[!is.na(Marine_area)] <- 1
plot(Marine_area) #check to see it was done correctly
#convert that to a vector polygon file
Marine_area <- as.polygons(Marine_area)
Marine_area <- st_as_sf(Marine_area)
#clip PU layer to Marine Area only
PU_Layer = st_intersection(PU_Layer,Marine_area)
PU_Layer <- st_sf(PU_Layer)
#add an area column to planning unit layer
#You will see that some cells are tiny bc of where
#clipping intersected hexagons w coast and edges tc
PU_Layer$area <- st_area(PU_Layer)
#Discarding cells that have less than half their original area
PU_Layer$area <- as.numeric(PU_Layer$area)
PU_Layer <- PU_Layer[PU_Layer$area>(cellarea/2),]
#look at it again
plot(env_vars[[4]])
plot(PU_Layer,add=T)
#Open a raster of fishing pressure, defined by swept area ratio
# SAR is the area fished / area of cell
SAR <- rast("SAR.tif")
plot(SAR,add=T)
#Extract the mean SAR in each cell
PU_Layer$SAR <- fast_extract(SAR,PU_Layer)
#where the Swept Area Ratio is NA, set to 0
PU_Layer$SAR[is.na(PU_Layer$SAR)] <- 0
#Scale that fishing effort between 0 and 1
#make a new function to do that
normalize <- function(x, na.rm = TRUE ){
  return((x-min(x))/ (max(x)-min(x)))
}
#apply the normalize function to the SAR column
PU_Layer$SAR <- normalize(PU_Layer$SAR)
#now adding 1 so there are no 0 vals
PU_Layer$SAR <- PU_Layer$SAR + 1
#normalizing cell area, and again add 1
PU_Layer$area <- normalize(PU_Layer$area)
PU_Layer$area <- PU_Layer$area + 1
#create a cost parameter, found to nearest whole integer
PU_Layer$cost <- as.integer(
  ceiling(10*(PU_Layer$SAR * PU_Layer$area))
  )
#read in Presence-Absence rasters from last week
PA_pred_rast <- rast(list.files(
  "/Users/kevinbuck/Desktop/AdvancedGIS/Lab8MarineSpatialPlanning/WorkingData/PA_Rasters",
  full.names = T, pattern = ".tif"))
#Define a systematic conservation problem
problem1 <- 
  #problem inputs are the Planning Unit layer and conservation features
  problem(PU_Layer,PA_pred_rast, cost_column = "cost") %>% 
  #minimize cost of the solution while achieving targets
  add_min_set_objective() %>%
  #We dont want the solution to be to spatially broken up (tiny marine protected areas)
  add_boundary_penalties(penalty = 0.00025) %>% 
  #what are the conservation targets? We will do for 30% of the extents modeled
  add_relative_targets(0.3) %>% 
  add_binary_decisions()

#check the problem setup
print(problem1)

#solve the problem
solution1 <- solve(problem1, force = T)

#look at the solution plotted out
par(mfrow=c(1,1))
plot(
  st_as_sf(solution1[, "solution_1"]), main = "Prioritization",
  pal=cm.colors(2,alpha=1) #pink is gonna be the marine reserves to save
)
#write out the shapefile
st_write(solution1, "Solution1.shp")


