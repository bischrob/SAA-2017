#################################################################################
## Name:Random Least Cost Paths
## Author: Robert J. Bischoff
## Date: last updated - 3/26/2019
## Purpose: This script requires a starting point and will generate least cost
## paths (LCP) from the starting point to randomly generated points within a 
## designated extent. The purpose is to determine whether LCP is beneficial for 
## linear regression. A DEM will be downloaded from the NED. This script only 
## works within the United States. The results provide the euclidean (straight-line)
## distance, the length of the LCPs, the average walking speed, the time to 
## travel the path using Tobler's Hiking Function, the time to walk the 
## euclidean distance if flat, and the time to walk the LCP if flat.
## Licence: MIT
##    (https://opensource.org/licenses/MIT)
#################################################################################

# load packages using pacman function
# install.packages('pacman') # uncomment if pacman is not installed
pacman::p_load(compiler,tidyverse,raster,gdistance,rgdal,rgeos,FedData,here,sf)

# Record run time
ptm <- proc.time()

# Load compiler to speed up code.
enableJIT(3)

# The section includes all variables that should be updated before running
# this script.

dir <- here("Random Least Cost Paths for SAA 2017 Project.R")
  
projName <- "Random Points"    # name of project
ptStart <- "StartingPoint"     # starting point
ptsRandom <- "RandomPoints"    # change the name of the shapefile with the
                              # randomly generated points
LCP <- "LCPs"                 # change the name of the LCPs
NED <- "Random Points"         # name of the saved NED
ptStartCoords <-  c(-112.065235, 33.332575) # coordinates for the starting point

# This code will determine the UTM zone based on the longitudinal coordinate
# and will assist in determining the appropriate coordinate reference system
# code from - http://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude
long2UTM <- function(long) {  # do not modify this function
  (floor((long + 180)/6) %% 60) + 1
  
}
long2UTM(ptStartCoords[1]) # this returns the UTM zone
crs1 <- 26912
areaSize <- 5000    # this is the size of the area in M to use in this script.
DEMAggregateN <- 6    # this number determines the cell size of the raster
                      # a value of 1 will keep the original 30 M DEM, a value
                      # of 2 will result in a 60 M cell, etc. - it is 
                      # best to increase this number for testing the script
                      # and on computers not designed for high power computing
ptsRandomN <- 10   # this determines the number of random points to generate
                      # LCPs to.


##############################################################################
# This section contains the script generating the LCPs and does not need to be
# modified unless reusing saved data or if other modifications are desired.
# To use saved files, place a "#" in front of the unneeded code and remove the
# "#" in front of the code to load files.

# use longitude / latitude coordinates to get raster
starting.point <- st_sf(data = tibble(id = 1),
                          geometry = 
                          st_sfc(st_point(ptStartCoords), crs = 4326))
starting.point <- st_transform(starting.point,crs1)
xMinus <- st_coordinates(starting.point)[,1] - areaSize  
xPlus <- st_coordinates(starting.point)[,1] + areaSize
yMinus <- st_coordinates(starting.point)[,2] - areaSize
yPlus <- st_coordinates(starting.point)[,2] + areaSize
coords <- matrix(c(xMinus,yMinus,
                   xMinus,yPlus,
                   xPlus,yPlus,
                   xPlus,yMinus,
                   xMinus,yMinus), ncol = 2, byrow = T)
poly <- st_sf(st_sfc(st_polygon(list(coords))),
              crs = crs1)
polySP <- as_Spatial(poly)

# writeOGR(starting.point, dsn = ".", layer = ptStart,
                   # driver = "ESRI Shapefile", overwrite = T)

# Get raster using FedData package; res is "1" for 1 arcsec or "13" for
# 1/3 arcsec
NED <- get_ned(template = polySP, label = projName, res = "1")
NED <- projectRaster(NED, crs=polySP@proj4string) # raster must be projected for
                                     # this script
  
# Save projected raster
# writeRaster(NED, filename=NED.name, format="GTiff", overwrite=TRUE)

# Or load data from file 
# NED <- raster(paste0(NED.name,".tif"))

# Generate LCP
# Downsize raster
NED <- aggregate(NED, fact = DEMAggregateN, fun = mean)

# Create function to generate elevation difference.
altDiff <- function(x){x[2] - x[1]}

# Create transition layer. The transition layer measures the difficulty
# of transitioning from one cell to another.
hd <- transition(NED, altDiff, 8, symm=FALSE) # negative values error
                                              # is expected
slope <- geoCorrection(hd) # Correct for slope.
adj <- adjacent(slope, cells=1:ncell(slope), pairs=TRUE, directions=8)
                                  # restricts calculations to adjacent cells
speed <- slope
# Tobler's Hiking Function
speed[adj] <- (6000 * exp(-3.5 * abs(slope[adj] + 0.05))) # 6000 is the value
                                  # required for a walking speed of 5k/h on
                                  # flat terrain
                   
# The following corrects for diagonal distances between cells
Conductance <- geoCorrection(speed)

# Generate random points
# Set the minimum and maximum values for the locations and ensures the values
# are not within a certain extent of the raster edges (1km default).
x.min <- NED@extent@xmin + 1000
x.max <- NED@extent@xmax - 1000
y.min <- NED@extent@ymin + 1000
y.max <- NED@extent@ymax - 1000
 
# Generate random x and y coordinates and pair them in a matrix.
x <- sample(x.min:x.max, ptsRandomN, replace = T)
y <- sample(y.min:y.max, ptsRandomN, replace = T)
xy.df <- data.frame(cbind(1:length(x),x,y))
xy.spdf <- SpatialPointsDataFrame(xy.df[,2:3], xy.df, proj4string = NED@crs)
# writeOGR(xy.spdf, dsn = ".", layer = ptsRandom, driver = "ESRI Shapefile",
         # overwrite = T)

# Or use points loaded from a shapefile
# xy.spdf <- readOGR(".", ptsRandom)

# Determine Least Cost path and calculate walking time and distances

# Load starting point from shapefile.
# starting.point <- readOGR('.', ptStart)
ending.point <- xy.spdf
# create dataframe for distance information
paths.df <- data.frame("ID" = 1:nrow(ending.point@data),
                       "Dist.km" = 1:length(ending.point),
                       "Path.km" = 1:length(ending.point),
                       "Speed.kh" = 1:length(ending.point),
                       "Time.hrs" = 1:length(ending.point))
a <- starting.point # origin of LCP
a.l <- nrow(ending.point) # number of loops

# Create 1st path
b <- SpatialPoints(matrix(ending.point@coords[1,], ncol = 2),
                   proj4string = NED@crs) # destination for LCP
if(identical(a,b)) # does not run script if start and end points are the same
{} else {
  sl <- shortestPath(Conductance, a, b, output = "SpatialLines") # create LCP
  paths <- sl
  time <- costDistance(Conductance, a,b) # Time from point a to b in hours
  e.dist <- gDistance(a,b)/1000 # euclidean or straight-line distance in km
  path.l <- gLength(sl)/1000 # length of path in km
  speed <- (gLength(sl)/1000) / time # average speed of travel in k/h
  paths.df[1,2:5] <- cbind(e.dist,path.l,speed,time) # add results to data frame
}

# Starts loop running above code for all destination points
for(z in 2:a.l){ 
  
  # Track time for each loop
  time.tracker <- Sys.time()
  
   b <- SpatialPoints(matrix(ending.point@coords[z,], ncol = 2),
                     proj4string = NED@crs)
  
  if(identical(a,b)) {} else {
    sl <- shortestPath(Conductance, a, b, output = "SpatialLines")
    id <- z
    sl@lines[[1]]@ID <- as.character(id)
    paths <- rbind(paths, sl)     # add LCP to all previous paths
    time <- costDistance(Conductance, a,b) # this result is in hours
    e.dist <- gDistance(a,b)/1000 # this result is in kilometers
    path.l <- gLength(sl)/1000 # this result is in kilometers
    speed <- (gLength(sl)/1000) / time # this result is kilometers per hour
    paths.df[z,2:5] <- cbind(e.dist,path.l,speed,time)
    
    # Print elapsed time per loop.
    etime <- Sys.time() - time.tracker
    print(paste0("Elapsed time for ", z, " in ",a.l, ":"))
    print(etime)
    
  }} # end of loop

# Convert distances into time
paths.df$Dist.hrs <- paths.df$Dist.km / (6 * exp(-3.5 *
                                                 abs(0 + 0.05)))
paths.df$Path.hrs <- paths.df$Path.km / (6 * exp(-3.5 *
                                                       abs(0 + 0.05)))

# Create spatial lines data frame  
paths <- SpatialLinesDataFrame(paths, data = paths.df)

# Plot results
plot(NED)
points(starting.point, col = "red", pch = 19)
points(ending.point, col = "purple", pch = 19)
lines(paths)

# Create shapefile
# writeOGR(paths, dsn = ".", layer = LCP,
#          driver = "ESRI Shapefile", overwrite = T)
# Create csv containing the attributes
write.csv(paths.df,paste0(LCP,".csv"), row.names = F)

# Print run time
proc.time() - ptm
print(proc.time() - ptm)
