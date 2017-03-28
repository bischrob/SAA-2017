#################################################################################
## Name:Random Least Cost Paths
## Author: Robert J. Bischoff
## Date: last updated - 3/27/2017
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
## Packages Used: "raster", "gdistance", "rgdal", "rgeos", "FedData"
#################################################################################

# Record run time
ptm <- proc.time()

# Load compiler to speed up code.
require(compiler)
enableJIT(3)

# Create a variable storing packages used, installs the package if missing, and
# loads packages.
my.packages <- c("raster", "gdistance", "rgdal",
                 "rgeos", "FedData")
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}  
lapply(my.packages, usePackage)

# The section includes all variables that should be updated before running
# this script.

wd <- "E:/SAA 2017 Poster/Random Points" # this is the directory where files will
                                        # be stored/accessed. Note "/" is used
                                        # and not "\" as is typical for windows
proj.name <- "Random Points"        # name of project - used for naming the NED
                                # the proj.name must be changed for every new
                                # point or it will cause an error
                               # download folder
shp1.name <- "StartingPoint"   # change the name of the shapefile for the 
                               # starting point
shp2.name <- "RandomPoints"    # change the name of the shapefile with the
                               # randomly generated points
shp3.name <- "LCPs"            # change the name of the LCPs
NED.name <- "Random Points"         # name of the saved NED
starting.pt <-  c(-109.2071, 37.4262) # coordinates for the starting point
                                        # in long/lat WGS84.
# This code will determine the UTM zone based on the longitudinal coordinate
# and will assist in determining the appropriate coordinate reference system
# code from - http://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude
long2UTM <- function(long) {  # do not modify this function
  (floor((long + 180)/6) %% 60) + 1
  
}
long2UTM(starting.pt[1]) # this returns the UTM zone
crs1 <- crs("+init=epsg:26914") # go to spatialreference.org to determine the
# epsg for the UTM zone the point is located in
area.size <- 10000    # this is the size of the area in M to use in this script.
dem.aggregate <- 6    # this number determines the cell size of the raster
                      # a value of 1 will keep the original 30 M DEM, a value
                      # of 2 will result in a 60 M cell, etc. - it is 
                      # best to increase this number for testing the script
                      # and on computers not designed for high power computing
random.points <- 10   # this determines the number of random points to generate
                      # LCPs to.


##############################################################################
# This section contains the script generating the LCPs and does not need to be
# modified unless reusing saved data or if other modifications are desired.
# To use saved files, place a "#" in front of the unneeded code and remove the
# "#" in front of the code to load files.

# Set working directory
setwd(wd)

# use longitude / latitude coordinates to get raster
starting.point <- SpatialPoints(matrix(starting.pt, ncol = 2),
                                CRS("+init=epsg:4326"))
starting.point <- spTransform(starting.point, crs1)
starting.point <- SpatialPointsDataFrame(matrix(starting.point@coords[1:2],
                       ncol = 2), data.frame("ID" = 1), proj4string = crs1)
starting.point@bbox[1] <- starting.point@bbox[1] - area.size  
starting.point@bbox[2] <- starting.point@bbox[2] - area.size
starting.point@bbox[3] <- starting.point@bbox[3] + area.size
starting.point@bbox[4] <- starting.point@bbox[4] + area.size
writeOGR(starting.point, dsn = ".", layer = shp1.name,
                   driver = "ESRI Shapefile", overwrite = T)

# Get raster using FedData package; res is "1" for 1 arcsec or "13" for
# 1/3 arcsec
NED <- get_ned(template = starting.point, label = proj.name, res = "1")
NED <- projectRaster(NED, crs=crs1) # raster must be projected for
                                     # this script
  
# Save projected raster
writeRaster(NED, filename=NED.name, format="GTiff", overwrite=TRUE)

# Or load data from file 
# NED <- raster(paste0(NED.name,".tif"))

# Generate LCP
# Downsize raster
NED <- aggregate(NED, fact = dem.aggregate, fun = mean)

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
                   
# Remove unnecessary variables
rm(slope)
rm(hd)
rm(adj)

# The following corrects for diagonal distances between cells
Conductance <- geoCorrection(speed)

# Remove variable
rm(speed)

# save Conductance layer for future use. 
saveRDS(Conductance, "Conductance.rds") # change name if more than one
                                        # conductance layer will be used

# Load previously generated conductance
# Conductance <- readRDS('Conductances.rds') # change name if necessary

# Generate random points
# Set the minimum and maximum values for the locations and ensures the values
# are not within a certain extent of the raster edges (1km default).
x.min <- NED@extent@xmin + 1000
x.max <- NED@extent@xmax - 1000
y.min <- NED@extent@ymin + 1000
y.max <- NED@extent@ymax - 1000
 
# Generate random x and y coordinates and pair them in a matrix.
x <- sample(x.min:x.max, random.points, replace = T)
y <- sample(y.min:y.max, random.points, replace = T)
xy.df <- data.frame(cbind(1:length(x),x,y))
xy.spdf <- SpatialPointsDataFrame(xy.df[,2:3], xy.df, proj4string = NED@crs)
writeOGR(xy.spdf, dsn = ".", layer = shp2.name, driver = "ESRI Shapefile",
         overwrite = T)

# Or use points loaded from a shapefile
# xy.spdf <- readOGR(".", shp2.name)

# Determine Least Cost path and calculate walking time and distances

# Load starting point from shapefile.
# starting.point <- readOGR('.', shp1.name)
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
writeOGR(paths, dsn = ".", layer = shp3.name,
         driver = "ESRI Shapefile", overwrite = T)
# Create csv containing the attributes
write.csv(paths.df,paste0(shp3.name,".csv"), row.names = F)

# Print run time
proc.time() - ptm
print(proc.time() - ptm)
