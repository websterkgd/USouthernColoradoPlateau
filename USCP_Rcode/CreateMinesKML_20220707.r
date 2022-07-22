#Creates a kml file of mine features

rm(list=ls())
getwd()

#loading some packages
require(vegan)
require(ggplot2)
require(devtools)
require(dplyr)
require(stringr)
require(mapdata)
require(ggmap)
require(gstat)
require(maps)
require(raster)
require(plotKML) #require plotKML for readGPX
        
# Pulling in U mines  

#In database 12 pulling out the database file

Mines <- read.csv("ULD_albers.csv")

f.Mines <- Mines[which(Mines[,4] > 33 & Mines[,4] < 38.5 & 
                         Mines[,5] > -112.5 & Mines[,5] < -107),]


# prepare coordinates, data, and proj4string
m.coords <- f.Mines[ , c(5,4)]   # coordinates
m.data   <- f.Mines[ , 11]          # data
crs    <- CRS("+proj=longlat +datum=WGS84 +no_defs") # proj4string of coords

# make the SpatialPoints object
mspat <- SpatialPointsDataFrame(coords = m.coords,
                                data = as.data.frame(m.data),
                                proj4string = crs)

##The one line of code below writes our SpatialPointsDataFrame to a KML File
writeOGR(mspat, dsn="Mspat.kml", layer= "mspat", driver="KML")


