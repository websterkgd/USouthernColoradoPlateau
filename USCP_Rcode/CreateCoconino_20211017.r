#Trying to Pull out the Coconino Aquifer

rm(list=ls())
getwd()

#Calling count color from the count colors package
require("countcolors")

ACP <- jpeg::readJPEG('AquifersColoradoPlateau.JPEG')

#identifying cluster of interest
kmeans.clusters <- 
  colordistance::getKMeanColors("AquifersColoradoPlateau.JPEG", 
                                n = 12, plotting = FALSE)
colordistance::extractClusters(kmeans.clusters)

#pulling out upper bound coconino
upper.r.C <- c(0.55, 0.90, 1)

#pulling out lower  coconino
lower.r.C <- c(0.20, 0.55, 0.85)  

T.C <- countcolors::rectangularRange(ACP, 
                                     upper = upper.r.C, 
                                     lower = lower.r.C,
                                     target.color = "blue",
                                     plotting = TRUE)

##Pulling out coconino
#pull out pixels
Index <- T.C$pixel.idx

#plotting pixels in a somewhat correct orientation
plot(Index[,2],-Index[,1])

##pulling out Coconino from Southern NN
r <- Index[c(which(Index[,2] > 140 & Index[,2] < 273)),]

t <- r[c(which(r[,1] > 425)),]

plot(t[,2],-t[,1])

#convert to lat + long
m <- as.matrix(cbind(t[,2],-t[,1]))

plot(m[,1],m[,2])

ccW <- -111.848477 #westmost point southern CC Aq
ccE <- -108.702512 #eastmost point southern CC Aq
ccN <-  36.676738 #northmost point southern CC Aq
ccS <-  34.108818 #southmost point southern CC Aq

bw <- max(t[,1])-min(t[,1])
bn <- max(t[,2])-min(t[,2])

dE <- m[,1]*(ccE-ccW)/bw + ccW - min(m[,1])*(ccE-ccW)/bw

dN <- m[,2]*(ccN-ccS)/bn + ccS - min(m[,2])*(ccN-ccS)/bn

plot(dE, dN)

southCC = cbind(dE, dN)

### Pulling out outline
require(concaveman)
CCp <- concaveman(southCC, concavity = 1, length_threshold = 0) 
plot(CCp)

### Converting to a polygon + kml for future use
#require(rgdal)
require(sp) 
require(rgdal) #used for writing the shape file

#pulling in states and reassigning
require(mapdata)
require(ggmap)
require(ggplot2)
require(raster)

E_usa <- map_data("state")

CoCop <- E_usa
CoCop <- CoCop[1:length(CCp[,1]),]

CoCop$long <- CCp[,1]
CoCop$lat <- CCp[,2]
CoCop$group <- 1
CoCop$region <- "C-dC Aquifer"

p <- Polygon(CCp)
ps <- Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps)

proj4string(sps) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

data = data.frame(f=99.9)
spdf = SpatialPolygonsDataFrame(sps,data)
spdf

plot(spdf)

require(rgeos)
# Create the KML
#writeOGR(spdf, dsn ='CoCoAq.kml', layer='CCAq', driver="KML")

# p.CA <- ggplot() + 
#   geom_polygon(data = CoCop, aes(x=long, y = lat, group = group))
# 
# p.CA # may not plot but it is there in pdf preview
# 
# ggsave(plot=p.CA, "CoCoAq.tiff", device = "tiff")
# 
# rasterMap <- stack("CoCoAq.tiff")
# 
# # Get the GeoSpatial Components
# lat_long <- ggplot_build(p.CA)$layout$panel_params[[1]][c("x.range","y.range")] 
# 
# # Supply GeoSpatial  data to the StackedRaster 
# extent(rasterMap) <- c(lat_long$x.range,lat_long$y.range)
# projection(rasterMap) <- CRS("+proj=longlat +datum=WGS84")


# 
# ###Plot an outline
# #plot(CCp, col= 'white')
# #polygon(x = CCp[,1], y = CCp[,2]) #for plotting
# 
# writeOGR(sps, dsn="CoconinoAquifer.kml", layer= "Coconino", driver="KML")
