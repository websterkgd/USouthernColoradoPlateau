#Pulling out the Coconino Aquifer greater than 3 K TDS

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/UraniumContamination/Data_WCM')

#Calling count color from the count colors package
require("countcolors")

DSC2 <- jpeg::readJPEG('CoconinoDeChellyTDS_Cropped.JPEG')

#identifying cluster of interest
kmeans.clusters <- 
  colordistance::getKMeanColors("CoconinoDeChellyTDS_Cropped.JPEG", 
                                n = 15, plotting = TRUE)
colors <- colordistance::extractClusters(kmeans.clusters)
colors_o <- colors[order(colors[,4], decreasing = TRUE),]
print(colors_o)

### Converting to a polygon + kml for future use
require(sp) 
require(rgdal) #used for writing the shape file #also will be retired by 2023
require(ggmap)
require(ggplot2)
require(raster)
require(rgeos)

#pulling out 3000 - 25,000
upper.r.C <- c(0.69, 0.31, 0.60)

#pulling out 3000 - 25,000
lower.r.C <- c(0.61, 0.23, 0.52)  

MR <- countcolors::rectangularRange(DSC2, 
                                     upper = upper.r.C, 
                                     lower = lower.r.C,
                                     target.color = "green",
                                     plotting = TRUE)

##Pulling out solids areas
#pull out pixels
Ix <- MR$pixel.idx

#plotting pixels in a somewhat correct orientation
plot(Ix[,2],-Ix[,1])

# the pixels stay fixed

np <- 39.524254   # -1N = 39.524254,
sp <- 33.742784 # -306 S = 33.742784
nssp <-(sp-np)/(-305) # NS line:  

wp <- -111.955907 # 1W = -111.955907
ep <- -107.808380 # 180E = -107.808380
ewsp <- (ep-wp)/(179) # EW line

ccW <- ewsp*range(Ix[,2])[1] + (wp-ewsp) #westmost point 0 - 1000 CC Aq
ccE <- ewsp*range(Ix[,2])[2] + (wp-ewsp) #eastmost point 0 - 1000 CC Aq
ccN <- nssp*-range(Ix[,1])[1] + (np-nssp) #northtmost point 0 - 1000 CC Aq
ccS <- nssp*-range(Ix[,1])[2] + (np-nssp)  #southmost point 0 - 1000 CC Aq

bw <- max(Ix[,2])-min(Ix[,2])
bn <- max(-Ix[,1])-min(-Ix[,1])

dE <- Ix[,2]*(ccE-ccW)/bw + ccW - min(Ix[,2])*(ccE-ccW)/bw

dN <- -Ix[,1]*(ccN-ccS)/bn + ccS - min(-Ix[,1])*(ccN-ccS)/bn

plot(dE, dN)


SH = cbind(dE, dN) #SH Second Highest

### Pulling out outline
require(concaveman)
SHo <- concaveman(SH, concavity = 1, length_threshold = 0.4) 
plot(SHo)

###########   converting SHo to polygon
SHop <- Polygon(SHo)
SHops <- Polygons(list(SHop),1)
SHosps = SpatialPolygons(list(SHops))
plot(SHosps)

data <- data.frame(99.9)
row.names(data) <- '3-25k'
SHospsdf = SpatialPolygonsDataFrame(SHosps, data, match.ID = F)
SHospsdf
plot(SHospsdf)

#writing file
#writeOGR(SHospsdf, dsn ='CoCoAq_3-25.kml', layer='CCAq_3-25', driver="KML")

p.CA_3t25 <- ggplot() +
  geom_polygon(data = SHosps, aes(x=long, y = lat, group = group))

p.CA_3t25 # may not plot but it is there in pdf preview

ggsave(plot=p.CA_3t25, "CoCoAq_3t25.tiff", device = "tiff")

rasterMap <- stack("CoCoAq_3t25.tiff")

# Get the GeoSpatial Components
lat_long <- ggplot_build(p.CA_3t25)$layout$panel_params[[1]][c("x.range","y.range")]

# Supply GeoSpatial  data to the StackedRaster
extent(rasterMap) <- c(lat_long$x.range,lat_long$y.range)
projection(rasterMap) <- CRS("+proj=longlat +datum=WGS84")

###write file
#writeOGR(SHospsdf, dsn ='CoCoAq_3-25.kml', layer='CCAq_3-25', driver="KML")

#####pulling out > 25000

c_h <- c(0.79, 0.31, 0.60)

upper.h.c <- c_h + 0.03

#pulling out 1 - 300
lower.h.c <- c_h - 0.08

h.C <- countcolors::rectangularRange(DSC2, 
                                     upper = upper.h.c, 
                                     lower = lower.h.c,
                                     target.color = "green",
                                     plotting = TRUE)
#converting area to points
h.C.I <- h.C$pixel.idx


#plotting pixels in a somewhat correct orientation
plot(h.C.I[,2],- h.C.I[,1])

# the pixels stay fixed

np <- 39.524254   # -1N = 39.524254,
sp <- 33.742784 # -306 S = 33.742784
nssp <-(sp-np)/(-305) # NS line:  

wp <- -111.955907 # 1W = -111.955907
ep <- -107.808380 # 180E = -107.808380
ewsp <- (ep-wp)/(179) # EW line

ccW <- ewsp*range(h.C.I[,2])[1] + (wp-ewsp) #westmost point 0 - 1000 CC Aq
ccE <- ewsp*range(h.C.I[,2])[2] + (wp-ewsp) #eastmost point 0 - 1000 CC Aq
ccN <- nssp*-range(h.C.I[,1])[1] + (np-nssp) #northtmost point 0 - 1000 CC Aq
ccS <- nssp*-range(h.C.I[,1])[2] + (np-nssp)  #southmost point 0 - 1000 CC Aq

bw <- max(h.C.I[,2])-min(h.C.I[,2])
bn <- max(-h.C.I[,1])-min(-h.C.I[,1])

dE <- h.C.I[,2]*(ccE-ccW)/bw + ccW - min(h.C.I[,2])*(ccE-ccW)/bw

dN <- -h.C.I[,1]*(ccN-ccS)/bn + ccS - min(-h.C.I[,1])*(ccN-ccS)/bn

plot(dE, dN)
ov = cbind(dE, dN) #Highest

### Pulling out outline
plot(ov)

oCC <- ov[which(ov[,2] > 34.6 & ov[,2] < 35.25 & ov[,1] < -109.2),]
plot(oCC)

oCCp <- concaveman(oCC, concavity = 1, length_threshold = 0.2) 
plot(oCCp)

### Converting to a polygon + kml for future use
oCCpp <- Polygon(oCCp)
oCCpps <- Polygons(list(oCCpp),1)
oCCspps = SpatialPolygons(list(oCCpps))
plot(oCCspps) ### Final polygon > 25000 TDS

data <- data.frame(99.9)
row.names(data) <- 'g25k'
oCCspdf = SpatialPolygonsDataFrame(oCCspps, data, match.ID = F)
oCCspdf
plot(oCCspdf)

#writing file
#writeOGR(oCCspdf, dsn ='CoCoAq_g25.kml', layer='CCAq_25', driver="KML")

p.CA_g25 <- ggplot() +
  geom_polygon(data = oCCspps, aes(x=long, y = lat, group = group))

p.CA_g25 # may not plot but it is there in pdf preview

ggsave(plot=p.CA_g25, "CoCoAq_g25.tiff", device = "tiff")

rasterMap <- stack("CoCoAq_g25.tiff")

# Get the GeoSpatial Components
lat_long <- ggplot_build(p.CA_g25)$layout$panel_params[[1]][c("x.range","y.range")]

# Supply GeoSpatial  data to the StackedRaster
extent(rasterMap) <- c(lat_long$x.range,lat_long$y.range)
projection(rasterMap) <- CRS("+proj=longlat +datum=WGS84")

