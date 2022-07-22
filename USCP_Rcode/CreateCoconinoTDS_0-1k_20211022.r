#Pulling out the Coconino Aquifer 0 - 1000 TDS

rm(list=ls())
getwd()

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

############## 

#pulling out 0 - 1000
#with 15 colors -> color 3

c_3 <-c(0.8519369, 0.7990071, 0.9171506)

upper.lm.C <- c_3 + 0.014

#pulling out
lower.lm.C <- c_3 - 0.014 

lm.C <- countcolors::rectangularRange(DSC2, 
                                    upper = upper.lm.C, 
                                    lower = lower.lm.C,
                                    target.color = "green",
                                    plotting = TRUE)
#converting area to points
lm.C.I <- lm.C$pixel.idx

#plotting pixels in a somewhat correct orientation
plot(lm.C.I[,2],- lm.C.I[,1])

range(lm.C.I[,1]) #
range(lm.C.I[,2]) #

# the pixels stay fixed

np <- 39.524254   # -1N = 39.524254,
sp <- 33.742784 # -306 S = 33.742784
nssp <-(sp-np)/(-305) # NS line:  

wp <- -111.955907 # 1W = -111.955907
ep <- -107.808380 # 180E = -107.808380
ewsp <- (ep-wp)/(179) # EW line

ccW <- ewsp*range(lm.C.I[,2])[1] + (wp-ewsp) #westmost point 0 - 1000 CC Aq
ccE <- ewsp*range(lm.C.I[,2])[2] + (wp-ewsp) #eastmost point 0 - 1000 CC Aq
ccN <- nssp*-range(lm.C.I[,1])[1] + (np-nssp) #northtmost point 0 - 1000 CC Aq
ccS <- nssp*-range(lm.C.I[,1])[2] + (np-nssp)  #southmost point 0 - 1000 CC Aq

bw <- max(lm.C.I[,2])-min(lm.C.I[,2])
bn <- max(-lm.C.I[,1])-min(-lm.C.I[,1])

dE <- lm.C.I[,2]*(ccE-ccW)/bw + ccW - min(lm.C.I[,2])*(ccE-ccW)/bw

dN <- -lm.C.I[,1]*(ccN-ccS)/bn + ccS - min(-lm.C.I[,1])*(ccN-ccS)/bn

plot(dE, dN)

L = cbind(dE, dN) #L Lowest     

##pulling out 0 - 1000

### Pulling out outline
#need to pull out shapes individually
require(concaveman)
LN <- L[which(L[,2] > 38.6),] #Northernmost region
plot(LN)

LNo <- concaveman(LN, concavity = 1, length_threshold = 0.15) 
plot(LNo)

### Converting to a polygon + kml for future use
LNop <- Polygon(LNo)
LNops <- Polygons(list(LNop),3)
LNosps = SpatialPolygons(list(LNops))
plot(LNosps) #Upper most polygon 0 - 1000 TDS

#need to pull out shapes individually
plot(L)
t <- L[c(which(L[,1] > -111 & L[,2] < 38.5)),] #Middle-High polygon 0 - 1000

plot(t)

LHM <- t[c(which(t[,1] < -109.4 & t[,2] > 36.6)),]

plot(LHM)

LHMo <- concaveman(LHM, concavity = 1, length_threshold = 0.2) 
plot(LHMo)

### Converting to a polygon + kml for future use
LHMop <- Polygon(LHMo)
LHMops <- Polygons(list(LHMop),3)
LHMosps = SpatialPolygons(list(LHMops))
plot(LHMosps) # Middle High polygon 0 - 1000 TDS

#need to pull out shapes individually
plot(L)

t <- L[c(which(L[,1] > -111.2 & L[,2] < 36.6)),] #Middle Low polygon 0 - 1000 TDS
plot(t)
 
LNM <- t[c(which(t[,1] < -108.7 & t[,2] > 35.2)),]
plot(LNM)
 
LNMo <- concaveman(LNM, concavity = 1, length_threshold = 0.2) 
plot(LNMo)
 
### Converting to a polygon + kml for future use
LNMop <- Polygon(LNMo)
LNMops <- Polygons(list(LNMop),1)
LNMosps = SpatialPolygons(list(LNMops))
plot(LNMosps) # Final middle low polygon 0 - 1000 TDS

#need to pull out shapes individually 
plot(L)
 
Ll <- L[-which(L[,1] > -111.4 & L[,1] < -108.7 
              & L[,2] > 35.2 ),] #Low Polygon 0 - 1000 TDS  

plot(Ll)
 
Llo <- concaveman(Ll, concavity = 1, length_threshold = 0.2) 
plot(Llo)
 
### Converting to a polygon 
Llop <- Polygon(Llo)
Llops <- Polygons(list(Llop),1)
Llosps = SpatialPolygons(list(Llops))
plot(Llosps) # Final low polygon 0 - 1000 TDS

#converting these polygons to a kml file

#trying to create list of polygons
zt1 <- Polygons(list(LNop, LHMop, LNMop, Llop),4)
zt1sp = SpatialPolygons(list(zt1))
plot(zt1sp)

proj4string(zt1sp) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

data <- data.frame(cbind(99.9,99.9,99.9,99.9))
row.names(data) <- '0-1k'
zt1spdf = SpatialPolygonsDataFrame(zt1sp, data, match.ID = F)
zt1spdf
plot(zt1spdf)

require(rgeos)
#Create the KML
writeOGR(zt1spdf, dsn ='CoCoAq_0-1.kml', layer='CCAq_0-1', driver="KML")

p.CA_0t1 <- ggplot() +
  geom_polygon(data = zt1sp, aes(x=long, y = lat, group = group))

p.CA_0t1 # may not plot but it is there in pdf preview

ggsave(plot=p.CA_0t1, "CoCoAq_0-1.tiff", device = "tiff")

rasterMap <- stack("CoCoAq_0-1.tiff")

# Get the GeoSpatial Components
lat_long <- ggplot_build(p.CA_0t1)$layout$panel_params[[1]][c("x.range","y.range")]

# Supply GeoSpatial  data to the StackedRaster
extent(rasterMap) <- c(lat_long$x.range,lat_long$y.range)
projection(rasterMap) <- CRS("+proj=longlat +datum=WGS84")

###write file
#writeOGR(zt1spdf, dsn="CoCoAq_0-1.kml", layer= "0-1kTDS", driver="KML")
