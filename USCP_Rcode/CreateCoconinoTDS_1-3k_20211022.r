#Pulling out the Coconino Aquifer for 1 - 3 TDS

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

############## 

#pulling out 1000 - 3000 
# with 15 colors = color 9

upper.m.C <- c(0.7479278, 0.7115530, 0.8726493)+0.028

#pulling out 1000 - 3000
lower.m.C <- c(0.7479278, 0.7115530, 0.8726493)-0.015 

m.C <- countcolors::rectangularRange(DSC2, 
                                      upper = as.numeric(upper.m.C), 
                                      lower = as.numeric(lower.m.C),
                                      target.color = "green",
                                      plotting = TRUE)
#converting area to points
m.C.I <- m.C$pixel.idx

#plotting pixels in a somewhat correct orientation
plot(m.C.I[,2],- m.C.I[,1])

#converting to lat and long
# the pixels stay fixed

np <- 39.524254   # -1N = 39.524254,
sp <- 33.742784 # -306 S = 33.742784
nssp <-(sp-np)/(-305) # NS line:  

wp <- -111.955907 # 1W = -111.955907
ep <- -107.808380 # 180E = -107.808380
ewsp <- (ep-wp)/(179) # EW line

ccW <- ewsp*range(m.C.I[,2])[1] + (wp-ewsp) #westmost point 0 - 1000 CC Aq
ccE <- ewsp*range(m.C.I[,2])[2] + (wp-ewsp) #eastmost point 0 - 1000 CC Aq
ccN <- nssp*-range(m.C.I[,1])[1] + (np-nssp) #northtmost point 0 - 1000 CC Aq
ccS <- nssp*-range(m.C.I[,1])[2] + (np-nssp)  #southmost point 0 - 1000 CC Aq

bw <- max(m.C.I[,2])-min(m.C.I[,2])
bn <- max(-m.C.I[,1])-min(-m.C.I[,1])

dE <- m.C.I[,2]*(ccE-ccW)/bw + ccW - min(m.C.I[,2])*(ccE-ccW)/bw

dN <- -m.C.I[,1]*(ccN-ccS)/bn + ccS - min(-m.C.I[,1])*(ccN-ccS)/bn

plot(dE, dN)

SL = cbind(dE, dN) #SL Second Lowest

#need to pull out shapes individually
require(concaveman)
plot(SL) 

t <- SL[c(which(SL[,2] > 38.4)),] 
 
plot(t)
 
SLh <- t[c(which(t[,1] > -111.0 & t[,1] < -110.2)),]
 
plot(SLh)

SLhh <- SLh[which(SLh[,2] < 38.9),]

plot(SLhh)
 
SLhho <- concaveman(SLhh, concavity = 0.0, length_threshold = 0.1) 
plot(SLhho)

### Converting to a polygon + kml for future use
SLhop <- Polygon(SLhho)
SLhops <- Polygons(list(SLhop),1)
SLhosps = SpatialPolygons(list(SLhops))
plot(SLhosps) # Upper polygon 1000 - 3000 TDS

#Second shape 1000 - 3000 TDS range
plot(SL)

t <- SL[c(which(SL[,2] > 36.1 & SL[,2] < 38.7)),] #Middle high Low
plot(t)

SLt <- t[-c(which(t[,1] < -110.8)),]
plot(SLt)

SLth <- SLt[-c(which(SLt[,2] > 38.4 & SLt[,1] < -110.2)),]
plot(SLth)

SLtho <- concaveman(SLth, concavity = 1, length_threshold = 0.4) 
plot(SLtho)

### Converting to a polygon + kml for future use 
SLhhoop <- Polygon(SLtho)
SLhhoops <- Polygons(list(SLhhoop),1)
SLhhoosps = SpatialPolygons(list(SLhhoops))
plot(SLhhoosps) #Final Middle Polygon 1000 - 3000 TDS #Close enough

# Southernmost polygon 1000 - 3000 TDS 
plot(SL)

t <- SL[c(which(SL[,2] < 35.9)),] 

plot(t)

to <- concaveman(t, concavity = .2, length_threshold = 0.3) 
plot(to)

###Pretty big errors in plotting here

tr <- to[c(5:90),]
plot(tr)
tt <- to[c(151:177),]
plot(rbind(tr,tt))
ts <- to[c(1:2),]
plot(rbind(tr,tt,ts))
tv <- to[c(143:148),]
plot(rbind(tr,tt,ts,tv))
tq <- to[c(100:140),] 
plot(rbind(tr,tt,ts,tv,tq))

tou <-rbind(tr,tt,ts,tv,tq)

# ### Converting to a polygon + kml for future use
toup <- Polygon(tou)
toups <- Polygons(list(toup),1)
tousps = SpatialPolygons(list(toups))
plot(tousps) ### Final polygon Southernmost 1000 - 3000 TDS #close enough

###### Creating list of polygons
ot3 <- Polygons(list(SLhop, SLhhoop, toup),3)
ot3sp = SpatialPolygons(list(ot3))
plot(ot3sp)

proj4string(ot3sp) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

data <- data.frame(cbind(99.9,99.9,99.9))
row.names(data) <- '1-3k'
ot3spdf = SpatialPolygonsDataFrame(ot3sp, data, match.ID = F)
ot3spdf
plot(ot3spdf)

require(rgeos)
#Create the KML
writeOGR(ot3spdf, dsn ='CoCoAq_1-3.kml', layer='CCAq_1-3', driver="KML")

p.CA_1t3 <- ggplot() +
  geom_polygon(data = ot3sp, aes(x=long, y = lat, group = group))

p.CA_1t3 # may not plot but it is there in pdf preview

ggsave(plot=p.CA_1t3, "CoCoAq_1-3.tiff", device = "tiff")

rasterMap <- stack("CoCoAq_1-3.tiff")

# Get the GeoSpatial Components
lat_long <- ggplot_build(p.CA_1t3)$layout$panel_params[[1]][c("x.range","y.range")]

# Supply GeoSpatial  data to the StackedRaster
extent(rasterMap) <- c(lat_long$x.range,lat_long$y.range)
projection(rasterMap) <- CRS("+proj=longlat +datum=WGS84")

#writeOGR(ot3spdf, dsn ='CoCoAq_1-3.kml', layer='CCAq_1-3', driver="KML")
