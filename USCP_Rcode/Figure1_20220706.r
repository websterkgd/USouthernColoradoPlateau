#Code to Generate Figure1

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/UraniumContamination/Data_WCM')

#loading some packages
require(vegan)
require(ggplot2)
require(devtools)
require(dplyr)
require(stringr)
require(mapdata)
require(ggmap)
require(gstat)
#require(tmap)
require(maps)
require(raster)
require(rgdal)

###Pulling in Map Data
usa <- map_data("usa")
us.xy <- cbind(usa$long[1:6850], usa$lat[1:6850])

library(sp)
p <- Polygon(us.xy, hole = T)
ps <- Polygons(list(p),1)
US <- SpatialPolygons(list(ps))
proj4string(US) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#function to import KMLs #requires rgdal
qimpt <- function(f) {
  lyr <- ogrListLayers(f)
  mykml <- lapply(lyr, function(i) readOGR(f, i))
  names(mykml) <- lyr
  return(mykml)
}

### pulling in Colorado Plateau region
CP <- qimpt("ColoradoPlateau.KML") #
CPf <- as.data.frame(cbind(
  CP$`Colorado Plateau`@polygons[[1]]@Polygons[[1]]@coords, 
  rep(1, length(CP$`Colorado Plateau`@polygons[[1]]@Polygons[[1]]@coords[,1]))))
colnames(CPf) <- c("long","lat","group") # ggplot friendly

#pulling in NN
NN <- qimpt("NavajoNation.KML")

NNl <- list(NN$NavajoNation.kmz@polygons[[1]]@Polygons[[1]])

for (i in 2:length(NN$NavajoNation.kmz@polygons[[1]]@Polygons)){
  NNl <- append(NNl, NN$NavajoNation.kmz@polygons[[1]]@Polygons[[i]])
}

NNlp <- Polygons(NNl, length(NNl))

NNlsp = SpatialPolygons(list(NNlp))

proj4string(NNlsp) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

data <- data.frame(rbind(rep(99.9, length(NNl))))
row.names(data) <- 'NN'
NNlspdf = SpatialPolygonsDataFrame(NNlsp, data, match.ID = F)

require(jpeg) #read_jpg
AZ_CP <- readJPEG('F1_PanelAwoMarkings.JPG', native = TRUE) 
# boundaries = -114.149388, -106.356675, 32.998294, 41.152667

#Creating map of the area
E_usa <- map_data("state", 
                  region = c("Arizona", "Colorado", "New Mexico", "Utah"))

#filtering E_usa to region of imagery

# E_usaf <- E_usa[which(E_usa[,1] < -106.356675 & E_usa[,1] > -114.149388
#                       & E_usa[,2] > 32.9982945 & E_usa[,2] < 41.152667),]
# #works but I know have false boundaries  #using lines for states instead

#creating inset
iset <- as.data.frame(rbind(c(-112, 37.5, 1),
                            c(-109, 37.5, 1),
                            c(-109, 34.8, 1),
                            c(-112, 34.8, 1)))
colnames(iset) <- c('long', 'lat', 'group')

require(grid) # custom annotation - plots background #for rasterGrob

#plot of panel A # looks nice by itself as a 5 X 5, issue with up arrow in pdf
# Can use a line segment for the up arrow
p.Pa <- ggplot() + 
  annotation_custom(rasterGrob(AZ_CP,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -114.149388, -106.356675, 32.998294, 41.152667) +
  geom_segment(
    aes(x = -114, y = 37, xend = -106.356675, yend = 37), 
    color = "light gray") + # N-S boundary Colo-NM
  geom_segment(
    aes(x = -111.046406, y = 41, xend = -106.356675, yend = 41), 
    color = "light gray") + # Northern most boundary
  geom_segment(
    aes(x = -114.047388, y = 41.152667, xend = -114.047388, yend = 36.177605), 
    color = "light gray") + #furthest W boundary
  geom_segment(
    aes(x = -109.043043, y = 41, xend = -109.043043, yend = 32.998294), 
    color = "light gray") +
  geom_segment(
    aes(x = -111.046406, y = 41, xend = -111.046406, yend = 41.152667), 
    color = "dark gray") +
  geom_segment(
    aes(x = -114.149388, y = 36.177605, xend = -114, yend = 36.115944), 
    color = "dark gray") + #squiggly boundary AZ
  geom_polygon(data=NNlspdf, aes(x=long, y = lat, group = group),
               fill = "blue", color = "blue", alpha = 0.3) +
  geom_polygon(data = CPf, aes(x=long, y = lat, group = group), 
               fill = alpha("gray", 0.4), color = "dark gray")+
  annotate(geom="text", x=-112, y=39, label="Utah",
           color="white", size=5)+
  annotate(geom="text", x=-107.5, y=39, label="Colorado",
           color="white", size=5)+
  annotate(geom="text", x=-112.5, y=34.5, label="Arizona",
           color="white", size=5)+
  annotate(geom="text", x=-107.5, y=36, label="New Mexico",
           color="white", size=5)+
  geom_segment(
    aes(x = -113.5, y = 33.7, xend = -113.5, yend = 34.2), 
    color = "white", size = 1, 
    arrow = arrow(length = unit(0.3, "cm"), type = "open", ends = "last"))+
  annotate(geom="text", x=-113.5, y=33.5, label="N",
           color="white", size=7)+
  geom_segment(
    aes(x = -112, y = 33.8, xend = -110.341, yend = 33.8), 
    color = "white", size = 2) +
  annotate(geom="text", x=-112, y=33.5, label="0",
           color="white", size=6)+
  annotate(geom="text", x=-110.341, y=33.5, label="100 km",
           color="white", size=6)+
  annotate(geom="text", x=-111.8, y=37.2, label="B",
           color="white", size=7)+
  annotate(geom="text", x=-113.9, y=40.8, label="A",
           color="white", size=7)+ 
  geom_polygon(data = iset, aes(x=long, y = lat, group = group), 
               fill =alpha("gray", 0), color = "black")+
  #coord_fixed(1.3) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
       axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plot(p.Pa)#  
#some errors plottingNavajo Nation so may not keep it
#Need to plot boundary of Utah to top of figure
#need to plot boundary of AZ on west side

#Loading some of the B panel Reqs
IS_CP <- readJPEG('F1_PanelBwoMarkings.JPG', native = TRUE) 
# boundaries = -109, 37.5, -112, 34.8

mykml <- qimpt("azgeol.KML") # brining in AZ geological units

azkml <- mykml$`Geologic units of Arizona`

#pulling in utah
utkml <- readOGR("utgeol.KML", 
                 "Geologic units of Utah",require_geomType="wkbPolygon")

ChF <-  azkml@polygons[[1]] #Chinle Formation

SF <-  azkml@polygons[[49]] #Shinarump member

ChF_U <- utkml@polygons[[31]] #the chinle formation in Utah

#converting these to a dataframe

pl <- list(ChF@Polygons[[1]])

for (i in 2:length(ChF@Polygons)){
  pl <- append(pl, ChF@Polygons[[i]])
}

for (i in 1:length(ChF_U@Polygons)){
  pl <- append(pl, ChF_U@Polygons[[i]])
}

for (i in 1:length(SF@Polygons)){
  pl <- append(pl, SF@Polygons[[i]])
}

ChL <- Polygons(pl, length(pl))

ChLsp = SpatialPolygons(list(ChL))

proj4string(ChLsp) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

data <- data.frame(rbind(rep(99.9, length(pl))))
row.names(data) <- 'Chinle'
ChLspdf = SpatialPolygonsDataFrame(ChLsp, data, match.ID = F)

ChLf <- ChLspdf

#for morrison
MF <- azkml@polygons[[36]] #Morrison formation

utkml <- readOGR("utgeol.KML", 
                 "Geologic units of Utah",require_geomType="wkbPolygon")

MF_ut <- utkml@polygons[[123]] #the Morrison formation in Utah

M3 <- list(MF@Polygons[[1]])

for (i in 2:length(MF@Polygons)){
  M3 <- append(M3, MF@Polygons[[i]])
}

for (i in 1:length(MF_ut@Polygons)){
  M3 <- append(M3, MF_ut@Polygons[[i]])
}

ML <- Polygons(M3, length(M3))

M3sp = SpatialPolygons(list(ML))

proj4string(M3sp) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

data <- data.frame(rbind(rep(99.9, length(M3))))
row.names(data) <- 'Morrison'
M3spdf = SpatialPolygonsDataFrame(M3sp, data, match.ID = F)
M3spdf

CS <- qimpt("CoCoAqSurface.KML") # brining in AZ geological units

CCAa <- CS$CCAq

######## Panel B
p.Pb <- ggplot() +
  annotation_custom(rasterGrob(IS_CP,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -112, -109, 34.8, 37.5)+
  annotate(geom="text", x=-111.9, y=37.4, label="B",
           color="white", size=7)+
  geom_polygon(data = iset, aes(x=long, y = lat, group = group), 
               fill =alpha("gray", 0), color = alpha("gray", 0))+ #facilitates plotting
  geom_segment(
    aes(x = -109.043043, y = 34.8, xend = -109.043043, yend = 37.5), 
    color = "light gray") +
  geom_segment(
    aes(x = -112, y = 37, xend = -109, yend = 37), 
    color = "light gray") +
  geom_polygon(data=ChLspdf, aes(x=long, y = lat, group = group),
               fill = "blue", color = "blue") +
  geom_polygon(data=M3spdf, aes(x=long, y = lat, group = group),
               fill = "red", color = "red") +
  geom_polygon(data=CCAa, aes(x=long, y = lat, group = group),
               fill = alpha("gray", 0.5), color = "gray") +
  geom_segment(
    aes(x = -111.810356, y = 35.2, xend = -111.254868, yend = 35.2), 
    color = "white", size = 2) +
  annotate(geom="text", x=-111.81, y=35.1, label="0",
           color="white", size=6)+
  annotate(geom="text", x=-111.20, y=35.1, label="50 km",
           color="white", size=6)+
  geom_segment(
    aes(x = -111.7, y = 35.5, xend = -111.7, yend = 35.7), 
    color = "white", size = 1, 
    arrow = arrow(length = unit(0.3, "cm"), type = "open", ends = "last"))+
  annotate(geom="text", x=-111.7, y=35.4, label="N",
           color="white", size=7)+
  coord_cartesian(ylim = c(34.95, 37.35), xlim =c(-111.85,-109.15), clip = "on")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) # 

plot(p.Pb) # Plotting 

require(cowplot) #for multipanel

plot_grid(p.Pa, p.Pb, ncol = 2)
#export as 5 X 10
