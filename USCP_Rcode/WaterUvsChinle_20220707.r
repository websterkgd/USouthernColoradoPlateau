# U vs Chinle

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
require(maps)
require(raster)

###Pulling in Map Data
usa <- map_data("usa")
us.xy <- cbind(usa$long[1:6850], usa$lat[1:6850])

library(sp)
p <- Polygon(us.xy, hole = T)
ps <- Polygons(list(p),1)
US <- SpatialPolygons(list(ps))
proj4string(US) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#Pulling in data
header <- read.csv('CredoEtal2019_238U.csv', nrow=1)
U <- read.csv('CredoEtal2019_238U.csv', skip = 2, col.names=header)

#Renaming the Uranium concentration
Uconc <- U$Overall.Reported.Average..ug.L. #average reported U conc. 
Uconc[which(Uconc == "BD")] <- 0.047  #detection limit = 0.047 ug/g
Uconc <- as.numeric(Uconc)

#pulling in the chinle formation
require(rgdal)    #retired by end of 2023 

lyr <- ogrListLayers("azgeol.KML")
mykml <- lapply(lyr, function(i) readOGR("azgeol.KML", i))
names(mykml) <- lyr

azkml <- mykml$`Geologic units of Arizona`

#pulling in utah
utkml <- readOGR("utgeol.KML", 
                "Geologic units of Utah",require_geomType="wkbPolygon")

#azkml@data[c(1,49),] #pulls out names and descriptions of chinle +shinarump
ChF <- azkml@polygons[[1]] #Chinle formation

SF <-  azkml@polygons[[49]] #Shinarump member

ChF_U <- utkml@polygons[[31]] #the chinle formation in Utah

#code to pull out single "labpt" numbers out of data 

#code to plot Chinle formation
c <- as.list(NA)
d <- as.list(NA)

for (i in 1:length(ChF@Polygons)){
  c[[i]] <- ChF@Polygons[[i]]@coords[,1]
  d[[i]] <- ChF@Polygons[[i]]@coords[,2]
}

CreatePolygon <- function(X,Y) {
  
  #Creates a polygon from data
  
  p <- geom_polygon(data= NULL, aes(x = X, y = Y), fill ="blue", alpha = 0.2)
  return(p)
}

mChF <- mapply(CreatePolygon, c, d) #polygons of Chinle in ggplot format

p.tr <- ggplot() + mChF  # plots
p.tr

###quick visual plot
p.UvChF <- p.tr + 
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  geom_point(aes(x = U[,4], y = U[,3]), 
             color = 'black', alpha = 0.5,
             pch = 16, 
             cex = sqrt((4/pi)*0.5*as.numeric(Uconc)))
p.UvChF # a few points- will try to look at U as is + slightly fuzz boundary

#Converting ChF to Spatial Polygons Data Frame

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
ChLspdf
plot(ChLspdf) #Both Utah and AZ

#creating spatial Y
Uc <- cbind(U[,c(4,3)],Uconc)
Uc <- Uc[complete.cases(Uc),]

# prepare coordinates, data, and proj4string
coords <- Uc[ , c(1,2)]   # coordinates
data   <- Uc[ , 3]          # data
crs    <- CRS("+proj=longlat +datum=WGS84 +no_defs") # proj4string of coords

# make the SpatialPoints object
uspat <- SpatialPoints(coords      = coords,
                       proj4string = crs)

#selecting Wells in Chinle
uinChF <- which(is.na(sp::over(uspat,ChLsp)) == FALSE) #returns u in 1 -3
uinChF <- as.data.frame(uinChF)

# prepare coordinates, data, and proj4string
finChF <- Uc[which(row.names(uspat@coords) %in% row.names(uinChF)),] # data

#want a better visual check
p.ChF <- ggplot() +
  geom_polygon(data=ChLspdf, aes(x=long, y = lat, group = group), fill = NA,
               color = "dark gray") +
  coord_fixed(1.3) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  geom_point(aes(x = finChF[,1], y = finChF[,2]),
             color = 'black', alpha = 0.5,
             pch = 16)
p.ChF # 

# #Quick function to add increase boundary
require(rgeos) #used for gUnion gUnion(polygons)  takes overlap of pgons
fuzzypolygon <- function(pgon, d) {
  pg.px <- pgon
  pg.mx <- pgon
  pg.py <- pgon
  pg.my <- pgon
  pg.pxpy <- pgon
  pg.mxpy <- pgon
  pg.mxmy <- pgon
  pg.pxmy <- pgon
  for(i in 1:length(pgon@polygons[[1]]@Polygons)) {
    pg.px@polygons[[1]]@Polygons[[i]]@coords[,1] <-
      pg.px@polygons[[1]]@Polygons[[i]]@coords[,1] + d

    pg.mx@polygons[[1]]@Polygons[[i]]@coords[,1] <-
      pg.mx@polygons[[1]]@Polygons[[i]]@coords[,1] - d

    pg.py@polygons[[1]]@Polygons[[i]]@coords[,2] <-
      pg.py@polygons[[1]]@Polygons[[i]]@coords[,2] + d

    pg.my@polygons[[1]]@Polygons[[i]]@coords[,2] <-
      pg.my@polygons[[1]]@Polygons[[i]]@coords[,2] - d

    pg.pxpy@polygons[[1]]@Polygons[[i]]@coords[,1] <-
      pg.pxpy@polygons[[1]]@Polygons[[i]]@coords[,1] + d
    pg.pxpy@polygons[[1]]@Polygons[[i]]@coords[,2] <-
      pg.pxpy@polygons[[1]]@Polygons[[i]]@coords[,2] + d

    pg.mxpy@polygons[[1]]@Polygons[[i]]@coords[,1] <-
      pg.mxpy@polygons[[1]]@Polygons[[i]]@coords[,1] - d
    pg.mxpy@polygons[[1]]@Polygons[[i]]@coords[,2] <-
      pg.mxpy@polygons[[1]]@Polygons[[i]]@coords[,2] + d

    pg.mxmy@polygons[[1]]@Polygons[[i]]@coords[,1] <-
      pg.mxmy@polygons[[1]]@Polygons[[i]]@coords[,1] - d
    pg.mxmy@polygons[[1]]@Polygons[[i]]@coords[,2] <-
      pg.mxmy@polygons[[1]]@Polygons[[i]]@coords[,2] - d

    pg.pxmy@polygons[[1]]@Polygons[[i]]@coords[,1] <-
      pg.pxmy@polygons[[1]]@Polygons[[i]]@coords[,1] + d
    pg.pxmy@polygons[[1]]@Polygons[[i]]@coords[,2] <-
      pg.pxmy@polygons[[1]]@Polygons[[i]]@coords[,2] - d
  }
  f.p <- gUnion(pg.px, pg.mx)
  f.p <- gUnion(f.p, pg.py)
  f.p <- gUnion(f.p, pg.my)
  f.p <- gUnion(f.p, pg.pxpy)
  f.p <- gUnion(f.p, pg.mxpy)
  f.p <- gUnion(f.p, pg.mxmy)
  f.p <- gUnion(f.p, pg.pxmy)
  return(f.p)
}
 
#simplifying outline of the chinle since it is complicated
#fixes generation of invalid polygons

Ch_simp <- gSimplify(ChLspdf, tol = 0.0001) #requires rgeos

#creating fuzzy Chinle
d <- 0.01
f.ChF <- fuzzypolygon(Ch_simp,d) # works
 
#selecting Wells in Chinle
uinfChF <- which(is.na(sp::over(uspat,f.ChF)) == FALSE) #returns u in 1 -3
uinfChF <- as.data.frame(uinfChF)

# prepare coordinates, data, and proj4string
finfChF <- Uc[which(row.names(uspat@coords) %in% row.names(uinfChF)),] # data

#want a better visual check
p.fChF <- ggplot() +
  geom_polygon(data=ChLspdf, aes(x=long, y = lat, group = group), fill = NA,
               color = "dark gray") +
  coord_fixed(1.3) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  geom_point(aes(x = finfChF[,1], y = finfChF[,2]),
             color = 'black', alpha = 0.5,
             pch = 16)
p.fChF # 

#basic stats
fluChF <- log10(finfChF[,3])
m.fChf <-  median(fluChF) #0.☺ #mid value
plot(fluChF) #easy beeswarm # some higher points
boxplot(fluChF)
