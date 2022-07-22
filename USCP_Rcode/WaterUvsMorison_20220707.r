#U vs Morrison

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

#pulling in the morrison formation
require(rgdal)    #retired by end of 2023 

lyr <- ogrListLayers("azgeol.KML")
mykml <- lapply(lyr, function(i) readOGR("azgeol.KML", i))
names(mykml) <- lyr

azkml <- mykml$`Geologic units of Arizona`

MF <- azkml@polygons[[36]] #Morrison formation

#pulling in utah
utkml <- readOGR("utgeol.KML", 
                 "Geologic units of Utah",require_geomType="wkbPolygon")

MF_ut <- utkml@polygons[[123]] #the Morrison formation in Utah

#code to pull out single "labpt" numbers out of data 

#code to plot morrison formation
c <- as.list(NA)
d <- as.list(NA)

for (i in 1:length(MF@Polygons)){
  c[[i]] <- MF@Polygons[[i]]@coords[,1]
  d[[i]] <- MF@Polygons[[i]]@coords[,2]
}

CreatePolygon <- function(X,Y) {
  
  #Creates a polygon from data
  
  p <- geom_polygon(data= NULL, aes(x = X, y = Y), fill ="blue", alpha = 0.2)
  return(p)
}

mMF <- mapply(CreatePolygon, c, d) #polygons of Morrison in ggplot format

p.tr <- ggplot() + mMF  # plots
p.tr

###quick visual plot
p.UvMF <- p.tr + 
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
p.UvMF # visual check

#Converting MF to Spatial Polygons Data Frame

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
plot(M3spdf)

#Quick function to add increase boundary
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

#simplifying outline of the Morrison since it is complicated
#fixes generation of invalid polygons

M3_simp <- gSimplify(M3spdf, tol = 0.00001) #requires rgeos

#creating fuzzy Morrison
d <- 0.01
f.MF <- fuzzypolygon(M3_simp,d) #

#going to use f.MF - since it's a little broader

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

#selecting Wells in Morrison
uinMF <- which(is.na(sp::over(uspat,f.MF)) == FALSE) #returns u in 1 -3
uinMF <- as.data.frame(uinMF)

# prepare coordinates, data, and proj4string
finMF <- Uc[which(row.names(uspat@coords) %in% row.names(uinMF)),] # data

#quick visual check
p.fMF <- ggplot() + 
  geom_polygon(data=M3spdf, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = finMF[,1], y = finMF[,2]), 
             color = 'black', alpha = 0.5,
             pch = 16)
p.fMF # nice filter  - so method works. nice to see a visual check 
# no new points with utah added

#basic stats
fluMF <- log10(finMF[,3])
m.fMF <-  median(fluMF) #0.091 #pretty low - some points have higher values
plot(fluMF) #easy beeswarm 
boxplot(fluMF)
