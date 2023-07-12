#ApU in Mines in Chinle in AZ. 

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
        
# pulling in the important data

#reading in mines + rad
ApUmine <- read.table('mines_kriged_U.txt')
colnames(ApUmine) <- Rmine[1,]
ApUmine <- ApUmine[-1,]

for (i in 1:length(ApUmine)){
  ApUmine[,i] <- as.numeric(ApUmine[,i])
}

#pulling in fuzzy polygon function
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

#pulling in the U rich rocks
require(rgdal)    #retired by end of 2023 

lyr <- ogrListLayers("azgeol.KML")
mykml <- lapply(lyr, function(i) readOGR("azgeol.KML", i))
names(mykml) <- lyr

azkml <- mykml$`Geologic units of Arizona`

#azkml@data[c(1,49),] #pulls out names and descriptions of chinle +shinarump
ChF <- azkml@polygons[[1]] #Chinle formation
MF <- azkml@polygons[[36]] #Morrison formation
SF <-  azkml@polygons[[49]] #Shinarump member

#converting these objects to polygons 
RockSpaPoly <- function(formation, ...) {
  pl <- list(formation@Polygons[[1]])
  for (i in 2:length(formation@Polygons)){
    pl <- append(pl, formation@Polygons[[i]])
    
    for(j in list(...)) {
      for (k in 1:length(j@Polygons)){
        pl <- append(pl, j@Polygons[[k]])
      }
    }
  }
  Fp <- Polygons(pl, length(pl))
  Fpsp = SpatialPolygons(list(Fp))
  proj4string(Fpsp) = 
  CRS("+proj=longlat +datum=WGS84 +no_defs")
  return(Fpsp)
}

ChFsp <- RockSpaPoly(ChF, SF) #Chinle Polygon

#Chinle spatial polygons data frame
data <- data.frame(rbind(rep(NA, length(ChFsp@polygons[[1]]@Polygons))))
row.names(data) <- 'Chinle'
ChFspdf = SpatialPolygonsDataFrame(ChFsp, data, match.ID = F)
ChFspdf
plot(ChFspdf) 

#creating fuzzy Chinle and Morrison
Ch_simp <- gSimplify(ChFspdf, tol = 0.0001) #requires rgeos

f.ChF <- fuzzypolygon(Ch_simp,0.01) #fuzzy Chinle

#pulling spatial rad coordinates
SR <- cbind(ApUmine[,c(1,2,3)])
SR <- SR[complete.cases(SR),]

# filtering SR data to be closer to the study region
SRf <- SR[which(SR[,2] > 34.8 & SR[,2] < 37.5 &
                  SR[,1] > -112 & SR[,1] < -109),]

# prepare coordinates, data, and proj4string
coords <- SRf[ , c(1,2)]   # coordinates
data   <- SRf[ , 3]          # data
crs    <- CRS("+proj=longlat +datum=WGS84 +no_defs") 
# proj4string of coords

# make the SpatialPoints object
ApUspat <- SpatialPoints(coords      = coords,
                       proj4string = crs)

#Creating TrueFalse vectors for analysis
MinChF <- sp::over(ApUspat,f.ChF) # runs
MinChF[is.na(MinChF)] = 0
MinChF <- as.numeric(MinChF) #Boolean Chinle 

#Visual check #
p.MC <- ggplot() +
  geom_polygon(data=f.ChF, aes(x=long, y = lat, group = group), fill = NA,
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
  geom_point(aes(x = SRf[which(MinChF > 0),1], y = SRf[which(MinChF > 0),2]),
             color = 'black', alpha = 0.5,
             pch = 16)
p.MC #

### Analyses of ApU concentrations at mines in chinle vs mines other regions

# is ApU higher at mines in the Chinle
# Is ApU higher in chinle
w.ApUMC <- wilcox.test(SRf[which(MinChF > 0),3], 
                      SRf[which(MinChF == 0),3]) 
w.ApUMC
#p = 3*10^-7, w = 59020 # mines in Chinle higher than mines outside of Chinle

boxplot(SRf[which(MinChF > 0),3], 
        SRf[which(MinChF == 0),3])

# ApU higher in Chinle 

# Is the radiation from mines in the Chinle higher than the chinle? 

#changing WD
#setwd()

# pulling in xyz files
SRR <- read.table("shiprock_rad.xyz", header= FALSE, sep = "") #reads file
FgR <- read.table("flagstaff_rad.xyz", header= FALSE, sep = "") #reads file
GpR <- read.table("gallup_rad.xyz", header= FALSE, sep = "") #reads file
MCR <- read.table("marble_canyon_rad.xyz", header= FALSE, sep = "") #reads file

#compiling info. 
RC <- rbind(FgR[,c(6,7,13)],SRR[,c(6,7,13)],MCR[,c(6,7,13)],GpR[,c(6,7,13)])

#removing larger datasets
rm(FgR,GpR,MCR,SRR)

#filtering the aerial rad data to the region closer to the study region
RC <- RC[which(RC[,3] != -9999.9),]

RC <- RC[-which(RC[,3] < 0),]

RCf <- RC[which(RC[,1] > 34.8 & RC[,1] < 37.5 &
                  RC[,2] > -112 & RC[,2] < -109),]

#which points in RCf are in the Chinle?
# prepare coordinates, data, and proj4string
coords <- RCf[ , c(2,1)]   # coordinates
data   <- RCf[ , 3]          # data
crs    <- CRS("+proj=longlat +datum=WGS84 +no_defs") 
# proj4string of coords

# make the SpatialPoints object
RCspat <- SpatialPoints(coords      = coords,
                       proj4string = crs)

#Creating TrueFalse vectors for analysis
ApUinChF <- sp::over(RCspat,f.ChF) # runs #takes time
ApUinChF[is.na(ApUinChF)] = 0
ApUinChF <- as.numeric(ApUinChF) #Boolean Chinle 

# ### aligning data for quick visualation
ChApUsp <- RCf[which(ApUinChF > 0),3]
ApUinCh <- SRf[which(MinChF > 0),3]
 
c(median(ChApUsp),median(ApUinCh)) #2.5 vs 3.2

boxplot(ChApUsp, ApUinCh) # Hard to tell - what's going on here

# Is Rad higher at mines in chinle compared to other areas in the chinle
w.RiC <- wilcox.test(ChApUsp, ApUinCh) #p = 0.13, w = 11392751
w.RiC # W = 7138435 p = 2.2*10^-16
# not sure how much I buy this due to differences in samples size

for (i in 1:10){
  is <- sample(ChApUsp,length(ApUinCh),replace=FALSE)
  print(i)
  print(c(median(is), median(ApUinCh))) 
  print(wilcox.test(is, ApUinCh))
} # All results are significant  # 2.5 vs 3.15
#Med in Chinle is higher than median outside Chinle

# testing ApU in chinle  vs ApU not in chinle (not related to mines)
#ChApUsp <- ApUinChf[which(ApUinCh > 0),3]
RApUsp <- RCf[which(ApUinChF == 0),3]

w.rad <- wilcox.test(ChApUsp, RApUsp) # Ap U higher in chinle
w.rad
#W = 1.2*10^10 p <2.2*10^-16
length(ChApUsp) #82959 
length(RApUsp) #267934
median(ChApUsp) #2.5
median(RApUsp) #2.1

#quick visualization 
boxplot(ChApUsp, RApUsp) 

# rad from rock in Chinle appears higher -> ever so slightly

for (i in 1:10){
  ChS <- sample(ChApUsp, 10000, replace=FALSE)
  ApUS <- sample(RApUsp, 10000, replace=FALSE)
  print(c(median(ChS),median(ApUS)))
  print(i)
  print(wilcox.test(ChS, ApUS))
} # ApU in chinle is higher - subtly


#Mines not in chinle vs rocks total Ap U env.
w.MnCvR <- wilcox.test(SRf[which(MinChF == 0),3], RCf[,3]) 
w.MnCvR
#different, p 2.2*10^-16, W = 82690684
median(SRf[which(MinChF == 0),3])
median(RCf[,3])
# median mine not in chinle = 2.7, median total env  = 2.2 
#length(RnnCh) #489
length(RCf[,3]) #350893

#Mines not in chinle vs rocks not in chinle
c(median(SRf[which(MinChF == 0),3]), median(RApUsp))
# median mine not in chinle = 2.7, median other = 2.1
w.MnCvRn <- wilcox.test(SRf[which(MinChF == 0),3], RApUsp)
w.MnCvRn
#different, p 2.2*10^-16, W = 65165954  
length(RApUsp) #267934

#mines in chinle vs rocks not in chinle
c(median(ApUinCh), median(RApUsp))
# median mine in chinle = 3.2, median other = 2.1
w.MiCvRn <- wilcox.test(ApUinCh, RApUsp)
w.MiCvRn  
#different, p 2.2*10^-16, W = 47060328
length(ApUinCh) #263
