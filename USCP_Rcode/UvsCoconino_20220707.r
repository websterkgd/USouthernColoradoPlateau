#Trying to Pull out the Coconino Aquifer

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/UraniumContamination/Data_WCM')

#Pulling in U data
header <- read.csv('CredoEtal2019_238U.csv', nrow=1)
U <- read.csv('CredoEtal2019_238U.csv', skip = 2, col.names=header)

#pulling in Coconino AQ
require(rgdal)
require(ggplot2)

lyr <- ogrListLayers("CoCoAqSurface.KML")
mykml <- lapply(lyr, function(i) readOGR("CoCoAqSurface.KML", i))
names(mykml) <- lyr

CCAq <- mykml$CCAq

plot(CCAq)
plot(U[,4],U[,3])

#Renaming the Uranium concentration
#average reported U conc.
Uconc <- as.numeric(U$Overall.Reported.Average..ug.L.) 

Uconc[which(Uconc == "BD")] <- 0.047  #detection limit = 0.047 ug/g

#plot of points + CocoAq Surface; quick visualization
p.UCA <- ggplot() + 
  geom_polygon(data=CCAq, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = U[,4], y = U[,3]), 
             color = 'black', alpha = 0.5,
             pch = 16, 
             cex = sqrt((4/pi)*0.5*Uconc))
p.UCA

##### Plotting points against other portions of the aquifer

#Quick function to pull in kml files of interest  

qimpt <- function(f) {
  lyr <- ogrListLayers(f)
  mykml <- lapply(lyr, function(i) readOGR(f, i))
  names(mykml) <- lyr
  return(mykml)
}

### O  - 1 k TDS

im_CC01 <- qimpt("CoCoAq_0-1.KML")

CCA01 <- im_CC01$CCAq_0_1
 
p.UC01 <- ggplot() + 
  geom_polygon(data=CCA01, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = U[,4], y = U[,3]), 
             color = 'black', alpha = 0.5,
             pch = 16, 
             cex = sqrt((4/pi)*0.5*Uconc))
p.UC01 #

### 1  - 1 3 TDS

im_CC13 <- qimpt("CoCoAq_1-3.KML")

CCA13 <- im_CC13$CCAq_1_3

p.UC13 <- ggplot() + 
  geom_polygon(data=CCA13, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = U[,4], y = U[,3]), 
             color = 'black', alpha = 0.5,
             pch = 16, 
             cex = sqrt((4/pi)*0.5*Uconc))
p.UC13

### 3  -  25 TDS

im_CC325 <- qimpt("CoCoAq_3-25.KML")

CCA325 <- im_CC325$CCAq_3_25

p.UC325 <- ggplot() + 
  geom_polygon(data=CCA325, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = U[,4], y = U[,3]), 
             color = 'black', alpha = 0.5,
             pch = 16, 
             cex = sqrt((4/pi)*0.5*Uconc))
p.UC325 # highest point are near the aquifer outflow

### > 25 TDS

im_CCg25 <- qimpt("CoCoAq_g25.KML")

CCAg25 <- im_CCg25$CCAq_25

p.UCg25 <- ggplot() + 
  geom_polygon(data=CCAg25, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = U[,4], y = U[,3]), 
             color = 'black', alpha = 0.5,
             pch = 16, 
             cex = sqrt((4/pi)*0.5*Uconc))
p.UCg25 # no points

#Looking at point in the coconino

#creating spatial u points

Uc <- cbind(U[,c(4,3)],Uconc) 
Uc <- Uc[complete.cases(Uc),]

# prepare coordinates, data, and proj4string
coords <- Uc[ , c(1,2)]   # coordinates
data   <- Uc[ , 3]          # data
crs    <- CRS("+proj=longlat +datum=WGS84 +no_defs") # proj4string of coords

# make the SpatialPoints object
uspat <- SpatialPoints(coords      = coords, 
                       proj4string = crs)

uin01 <- na.omit(sp::over(uspat,CCA01)) #returns u in 01

# prepare coordinates, data, and proj4string
fin01 <- Uc[which(row.names(uspat@coords) %in% row.names(uin01)),] # data
crs    <- CRS("+proj=longlat +datum=WGS84 +no_defs") # proj4string of coords

fuin01 <- SpatialPoints(coords = fin01[,1:2], proj4string = crs)

### Visual check of points in 01
p.f01 <- ggplot() + 
  geom_polygon(data=CCA01, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = fin01[,1], y = fin01[,2]), 
             color = 'black', alpha = 0.5,
             pch = 16)
p.f01 # visual check

#Pulling average U for 0-1
flu01 <- log10(fin01[,3])
m.f01 <-  median(flu01) #0.34

#Pulling Data for 1 - 3

#I need to add in a little error around edges since not quite exact

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

d <- 0.025 # setting expansion for polygons

fuz.CCA13 <- fuzzypolygon(CCA13, d) #change by 0.02 in lat + long

uin13 <- which(is.na(sp::over(uspat,fuz.CCA13)) == FALSE) #returns u in 1 -3
uin13 <- as.data.frame(uin13)

# prepare coordinates, data, and proj4string
fin13 <- Uc[which(row.names(uspat@coords) %in% row.names(uin13)),] # data

p.f13 <- ggplot() + 
  geom_polygon(data=CCA13, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = fin13[,1], y = fin13[,2]), 
             color = 'black', alpha = 0.5,
             pch = 16)
p.f13 # visual check

#Pulling average U for 0-1
flu13 <- log10(fin13[,3])
m.f13 <-  median(flu13) #0.64

#Pulling Data for 3 - 25
fuz.CCA325 <- fuzzypolygon(CCA325, d) #change by 0.02 in lat + long

uin325 <- which(is.na(sp::over(uspat,fuz.CCA325)) == FALSE) #returns u in 1 -3
uin325 <- as.data.frame(uin325)

# prepare coordinates, data, and proj4string
fin325 <- Uc[which(row.names(uspat@coords) %in% row.names(uin325)),] # data

p.f325 <- ggplot() + 
  geom_polygon(data=CCA325, aes(x=long, y = lat, group = group), fill = NA, 
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
  geom_point(aes(x = fin325[,1], y = fin325[,2]), 
             color = 'black', alpha = 0.5,
             pch = 16)
p.f325 # 


#Pulling average U for 3 to 24
flu325 <- log10(fin325[,3]) #9 pts
l.flu325 <- c(rep("3 - 25k TDS",length(flu325)))
m.f325 <-  median(flu325) #1.53

#doing a little plotting

#adding labels to vectors
l01 <- cbind(flu01, (rep('0-1k', length(flu01))))
l13 <- cbind(flu13, (rep('1-3k', length(flu13))))
l325 <- cbind(flu325, (rep('3-25k', length(flu325))))

##Stack Vectors
T.U <- as.data.frame(rbind(l01,l13,l325))

boxplot(as.numeric(as.character(T.U[,1]))~as.character(T.U[,2]),
        data=T.U) 

#ANOVA
AOO.aov <- aov(as.numeric(as.character(T.U[,1]))
               ~ as.character(T.U[,2]), data = T.U)
summary(AOO.aov) #significant 7*10^-6
TukeyHSD(AOO.aov) #3 to 25 is higher than other two 

### Surface Analysis of Coconino and U 

uinCsr <- which(is.na(sp::over(uspat,CCAq)) == FALSE) #returns u in surface
uinCsr <- as.data.frame(uinCsr)

finCsr <- Uc[which(row.names(uspat@coords) %in% row.names(uinCsr)),] # data

unnCsr <- which(is.na(sp::over(uspat,CCAq)) == TRUE) #returns u outside surface
unnCsr <- as.data.frame(unnCsr)

fnnCsr <- Uc[which(row.names(uspat@coords) %in% row.names(unnCsr)),] 

wilcox.test(log10(finCsr[,3]),log10(fnnCsr[,3]))

            