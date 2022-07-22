#multi regression of U data

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
require(plotKML) #require plotKML for readGPX
        
# pulling in the important data

#reading in well elevations
usel <- read.table('UspatElev.txt', sep = "\t")

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

#Pulling in data
header <- read.csv('CredoEtal2019_238U.csv', nrow=1)
U <- read.csv('CredoEtal2019_238U.csv', skip = 2, col.names=header)

#Renaming the Uranium concentration
Uconc <- U$Overall.Reported.Average..ug.L. #average reported U conc. 
Uconc[which(Uconc == "BD")] <- 0.047  #detection limit = 0.047 ug/g
Uconc <- as.numeric(Uconc)

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

#pulling in Utah
utkml <- readOGR("utgeol.KML", 
                 "Geologic units of Utah",require_geomType="wkbPolygon")
ChF_U <- utkml@polygons[[31]] #the chinle formation in Utah
MF_ut <- utkml@polygons[[123]] #the Morrison formation in Utah

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
  CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  return(Fpsp)
}

ChFsp <- RockSpaPoly(ChF, SF, ChF_U) #Chinle Polygon
MFsp <- RockSpaPoly(MF, MF_ut) #Morison polygon 

#Chinle spatial polygons data frame
data <- data.frame(rbind(rep(NA, length(ChFsp@polygons[[1]]@Polygons))))
row.names(data) <- 'Chinle'
ChFspdf = SpatialPolygonsDataFrame(ChFsp, data, match.ID = F)
ChFspdf
plot(ChFspdf) #plots - looks good

#morrison spatial polygons data frame
data <- data.frame(rbind(rep(NA, length(MFsp@polygons[[1]]@Polygons))))
row.names(data) <- 'Morrison'
MFspdf = SpatialPolygonsDataFrame(MFsp, data, match.ID = F)
MFspdf
plot(MFspdf) #plots

#creating fuzzy Chinle and Morrison
Ch_simp <- gSimplify(ChFspdf, tol = 0.0001) #requires rgeos
M3_simp <- gSimplify(MFspdf, tol = 0.00001) #requires rgeos

f.ChF <- fuzzypolygon(Ch_simp,0.01) #fuzzy Chinle
f.MF <- fuzzypolygon(M3_simp,0.01) #fuzzy Morrison


#pulling in uspat
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

#Creating TrueFalse vectors for analysis
uinChF <- sp::over(uspat,f.ChF) == 1 
uinChF[is.na(uinChF)] = 0
uinChF <- as.numeric(uinChF) #Boolean Chinle

uinMF <- sp::over(uspat,f.MF) == 1 
uinMF[is.na(uinMF)] = 0
uinMF <- as.numeric(uinMF) #Boolean morrison

#pulling in elevation of wells
uEl <- as.numeric(usel[2:length(usel[,4]),4]) #

#pulling in Aquifer regions

#function to import aquifer KMLs
qimpt <- function(f) {
  lyr <- ogrListLayers(f)
  mykml <- lapply(lyr, function(i) readOGR(f, i))
  names(mykml) <- lyr
  return(mykml)
}

### pulling in Aquifer regions
im_CC01 <- qimpt("CoCoAq_0-1.KML") # 0-1 k TDS
im_CC13 <- qimpt("CoCoAq_1-3.KML") # 1-3 k TDS
im_CC325 <- qimpt("CoCoAq_3-25.KML") # 3-25 k TDS

#pulling in particular regions
CCA01 <- im_CC01$CCAq_0_1
CCA13 <- im_CC13$CCAq_1_3
CCA325 <- im_CC325$CCAq_3_25

#pulling out wells in lowest region
#return boolean in 0 - 1
uin01 <- as.numeric(is.na(sp::over(uspat,CCA01)[,1]) == FALSE) 

#creating fuzzy CCA13 + CCA325
f.CCA13 <- fuzzypolygon(CCA13, 0.025) #change by 0.02 in lat + long
f.CCA325 <- fuzzypolygon(CCA325, 0.025) #change by 0.02 in lat + long

#pulling out wells in middle region
#return boolean in 0 - 1
uin13 <- as.numeric(is.na(sp::over(uspat,f.CCA13)) == FALSE)
uin325 <- as.numeric(is.na(sp::over(uspat,f.CCA325)) == FALSE)

#running regressions
lU <- log10(data)

#against elevation
RE <- lm(lU ~ uEl)
summary(RE)
plot(lU ~ uEl)
abline(RE,col="red") # significant against elevation 
# p = 5*10^-7 r2 = 0.11 U is higher in lower elevation wells 
#consider plotting this with elevation on the Y

#Adding in mine data
Uce <- cbind(Uc,as.numeric(usel[2:length(usel[,4]),4]))

mines <- qimpt("Mspat.KML") # mine locations

minesElev <- read.table("Mspatelev.txt", sep = "\t", fill = TRUE) 

f.mines <- minesElev[2:length(minesElev[,1]),2:4]
colnames(f.mines) <- c("lat","long","alt")

upstream <- function(df1, df2, d) {
  mlist <- rep(NA, length(df1[,1]))
  
  for (i in 1:length(df1[,1])) {
    #filter fmines > with elevations > greater than well of interest
    fdf <- df2[which(df2[,3] >  df1[i,4]),]
    
    #find mines within a distance from a well
    ffdf <- fdf[which(((df1[i,2] - as.numeric(fdf[,1]))^2 + 
                         (df1[i,1] - as.numeric(fdf[,2]))^2)^(0.5) < d),]
    
    #add number of mines upstream to well list
    mlist[i] <- length(ffdf[,1])
  }
  return(mlist)
}

#check function #run at 2 km 0.016 and 0.008
mu2 <- upstream(Uce, f.mines, 0.016)
mu1 <- upstream(Uce, f.mines, 0.008)

#### Multiple Regression Models

#with additive factors
AR <- lm(lU ~ uEl+uin01+uin13+uin325+uinChF+uinMF+mu2)
summary(AR) 
#relationships with El, 13, 325, ChF; mu2; r^2 0.24; adj r2 = 0.21 

require(car)
#### looking for colinearities
iv <- as.data.frame(cbind(lU, uEl,uin01,uin13,uin325,uinChF,uinMF,mu2))

vif(AR)

# running model removing the morrison formation
AR_mr <- lm(lU ~ uEl+uin01+uin13+uin325+uinChF+mu2)
summary(AR_mr) # All sig except u in01; r2 = 0.23, ajd r2 = 0.21

vif(AR_mr) 

# running model removing 01 + morrison
AR_mr_01 <- lm(lU ~ uEl+uin13+uin325+uinChF+mu2)
summary(AR_mr_01) # All sig; r2 = 0.23, ajd r2 = 0.21

vif(AR_mr_01) # this is the best "model"
#the other model is important says whats not sig

# running model removing 01 + 13 + morrison formation
AR_mr_01_13 <- lm(lU ~ uEl+uin325+uinChF+mu2)
summary(AR_mr_01_13) # All sig; r2 = 0.21, ajd r2 = 0.20

vif(AR_mr_01_13) 

# running model removing 01 + 13 + morrison + Ch formation
AR_mr_01_13_ch <- lm(lU ~ uEl+uin325+mu2)
summary(AR_mr_01_13_ch) # All sig except u in01; r2 = 0.20, ajd r2 = 0.19

vif(AR_mr_01_13_ch) 

#multiplicative model
MR <- lm(lU ~ uEl*uin01*uin13*uin325*uinChF*uinMF*mu2)
summary(MR)
#mu2, 01*mu2, only sig factors,
#total model is very sig 
#r^2 = 0.44; adj r2 = 0.34 too much power

vif(MR, type = 'predictor') 
#the error in this model is the colinearity and interactions among data - 

#### Making nice figures

### Fig 4a  - Coconino Box plot
fin01 <- iv[which(iv[,3] == 1),] # data
fin13 <- iv[which(iv[,4] == 1),] # data
fin325 <- iv[which(iv[,5] == 1),] # data

flu01 <- fin01[,1]
flu13 <- fin13[,1]
flu325 <- fin325[,1]

l01 <- cbind(flu01, (rep('0-1k', length(flu01))))
l13 <- cbind(flu13, (rep('1-3k', length(flu13))))
l325 <- cbind(flu325, (rep('3-25k', length(flu325))))

##Stack Vectors
T.U <- as.data.frame(rbind(l01,l13,l325))

boxplot(as.numeric(as.character(T.U[,1]))~as.character(T.U[,2]),
        data=T.U) 

p.4a <- ggplot() +
  geom_boxplot(aes(T.U[,2], as.numeric(T.U[,1])),outlier.shape = NA)+
  geom_jitter(aes(T.U[,2], as.numeric(T.U[,1])),width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
               labels = c(expression(0.01,0.1,1,10,100,1000)),
                          limits = c(-2, 3))+
  labs(y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), 
       x="TDS in the Coconino-De Chelly Aquifer (mg/L)")+
  annotate(geom="text", x="3-25k", y=2.9, label="*",
          color="red", size=7)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  theme(panel.background = element_blank()) 

plot(p.4a) # Plotting 

##### Elevation vs U  = plot 4b
p.4b <- ggplot()+
  geom_point(data = iv, aes(x=iv[,1], y=iv[,2]),
             color = "black", alpha = 0.5, pch = 16, size = 3)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  geom_segment(
    aes(x = -0.585, y = 2500, xend = 1.32, yend = 1000),
    color = "blue", size = 1) +
  scale_x_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  annotate(geom="text", x=1.8, y=2500, 
           label=expression(paste('[U] = 10^(-0.0013*'*italic(z)~'+ 2.59)')),
           color="black", size=3)+
  annotate(geom="text", x=1.8, y=2400, 
           label=expression(paste(r^{2},'= 0.11')),
           color="black", size=3)+
  labs(y="Elevation (m)", 
       x=expression(paste('[U] ('*italic("\u00B5")*'g/L)')))+
  theme(panel.background = element_blank())

plot(p.4b)

#### Panel 4c  U vs Chinle
InCh <- iv[which(iv[,6] == 1),1] # data
NnCh <- iv[which(iv[,7] == 0 & iv[,6] == 0),1] # data

nInCh <- cbind(InCh, (rep('In Chinle Fm.', length(InCh))))
nNnCh <- cbind(NnCh, (rep('Not in Chinle Fm.', length(NnCh))))

Ch.U <- as.data.frame(rbind(nInCh,nNnCh))

wCh <- wilcox.test(NnCh, InCh) #W = 2521, p = 0.002

p.4c <- ggplot() +
  geom_boxplot(aes(Ch.U[,2], as.numeric(Ch.U[,1])),outlier.shape = NA)+
  geom_jitter(aes(Ch.U[,2], as.numeric(Ch.U[,1])),width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  labs(y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), x = "")+
  annotate(geom="text", x="In Chinle Fm.", y=2.9, label="*",
           color="red", size=7)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  theme(panel.background = element_blank()) 

plot(p.4c) # Plotting 

#### Panel 4d  U vs Morrison
InMF <- iv[which(iv[,7] == 1),1] # data
NnMF <- iv[which(iv[,7] == 0 & iv[,6] == 0),1] # data

nInMF <- cbind(InMF, (rep('In Morrison Fm.', length(InMF))))
nNnMF <- cbind(NnMF, (rep('Not in Morrison Fm.', length(NnMF))))

MF.U <- as.data.frame(rbind(nInMF,nNnMF))

wMF <- wilcox.test(NnMF, InMF) #W = 1032, p-value = 0.6857

p.4d <- ggplot() +
  geom_boxplot(aes(MF.U[,2], as.numeric(MF.U[,1])),outlier.shape = NA)+
  geom_jitter(aes(MF.U[,2], as.numeric(MF.U[,1])),width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  labs(y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), x = "")+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  theme(panel.background = element_blank()) 

plot(p.4d) # Plotting 

require("cowplot")

#export as 8 x 8
plot_grid(p.4a, p.4b, p.4c, p.4d,
          labels = c("A","B","C","D"),
          ncol = 2)

#### Panel 5  U vs Mine
InM2 <- iv[which(iv[,8] > 0),1] # data
NnM2 <- iv[which(iv[,8] == 0),1] # data

nInM2 <- cbind(InM2, (rep('Mine within 2 km', length(InM2))))
nNnM2 <- cbind(NnM2, (rep('Mine not within 2 km', length(NnM2))))

M2.U <- as.data.frame(rbind(nInM2,nNnM2))

p.5 <- ggplot() +
  geom_boxplot(aes(M2.U[,2], as.numeric(M2.U[,1])),outlier.shape = NA)+
  geom_jitter(aes(M2.U[,2], as.numeric(M2.U[,1])),width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  annotate(geom="text", x="Mine within 2 km", y=2.9, label="*",
           color="red", size=7)+
  labs(y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), x = "")+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  theme(panel.background = element_blank()) 

plot(p.5) # Plotting export as 5 X 5


#### Output a concise csv of the data of interest here - 

# doi <- as.data.frame(
#   cbind(lU, Uc[,1:2], uEl,uin01,uin13,uin325,uinChF,uinMF,
#         mu1, mu2))
# 
# write.csv(doi, "UCaus.csv")
