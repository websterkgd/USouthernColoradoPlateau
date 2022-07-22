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

#reading in elevation of Wells
usel <- read.table('UspatElev.txt', sep = "\t")

#Pulling in data
header <- read.csv('CredoEtal2019_238U.csv', nrow=1)
U <- read.csv('CredoEtal2019_238U.csv', skip = 2, col.names=header)

#Renaming the Uranium concentration
Uconc <- U$Overall.Reported.Average..ug.L. #average reported U conc. 
Uconc[which(Uconc == "BD")] <- 0.047  #detection limit = 0.047 ug/g
Uconc <- as.numeric(Uconc)

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

#read in kml of mines

#function to import KMLs

require(rgdal)

qimpt <- function(f) {
  lyr <- ogrListLayers(f)
  mykml <- lapply(lyr, function(i) readOGR(f, i))
  names(mykml) <- lyr
  return(mykml)
}

### pulling in Mines
mines <- qimpt("Mspat.KML") # mine locations

minesElev <- read.table("Mspatelev.txt", sep = "\t", fill = TRUE) 

f.mines <- minesElev[2:length(minesElev[,1]),2:4]
colnames(f.mines) <- c("lat","long","alt")

# Now which wells are downstream from a mine

Uce <- cbind(Uc,as.numeric(usel[2:length(usel[,4]),4]))

nmuw <- rep(NA, length(Uce[,1]))

for (i in 1:length(Uce[,1])) {
  #filter fmines > with elevations > greater than well of interest
  fdf <- f.mines[which(f.mines[,3] >  Uce[i,4]),]
  
  #find mines within 0.040 deg (5 km) from well
  dmw <- 0.040
  ffdf <- fdf[which(((Uce[i,2] - as.numeric(fdf[,1]))^2 + 
                       (Uce[i,1] - as.numeric(fdf[,2]))^2)^(0.5) < dmw),]
  
  #add number of mines upstream to well list
  nmuw[i] <- length(ffdf[,1])
}

plot(nmuw, log10(Uce[,3]))

# I would like to write a function to do this because I want to
# look at more distances

#function is specific to the dataframes used here,
#could be made more general

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

#check function #run at 0.04 degrees (5 km)
mu5 <- upstream(Uce, f.mines, 0.04) # function is correct
mu4 <- upstream(Uce, f.mines, 0.032)
mu3 <- upstream(Uce, f.mines, 0.024)
mu2 <- upstream(Uce, f.mines, 0.016)
mu1 <- upstream(Uce, f.mines, 0.008)

#### exploratory linear models

#1 km
nmu1 <- length(which(mu1 > 0)) #13
plot(log10(Uce[,3]) ~ mu1)
summary(lm(log10(Uce[,3]) ~ mu1)) # p = 0.003 max = 2 r2 = 0.04

#2 km
nmu2 <- length(which(mu2 > 0)) #30
plot(log10(Uce[,3]) ~ mu2)
summary(lm(log10(Uce[,3]) ~ mu2)) # p = 0.0004 max = 6 r2 = 0.05

#3 km
nmu3 <- length(which(mu3 > 0)) #41
plot(log10(Uce[,3]) ~ mu3)
summary(lm(log10(Uce[,3]) ~ mu3)) # p = 0.0002 max = 10

#4 km
nmu4 <- length(which(mu4 > 0)) # 50
plot(log10(Uce[,3]) ~ mu4)
summary(lm(log10(Uce[,3]) ~ mu4)) # p = 0.0003 max = 20

#5 km
nmu5 <- length(which(mu5 > 0)) # 56
plot(log10(Uce[,3]) ~ mu5)
summary(lm(log10(Uce[,3]) ~ mu5)) # p = 0.0001 max ~ 22 r2 = 0.067

#examining at 2 km
ef <- cbind(log10(Uce[,3]), mu2)
efM <- ef[which(ef[,2] > 0), 1]
efnM <- ef[which(ef[,2] == 0), 1]

mtest <- wilcox.test(efM, efnM)
mtest
# p = 5*10^-5 W = 3978 U is higher downstream from a mine

# 1 km to check
ef_1 <- cbind(log10(Uce[,3]), mu1)
efM_1 <- ef[which(ef_1[,2] > 0), 1]
efnM_1 <- ef[which(ef_1[,2] == 0), 1]

mtest_1 <- wilcox.test(efM_1, efnM_1)
mtest_1
# p = 0.003 W = 1920 U is higher downstream from a mine

## 3 km
ef_3 <- cbind(log10(Uce[,3]), mu3)
efM_3 <- ef[which(ef_3[,2] > 0), 1]
efnM_3 <- ef[which(ef_3[,2] == 0), 1]

mtest_3 <- wilcox.test(efM_3, efnM_3)
mtest_3 # W = 5180 p = 1*10^-6

## 4 km
ef_4 <- cbind(log10(Uce[,3]), mu4)
efM_4 <- ef[which(ef_4[,2] > 0), 1]
efnM_4 <- ef[which(ef_4[,2] == 0), 1]

mtest_4 <- wilcox.test(efM_4, efnM_4)
mtest_4 # W = 5180 p = 1*10^-6

## 5 km
ef_5 <- cbind(log10(Uce[,3]), mu5)
efM_5 <- ef[which(ef_5[,2] > 0), 1]
efnM_5 <- ef[which(ef_5[,2] == 0), 1]

mtest_5 <- wilcox.test(efM_5, efnM_5)
mtest_5 # W = 6221 p = 2*10^-6
