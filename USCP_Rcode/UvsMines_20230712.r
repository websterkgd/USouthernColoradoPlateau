#multi regression of U data

rm(list=ls())
getwd()
#setwd('')

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
        
# pulling in data
Upred <- read.csv("Upred_matrix.csv")

##### Creating spatial Uranium concentrations object

#pulling in uspat
#creating spatial Y
Uc <- cbind(Upred[,c(4,3)],Upred[,2])
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

# qimpt <- function(f) {
#   lyr <- ogrListLayers(f)
#   mykml <- lapply(lyr, function(i) readOGR(f, i))
#   names(mykml) <- lyr
#   return(mykml)
# }
# 
# ### pulling in Mines
# mines <- qimpt("Mspat.KML") # mine locations   

## turning off above code

minesElev <- read.table("Mspatelev.txt", sep = "\t", fill = TRUE) 

f.mines <- minesElev[2:length(minesElev[,1]),2:4]
colnames(f.mines) <- c("lat","long","alt")

#want to ask if U is high in vicinity of mine 

#any direction, then downstream, then upstream, then Chinle vs nn Chinle

#Need distance from a water source to a mine

dm <- rep(NA, length(Uc[,1]))
for (i in 1:length(Uc[,1])){
  print(i)
  dwm <- ((((as.numeric(f.mines[,2])-as.numeric(Uc[i,2]))*90)^2) + 
            (((as.numeric(f.mines[,1])-as.numeric(Uc[i,1]))*111)^2))^(1/2)
  dm[i] <- sort(dwm)[1]
}
# 111 km per degree lat; 90 km per degree lon

#quick plotting
plot(dm, Upred$lU)
abline(median(Upred$lU), 0)
abline(lm(Upred$lU~dm)$coefficients[1],lm(Upred$lU~dm)$coefficients[2], 
       col= "blue")

summary(lm(Upred$lU~dm)) # r^2 0.14; p = 2.2*10^-8; m = -0.03

# removing lowest lU
plot(dm[-which(Upred$lU == min(Upred$lU))],
     Upred$lU[-which(Upred$lU == min(Upred$lU))])
abline(median(Upred$lU[-which(Upred$lU == min(Upred$lU))]), 0)
abline(lm(Upred$lU[-which(Upred$lU == min(Upred$lU))]~
            dm[-which(Upred$lU == min(Upred$lU))])$coefficients[1],
       lm(Upred$lU[-which(Upred$lU == min(Upred$lU))]~
           dm[-which(Upred$lU == min(Upred$lU))])$coefficients[2], 
       col= "blue")

summary(lm(Upred$lU[-which(Upred$lU == min(Upred$lU))]~
             dm[-which(Upred$lU == min(Upred$lU))])) 
# r^2 0.11; p = 1.8*10^-6; m = -0.02

#let's look at this pattern in the chinle
UpD <- cbind(Upred, dm)

inCH <- UpD[which(UpD$uinChF == 1),]

plot(inCH$dm, inCH$lU)
abline(median(inCH$lU), 0)
abline(lm(inCH$lU~inCH$dm)$coefficients[1],
       lm(inCH$lU~inCH$dm)$coefficients[2], 
       col= "blue")

summary(lm(inCH$lU~inCH$dm)) #r2 = 0.21; p = 0.001; m = -0.05

#outside of the Chinle
nnCH <- UpD[which(UpD$uinChF == 0),]

plot(nnCH$dm, nnCH$lU)
abline(median(nnCH$lU), 0)
abline(lm(nnCH$lU~nnCH$dm)$coefficients[1],
       lm(nnCH$lU~nnCH$dm)$coefficients[2], 
       col= "blue")

summary(lm(nnCH$lU~nnCH$dm)) #r2 = 0.11; p = 1.5*10^-5; m = -0.03

summary(lm(nnCH$lU[-which(nnCH$lU == min(nnCH$lU))]~
             nnCH$dm[-which(nnCH$lU == min(nnCH$lU))])) 
#r2 = 0.09; p = 4*10^-4; m = -0.018

### median distance to the nearest mines 
mdm <- rep(NA, length(Uc[,1]))
for (i in 1:length(Uc[,1])){
  print(i)
  dwm <- ((((as.numeric(f.mines[,2])-as.numeric(Uc[i,2]))*90)^2) + 
            (((as.numeric(f.mines[,1])-as.numeric(Uc[i,1]))*111)^2))^(1/2)
  mdm[i] <- median(sort(dwm)[1:1])
}
# 90 km to 1 deg longitude,   111 km to 1 deg latitude


#quick plotting
plot(mdm, Upred$lU)
abline(median(Upred$lU), 0)
abline(lm(Upred$lU~mdm)$coefficients[1],lm(Upred$lU~mdm)$coefficients[2], 
       col= "blue")

summary(lm(Upred$lU~mdm)) 
# r^2 0.14; p 2*10^-8; m = -0.03 # for nearest mines
# r^2 0.19; p 2*10^-11; m = -0.03 # for nearest 5 mines
# r^2 0.18; p 2*10^-10; m = -0.03 # for nearest 10 mines
# r^2 0.17; p 3*10^-10; m = -0.03 # for nearest 20 mines
# r^2 0.16; p 1*10^-9; m = -0.02 # for nearest 40 mines
# r^2 0.13; p 2*10^-7; m = -0.02 # for nearest 80 mines

#writing a function to save these vectors
med.dist.n <- function(df1, df2, n) {
  mlist <- rep(NA, length(df1[,1]))
  
  for (i in 1:length(df1[,1])) {
    print(i)
    
    #find distance to a mine
    dm <- ((((as.numeric(df2[,2])-as.numeric(df1[i,3]))*90)^2) + 
              (((as.numeric(df2[,1])-as.numeric(df1[i,4]))*111)^2))^(1/2)
    
    #median dist of n closest mines to well 
    mlist[i] <- median(sort(dm)[1:n])
  }
  
  return(mlist)
}

md1 <- med.dist.n(Upred, f.mines, 1)
md5 <- med.dist.n(Upred, f.mines, 5)


#view this pattern in the chinle
UpmD <- cbind(Upred, mdm)

inCHm <- UpmD[which(UpmD$uinChF == 1),]

plot(inCHm$mdm, inCHm$lU)
abline(median(inCHm$lU), 0)
abline(lm(inCHm$lU~inCHm$mdm)$coefficients[1],
       lm(inCHm$lU~inCHm$mdm)$coefficients[2], 
       col= "blue")

summary(lm(inCHm$lU~inCHm$mdm)) 
#m = -0.048,  r^2 = 0.26; p = 2*10^-4; nearest 5 mines 
#m = -0.037,  r^2 = 0.20; p = 1*10^-3; nearest 10 mines
#m = -0.035,  r^2 = 0.21; p = 8*10^-3; nearest 20 mines
#m = -0.027,  r^2 = 0.17; p = 3*10^-3; nearest 40 mines
#m = -0.021,  r^2 = 0.14; p = 8*10^-3; nearest 80 mines

### looking outside of chinle - don't want to write a function

mdm <- rep(NA, length(nnCH[,1]))
for (i in 1:length(nnCH[,1])){
  print(i)
  dwm <- ((((as.numeric(f.mines[,2])-as.numeric(nnCH[i,3]))*90)^2) + 
            (((as.numeric(f.mines[,1])-as.numeric(nnCH[i,4]))*111)^2))^(1/2)
  mdm[i] <- median(sort(dwm)[1:80])
}

plot(mdm, nnCH$lU)
abline(median(nnCH$lU), 0)
abline(lm(nnCH$lU~mdm)$coefficients[1],
       lm(nnCH$lU~mdm)$coefficients[2], 
       col= "blue")

summary(lm(nnCH$lU~mdm)) 
#m = -0.025,  r^2 = 0.11; p = 1.5*10^-5; nearest mine 
#m = -0.026,  r^2 = 0.16; p = 1.3*10^-7; median nearest 5 mines 
#m = -0.023,  r^2 = 0.15; p = 4.3*10^-7; median nearest 10 mines 
#m = -0.022,  r^2 = 0.14; p = 9.7*10^-7; median nearest 20 mines 
#m = -0.021,  r^2 = 0.13; p = 2.1*10^-6; median nearest 40 mines 
#m = -0.017,  r^2 = 0.09; p = 8.2*10^-5; median nearest 80 mines 


#### Now which wells are downstream from a mine
nmuw <- rep(NA, length(Upred[,1]))

for (i in 1:length(Upred[,1])) {
  print(i)
  #filter fmines > with elevations > greater than well of interest
  fdf <- f.mines[which(f.mines[,3] >  Upred[i,5]),]
  
  #find distance to all mines downstream
  ffdf <- ((((as.numeric(fdf[,2])-as.numeric(Uc[i,2]))*90)^2) + 
             ((as.numeric(fdf[,1])-as.numeric(Uc[i,1]))*111)^2)^(1/2)
  #90 km per degree lon; 111 km per degree lat
  
  #add nearest upstream mine
  nmuw[i] <- median(sort(ffdf)[1:5])
}   

plot(nmuw, Upred[,2])
abline(median(Upred[,2]),0)

summary(lm(Upred[,2]~nmuw))
abline(lm(Upred[,2]~nmuw)$coefficients[1],
       lm(Upred[,2]~nmuw)$coefficients[2], 
       col= "blue")
#m = -0.020,  r^2 = 0.21; p = 1*10^-12;  Nearest mine downstream
#m = -0.018,  r^2 = 0.21; p = 5*10^-12;  Median nearest 5 mines downstream
#m = -0.015,  r^2 = 0.17; p = 4*10^-10;  Median nearest 10 mines downstream
#m = -0.012,  r^2 = 0.15; p = 4*10^-9;  Median nearest 20 mines downstream
#m = -0.009,  r^2 = 0.10; p = 4*10^-6;  Median nearest 40 mines downstream
#m = -0.007,  r^2 = 0.07; p = 8*10^-5;  Median nearest 80 mines downstream

# writing a function to save these vectors
ddmines <- function(df1, df2, n) {
  
  nmuw <- rep(NA, length(df1[,1]))
  
  for (i in 1:length(df1[,1])) {
    print(i)
    #filter fmines > with elevations > greater than well of interest
    fdf <- df2[which(df2[,3] >  df1[i,5]),]
    
    #find distance to all mines downstream
    ffdf <- ((((as.numeric(fdf[,2])-as.numeric(df1[i,3]))*90)^2) + 
               ((as.numeric(fdf[,1])-as.numeric(df1[i,4]))*111)^2)^(1/2)
    #90 km per degree lon; 111 km per degree lat
    
    #add nearest upstream mine
    nmuw[i] <- median(sort(ffdf)[1:n])
    }
  
  return(nmuw)
}

dd1 <- ddmines(Upred, f.mines, 1)
dd5 <- ddmines(Upred, f.mines, 5)
dd10 <- ddmines(Upred, f.mines, 10)
dd20 <- ddmines(Upred, f.mines, 20)
dd40 <- ddmines(Upred, f.mines, 40)
dd80 <- ddmines(Upred, f.mines, 80)

##looking at this in the chinle

#### Now which wells are downstream from a mine
nmuw <- rep(NA, length(inCH[,1]))

for (i in 1:length(inCH[,1])) {
  print(i)
  #filter fmines > with elevations > greater than well of interest
  fdf <- f.mines[which(f.mines[,3] >  inCH[i,5]),]
  
  #find distance to all mines downstream
  ffdf <- ((((as.numeric(fdf[,2])-as.numeric(inCH[i,3]))*90)^2) + 
             ((as.numeric(fdf[,1])-as.numeric(inCH[i,4]))*111)^2)^(1/2)
  #90 km per degree lon; 111 km per degree lat
  
  #add nearest upstream mine
  nmuw[i] <- median(sort(ffdf)[1:5])
}   

plot(nmuw, inCH[,2])
abline(median(inCH[,2]),0)
abline(lm(inCH[,2]~nmuw)$coefficients[1],
       lm(inCH[,2]~nmuw)$coefficients[2], 
       col= "blue")
summary(lm(inCH[,2]~nmuw))  # these stats are wrong - redo them
#m = -0.04,  r^2 = 0.21; p = 0.001; Nearest upstream mine in the Chinle 
#m = -0.02,  r^2 = 0.16; p = 0.005; Median Nearest 5 upstream mines in the Chinle 

## Looking outside of chinle
nmuw <- rep(NA, length(nnCH[,1]))

for (i in 1:length(nnCH[,1])) {
  print(i)
  #filter fmines > with elevations > greater than well of interest
  fdf <- f.mines[which(f.mines[,3] >  nnCH[i,5]),]
  
  #find distance to all mines downstream
  ffdf <- ((((as.numeric(fdf[,2])-as.numeric(nnCH[i,3]))*90)^2) + 
             ((as.numeric(fdf[,1])-as.numeric(nnCH[i,4]))*111)^2)^(1/2)
  #90 km per degree lon; 111 km per degree lat
  
  #add nearest upstream mine
  nmuw[i] <- median(sort(ffdf)[1:1]) 
}   

plot(nmuw, nnCH[,2])
abline(median(nnCH[,2]),0)
abline(lm(nnCH[,2]~nmuw)$coefficients[1],
       lm(nnCH[,2]~nmuw)$coefficients[2], 
       col= "blue")
summary(lm(nnCH[,2]~nmuw))

#m = -0.018,  r^2 = 0.20; p = 3*10^-9; Nearest upstream mine outside Chinle 
#m = -0.016,  r^2 = 0.19; p = 7*10^-9; Nearest 5 upstream mines outside Chinle 
#m = -0.014,  r^2 = 0.17; p = 7*10^-8; Nearest 10 upstream mines outside Chinle 
#m = -0.011,  r^2 = 0.15; p = 4*10^-7; Nearest 20 upstream mines outside Chinle 
#m = -0.008,  r^2 = 0.10; p = 6*10^-5; Nearest 40 upstream mines outside Chinle 
#m = -0.007,  r^2 = 0.07; p = 0.0005; Nearest 80 upstream mines outside Chinle 


# I would like to write a function to get the number of nearest mines 
# in all directions

#function is specific to the dataframes used here,
#could be made more general

allmines <- function(df1, df2, d) {
  mlist <- rep(NA, length(df1[,1]))
  
  for (i in 1:length(df1[,1])) {
    print(i)
      
    #find mines within a distance from a well
    dwm <- ((((as.numeric(df2[,2])-as.numeric(df1[i,3]))*90)^2) + 
              (((as.numeric(df2[,1])-as.numeric(df1[i,4]))*111)^2))^(1/2)
    f.dwm <- dwm[which(dwm < d)]
    
    #add number of mines within d kmof a well
    mlist[i] <- length(f.dwm)
  }
  
  return(mlist)
}

#check function #run at 0.04 degrees (5 km)
am5 <- allmines(Upred, f.mines, 5) 
am4 <- allmines(Upred, f.mines, 4)
am3 <- allmines(Upred, f.mines, 3)
am2 <- allmines(Upred, f.mines, 2) # max here is eleven mines
am1 <- allmines(Upred, f.mines, 1) # this has a max of 4 mines
am6 <- allmines(Upred, f.mines, 6)
am7 <- allmines(Upred, f.mines, 7)
am_p5 <- allmines(Upred, f.mines, 7)

#### exploratory linear models

#1 km
nam1 <- length(which(am1 > 0)) #24
plot(Upred$lU ~ am1)
summary(lm(Upred$lU ~ am1)) # m = 0.31 p = 0.001 max = 4 r2 = 0.05

#2 km
nam2 <- length(which(am2 > 0)) #44
plot(Upred$lU ~ am2)
summary(lm(Upred$lU ~ am2)) # m = 0.14 p = 2.1*10^-5 max = 11 r2 = 0.08

#3 km
nam3 <- length(which(am3 > 0)) #57
plot(Upred$lU ~ am3)
summary(lm(Upred$lU ~ am3)) # m = 0.08 p = 2.4*10^-5 max = 22 r2 = 0.08

#4 km
nam4 <- length(which(am4 > 0)) #68
plot(Upred$lU ~ am4)
summary(lm(Upred$lU ~ am4)) # m = 0.06 p = 1.2*10^-6 max = 26 r2 = 0.11

#5 km
nam5 <- length(which(am5 > 0)) #78
plot(Upred$lU ~ am5)
summary(lm(Upred$lU ~ am5)) # m = 0.04 p = 1.9*10^-7 max = 35 r2 = 0.12


#####################################

#create a matrix of important vectors 

m.table <- cbind(am1, am2, md1, md5, dd1, dd5)


# Write this table to a matrix 
## write.csv(m.table, "MineDistances.csv")

