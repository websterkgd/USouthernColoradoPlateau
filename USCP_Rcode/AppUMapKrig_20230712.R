#pulling in radiation values into the analysis

rm(list=ls())
getwd()

require(ggplot2)

# pulling in xyz files
SRR <- read.table("shiprock_rad.xyz", header= FALSE, sep = "") #reads file
FgR <- read.table("flagstaff_rad.xyz", header= FALSE, sep = "") #reads file
GpR <- read.table("gallup_rad.xyz", header= FALSE, sep = "") #reads file
MCR <- read.table("marble_canyon_rad.xyz", header= FALSE, sep = "") #reads file

#From Meta File

# 1  line        flight line number
# 2  fid         fiducial number (integer)
# 3  time        time (hhmmss)
# 4  day         Julian day flown (integer)
# 5  year        year flown (integer)
# 6  latitude    latitude (decimal degrees)
# 7  longitude   longitude (decimal degrees)
# 8  radalt      radar altimeter reading above ground (meters)
# 9  resmag      residual magnetic value (nT)
# 10 geology     surficial geology beneath flight line (coded)
# 11 qual        quality flags of the radiometrics (integer)
# 12 app_K       apparent Potassium as Potassium 40 (percent potassium)
# 13 app_U       apparent Uranium as Bismuth 214 (parts per million 
#                                                 equivalent uranium)
# 14 app_Th      apparent Thorium as Thallium 208 (parts per million 
#                                                  equivalent thorium)
# 15 U_Th_ratio  ratio of Uranium and Thorium
# 16 U_K_ratio   ratio of Uranium and Potassium
# 17 Th_K_ratio  ratio of Thorium and Potassium
# 18 total_count total count of radioactivity (counts/second)
# 19 atmos_BI214 atmospheric Uranium as Bi214 (counts/second)
# 20 air_temp    air temperature (degrees Celsius)
# 21 air_press   air pressure (mmHg)

#compiling info. 
RC <- rbind(FgR[,c(6,7,13)],SRR[,c(6,7,13)],MCR[,c(6,7,13)],GpR[,c(6,7,13)])

#quick ish plot
plot(SRR$V7,SRR$V6, cex = 0.02*(SRR$V18)^(1/2))

#removing larger datasets
rm(FgR,GpR,MCR,SRR)

###setwd()

MSE <- read.table('Mspatelev.txt', header = TRUE, sep ="\t", fill = TRUE)

Ucaus <- read.csv('Ucaus.csv')
Ucaus <- Ucaus[,-1] 

#quick visual check here
plot(RC[,2],RC[,1], cex = 0.1)
points(MSE$longitude, MSE$latitude, pch = 4)
points(Ucaus[,2], Ucaus[,3], pch = 5)

#filtering out negative values in SRR
RC <- RC[which(RC[,3] >= 0),]

# boxplot(RC[which(RC[,3] >= 0),3])
# 
# length(RC[which(RC[,3] < 0),3])

m_SR <- MSE[which(MSE[,3] >= min(RC[,2]) & 
                    MSE[,3] <= max(RC[,2]) &
                    MSE[,2] >= min(RC[,1]) & 
                    MSE[,2] <= max(RC[,1])), ]

w_IR <- Ucaus[which(Ucaus[,2] >= min(RC[,2]) & 
                    Ucaus[,2] <= max(RC[,2]) &
                    Ucaus[,3] >= min(RC[,1]) & 
                    Ucaus[,3] <= max(RC[,1])), ]

#quick visual check
plot(RC[,2],RC[,1], cex = 0.1)
points(m_SR$longitude, m_SR$latitude, pch = 4)
points(w_IR[,2], w_IR[,3], pch = 5)

#linear interpolation function  

### breaks at i = 305, 366, 932, 942. 948  for reasons unknown

radinterp1 <- function(iLat, iLon, kLat, kLon, kRad) {
  rad.mine <- rep(NA, length(iLat))
  for (i in c(1:304,306:365,367:931,933:941,943:947,949:length(iLat))){
    lon_dif <- abs(rep(iLon[i], length(kLon)) - kLon)
    rmat <- cbind(kLat, kLon, kRad)
    mlon <- rmat[which(lon_dif == min(lon_dif)),]
    alat_dif <- abs(rep(iLat[i], length(mlon[, 1])) - mlon[, 1])
    lat_dif <- rep(iLat[i], length(mlon[, 1])) - mlon[, 1]
    
    print(i)
    
    if ((length(alat_dif) > 1) & (sum(alat_dif + lat_dif) != 0) &
        (iLat[i] - mlon[which(alat_dif == min(alat_dif)),1]) < 0){
      h <- mlon[which(alat_dif == min(alat_dif)),c(1,3)]
      lat_dif[lat_dif < 0] <- 10
      l <- mlon[which(abs(lat_dif) == min(abs(lat_dif))),c(1,3)]
      rad.mine[i] <- approx(x = c(as.numeric(h[1]),as.numeric(l[1])), 
                            y = c(as.numeric(h[2]),as.numeric(l[2])),
                            method = "linear", xout = iLat[i])[2]
      # this section works
    } else if ((length(alat_dif) > 1) & (sum(alat_dif + lat_dif) == 0)) {
      rad.mine[i] <- rmat[which(lon_dif == min(lon_dif)),3]
    } else if ((length(alat_dif) > 1) & (sum(alat_dif - lat_dif) != 0) &
               (iLat[i] - mlon[which(alat_dif == min(alat_dif)),1]) > 0){
      h <- mlon[which(alat_dif == min(alat_dif)),c(1,3)]
      lat_dif[lat_dif > 0] <- 10
      l <- mlon[which(abs(lat_dif) == min(abs(lat_dif))),c(1,3)]
      rad.mine[i] <- approx(x = c(as.numeric(h[1]),as.numeric(l[1])), 
                            y = c(as.numeric(h[2]),as.numeric(l[2])),
                            method = "linear", xout = iLat[i])[2]
    } 
    else {
      rad.mine[i] <- rmat[which(lon_dif == min(lon_dif)),3] # fills out
    }
  }
  rad.mine <- as.numeric(rad.mine)
  return(rad.mine) #final return
}

mine.rad <- radinterp1(m_SR[,2],m_SR[,3],RC[,1],RC[,2],RC[,3])


#linear interpolation at waters

radinterp <- function(iLat, iLon, kLat, kLon, kRad) {
  rad.mine <- rep(NA, length(iLat))
  for (i in 1:length(iLat)){
    lon_dif <- abs(rep(iLon[i], length(kLon)) - kLon)
    rmat <- cbind(kLat, kLon, kRad)
    mlon <- rmat[which(lon_dif == min(lon_dif)),]
    alat_dif <- abs(rep(iLat[i], length(mlon[, 1])) - mlon[, 1])
    lat_dif <- rep(iLat[i], length(mlon[, 1])) - mlon[, 1]
    
    print(i)
    
    if ((length(alat_dif) > 1) & (sum(alat_dif + lat_dif) != 0) &
        (iLat[i] - mlon[which(alat_dif == min(alat_dif)),1]) < 0){
      h <- mlon[which(alat_dif == min(alat_dif)),c(1,3)]
      lat_dif[lat_dif < 0] <- 10
      l <- mlon[which(abs(lat_dif) == min(abs(lat_dif))),c(1,3)]
      rad.mine[i] <- approx(x = c(as.numeric(h[1]),as.numeric(l[1])), 
                            y = c(as.numeric(h[2]),as.numeric(l[2])),
                            method = "linear", xout = iLat[i])[2]
      # this section works
    } else if ((length(alat_dif) > 1) & (sum(alat_dif + lat_dif) == 0)) {
      rad.mine[i] <- rmat[which(lon_dif == min(lon_dif)),3]
    } else if ((length(alat_dif) > 1) & (sum(alat_dif - lat_dif) != 0) &
               (iLat[i] - mlon[which(alat_dif == min(alat_dif)),1]) > 0){
      h <- mlon[which(alat_dif == min(alat_dif)),c(1,3)]
      lat_dif[lat_dif > 0] <- 10
      l <- mlon[which(abs(lat_dif) == min(abs(lat_dif))),c(1,3)]
      rad.mine[i] <- approx(x = c(as.numeric(h[1]),as.numeric(l[1])), 
                            y = c(as.numeric(h[2]),as.numeric(l[2])),
                            method = "linear", xout = iLat[i])[2]
    } 
    else {
      rad.mine[i] <- rmat[which(lon_dif == min(lon_dif)),3] # fills out
    }
  }
  rad.mine <- as.numeric(rad.mine)
  return(rad.mine) #final return
}

waters.rad <- radinterp(w_IR[,3],w_IR[,2],RC[,1],RC[,2],RC[,3])


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue','green','yellow','orange','red'))
Col2 <- rbPal(30)[as.numeric(cut(mine.rad,breaks = 30))]
Col1 <- rbPal(30)[as.numeric(cut(waters.rad,breaks = 30))]

# Quick plotting for a visual check
ggplot() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  geom_point(aes(x =m_SR$longitude, y = m_SR$latitude), 
             color = Col2, pch = 15, 
             size = 1)+
  # geom_point(aes(x = RC[,2],
  #                y = RC[,1]),
  #            color = "black", pch = 2,
  #            size = 0.1)+
  geom_point(aes(x = w_IR[,2], 
                 y = w_IR[,3]),
             color = Col1, pch = 18, 
             size = 4)

# #Writing vectors and data
# 
# wwr <- cbind(w_IR, waters.rad)
# 
# write.csv(wwr, "WatersWApU.csv")
# 
# ram <- cbind(m_SR[2:4], mine.rad)
# 
# write.csv(ram, "MinesWApU.csv")

#### Interpolation by kriging
require(sp)
require(gstat)
require(dplyr) # for "glimpse"
require(scales) # for "comma"
require(magrittr)

#making g stat RC to calculate the variogram
n_RC <- RC

colnames(n_RC) <- cbind("y","x","lU")
coordinates(n_RC) <- ~x+y

# lU.vgm <- variogram(lU~1, n_RC) # ~ 443 k rows - 
# Takes close to 30 min to run

# 
# write.table(lU.vgm, file =  "U_variogram.txt", sep = "\t", dec =".",
#             col.names = TRUE, row.names = FALSE)
## writing the variogram to a table to load in 
# 
# lU.vgm <- read.table('U_variogram.txt', sep = "\t", dec =".")

# colnames(lU.vgm) <- lU.vgm[1,]
# lU.vgm <- lU.vgm[-1,]
# 
# lU.vgm[1,]

# Variogram fit
lU.fit <- fit.variogram(lU.vgm, model = vgm(600, "Exp", 350))

plot(lU.vgm, lU.fit) 
# calculates sample variogram values

# load spatial domain to interpolate over
#data("meuse.grid")

U_new <- Ucaus 

colnames(U_new) <- cbind("logU","x","y",  "uEl", "uin01", "uin13", "uin325",
                         "uinChF", "uinMF", "mu1", "mu2")      

locs <- U_new[,2:3] 

coordinates(locs) <- ~ x + y # step 3 above

lU.kriged <- krige(n_RC$lU ~ 1, n_RC, locs, model=lU.fit, maxdist = 0.1)
# takes 4 min

wwapu <- read.csv("WatersWApU.csv")

#assessing interpolation vs kriging data
apU_k <- cbind(lU.kriged@coords, lU.kriged@data)

apU_kf <- apU_k[which(apU_k$y %in% wwapu$Average.Latitude),] #filters

#visual check
plot(wwapu$waters.rad ~ apU_kf$var1.pred) # some differences 

summary(lm(wwapu$waters.rad ~ apU_kf$var1.pred)) # r2 = 0.29
# want to pull vector

# write.table(apU_kf, "apU_Krige.txt", sep = "\t", col.names = TRUE, 
#             row.names = FALSE)


##kriging for mines

#calling in mine locs

mlocs <- m_SR[,3:2] 
colnames(mlocs) <- cbind("x", "y")

coordinates(mlocs) <- ~ x + y # step 3 above

# mines.kriged <- krige(n_RC$lU ~ 1, n_RC, mlocs, model=lU.fit, maxdist = 0.1)
# #started at 12:02 pm - finished at 12:20 pm

plot(mine.rad ~ mines.kriged$var1.pred) # 

summary(lm(mine.rad ~ mines.kriged$var1.pred)) # r2 = 0.16

# Writing vector to table
mines.k.table <- cbind(mines.kriged@coords, mines.kriged@data)

# write.table(mines.k.table, "mines_kriged.txt", sep = "\t", col.names = TRUE,
#             row.names = FALSE)
