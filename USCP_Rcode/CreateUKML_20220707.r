#Creates a kml file of water source locations on NN

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/UraniumContamination/Data_WCM')

#loading packages needed to create a kml file
require(rgdal) #used to create km l
require(sp) # used for ????
require(ggplot2)
#require(plotKML)
require(raster)
#require(mapview)

#Pulling in data
header <- read.csv('CredoEtal2019_238U.csv', nrow=1)
U <- read.csv('CredoEtal2019_238U.csv', skip = 2, col.names=header)

#Renaming the Uranium concentration
Uconc <- U$Overall.Reported.Average..ug.L. #average reported U conc. 

Uconc[which(Uconc == "BD")] <- 0.047  #detection limit = 0.047 ug/g
 
Uspatial <- ggplot() + geom_point(aes(x = U[,4], y = U[,3]), 
                       color = 'black', alpha = 0.5,
                       pch = 16, cex = sqrt((4/pi)*0.5*as.numeric(Uconc))) +
  theme(axis.ticks=element_blank(), 
        axis.text.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.position = 'none') +
  labs(x=NULL, y=NULL)

plot(Uspatial)
                       
# Save the plot
ggsave(plot=Uspatial, "Uspatial.tiff", device = "tiff")
 
Uc <- cbind(U[,c(4,3)],as.numeric(Uconc))
Uc <- Uc[complete.cases(Uc),]

# prepare coordinates, data, and proj4string
coords <- Uc[ , c(1,2)]   # coordinates
data   <- Uc[ , 3]          # data
crs    <- CRS("+init=epsg:28992") # proj4string of coords

# make the SpatialPointsDataFrame object
Uspat <- SpatialPointsDataFrame(coords      = coords,
                               data        = as.data.frame(data), 
                               proj4string = crs)

##The one line of code below writes our SpatialPointsDataFrame to a KML File
writeOGR(Uspat, dsn="Uspat.kml", layer= "Usp", driver="KML")
