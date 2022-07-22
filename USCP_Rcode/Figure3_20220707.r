#Figure 3. U concentrations on the colorado plateau

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
require(grid)
require(jpeg)

#Loading in some of the background imagery
IS_CP <- readJPEG('F1_PanelBwoMarkings.JPG', native = TRUE) 
# boundaries = -109, 37.5, -112, 34.8


#Pulling in data
header <- read.csv('CredoEtal2019_238U.csv', nrow=1)
U <- read.csv('CredoEtal2019_238U.csv', skip = 2, col.names=header)

#Renaming the Uranium concentration
Uconc <- U$Overall.Reported.Average..ug.L. #average reported U conc.

Uconc[which(Uconc == "BD")] <- 0.047  #detection limit = 0.047 ug/g

#creating iset 
iset <- as.data.frame(rbind(c(-112, 34.8, 1),
                            c(-109, 34.8, 1),
                            c(-109, 37.5, 1),
                            c(-112, 37.5, 1)))
colnames(iset) <- c('long', 'lat', 'group')


p.3 <- ggplot() +  
  annotation_custom(rasterGrob(IS_CP,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -112, -109, 34.8, 37.5)+
  geom_polygon(data = iset, aes(x=long, y = lat, group = group),
               fill =alpha("gray", 0), color = alpha("gray", 0))+ #facilitates plotting
  geom_segment(
    aes(x = -109.043043, y = 34.8, xend = -109.043043, yend = 37.5), 
    color = "light gray") +
  geom_segment(
    aes(x = -112, y = 37, xend = -109, yend = 37), 
    color = "light gray") +
  geom_point(aes(x = U[,4], y = U[,3]),
             fill="light gray", color="black", alpha = 0.6, pch=21, 
             size=1.6*(sqrt((4/pi)*0.5*as.numeric(Uconc))))+
  geom_point(aes(x = -109.5525, y = 36.1544),
             color = 'dark gray', fill= "red", pch = 23, cex = 3)+
  annotate("text", x=-109.3, y=36.1, label="Chinle",
           color ="red", parse=TRUE, size = 5)+
  geom_point(aes(x = -111.2399, y = 36.1350),
             color = 'dark gray', fill = "red", pch = 23, cex = 3)+
  annotate("text", x=-110.9, y=36.1, label="Tuba*' City'",
           color ="red", parse=TRUE, size = 5)+
  geom_polygon(data = NULL,
               aes(x=c(-112, -109,-109,-112), y = c(34, 34, 34.8, 34.8)),
               fill = "white", color = "black")+
  geom_point(aes(x = -111.8, y = 34.55), 
             color = 'black', alpha = 0.5, fill = "light gray",
             pch = 21, 
             cex = 1.6*sqrt((4/pi)*0.5*as.numeric(0.5)))+
  geom_point(aes(x = -111.8, y = 34.35), 
             color = 'black', alpha = 0.6, fill = "light gray",
             pch = 21, 
             cex = 1.6*sqrt((4/pi)*0.5*as.numeric(10)))+
  geom_point(aes(x = -111.8, y = 34.15), 
             color = 'black', alpha = 0.6, fill = "light gray",
             pch = 21, 
             cex = 1.6*sqrt((4/pi)*0.5*as.numeric(50)))+
  geom_point(aes(x = -110.2, y = 34.4), 
             color = 'black', alpha = 0.6, fill = "light gray",
             pch = 21, 
             cex = 1.6*sqrt((4/pi)*0.5*as.numeric(500)))+
  annotate("text", x=-111.45, y=34.7, label="Uranium*' Concentration'", 
           parse=TRUE)+
  annotate("text", x=-111.5, y=34.55, 
           label=expression(paste('0.5 '*italic("\u00B5")*'g/L')), parse=TRUE)+ 
  annotate("text", x=-111.5, y=34.35, 
           label=expression(paste('5 '*italic("\u00B5")*'g/L')), parse=TRUE)+
  annotate("text", x=-111.5, y=34.15, 
           label=expression(paste('50 '*italic("\u00B5")*'g/L')), parse=TRUE)+
  annotate("text", x=-109.65, y=34.4, 
           label=expression(paste('500 '*italic("\u00B5")*'g/L')), parse=TRUE)+
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

plot(p.3) # Plotting 




