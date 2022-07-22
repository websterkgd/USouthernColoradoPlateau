#Code to Generate Figure2

rm(list=ls())
getwd()
setwd('C:/Current/Dine/Research/UraniumContamination/Data_WCM')

#loading some packages
require(ggplot2)
require(jpeg) #read_jpg
require(grid) #rastergrob

CCA_FP <- readJPEG('CoconinoDeChellyTDS_ForPub.JPEG', native = TRUE) 
iset <- as.data.frame(rbind(c(0, 0, 1),
                            c(1, 0, 1),
                            c(1, 1, 1),
                            c(0, 1, 1)))
colnames(iset) <- c('long', 'lat', 'group')

p.2a <- ggplot() + 
  annotation_custom(rasterGrob(CCA_FP,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    0, 1, 0, 1) +
  geom_polygon(data = iset, aes(x=long, y = lat, group = group), 
               fill =alpha("gray", 0), color = alpha("black", 0))+ #facilitates plotting
  annotate(geom="text", x=0.05, y=0.95, label="A",
           color="black", size=7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plot(p.2a)

##Figure panel b
CCA_01G <- readJPEG('CCA_01Green.JPEG', native = TRUE) 
iset2 <- as.data.frame(rbind(c(0, 0, 1),
                            c(0.5, 0, 1),
                            c(0.5, 1, 1),
                            c(0, 1, 1)))
colnames(iset2) <- c('long', 'lat', 'group')

p.2b <- ggplot() + 
  annotation_custom(rasterGrob(CCA_01G,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    0, .5, 0, 1) +
  geom_polygon(data = iset, aes(x=long, y = lat, group = group), 
               fill =alpha("gray", 0), color = alpha("black", 0))+ #facilitates plotting
  annotate(geom="text", x=0.03, y=0.95, label="B",
           color="black", size=7) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

plot(p.2b)

###### Figure 2c
DSC2 <- jpeg::readJPEG('CoconinoDeChellyTDS_Cropped.JPEG')

#identifying cluster of interest
kmeans.clusters <- 
  colordistance::getKMeanColors("CoconinoDeChellyTDS_Cropped.JPEG", 
                                n = 15, plotting = TRUE)
colors <- colordistance::extractClusters(kmeans.clusters)
colors_o <- colors[order(colors[,4], decreasing = TRUE),]
print(colors_o)

### Converting to a polygon + kml for future use
require(sp) 
require(rgdal) #used for writing the shape file #also will be retired by 2023
require(ggmap)
require(ggplot2)
require(raster)
require(concaveman)

############## 

#pulling out 0 - 1000
#with 15 colors -> color 3

c_3 <-c(0.8519369, 0.7990071, 0.9171506)

upper.lm.C <- c_3 + 0.014

#pulling out
lower.lm.C <- c_3 - 0.014 

lm.C <- countcolors::rectangularRange(DSC2, 
                                      upper = upper.lm.C, 
                                      lower = lower.lm.C,
                                      target.color = "green",
                                      plotting = TRUE)
#converting area to points
lm.C.I <- lm.C$pixel.idx

#plotting pixels in a somewhat correct orientation
plot(lm.C.I[,2],- lm.C.I[,1])

range(lm.C.I[,1]) #
range(lm.C.I[,2]) #

# the pixels stay fixed

np <- 39.524254   # -1N = 39.524254,
sp <- 33.742784 # -306 S = 33.742784
nssp <-(sp-np)/(-305) # NS line:  

wp <- -111.955907 # 1W = -111.955907
ep <- -107.808380 # 180E = -107.808380
ewsp <- (ep-wp)/(179) # EW line

ccW <- ewsp*range(lm.C.I[,2])[1] + (wp-ewsp) #westmost point 0 - 1000 CC Aq
ccE <- ewsp*range(lm.C.I[,2])[2] + (wp-ewsp) #eastmost point 0 - 1000 CC Aq
ccN <- nssp*-range(lm.C.I[,1])[1] + (np-nssp) #northtmost point 0 - 1000 CC Aq
ccS <- nssp*-range(lm.C.I[,1])[2] + (np-nssp)  #southmost point 0 - 1000 CC Aq

bw <- max(lm.C.I[,2])-min(lm.C.I[,2])
bn <- max(-lm.C.I[,1])-min(-lm.C.I[,1])

dE <- lm.C.I[,2]*(ccE-ccW)/bw + ccW - min(lm.C.I[,2])*(ccE-ccW)/bw

dN <- -lm.C.I[,1]*(ccN-ccS)/bn + ccS - min(-lm.C.I[,1])*(ccN-ccS)/bn

plot(dE, dN)

L = cbind(dE, dN) #L Lowest     

##pulling out 0 - 1000

### Pulling out outline

#need to pull out shapes individually
plot(L)

t <- L[c(which(L[,1] > -111.2 & L[,2] < 36.6)),] #Middle Low polygon 0 - 1000 TDS
plot(t)

LNM <- t[c(which(t[,1] < -108.7 & t[,2] > 35.2)),]
plot(LNM)

LNMo <- concaveman(LNM, concavity = 1, length_threshold = 0.2) 
plot(LNMo)

#### Plot 2c
p.2c <- ggplot() +  
  geom_point(aes(x = LNMo[,1], y = LNMo[,2]), 
             color = 'black', alpha = 1,
             pch = 16)+
  ggtitle("Black Mesa Basin 0 - 1 k TDS")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(axis.line = element_line(color = "black", linetype = "solid"))+
  annotate("text", x = -112.5, y = 37.2, label = "C", size = 8) +
  coord_cartesian(ylim = c(35, 37), xlim =c(-112,-108), clip = "off")
    
plot(p.2c)

#plot 2d
p.2d <- ggplot()+ 
  geom_polygon(data= NULL, aes(x = LNMo[,1], y = LNMo[,2]), 
               fill ="blue", alpha = 0.2)+
  ggtitle("Polygon of Region")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(axis.line = element_line(color = "black", linetype = "solid"))+
  annotate("text", x = -112.5, y = 37.2, label = "D", size = 8) +
  coord_cartesian(ylim = c(35, 37), xlim =c(-112,-108), clip = "off")

plot(p.2d)

require(cowplot) #for multipanel
plot_grid(p.2a, p.2b, p.2c, p.2d, ncol = 2)
#export as 7 x 7




