#Multipanel2 generation and analyses

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
        

#reading data
Ucaus <- read.csv('Ucaus.csv')
Ucaus <- Ucaus[,-1] 

mtable <- read.csv("MineDistances.csv",row.names =1)

Ucaus <- cbind(Ucaus, mtable)

#Only looking at points in the Chinle
UinCH <- Ucaus[which(Ucaus[,8] > 0),] #length = 49
CH_325 <- UinCH[which(UinCH[,7] > 0),] # length = 4

plot(UinCH[,1], UinCH[,4]) #elevation looks to be important in Chinle

s.UvEnCH <- summary(lm(UinCH[,1] ~ UinCH[,4])) #p = 2*10^-5 r2 = 0.31

w.MnCh <- wilcox.test(UinCH[which(UinCH[,11] > 0),1], 
            UinCH[which(UinCH[,11] == 0),1]) # p = 0.06 
# presence of a mine in chinle not predictive of U

#Looking at points outside of the Chinle
UnnCH <- Ucaus[which(Ucaus[,8] == 0),]

plot(UnnCH[,1], UnnCH[,4]) # May not be sig

s.UvENnCH <- summary(lm(UnnCH[,1] ~ UnnCH[,4])) #sig but r2 is 0.03 p = 0.008

# presence of a mine outside of the Chinle
w.MNnCh <- wilcox.test(UnnCH[which(UnnCH[,11] > 0),1], 
            UnnCH[which(UnnCH[,11] == 0),1]) #p = 0.001, w = 2143 
# this is an important factor 

#Analyzing the presence of a mine in elevation bands outside of Chinle
w.e1 <- wilcox.test(UnnCH[which(UnnCH[,11] > 0 & UnnCH[,4] > 1350 
                        & UnnCH[,4] < 1600),1], 
            UnnCH[which(UnnCH[,11] == 0 & UnnCH[,4] > 1350 
                        & UnnCH[,4] < 1600),1]) #p = 0.52

w.e2 <- wilcox.test(UnnCH[which(UnnCH[,11] > 0 & UnnCH[,4] > 1600 
                        & UnnCH[,4] < 1801),1], 
            UnnCH[which(UnnCH[,11] == 0 & UnnCH[,4] > 1600 
                        & UnnCH[,4] < 1801),1]) #p = 0.04

w.e3 <- wilcox.test(UnnCH[which(UnnCH[,11] > 0 & UnnCH[,4] > 1800 
                        & UnnCH[,4] < 2001),1], 
            UnnCH[which(UnnCH[,11] == 0 & UnnCH[,4] > 1800 
                        & UnnCH[,4] < 2001),1]) # p = 0.02
 
w.e4 <- wilcox.test(UnnCH[which(UnnCH[,11] > 0 & UnnCH[,4] > 2000 
                        & UnnCH[,4] < 2300),1], 
            UnnCH[which(UnnCH[,11] == 0 & UnnCH[,4] > 2000 
                        & UnnCH[,4] < 2300),1]) # p = 0.005

## p.6a
p.6a <- ggplot()+
  geom_segment(aes(y = 1.477, x = 2200, yend = 1.477, xend = 1200),
               color = "gray", alpha = 0.5, size = 1, linetype = 2) +
  annotate(geom="text", y=1.7, x=2100, 
           label=expression(paste('30 '*italic("\u00B5")*'g/L')),
           color="gray", size=4)+
  geom_point(data = UinCH, aes(y=UinCH[,1], x=UinCH[,4]),
             color = "black", alpha = 0.5, pch = 16, size = 3)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  geom_segment(
    aes(y = -0.0415, x = 2200, yend = 1.3585, xend = 1200),
    color = "blue", size = 1) +
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_x_continuous(breaks=c(1100,1400,1700,2000,2300),
                     limits = c(1100, 2300))+
  annotate(geom="text", y=-1.4, x=1700, 
           label=expression(paste('[U] = 10^(-0.0014*'*italic(z)~'+ 3.03)')),
           color="black", size=5)+
  annotate(geom="text", y=-1.7, x=1700, 
           label=expression(paste(R^{2},'= 0.32')),
           color="black", size=5)+
  labs(title = 'Within the Chinle Fm.',
       y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')),
       x="Elevation (m)")+
  theme(panel.background = element_blank())

plot(p.6a)

#### Plot 6b
p.6b <-  ggplot()+
  geom_segment(aes(y = 1.477, x = 0, yend = 1.477, xend = 40),
               color = "gray", alpha = 0.5, size = 1, linetype = 2) +
  annotate(geom="text", y=1.7, x=30, 
           label=expression(paste('30 '*italic("\u00B5")*'g/L')),
           color="gray", size=4)+
  geom_point(data = UinCH, aes(y=UinCH$lU, x=UinCH$dd1),
                   color = "black", alpha = 0.5, pch = 16, size = 4)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  geom_segment(
    aes(y = 1.14, x = 0, yend = -0.26, xend = 40),
    color = "blue", alpha = 0.5, size = 1) +
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_x_continuous(breaks=c(0,10,20,30,40),
                     limits = c(0, 40))+
  annotate(geom="text", y=-1.4, x=20, 
           label=expression(paste('[U] = 10^(-0.035*'*italic(d)~'+ 1.14)')),
           color="black", size=5)+
  annotate(geom="text", y=-1.7, x=20, 
           label=expression(paste(R^{2},'= 0.19')),
           color="black", size=5)+
  labs(title = 'Within the Chinle Fm.',
       y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), 
       x="Distance to the nearest mine upstream (km)")+
  theme(panel.background = element_blank())

plot(p.6b)

#### Panel 6c
p.6c <- ggplot()+
  geom_segment(aes(y = 1.477, x = 1200, yend = 1.477, xend = 2300),
               color = "gray", alpha = 0.5, size = 1, linetype = 2) +
  annotate(geom="text", y=1.7, x=1300, 
           label=expression(paste('30 '*italic("\u00B5")*'g/L')),
           color="gray", size=4)+
  geom_point(data = UnnCH, aes(x=UnnCH[,4], y=UnnCH[,1]),
             color = "black", alpha = 0.5, pch = 16, size = 3)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_x_continuous(breaks = c(1100,1400,1700,2000,2300),
                     limits = c(1100, 2300))+
  labs(title = 'Outside of the Chinle Fm.',
       x="Elevation (m)", 
       y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')))+
  theme(panel.background = element_blank())

plot(p.6c)

#Plot 6d
p.6d <- ggplot() +
  geom_segment(aes(y = 1.477, x = 0, yend = 1.477, xend = 90),
               color = "gray", alpha = 0.5, size = 1, linetype = 2) +
  annotate(geom="text", y=1.7, x=20, 
           label=expression(paste('30 '*italic("\u00B5")*'g/L')),
           color="gray", size=4)+
  geom_point(data = UnnCH, aes(y=lU, x=dd1),
             color = "black", alpha = 0.5, pch = 16, size = 4)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  geom_segment(
    aes(y = 0.62, x = 0, yend = -1, xend = 90),
    color = "blue", alpha = 0.5, size = 1) +
  scale_y_continuous(breaks = c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90),
                     limits = c(0, 90))+
  annotate(geom="text", y=2.6, x=60, 
           label=expression(paste('[U] = 10^(-0.018*'*italic(d)~'+ 0.62)')),
           color="black", size=5)+
  annotate(geom="text", y=2.3, x=60, 
           label=expression(paste(R^{2},'= 0.20')),
           color="black", size=5)+
  labs(title = 'Outside of the Chinle Fm.',
       y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), 
       x="Distance to the nearest mine upstream (km)")+
  theme(panel.background = element_blank())

plot(p.6d)

require("cowplot")

#export as 8 x 8
plot_grid(p.6a, p.6b, p.6c, p.6d,
          labels = c("A","B","C","D"),
          ncol = 2)

#### Plot Mines Total
p.mt <-  ggplot()+
  geom_segment(aes(y = 1.477, x = 0, yend = 1.477, xend = 90),
               color = "gray", alpha = 0.5, size = 1, linetype = 2) +
  annotate(geom="text", y=1.7, x=50, 
           label=expression(paste('30 '*italic("\u00B5")*'g/L')),
           color="gray", size=4)+
  geom_point(data = Ucaus, aes(y=lU, x=dd1),
             color = "black", alpha = 0.5, pch = 16, size = 4)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  geom_segment(
    aes(y = 0.74, x = 0, yend = -1.06, xend = 90),
    color = "blue", alpha = 0.5, size = 1) +
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90),
                     limits = c(0, 90))+
  annotate(geom="text", y=2.7, x=60, 
           label=expression(paste('[U] = 10^(-0.020*'*italic(d)~'+ 0.74)')),
           color="black", size=5)+
  annotate(geom="text", y=2.3, x=60, 
           label=expression(paste(R^{2},'= 0.22')),
           color="black", size=5)+
  labs(y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), 
       x="Distance to the nearest mine upstream (km)")+
  theme(panel.background = element_blank())

plot(p.mt) # export as 5x5 

