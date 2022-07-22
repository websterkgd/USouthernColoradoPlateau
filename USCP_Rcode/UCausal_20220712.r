#Loose Causal Analysis of U data

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
        

#reading data
Ucaus <- read.csv('Ucaus.csv')
Ucaus <- Ucaus[,-1] 

#Filtering the Chinle FM

UinCH <- Ucaus[which(Ucaus[,8] > 0),] #length = 49
CH_325 <- UinCH[which(UinCH[,7] > 0),] # length = 4

plot(UinCH[,1], UinCH[,4]) #elevation related to U in Chinle

s.UvEnCH <- summary(lm(UinCH[,1] ~ UinCH[,4])) #p = 2*10^-5 r2 = 0.31

w.MnCh <- wilcox.test(UinCH[which(UinCH[,11] > 0),1], 
            UinCH[which(UinCH[,11] == 0),1]) # p = 0.06 
# presence of a mine in chinle not predictive of U

#Looking at points outside of the Chinle
UnnCH <- Ucaus[which(Ucaus[,8] == 0),]

plot(UnnCH[,1], UnnCH[,4]) # 

s.UvENnCH <- summary(lm(UnnCH[,1] ~ UnnCH[,4])) #sig but r2 is 0.03 p = 0.008

# presence of a mine outside of the Chinle
w.MNnCh <- wilcox.test(UnnCH[which(UnnCH[,11] > 0),1], 
            UnnCH[which(UnnCH[,11] == 0),1]) #p = 0.001, w = 2143 

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

#### making nice Figures

## p.6a
p.6a <- ggplot()+
  geom_point(data = UinCH, aes(x=UinCH[,1], y=UinCH[,4]),
             color = "black", alpha = 0.5, pch = 16, size = 3)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  geom_segment(
    aes(x = -0.0415, y = 2200, xend = 1.3585, yend = 1200),
    color = "blue", size = 1) +
  scale_x_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_y_continuous(breaks=c(1100,1400,1700,2000,2300),
                     limits = c(1100, 2300))+
  annotate(geom="text", x=1.8, y=2200, 
           label=expression(paste('[U] = 10^(-0.0014*'*italic(z)~'+ 3.03)')),
           color="black", size=3)+
  annotate(geom="text", x=1.8, y=2100, 
           label=expression(paste(r^{2},'= 0.32')),
           color="black", size=3)+
  labs(title = 'Within the Chinle Fm.',
       y="Elevation (m)", 
       x=expression(paste('[U] ('*italic("\u00B5")*'g/L)')))+
  theme(panel.background = element_blank())

plot(p.6a)

#### Plot 6b
MiCH <- UinCH[which(UinCH[,11] > 0),1] # data
MnCH <- UinCH[which(UinCH[,11] == 0),1] # data

lMiCH <- cbind(MiCH, (rep('Mine within 2 km', length(MiCH))))
lMnCH <- cbind(MnCH, (rep('Mine not within 2 km', length(MnCH))))

M2.CH <- as.data.frame(rbind(lMiCH,lMnCH))


p.6b <- ggplot() +
  geom_boxplot(aes(x = M2.CH[,2], y = as.numeric(M2.CH[,1])), 
               outlier.shape = NA)+
  geom_jitter(aes(x = M2.CH[,2], y = as.numeric(M2.CH[,1])), 
              width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-1,0,1,2,3),
                     labels = c(expression(0.1,1,10,100,1000)),
                     limits = c(-1, 3))+
  labs(title = 'Within the Chinle Fm.', 
       y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), 
       x = "")+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  theme(panel.background = element_blank()) 

plot(p.6b)

#### Panel 6c
p.6c <- ggplot()+
  geom_point(data = UnnCH, aes(x=UnnCH[,1], y=UnnCH[,4]),
             color = "black", alpha = 0.5, pch = 16, size = 3)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  scale_x_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_y_continuous(breaks=c(1100,1400,1700,2000,2300),
                     limits = c(1100, 2300))+
  labs(title = 'Outside of the Chinle Fm.',
       y="Elevation (m)", 
       x=expression(paste('[U] ('*italic("\u00B5")*'g/L)')))+
  theme(panel.background = element_blank())

plot(p.6c)

#Plot 6d
iMCH <- UnnCH[which(UnnCH[,11] > 0),1] # data
nMCH <- UnnCH[which(UnnCH[,11] == 0),1] # data

liMCH <- cbind(iMCH, (rep('Mine within 2 km', length(iMCH))))
lnMCH <- cbind(nMCH, (rep('Mine not within 2 km', length(nMCH))))

M2.N <- as.data.frame(rbind(liMCH,lnMCH))

p.6d <- ggplot() +
  geom_boxplot(aes(x = M2.N[,2], y = as.numeric(M2.N[,1])), 
               outlier.shape = NA)+
  geom_jitter(aes(x = M2.N[,2], y = as.numeric(M2.N[,1])), 
              width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  labs(title = 'Outside of the Chinle Fm.', 
       y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), x = "")+
  annotate(geom="text", x="Mine within 2 km", y=2.9, label="*",
           color="red", size=7)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  theme(panel.background = element_blank()) 

plot(p.6d)

require("cowplot")

#export as 8 x 8
plot_grid(p.6a, p.6b, p.6c, p.6d,
          labels = c("A","B","C","D"),
          ncol = 2)

