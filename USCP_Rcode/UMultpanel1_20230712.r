#multipanel figure 1

rm(list=ls())
getwd()
setwd('')

UCaus <- read.csv("UCaus.csv")

require(ggplot2)

##Stack Vectors

m_sal <- rbind(cbind(UCaus[which(UCaus$uin01 == 1),2], 
                     rep("0-1000",length(UCaus[which(UCaus$uin01 == 1),1]))),
               (cbind(UCaus[which(UCaus$uin13 == 1),2], 
                      rep("1000-3000",length(UCaus[which(UCaus$uin13 == 1),1])))),
               (cbind(UCaus[which(UCaus$uin325 == 1),2], 
                      rep("3000-25000",length(UCaus[which(UCaus$uin325 == 1),1])))))
               
m_sal <- as.data.frame(m_sal)

p.4a <- ggplot() +
  geom_hline(yintercept = 1.477, color = "gray", alpha = 0.5, size = 1, 
             linetype = 2)+
  annotate("text", x = 1.5, y = 1.7, 
           label = expression(paste('30 '*italic("\u00B5")*'g/L')),
           col = "gray", size = 4)+
  geom_boxplot(aes(m_sal$V2, as.numeric(m_sal$V1)), outlier.shape = NA)+
  geom_jitter(aes(m_sal$V2, as.numeric(m_sal$V1)),width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
               labels = c(expression(0.01,0.1,1,10,100,1000)),
                          limits = c(-2, 3))+
  scale_x_discrete(breaks=c("0-1000","1000-3000","3000-25000"))+
  labs(y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), 
       x="TDS in the Coconino-De Chelly Aquifer (mg/L)")+
  annotate(geom="text", x="3000-25000", y=2.9, label="*",
          color="red", size=7)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  theme(panel.background = element_blank()) 

#plot(p.4a) # Plotting 

##### Elevation vs U  = plot 4b
p.4b <- ggplot()+
  geom_segment(aes(y = 1.477, x = 2500, yend = 1.477, xend = 1000),
               color = "gray", alpha = 0.5, size = 1, linetype = 2) +
  annotate(geom="text", y=1.7, x=1100, 
           label=expression(paste('30 '*italic("\u00B5")*'g/L')),
           color="gray", size=4)+
  geom_point(data = UCaus, aes(x=uEl, y=lU),
             color = "black", alpha = 0.5, pch = 16, size = 3)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"))+
  geom_segment(
    aes(y = -0.585, x = 2500, yend = 1.32, xend = 1000),
    color = "blue", size = 1) +
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_x_continuous(breaks= c(1000,1500,2000,2500),
                     labels = c(1000,1500,2000,2500),
                     limits = c(900, 2600))+
  annotate(geom="text", y=-1.65, x=1400, 
           label=expression(paste('[U] = 10^(-0.0013*'*italic(z)~'+ 2.59)')),
           color="black", size=4)+
  annotate(geom="text", y=-1.9, x=1400, 
           label=expression(paste(R^{2},'= 0.11')),
           color="black", size=4)+
  labs(x="Elevation (m)", 
       y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')))+
  theme(panel.background = element_blank())

#plot(p.4b)

#### Panel 4c  U vs Chinle
InCh <- UCaus[which(UCaus$uinChF == 1),2] # data
NnCh <- UCaus[which(UCaus$uinChF == 0 & UCaus$uinMF == 0),2] # data

nInCh <- cbind(InCh, (rep('In Chinle Fm.', length(InCh))))
nNnCh <- cbind(NnCh, (rep('Not in Chinle Fm.', length(NnCh))))

Ch.U <- as.data.frame(rbind(nInCh,nNnCh))

wCh <- wilcox.test(NnCh, InCh) 
wCh #W = 2521, p = 0.002

p.4c <- ggplot() +
  geom_hline(yintercept = 1.477, color = "gray", alpha = 0.5, size = 1, 
             linetype = 2)+
  annotate("text", x= 1.5,  y = 1.7,
           label = expression(paste('30 '*italic("\u00B5")*'g/L')),
           col = "gray", size = 4)+
  geom_boxplot(aes(Ch.U$V2, as.numeric(Ch.U$InCh)), outlier.shape = NA)+
  geom_jitter(aes(Ch.U$V2, as.numeric(Ch.U$InCh)),width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_x_discrete(breaks=c("In Chinle Fm.","Not in Chinle Fm."))+
  labs(y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), x = "")+
  annotate(geom="text", x="In Chinle Fm.", y=2.9, label="*",
           color="red", size=7)+
  theme(axis.line = element_line(colour = "black", linetype = "solid"),
        axis.text.x = element_text(size = 12))+
  theme(panel.background = element_blank()) 

#plot(p.4c) # Plotting 

#### Panel 4d  U vs Morrison
InMF <- UCaus[which(UCaus$uinMF == 1),2] 
NnMF <- UCaus[which(UCaus$uinMF == 0 & UCaus$uinChF == 0),2] # data

nInMF <- cbind(InMF, (rep('In Morrison Fm.', length(InMF))))
nNnMF <- cbind(NnMF, (rep('Not in Morrison Fm.', length(NnMF))))

MF.U <- as.data.frame(rbind(nInMF,nNnMF))

wMF <- wilcox.test(NnMF, InMF)
wMF #W = 1032, p-value = 0.6857

p.4d <- ggplot() +
  geom_hline(yintercept = 1.477, color = "gray", alpha = 0.5, size = 1, 
             linetype = 2)+
  annotate("text", x= 1.5,  y = 1.7,
           label = expression(paste('30 '*italic("\u00B5")*'g/L')),
           col = "gray", size = 4)+
  geom_boxplot(aes(MF.U[,2], as.numeric(MF.U[,1])),outlier.shape = NA)+
  geom_jitter(aes(MF.U[,2], as.numeric(MF.U[,1])),width = 0.15, alpha = 0.5)+
  scale_y_continuous(breaks=c(-2,-1,0,1,2,3),
                     labels = c(expression(0.01,0.1,1,10,100,1000)),
                     limits = c(-2, 3))+
  scale_x_discrete(breaks=c("In Morrison Fm.","Not in Morrison Fm."))+
  labs(y=expression(paste('[U] ('*italic("\u00B5")*'g/L)')), x = "")+
  theme(axis.line = element_line(colour = "black", linetype = "solid"),
        axis.text.x = element_text(size = 12))+
  theme(panel.background = element_blank()) 

#plot(p.4d) # Plotting 

require("cowplot")

#export as 8 x 8
plot_grid(p.4a, p.4b, p.4c, p.4d,
          labels = c("A","B","C","D"),
          ncol = 2)
