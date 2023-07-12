#Exploration of radiation values in "main" data set

rm(list=ls())
getwd()

require(ggplot2)

# pulling in files with kriged U values
WWR <- read.table("apU_Krige.txt", sep = "\t") 
colnames(WWR) <- WWR[1,]
WWR <- WWR[-1,]

for (i in 1:length(WWR)){
  WWR[,i] <- as.numeric(WWR[,i])
}

#pulling in mines with U data
ApUmine <- read.table('mines_kriged_U.txt')
colnames(ApUmine) <- Rmine[1,]
ApUmine <- ApUmine[-1,]

for (i in 1:length(ApUmine)){
  ApUmine[,i] <- as.numeric(ApUmine[,i])
}

# read in major data set

Upred  <- read.csv("Upred_matrix.csv")

mtable <- read.csv("MineDistances.csv",row.names =1)

Upred <- cbind(Upred,mtable$dd1)

### aligning data
Upredf <- Upred[which(Upred$Average.Latitude %in% WWR$y),]

Upredf$Average.Latitude - WWR$y #checks for difference in latitude

Upredf <- cbind(Upredf,WWR$var1.pred)

# Explorations examining waters within in the chinle formation
InCH <- Upredf[which(Upredf$uinChF > 0),] #length = 47

plot(InCH$`WWR$var1.pred`, InCH$uEl) 
summary(lm(InCH$uEl ~ InCH$`WWR$var1.pred`)) # p = 0.3

#looking at mines upstream
w.AUCh <- wilcox.test(InCH[,19], 
                      InCH[which(InCH[,11] == 0),19]) 
w.AUCh
# p = 0.83 #w = 595.5 
# ApU is not to higher at sites with a mine upstream  compared to whole dataset

#### Is U higher at water sources in the chinle?  
median(na.omit(Upredf[which(Upredf[,9] > 0),19])) #2.65
median(na.omit(Upredf[which(Upredf[,9] == 0),19])) #2.03
w.AUnC <- wilcox.test(Upredf[which(Upredf[,9] > 0),19], 
                      Upredf[which(Upredf[,9] == 0),19]) 
w.AUnC
#p = 7*10^-5, w = 3516 
# ApU appears to be higher at water sources in the chinle 

#quick visualization
boxplot(Upredf[which(Upredf[,9] == 0),18], Upredf[which(Upredf[,9] > 0),18])

#Is ApU higher in the whole dataset at locations where a 
#mine is present upstream compared to areas where not
w.AUaM <- wilcox.test(Upredf[which(Upredf[,11] > 0),18], 
                      Upredf[which(Upredf[,11] == 0),18]) 
w.AUaM
#p = 0.005, w = 1357, sig. Ap U higher where a mine is located upstream

#quick plot
boxplot(Upredf[which(Upredf[,11] == 0),18], Upredf[which(Upredf[,11] > 0),18])

#relationships between lU and ApU?
cor.test(x=Upredf$`WWR$var1.pred`, y=Upredf$lU, method = 'spearman')
#p 4.6 * 10^-5 rho = 0.3

summary(lm(Upredf$lU~Upredf$`WWR$var1.pred`))
#p 3.4 * 10 ^-5, m = 0.23, r^2 = 0.09

#quick plot rad in morrison
boxplot(Upredf[which(Upredf[,10] == 0),18], Upredf[which(Upredf[,10] > 0),18])

w.AU_M <- wilcox.test(Upredf[which(Upredf[,10] > 0),18], 
                      Upredf[which(Upredf[,10] == 0),18]) 
w.AU_M # p = 0.15 not higher in morrison

#ApU vs distance to nearest mine?
plot(Upredf[,18], Upredf[,19])

summary(lm(Upredf[,19] ~ Upredf[,18]))
# m -0.19, p 1*10^-6, r2 = 0.13 

m_sal <- Upredf$uin01+Upredf$uin13*2+Upredf$uin325*3

#multiple regression analysis with rad (additive)
AR_mrr <- lm(Upredf$lU ~ Upredf$uEl+m_sal+Upredf$`mtable$dd1`+Upredf$uinChF+
               Upredf$`WWR$var1.pred`)
summary(AR_mrr)
#r^2 = 32

#dropping chinle
AR_mr <- lm(Upredf$lU ~ Upredf$uEl+m_sal+Upredf$`mtable$dd1`+
              Upredf$`WWR$var1.pred`)
summary(AR_mr) # 
#r^2 = 31


### looking within the Chinle
Ch_mr <- lm(InCH$lU ~ InCH$uEl+m_sal[which(Upredf$uinChF == 1)]+
              InCH$`mtable$dd1`+ InCH$`WWR$var1.pred`)
summary(Ch_mr) #r^2 = 0.44 ## 

#dropping dd1
Ch_mr_d <- lm(InCH$lU ~ InCH$uEl+m_sal[which(Upredf$uinChF == 1)]+
              InCH$`WWR$var1.pred`)
summary(Ch_mr_d) #r^2 = 0.43 

Ch_mr_du <- lm(InCH$lU ~ InCH$uEl+m_sal[which(Upredf$uinChF == 1)])
summary(Ch_mr_du) #r^2 = 0.47    

