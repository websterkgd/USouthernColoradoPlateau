#Analysis of data using natural experiment approach

rm(list=ls())
getwd()

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

#bring in multivariate matrix        
upred <- read.csv("Upred_matrix.csv")
uprd <- upred[,-1]
rm(list="upred")

#creating a matrix 
up_saldrop <- uprd[,-c(5:7, 13:15)]

#bringing in new mining vectors
mtable <- read.csv("MineDistances.csv",row.names =1)

up_sd_nm <- cbind(up_saldrop, mtable)

up_sd_nm <- na.omit(up_sd_nm)

##linear models 
AR <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$uinMF+
           up_sd_nm$mu1+up_sd_nm$mu2+up_sd_nm$u_sal+up_sd_nm$pH_m+
           up_sd_nm$am1+up_sd_nm$am2+up_sd_nm$md1+up_sd_nm$md5+up_sd_nm$dd1+
           up_sd_nm$md5+up_sd_nm$dd5)
summary(AR) 
#relationships with: El 0.01;  mu2 0.05; u_sal 1*10^-4; am2 0.02;   
#r^2 0.38; adj r2 = 0.33; p = 7*10^-13

require(car)
#### looking for colinearities
iv <- as.data.frame(cbind(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$uinMF+
                            up_sd_nm$mu1+up_sd_nm$mu2+up_sd_nm$u_sal+up_sd_nm$pH_m+
                            up_sd_nm$am1+up_sd_nm$am2+up_sd_nm$md1+up_sd_nm$md5+up_sd_nm$dd1+
                            up_sd_nm$md5+up_sd_nm$dd5)) # nor sure what running this is doing

vif(AR) # colinearity: mu1 3; mu2 9; am1 3; am2 9; md1 7; md5 10; dd1 14; dd5 17 

# what to take away is that the mining vectors fairly colinear

# will filter by removing the mining vectors that are least predictive individually
# Filtering out md1, dd5, am1, 
AR_mi_r <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$uinMF+
           up_sd_nm$mu1+up_sd_nm$mu2+up_sd_nm$u_sal+up_sd_nm$pH_m+up_sd_nm$am2+
             up_sd_nm$md5+up_sd_nm$dd1+up_sd_nm$md5)

summary(AR_mi_r)
#relationships with: El 0.009;  mu2 0.05; u_sal 7*10^-5; am2 0.02; dd1 0.004   
#r^2 0.37; adj r2 = 0.33; p = 7*10^-13

vif(AR_mi_r) #mu2 9; am2; md5 2; dd1 2

# running model removing the morrison formation
AR_mi_mr <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$mu1+
                 up_sd_nm$mu2+up_sd_nm$u_sal+up_sd_nm$pH_m+up_sd_nm$am2+
                 up_sd_nm$md5+up_sd_nm$dd1+up_sd_nm$md5)
summary(AR_mi_mr) 
# only El 9*10^-3; sal 7*10^-5; dd1 0.004;mu2 0.05; am2 0.02

vif(AR_mi_mr) #mu2 9; am2 8; md5 2; dd1 2
# even though mu2 is significant will drop since it has high colinearity

# running model removing the morrison formation + Mu2
AR_mr_2 <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$mu1+
                up_sd_nm$u_sal+up_sd_nm$pH_m+up_sd_nm$am2+
                up_sd_nm$md5+up_sd_nm$dd1)
summary(AR_mr_2) 
# El 0.02; sal 7*10^-5, dd1 0.004 r^2 0.35

vif(AR_mr_2) #md5 2; dd1 2

# running model removing the morrison formation + Mu2 + ph
AR_mr_2_pH <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$mu1+
                  up_sd_nm$u_sal+up_sd_nm$am2+up_sd_nm$md5+
                   up_sd_nm$dd1)
summary(AR_mr_2_pH) 
# El 0.03, sal 3*10^-5,  dd1 0.003; r^2 0.35

vif(AR_mr_2_pH) #md5 2; dd1 2

# running model removing the morrison formation + Mu2 + pH + md5
AR_mr_2_md5 <-  lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$mu1+
                   up_sd_nm$u_sal+up_sd_nm$am2+up_sd_nm$dd1)
summary(AR_mr_2_md5) 
# El 0.02, sal 2*10^-5, dd1 7*10^-7; am2 0.07; r^2 0.34

vif(AR_mr_2_md5) #nothing overly colinear 

# dropping mu1
AR_mr_2_1 <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$u_sal+
                  up_sd_nm$am2+up_sd_nm$dd1)
summary(AR_mr_2_1) 
# El 0.02, sal 2*10^-5, dd1 5*10^-7; am2 0.07; r^2 0.34

#dropping CH
AR_mr_2_CH <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$u_sal+up_sd_nm$am2+
                   up_sd_nm$dd1)
summary(AR_mr_2_CH) # El 0.02, sal 2*10^-5, dd1 5*10^-7; am2 0.005; r^2 0.33

# a problem with this is that even though am2 and dd1 are not able to be 
# measured as colinear they are. They measure mining. 
# so will re run model with only dd1

AR_1m <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$uinMF+
           up_sd_nm$u_sal+up_sd_nm$pH_m+up_sd_nm$dd1)
summary(AR_1m) 
# El 0.02, sal 3*10^-4, dd1 1*10^-7; r^2 0.32 

vif(AR_1m)# nothing super colinear

# chinle not important, but
ddCH <- wilcox.test(up_sd_nm$dd1[which(up_sd_nm$uinChF == 0)], 
            up_sd_nm$dd1[which(up_sd_nm$uinChF == 1)])

#distance to a mine upstream is smaller in the Chinle

# dropping pH and Morrion
AR_1m_mp <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$u_sal+
                 up_sd_nm$dd1)
summary(AR_1m_mp) 
# El 0.01, CH 0.07;  sal 1*10^-4, dd1 6*10^-9; r^2 0.32 

# What happens if Chinle Formation is dropped
AR_1m_mpc <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$u_sal+up_sd_nm$dd1)
summary(AR_1m_mpc) 
# El 0.006;  al 6*10^-5, dd1 1*10^-9; r^2 0.31 

##### Want to run the model with a more categorical sal vector
uprnna <- na.omit(uprd[,-c(13:15)])

mc_sal <- uprnna$uin01+uprnna$uin13*2+uprnna$uin325*3

AR_ns <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+up_sd_nm$uinMF+
              mc_sal+up_sd_nm$pH_m+up_sd_nm$dd1) 
summary(AR_ns)
# El 0.004, sal 0.008, Ch 0.07; dd1 = 3*10^-7 r^2 0.30

#dropping mf and ph
AR_ns_2 <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+mc_sal+up_sd_nm$dd1) 
summary(AR_ns_2)
# El 0.004, sal 0.007, Ch 0.05; dd1 = 2*10^-8 r^2 0.30

#messing around with sensitivity
v <- rep(1, length(uprnna$u_sal))
v[which(uprnna$u_sal == 14000)] <- 5000/14000

a_sal <- v*uprnna$u_sal

# dropping mf and pH
AR_as <- lm(up_sd_nm[,1] ~ up_sd_nm$uEl+up_sd_nm$uinChF+a_sal+up_sd_nm$dd1)
summary(AR_as) # Chinle sig at 4000; not sig at 7000 or 5000; 
# for 4000
# El 0.008, sal 0.001, Ch 0.05; dd1 = 1*10^-8; r^2 0.31

###is pH related to the chinle formation?
pHinCF <- na.omit(uprd$pH_m[which(uprd$uinChF == 1)])
pHnnCF <- na.omit(uprd$pH_m[which(uprd$uinChF == 0)])

median(pHinCF) #ph = 7.5
median(pHnnCF) #ph = 6.8 

mean(pHinCF) # 7.5
mean(pHnnCF) # 6.9

sd(pHinCF) #0.67
sd(pHnnCF) # 0.87

#quick plot
plot(pHnnCF)
points(pHinCF, pch = 9)  

#quick boxplot
boxplot(pHinCF, pHnnCF) #

test_pHCF <- wilcox.test(pHinCF,pHnnCF) 
test_pHCF
#p = 1*10^-5

t.test(pHinCF,pHnnCF) #p = 3.2 * 10^-5

###   pH vs lu
plot(uprd$lU ~ uprd$pH_m)
lm(uprd$lU ~ uprd$pH_m)
summary(lm(uprd$lU ~ uprd$pH_m))
abline(summary(lm(uprd$lU ~ uprd$pH_m))$coefficients[1],
      summary(lm(uprd$lU ~ uprd$pH_m))$coefficients[2],col='green')
# p = 5.5 * 10 ^-3; r^2 = 0.06

### what about dd1 - how does it relate to chinle?
ddinCF <- na.omit(up_sd_nm$dd1[which(up_sd_nm$uinChF == 1)])
ddnnCF <- na.omit(up_sd_nm$dd1[which(up_sd_nm$uinChF == 0)])

median(ddinCF) #ph = 7.7
median(ddnnCF) #ph = 9.3 - higher outside Chinle

mean(ddinCF) # 9.9
mean(ddnnCF) # 20.3

sd(ddinCF) # 9
sd(ddnnCF) # 22

#wquick boxplot
boxplot(ddinCF, ddnnCF) #

test_ddCF <- wilcox.test(ddinCF,ddnnCF) 
test_ddCF
# yes these groups are different p = 0.03

plot(up_sd_nm$dd1[which(up_sd_nm$uinChF == 1)],
     up_sd_nm$lU[which(up_sd_nm$uinChF == 1)]) # 

summary(lm(up_sd_nm$lU[which(up_sd_nm$uinChF == 1)]~
        up_sd_nm$dd1[which(up_sd_nm$uinChF == 1)])) 
           
