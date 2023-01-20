#Submitted by Melissa Head (melissa.head@noaa.gov) 

##Very basic glm maturity analysis###

data <- #######put data set in here#######
maturityglm <- glm (maturity ~ 1 + length, data <-data.frame(length = data$Length, maturity <- data$Maturity),
                   family = binomial(link ="logit"))
###data$Length and data$Maturity are just calling length and maturity columns in your dataset. So just make sure you column names match these or change naming in script##
####You can change Length to Age####
summary (maturityglm) ###give you A, B, SA, SB, and n ###
cor(data$Length, data$Maturity) ###gives you r##

A=
B= 
sA<
sB<-
r <-
n <- 

deltamethod <- ((sA^2)/(B^2))- ((2*A*sA*sB*r)/(B^3))+ (((A^2)*(sB^2))/(B^4))
deltamethod


1.96*(sqrt(deltamethod)/sqrt(n)) ###Gives you 95% CI##

-A/B ###Give you L50###

