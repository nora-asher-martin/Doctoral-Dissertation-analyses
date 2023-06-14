# GLM related packages
library(AER)
library(MASS)
library(pscl)
library(AEDForecasting)
library(CARS)
library(corrplot)
library(mctest)
library(glmmTMB)
library(bbmle) ## for AICtab
library(ggplot2)
library(TMB)
library(emmeans)

#import data
osylv <- read.csv("osylv2019_raw2.csv")

#store your factors
population <- osylv[["population"]]
elevation <- osylv[["elevation.meters"]]
ants <- osylv[["ants"]]
mites <- osylv[["mites"]]
beetles <- osylv[["beetle"]]
larvae <- osylv[["larvae"]]
other <- osylv[["other"]]
sex <- osylv[["sex"]]
frog <- osylv[["id_simple"]]
total <- osylv[["total"]]

#after looking at our data, we know that some of the prey categories are zero-inflated, but not all of them. Each model has one Y variable - a specific prey category, and considers population, the total # of prey items consumed as an offset, and frog ID as a random effect. Sex has been removed from this model because it does not predict prey consumption, ziformula=~1 signifies that these observations are considered zero-inflated by the model.
#####
#ant consumption - NOT zero-inflated
zipA <- glmmTMB(ants ~ population + offset(log(total)) + (1 | frog), family = "poisson", data = osylv)
summary(zipA)
#####
#mite consumption - NOT zero-inflated
zipM <- glmmTMB(mites ~ population + offset(log(total)) + (1 | frog), family = "poisson", data = osylv)
summary(zipM)
#####
#larvae consumption - zero-inflated dataset
zipL <- glmmTMB(larvae ~ population + offset(log(total)) + (1 | frog), family = "poisson",   ziformula=~1, data = osylv)
summary(zipL)
#####
#beetle consumption - zero-inflated dataset
zipB <- glmmTMB(beetles ~ population + offset(log(total)) + (1 | frog), family = "poisson",   ziformula=~1, data = osylv)
summary(zipB)
####
#other consumption - zero-inflated dataset, special case -negative binomial distribution, much better fit!
zipO <- glmmTMB(other ~ population + offset(log(total)) + (1 | frog), family = "nbinom2",   ziformula=~1, data = osylv)
summary(zipO)

#posthoc tests - pairwise comparison between populations.
emmeans(zipA, pairwise~population, adjust="tukey")
emmeans(zipM, pairwise~population, adjust="tukey")
emmeans(zipL, pairwise~population, adjust="tukey")
emmeans(zipB, pairwise~population, adjust="tukey")
emmeans(zipO, pairwise~population, adjust="tukey")

#plot residuals
#ants
zipA.resid <- resid(zipA)
plot(ants, zipA.resid)
#mites
zipM.resid <- resid(zipM)
plot(mites, zipM.resid)
#larvae
zipL.resid <- resid(zipL)
plot(larvae, zipL.resid)
#beetles
zipB.resid <- resid(zipB)
plot(beetles, zipB.resid)
#other
zipO.resid <- resid(zipO)
plot(other, zipO.resid)