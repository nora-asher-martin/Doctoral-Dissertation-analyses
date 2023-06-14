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
#other consumption - zero-inflated dataset, special case - negative binomial distribution, much better fit!
zipO <- glmmTMB(other ~ population + offset(log(total)) + (1 | frog), family = "nbinom2",   ziformula=~1, data = osylv)
summary(zipO)

#significant -- PAIRWISE POSTHOC.
emmeans(zipA, pairwise~population, adjust="hochberg")



#these models do not work for predicting larvae, beetle or other consumption . it does work for ants,
#so to compare the other groups we will do a k-wallis.
kruskal.test(ants ~ population, data = osylv)
kruskal.test(mites ~ population, data = osylv)
kruskal.test(larvae ~ population, data = osylv)
kruskal.test(beetles ~ population, data = osylv)
kruskal.test(other ~ population, data = osylv)

#posthoc tests - pairwise comparison between populations per prey group. only conduct this on groups with significant KW results. DO NOT USE THIS FOR ANTS.
library(FSA)
dunnTest(ants ~ population, data = osylv, method="bh")
dunnTest(other ~ population, data = osylv, method="bh")

#plot your residuals to assess model fit.  
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

#boxplot of total # of prey consumed
totals1 <- read.csv(file.choose(), header=TRUE)
totals <- read.csv("total_boxplot.csv")
boxplot(total~population,data=totals, main="Total # consumed prey items",
        xlab="Population", ylab="# prey items consumed")


################# PLEASE IGNORE EVERYTHING BELOW HERE ###################
library(vegan)
#############
#taken from the sex differences doc, I'm trying the same analysis as above but with hellinger transformed data
#for nmds, transform it using "Hellinger's transformation" of the vegan package
osylv2 <- decostand(osylv[, 4:9], "hellinger")
write.csv(osylv2,"~\\hellinger_osylv2019.csv", row.names = FALSE)

#store your factors (IF you use the hellingers, which means you're turning your values from absolute to relative.)
population <- osylv[["population"]]
ants <- osylv2[["ants"]]
mites <- osylv2[["mites"]]
beetles <- osylv2[["beetle"]]
larvae <- osylv2[["larvae"]]
other <- osylv2[["other"]]
sex <- osylv[["sex"]]
frog <- osylv[["id_simple"]]
total <- osylv[["total"]]

#for glmm of sex, log transform the numeric variables of your data
library(dplyr)
##still working on the log transformations
#transform(mites, method = "log+1")
#transform(beetles, method = "log+1")
#transform(larvae, method = "log+1")
#transform(other, method = "log+1")
#osylv_log


#for this analysis, we are interested in sex differences across populations. sample sizes are too small or skewed towards one sex to reliably test sex differences within them 
#here, we examine sex differences across prey categories (total #arthropods consumed), using frog as a random effect.
zipsex <- glmmTMB(total ~ sex + (1 | frog), family = "poisson", data = osylv)
summary(zipsex)
#here, we examine sex differences per prey category, correcting for total consumed, and using frog as a random effect.
#####
#ant consumption - NOT zero-inflated
zipA <- glmmTMB(ants ~ population + offset(log(total)) + (1 | frog), family = "poisson", ziformula=~0, data = osylv)
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
#other consumption - zero-inflated dataset, special case - negative binomial distribution, much better fit!
zipO <- glmmTMB(other ~ population + offset(log(total)) + (1 | frog), family = "nbinom2",   ziformula=~1, data = osylv)
summary(zipO)

#posthoc tests - pairwise comparison of per category and total consumption between sexes
emmeans(zipsex, pairwise~sex, adjust="hochberg")v
emmeans(zipA, pairwise~population, adjust="hochberg")
emmeans(zipM, pairwise~population, adjust="hochberg")
emmeans(zipL, pairwise~population, adjust="hochberg")
emmeans(zipB, pairwise~population, adjust="hochberg")
emmeans(zipO, pairwise~population, adjust="hochberg")


#plot your residuals to assess model fit. 
#sex totals
zipsex.resid <- resid(zipsex)
plot (total, zipsex.resid)
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
