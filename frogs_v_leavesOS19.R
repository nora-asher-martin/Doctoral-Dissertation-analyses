#analyzing leaf litter and frog stomach contents for differences in proportional representation of ant genera

setwd("~/Documents/dissertation_analyses_local/Osylv_2019_local")

#Load these packages
library(ade4)
library(vegan)
library(dplyr)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)
library(factoextra)
library(plotly)
library(dendextend)
library(ggrepel)
library(ggplot2)


#before we do our whole analysis, here is the code to make pie charts and stacked bar charts.

# import  data, which are the combined raw #s of ant genera shared between frogs and leaf litter.
data <- read.csv("ants_pie.csv", header = TRUE)
frog <- data[1:16,1:3]
frog_genus <- data[1:16,2]
frog_freq <- data[1:16,3]

leaf <- data[17:32,1:3]
leaf_genus <- data[17:32,2]
leaf_freq <- data[17:32,3]

#piecharts
frogpie <- ggplot(frog, aes(x="", y=frog_freq, fill=frog_genus)) +
  geom_bar(stat="identity", width=2, color="#f6f6f6") +
  coord_polar("y", start=0) + theme_void() 

frogpie #view chart

leafpie <- ggplot(leaf, aes(x="", y=leaf_freq, fill=leaf_genus)) +
  geom_bar(stat="identity", width=2, color="#f6f6f6") +
  coord_polar("y", start=0) + theme_void() 

leafpie #view chart

#barcharts
frogleafbar <- ggplot(data, aes(fill=genus, y=raw_freq, x=group)) + 
  geom_bar(position="fill", stat="identity") + theme_void() 

frogleafbar




#########
#analysis code begins here

#import file
all_ants <- read.csv("frogants_shared3.csv")

#check the file
all_ants

#subset variables. pop and group are your metadata and allants2 aree your numeric data.
pop <- all_ants[["population"]]
group <- all_ants[["group"]]
all_ants2 <- all_ants[1:169,4:19]

#check numeric variable
all_ants2

#to compare ants in stomach contents by population
frog_pop <- pop[1:58]
frog_ants <- all_ants[1:58, 4:19]
frog_phe <- frog_ants[["Pheidole"]]
frog_was <- frog_ants[["Wasmannia"]]
frog_sol <- frog_ants[["Solenopsis"]]

#to compare ants in leaf litter by population
leaf_ants <- all_ants[59:169, 4:19]
leaf_phe <- leaf_ants[["Pheidole"]]
leaf_was <- leaf_ants[["Wasmannia"]]
leaf_sol <- leaf_ants[["Solenopsis"]]
leaf_pop <- leaf_ants[["population"]]

##########################################
#BY GROUP (STOMACHS VS LEAF LITTER), SHARED ANT GENERA ONLY. ANALYSIS DOESNT WORK WITH ALL ANTS.
allants <- read.csv("frogants_shared3_no9503_v2.csv", header = TRUE)

#subset data for point colors
ants <- allants[1:168, 4:19]
group <- allants[1:168, 2]

#store as factors so you can color points later
groups<- factor(c('frog', 'leaf'))

#do the ordination           
ord<-metaMDS(ants, perm=9999, distance = "bray", k = 3, autotransform = FALSE)
#Always run the ordination a second time starting at the previous best solution 
ord2<-metaMDS(ants,previous.best = TRUE)
ord2

# Identifies colors for group assignments
colvec <- c("#85C8CF","#AAD8B5")
#make your plot -- YES THIS CODE WORKS
plot(ord2, type = "n") #displays empty ordination space
points(ord2, display = "sites", pch=16, col = c("#85C8CF","#AAD8B5") [as.numeric(groups)]) 
ordihull(ord2,group,col=c('#85C8CF','#AAD8B5', lwd = 3))
# displays site points where symbols (pch) are different management options and colour (col) are different land uses
#Permanova - you can write the equation to include any additional parameters which might explain observed patterns in your community/data table
x<-adonis2(ants~group, permutations = 9999,method="bray")
x
#time for the posthoc test
library(devtools)
library(pairwiseAdonis)
pair.mod <-pairwise.adonis(all_ants,group, p.adjust = "BH")
pair.mod
#export as a csv
write.csv(pair.mod,"~/Documents/dissertation_analyses_local/Osylv_2019_local\\posthocOSallants.csv")

#Now it's time to run GLMMs. like your other diet tests, use prey total as an offset and frog id as a random effect. if this doesn't work you need to use K-Wallis.

#import file
all_ants2 <- read.csv("frog_ants_psw1.csv")

#check the file
all_ants2

#subset variables.
frog_phe <- all_ants2[["pheidole"]]
frog_was <- all_ants2[["wasmannia"]]
frog_sol <- all_ants2[["solenopsis"]]
frog_pop <- all_ants2[["population"]]
frog_id <- all_ants2[["frog_id"]]
total <- all_ants2[["total"]]

library(glmmTMB)
#pheidole consumption - NOT zero-inflated
zip_fphe <- glmmTMB(frog_phe ~ frog_pop + offset(log(total)) + (1 | frog_id), family = "poisson", ziformula=~0, data = all_ants2)
summary(zip_fphe)
#####
#wasmannia consumption - NOT zero-inflated
zip_fwas <- glmmTMB(frog_was ~ frog_pop + offset(log(total)) + (1 | frog_id), family = "poisson", ziformula=~0, data = all_ants2)
summary(zip_fwas)
#####
#solenopsis consumption - NOT zero-inflated
zip_fsol <- glmmTMB(frog_sol ~ frog_pop + offset(log(total)) + (1 | frog_id), family = "poisson",   ziformula=~0, data = all_ants2)
summary(zip_fsol)
#####

#look at residuals to assess model fit.
#pheidole in frogs by pop
zipA.resid <- resid(zip_fphe)
plot(frog_phe, zipA.resid)
#wasmannia in frogs by pop
zipM.resid <- resid(zip_fwas)
plot(frog_was, zipM.resid)
#solenopsis in frogs by pop
zipL.resid <- resid(zip_fsol)
plot(frog_sol, zipL.resid)

#BY POPULATION. LEAF LITTER ONLY.
#import file
all_ants<- read.csv("leafantsonly.csv")

#subset data
leafants <- all_ants[1:127,4:49]
popleaf <- all_ants[1:127,3]

#check the data
head(leafants)

#subset data for point colors
c_colon_leaf <- leafants[1:20,1:46]
ceiba_leaf <- leafants[21:36,1:46]
la_mana_leaf <- leafants[37:73,1:46]
p_quito_leaf <- leafants[74:91,1:46]
s_domingo_leaf <- leafants[92:127,1:46]

#store as factors so you can color points later
leafpops<- factor(c('c_colon_leaf', 'ceiba_leaf', 'la_mana_leaf', 'p_quito_leaf', 's_domingo_leaf'))

#do the ordination           
ord<-metaMDS(leafants, perm=9999, distance = "bray", k = 2, autotransform = FALSE)
#Always run the ordination a second time starting at the previous best solution 
ord2<-metaMDS(leafants,previous.best = TRUE)
ord2

# Identifies colors for group assignments
colvec <- c("#F3B8C1","#A9939A","#E7DFE9","#B2B7D5","#596780")  
#make your plot -- YES THIS CODE WORKS
plot(ord2, type = "n") #displays empty ordination space
points(ord2, display = "sites", pch=16, col = c("#F3B8C1","#A9939A","#E7DFE9", "#B2B7D5", "#596780") [as.numeric(popleaf)]) # displays site points where symbols (pch) are different management options and colour (col) are different land uses
#ordihull places outlines around your assigned groups
ordihull(ord2,popleaf,col=c('#F3B8C1','#A9939A','#E7DFE9','#B2B7D5','#596780'),lwd=3)
#Permanova - you can write the equation to include any additional parameters which might explain observed patterns in your community/data table
x<-adonis2(leafants~popleaf, permutations = 9999,method="bray")
x
#time for the posthoc test
library(devtools)
library(pairwiseAdonis)
pair.mod <-pairwise.adonis(leafants,popleaf, p.adjust = "BH")
pair.mod
write.csv(pair.mod,"~/Documents/dissertation_analyses_local/Osylv_2019_local\\leafants_bypop.csv")
#####
#BY POPULATION. FROGS ONLY.
#import file
all_ants <- read.csv("frogantsonly_no9503_2.csv")

#subset data
frogants <- all_ants[1:57,4:19]
pop <- all_ants[1:57,3]

#subset data for point colors
c_colon_leaf <- frogants[1:10, 1:16]
ceiba_leaf <- frogants[11:20, 1:16]
la_mana_leaf <- frogants[21:37, 1:16]
p_quito_leaf <- frogants[38:48, 1:16]
s_domingo_leaf <- frogants[49:57, 1:16]

#store as factors so you can color points later
frogpops<- factor(c('c_colon_leaf', 'ceiba_leaf', 'la_mana_leaf', 'p_quito_leaf', 's_domingo_leaf'))

#do the ordination           
ord<-metaMDS(frogants, perm=9999, distance = "bray", k = 3, autotransform = FALSE)
#Always run the ordination a second time starting at the previous best solution 
ord2<-metaMDS(frogants,previous.best = TRUE)
ord2

# Identifies colors for group assignments
colvec <- c("#F3B8C1","#A9939A","#E7DFE9", "#B2B7D5", "#596780")  
#make your plot -- YES THIS CODE WORKS
plot(ord2, type = "n") #displays empty ordination space
points(ord2, display = "sites", pch=16, col = c("#F3B8C1","#A9939A","#E7DFE9", "#B2B7D5", "#596780") [as.numeric(pop)]) 
ordihull(ord2,pop,col=c("#F3B8C1","#A9939A","#E7DFE9", "#B2B7D5", "#596780", lwd = 3))# displays site points where symbols (pch) are different management options and colour (col) are different land uses
#Permanova - you can write the equation to include any additional parameters which might explain observed patterns in your community/data table
x<-adonis2(frogants~pop, permutations = 9999,method="bray")
x
#time for the posthoc test
library(devtools)
library(pairwiseAdonis)
pair.mod <-pairwise.adonis(frogants,pop, p.adjust = "BH")
pair.mod
write.csv(pair.mod,"~/Documents/dissertation_analyses_local/Osylv_2019_local\\frogants_bypop.csv")

#BY POPULATION -- FROGS VS LEAF LITTER WITHIN POPULATIONS. SHARED ANT GENERA ONLY.
#Permanova - you can write the equation to include any additional parameters which might explain observed patterns in your community/data table
#import file
all_ants <- read.csv("frogants_shared3_bypop.csv")

#check the file
all_ants

#subset variables.
c_colon <- all_ants[1:30, 4:19]
cc_group <- all_ants[1:30, 2]
#
ceiba <- all_ants[31:55, 4:19]
lc_group <- all_ants[31:55, 2]
#
la_mana <- all_ants[56:95, 4:19]
lm_group <- all_ants[56:95, 2]
#
p_quito <- all_ants[96:124, 4:19]
pq_group <- all_ants[96:124, 2]
#
s_domingo <- all_ants[125:169, 4:19]
sd_group <- all_ants[125:169, 2]

#permanova                      
x<-adonis2(c_colon~cc_group, permutations = 9999,method="bray")
x
x<-adonis2(ceiba~lc_group, permutations = 9999,method="bray")
x
x<-adonis2(la_mana~lm_group, permutations = 9999,method="bray")
x
x<-adonis2(p_quito~pq_group, permutations = 9999,method="bray")
x
x<-adonis2(s_domingo~sd_group, permutations = 9999,method="bray")
x

#COMPARISON BETWEEN 16 ANT GENERA FOUND IN ALL LEAF LITTER AND FROG SAMPLES. VISUALIZING BY POPULATION, we've done this re permanova. but now we want to visualize.

#import file
all_ants <- read.csv("allants_nmds.csv")

#subset data
all_ants3 <- all_ants[1:169,5:20]
pop <- all_ants["group_pop"]

#do the ordination
ord<-metaMDS(all_ants3, perm=9999, distance = "bray", k = 3, autotransform = FALSE)
#Always run the ordination a second time starting at the previous best solution 
ord2<-metaMDS(all_ants3,previous.best = TRUE)
ord2
#creates an open plot
plot(ord2, type="n",ylim=c(-1.1,1.1),xlim=range(-1:1.1))
#fill in the points on your plot
x<-ordipointlabel(ord2,display = "site",cex=0.85,pch=21,col="black",bg="yellow")
#ordihull places outlines around your assigned groups
ordihull(ord2,pop,col=3:16,lwd=3)

#1 wk 2 days before defense (May 18 2022), here is additional analysis requested by Lauren after handing the dissertation in to my committee. This is to go in the final dissertation and published manuscript.
#
#First we want to see if solenopsis/pheidole/wasmannia consumption differs consistently between populations.
#
#upload shared ants file only contining pheidole, solenopsis and wasmannia values.
psw_shared <- read.csv("frogants_shared4.csv")
group <- psw_shared[, 2]
sol <- psw_shared[, 5]
was <- psw_shared[, 6]
phe <- psw_shared[, 4]

#then do stats on count data between frog stomachs and leaf litter group. do wilcox, not kruskal. because you only have 2 groups. you do not have enough data to create a mixed model.

#whole comparison - population not considered
#then do stats on count data. do wilcox, not kruskal. because you only have 2 groups.
wilcox.test(sol ~ group, data = psw_shared)
wilcox.test(was ~ group, data = psw_shared)
wilcox.test(phe ~ group, data = psw_shared)

#Create variables for populations frog/leaf groups.
cc_group <- psw_shared[1:30, 2]
ceiba_group <- psw_shared[31:55, 2]
mana_group <- psw_shared[56:95, 2]
pq_group <- psw_shared[96:124, 2]
sd_group <- psw_shared[125:169, 2]

#solenopsis values within pops
sol_cc <- psw_shared[1:30, 5]
sol_ceiba <- psw_shared[31:55, 5]
sol_mana <- psw_shared[56:95, 5]
sol_pq <- psw_shared[96:124, 5]
sol_sd <- psw_shared[125:169, 5]
#wasmannia values within pops
was_cc <- psw_shared[1:30, 6]
was_ceiba <- psw_shared[31:55, 6]
was_mana <- psw_shared[56:95, 6]
was_pq <- psw_shared[96:124, 6]
was_sd <- psw_shared[125:169, 6]
#pheidole values within pops
phe_cc <- psw_shared[1:30, 4]
phe_ceiba <- psw_shared[31:55, 4]
phe_mana <- psw_shared[56:95, 4]
phe_pq <- psw_shared[96:124, 4]
phe_sd <- psw_shared[125:169, 4]

#now we test for group differences (frog stomachs vs leaf litter prevalance)
#Cristobal Colon
wilcox.test(sol_cc ~ cc_group, data = psw_shared)
wilcox.test(was_cc ~ cc_group, data = psw_shared)
wilcox.test(phe_cc ~ cc_group, data = psw_shared)
#Ceiba
wilcox.test(sol_ceiba ~ ceiba_group, data = psw_shared)
wilcox.test(was_ceiba ~ ceiba_group, data = psw_shared)
wilcox.test(phe_ceiba ~ ceiba_group, data = psw_shared)
#La Mana
wilcox.test(sol_mana ~ mana_group, data = psw_shared)
wilcox.test(was_mana ~ mana_group, data = psw_shared)
wilcox.test(phe_mana ~ mana_group, data = psw_shared)
#Puerto Quito
wilcox.test(sol_pq ~ pq_group, data = psw_shared)
wilcox.test(was_pq ~ pq_group, data = psw_shared)
wilcox.test(phe_pq ~ pq_group, data = psw_shared)
#Santo Domingo
wilcox.test(sol_sd ~ sd_group, data = psw_shared)
wilcox.test(was_sd ~ sd_group, data = psw_shared)
wilcox.test(phe_sd ~ sd_group, data = psw_shared)

#Do your p adjustments. put values in ascending order, separate adjustments for each genus.
solp <- c(5.832e-06, 2.014e-05, 0.002029, 0.002989, 0.01274)
phep <-c(0.001211, 0.00389, 0.06053, 0.1248, 0.76)
wasp <- c(8.18e-05, 0.0001553, 0.01102, 0.1184, 0.2602)
allp <- c(5.83E-06,
          2.01E-05,
          8.18E-05,
          0.0001553,
          0.001211,
          0.002029,
          0.002989,
          0.00389,
          0.01102,
          0.01274,
          0.06053,
          0.1184,
          0.1248,
          0.2602,
          0.76)

p.adjust(solp,method="BH")
p.adjust(phep,method="BH")
p.adjust(wasp,method="BH")

#in case they should all be adjusted together.
p.adjust(allp,method="BH")

#p.adjusting our overall stats, the ones not considering populations.
overallp <-c(9.671e-15, 0.0002742, 0.001485)
p.adjust(overallp,method="BH")


