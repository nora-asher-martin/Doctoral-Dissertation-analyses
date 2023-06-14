#Load these packages
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)
library(factoextra)
library(plotly)
library(dendextend)
library(ggrepel)
library(ggplot2)

setwd("~/Documents/dissertation_analyses_local/Osylv_2019_local")

#import csv data files
os_hi_raw <- read.csv("os_hi_counts.csv") 
os_hi_percent <- read.csv("oshi_percents.csv") #using for barchart, only a few values here bc they are population summaries

#store your factors
group <- os_hi_raw[["group"]]
ants <- os_hi_raw[["ants"]]
mites <- os_hi_raw[["mites"]]
beetles <- os_hi_raw[["beetles"]]
larvae <- os_hi_raw[["larvae"]]
other <- os_hi_raw[["other"]]
frog_id <- os_hi_raw[["frog_id"]]
total <- os_hi_raw[["total"]]

#start with a GLM.
#after looking at our data, we know that some of the prey categories are zero-inflated. we also know there is a substantial difference in the # of prey eaten between our two species: oophaga sylvatica and hyloxalus infraguttatus. so our "total" offset is even more crucial.

#####
library(glmmTMB)
#ant consumption - zero-inflated
zipA <- glmmTMB(ants ~ group + offset(log(total)) + (1 | frog_id), family = "poisson", ziformula=~1, data = os_hi_raw)
summary(zipA)
#####
#mite consumption - zero-inflated
zipM <- glmmTMB(mites ~ group + offset(log(total)) + (1 | frog_id), family = "poisson", ziformula=~1, data = os_hi_raw)
summary(zipM)
#####
#larvae consumption - zero-inflated
zipL <- glmmTMB(larvae ~ group + offset(log(total)) + (1 | frog_id), family = "poisson",   ziformula=~1, data = os_hi_raw)
summary(zipL)
#####
#beetle consumption - zero-inflated
zipB <- glmmTMB(beetles ~ group + offset(log(total)) + (1 | frog_id), family = "poisson",   ziformula=~1, data = os_hi_raw)
summary(zipB)
####
#other consumption - zero-inflated dataset
zipO <- glmmTMB(other ~ group + offset(log(total)) + (1 | frog_id), family = "poisson",   ziformula=~1, data = os_hi_raw)
summary(zipO)

#look at residuals to assess model fit.
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

#then do stats on count data. do wilcox, not kruskal. because you only have 2 groups.
wilcox.test(ants ~ group, data = os_hi_raw)
wilcox.test(mites ~ group, data = os_hi_raw)
wilcox.test(larvae ~ group, data = os_hi_raw)
wilcox.test(beetles ~ group, data = os_hi_raw)
wilcox.test(other ~ group, data = os_hi_raw)

#define variables to be used in box plot. 
os_hi_box <- read.csv("os_hi_boxplot.csv")
group <- os_hi_box[["group"]]
prey_type <- os_hi_box[["prey_type"]]
freq <- os_hi_box[["freq"]]

#make a boxplot to look at the spread of your raw data
library(reshape2)
os_hi_melt<- melt(os_hi_box,id.var=c('group','freq','prey_type'))
boxplot_all <- ggplot(os_hi_melt, aes(group, freq, fill=prey_type))
boxplot_all2 <- boxplot_all + geom_boxplot() + scale_fill_manual(values=c("#EB6084", "#EBDE31", "#5CD4EB", "#5DDB81", "#AAA8ED")) + scale_y_sqrt() + scale_y_continuous(trans = 'sqrt') + geom_point(position=position_jitterdodge(jitter.width = 0, jitter.height = 0)) + theme(axis.line = element_line(color = "black")+ theme(legend.key = element_rect(fill = "white")),
                                                                                                                                                                                                                           panel.grid.major = element_blank(),
                                                                                                                                                                                                                           panel.grid.minor = element_blank(),
                                                                                                                                                                                                                           panel.border = element_blank(),
                                                                                                                                                                                                                           panel.background = element_blank())

# view your plot!
boxplot_all2

#barchart

#create a stacked bar chart using the ggplot package, to examine the # of times frogs interacted with large or small prey items.
#But first, we have to manually reorder our variables so that they show up in the correct order on the bar chart. We will transform the table as follows:
os_hi_percent <- within(os_hi_percent, 
                      group <- factor(group, 
                                      levels=names(sort(table(group), 
                                                        decreasing=TRUE))))
##Now, you are ready to make your plot.
percent_plot <- ggplot(os_hi_percent,
                       aes(x = group, y = percent, fill = prey_type)) + 
  geom_bar(width=0.5, stat="identity") +  
  scale_fill_manual(values=c("#994455", "#004488","#eecc66", "#ee99aa", "#6699cc")) + 
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))

#view your plot
percent_plot


### ordinations & statistics ###
os_hi<- read.csv("oshi_nmds.csv")
#remove unnecessary variables
os_hi$frog_id <- NULL
os_hi$group <- NULL
meta<- read.csv("oshi_meta.csv")
#do the ordination, you want bray-curtis for distance.
ord <- metaMDS(os_hi, perm=9999, distance = "bray", k = 3, autotransform = FALSE)

#Always run the ordination a second time starting at the previous best solution 
ord2<-metaMDS(os_hi, previous.best = TRUE)
ord2 

#creates an open plot 
#use this to create your x and y limits, within your "plot" parentheses. "ylim=c(-1.1,1.1),xlim=range(-1:1.1)" change values as needed.
plot(ord2, type="n")

#fill in the points on your plot
x<-ordipointlabel(ord2,display = "sites",cex=0.5,pch=21,col="black",bg="black")
#ordihull places outlines, polygons around your assigned groups! If they overlap a lot they're not that different. These are easy to fix in Illustrator if you export as a PDF.


#Permanova - you can write the test to include any additional parameters which might explain observed patterns in your community/data table
permanova.output<-adonis2(os_hi~meta$group, permutations = 9999,method="bray")
summary(permanova.output)
permanova.output

#Posthoc permanova because our results were significant
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwiseAdonis)
pair.mod <-pairwise.adonis(et_hi,factors=meta$group[])
pair.mod

#Analysis of Similarity
an<-anosim(et_hi,meta$group, permutations = 9999, distance = "bray", strata = NULL)
an
############################################################################################33333

##Let's do a more specific analysis using Colin's suborder ID's of bugs in our "other" category.

### ordinations & statistics ###
colin<- read.csv("colin.csv")
#remove unnecessary variables
colin$population <- NULL

meta<- read.csv("colin_meta.csv")
#do the ordination, you want bray-curtis for distance.
ord <- metaMDS(colin, perm=9999, distance = "bray", k = 3, autotransform = FALSE)

#Always run the ordination a second time starting at the previous best solution 
ord2<-metaMDS(colin, previous.best = TRUE)
ord2 
#creates an open plot 
#use this to create your x and y limits, within your "plot" parentheses. "ylim=c(-1.1,1.1),xlim=range(-1:1.1)" change values as needed.
plot(ord2, type="n")

#fill in the points on your plot
x<-ordipointlabel(ord2,display = "sites",cex=0.5,pch=21,col="black",bg="black")
#ordihull places outlines, polygons around your assigned groups! If they overlap a lot they're not that different. These are easy to fix in Illustrator if you export as a PDF.
ordihull(
  ord2,
  meta$population,
  display = "sites",
  draw = c("polygon"),
  col = c("gray0", "pink4", "midnightblue", "green3", "darkgoldenrod"),
  border = c("gray0", "pink4", "midnightblue", "green3", "darkgoldenrod"),
  lty = c(1, 1, 1, 1),
  lwd = 2.5
)

#Permanova - you can write the test to include any additional parameters which might explain observed patterns in your community/data table
permanova.output<-adonis2(colin~meta$population, permutations = 9999,method="bray")
summary(permanova.output)
permanova.output

#Posthoc permanova because our results were significant
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwiseAdonis)
pair.mod <-pairwise.adonis(colin,factors=meta$population, p.adjust = "BH")
pair.mod




