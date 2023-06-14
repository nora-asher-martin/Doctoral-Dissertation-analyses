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

et_hi<- read.csv("ET_HI_dietdata.csv")
meta<- read.csv("et_hi_meta.csv")

et_hi$group <- NULL
et_hi$frog_id <- NULL
#make a boxplot to look at the spread of your data
library(reshape2)
et_hi_box <- read.csv("et_hi_boxplot.csv")
et_hi_melt<- melt(et_hi_box,id.var=c('group','number','prey_type'))
boxplot_all <- ggplot(et_hi_melt, aes(group, number, fill=prey_type))
boxplot_all2 <- boxplot_all + geom_boxplot() + scale_fill_manual(values=c("slateblue4", "slateblue1", "skyblue1", "seashell4", "seashell3")) + theme(axis.line = element_line(color = "black")+ theme(legend.key = element_rect(fill = "white")),
                                                                                                                                                     panel.grid.major = element_blank(),
                                                                                                                                                     panel.grid.minor = element_blank(),
                                                                                                                                                     panel.border = element_blank(),
                                                                                                                                                     panel.background = element_blank())

# view your plot!
boxplot_all2

#do the ordination, you want bray-curtis for distance.
ord <- metaMDS(et_hi, perm=9999, distance = "bray", k = 3, autotransform = FALSE)

#Always run the ordination a second time starting at the previous best solution 
ord2<-metaMDS(et_hi, previous.best = TRUE)
ord2 

#creates an open plot 
#use this to create your x and y limits, within your "plot" parentheses. "ylim=c(-1.1,1.1),xlim=range(-1:1.1)" change values as needed.
plot(ord2, type="n")

#fill in the points on your plot
x<-ordipointlabel(ord2,display = "sites",cex=0.5,pch=21,col="black",bg="black")
#ordihull places outlines, polygons around your assigned groups! If they overlap a lot they're not that different. These are easy to fix in Illustrator if you export as a PDF.
 

#Permanova - you can write the test to include any additional parameters which might explain observed patterns in your community/data table
permanova.output<-adonis2(et_hi~meta$group, permutations = 9999,method="bray")
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
