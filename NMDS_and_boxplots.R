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

##VISUALIZATIONS
#create a stacked bar chart using the ggplot package, to examine the # each broad bug categoriy.
#But first, we have to manually reorder our variables so that they show up in the correct order on the bar chart. We will transform the table as follows:
os19_percents <- read.csv("osylv_barchart.csv")

##Now, you are ready to make your plot.
percent_plot <- ggplot(os19_percents,
                       aes(x = population, y = percent, fill = prey_type)) + 
  geom_bar(width=0.5, stat="identity") +  
  scale_fill_manual(values=c("#EB6084", "#EBDE31", "#5CD4EB", "#5DDB81", "#AAA8ED")) + 
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))

#view your plot
percent_plot



#make a boxplot to look at the spread of your data
library(reshape2)
osylvbox <- read.csv("osylv2019_rawforbox.csv")
osylv_melt<- melt(osylvbox,id.var=c('population','number','prey_type'))
boxplot_all <- ggplot(osylv_melt, aes(population, number, fill=prey_type))
boxplot_all2 <- boxplot_all + geom_boxplot() + scale_fill_manual(values=c("#EB6084", "#EBDE31", "#5CD4EB", "#5DDB81", "#AAA8ED")) + scale_y_continuous(trans = 'sqrt') + geom_point(position=position_jitterdodge(jitter.width = 0, jitter.height = 0)) + theme(axis.line = element_line(color = "black")+ theme(legend.key = element_rect(fill = "white")),
                                                                                                                                       panel.grid.major = element_blank(),
                                                                                                                                       panel.grid.minor = element_blank(),
                                                                                                                                       panel.border = element_blank(),
                                                                                                                                       panel.background = element_blank())

# view your plot!
boxplot_all2

#there are so many ants that we lose resolution of our other arthropod categories. so here is a boxplot code that links to a file in your working directory that omits ants.
#NOTE: there was an outlier in this plot that was removed. (bigger than all other points). it is because this individual ate 111 mites, which is not a mistake. for the purposes of visualization only, this is removed. it is only missing from the "osylv2019_rawforbox_NOANTS.csv" file. the main data files remain complete.

#make a boxplot to look at the spread of your data
library(reshape2)
osylvbox <- read.csv("osylv2019_rawforbox_NOANTS.csv")
osylv_melt<- melt(osylvbox,id.var=c('population','number','prey_type'))
boxplot_all <- ggplot(osylv_melt, aes(population, number, fill=prey_type))
boxplot_all2 <- boxplot_all + geom_boxplot() + scale_fill_manual(values=c("#EBDE31", "#5CD4EB", "#5DDB81", "#AAA8ED"))+ scale_y_continuous(trans = 'log2') + geom_point(position=position_jitterdodge(jitter.width = 0, jitter.height = 0)) + theme(axis.line = element_line(color = "black")+ theme(legend.key = element_rect(fill = "white")),
                                                                                                                                                                                                                           panel.grid.major = element_blank(),
                                                                                                                                                                                                                           panel.grid.minor = element_blank(),
                                                                                                                                                                                                                           panel.border = element_blank(),
                                                                                                                                                                                                                           panel.background = element_blank())

# view your plot!
boxplot_all2

##OVERALL stats
meta <- read.csv("osylv_META.csv")
osylvN <- read.csv("osylv2019_raw2.csv")

#remove everything from the sheet that isn't your diet data, but REMOVE TOTALS.
osylvN$population <- NULL
osylvN$sex <- NULL
osylvN$elevation.meters <- NULL
osylvN$frog_id <- NULL
osylvN$id_simple <- NULL
osylvN$total <- NULL
#do the ordination, you want bray-curtis for distance.
ord <- metaMDS(osylvN, perm=9999, distance = "bray", k = 3, autotransform = FALSE)

#Always run the ordination a second time starting at the previous best solution 
ord2<-metaMDS(osylvN, previous.best = TRUE)
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
permanova.output<-adonis2(osylvN~meta$population, permutations = 9999,method="bray")
summary(permanova.output)
permanova.output

#Posthoc permanova because our results were significant
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwiseAdonis)
pair.mod <-pairwise.adonis(osylvN,factors=meta$population, p.adjust = "BH")
pair.mod

#Analysis of Similarity
an<-anosim(osylvN,meta$population, permutations = 9999, distance = "bray", strata = NULL)
an
