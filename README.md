# HbilisASFmodel
Raw data and R scripts for community structure analysis of the ASF microbiota.
---
title: "ASF community structure analysis"
author: "Joao Carlos Gomes Neto"
date: "2/22/2017"
output: html_document
---

```{r, echo=TRUE}
library(permute)
library(vegan)
library(lattice)
library(ggplot2)
library(ggfortify)
library(scales)
setwd("~/Documents/pathobiont paper/communitystructure")
```


Does the order of colonization with the pathobiont and ASF community alter the ASF composition?
all data comes from C3H/HeN mice 
```{r, echo=TRUE}
data2<-read.csv("data2cont.csv")
sum(is.na(data2)) #checking for NA
str(data2)
head(data2)
tail(data2)
data2a<-na.omit(data2)
str(data2a)
sum(is.na(data2a))
table(data2a$groups)
df2 <- data2a[c(5:12)]
str(df2)

#PCA analysis 
pca2<-prcomp(df2, scale= TRUE) 
summary(pca2)
pca2
plot(pca2, type='l', main="Scree plot")


p5<-autoplot(prcomp(df2, scale=TRUE), data = data2a, colour = 'asf', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (41.21%)", ylab="PC2 (19.22%)")
p6<-autoplot(prcomp(df2, scale=TRUE), data = data2a, colour = 'dss', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (41.21%)", ylab="PC2 (19.22%)")
p7<-autoplot(prcomp(df2, scale=TRUE), data = data2a, colour = 'pathobiont', loadings= TRUE, loadings.label = TRUE, frame= TRUE, frame.type= 't', xlab="PC1 (41.21%)", ylab="PC2 (19.22%)")
p8<-autoplot(prcomp(df2, scale=TRUE), data = data2a, colour = 'groups', loadings= FALSE, loadings.label = FALSE, frame = TRUE, frame.type= 't', xlab="PC1 (41.21%)", ylab="PC2 (19.22%)", label.size = 50, size = 3, label.n = 30)+
theme(axis.text=element_text(size=20, family="Times New Roman"),
        axis.title=element_text(size=20,face="bold", family="Times New Roman"),
        legend.text=element_text(size=16, family="Times New Roman"),
        legend.title=element_text(size=16,face="bold", family="Times New Roman"))


p8a
new('ggmultiplot', plots = list(p5, p6, p7, p8))
#p5
#p6
#p7
#p8

#dataframe from wide to long format for relative abundance plot
library(tidyr)
datalong2a <- gather(data2a, taxa, abundance, ASF_356:ASF_519, factor_key=TRUE)
head(datalong2a)
str(datalong2a)

# calculate proportions
library(plyr)
require(plyr)
datalong2b <- ddply(datalong2a, .(groups, id, taxa))
head(datalong2b)

datalong2c<-ddply(datalong2b, .(groups, id), transform, total=sum(abundance))
head(datalong2c)

datalong2c$proportion<-datalong2c$abundance/datalong2c$total
head(datalong2c, 20)
testnew<-datalong2c[1:8,]
testnew
sum(testnew$proportion)
testnew1<-datalong2c[9:16,]
testnew1
sum(testnew1$proportion)

#relative abundances for each ASF considering the total number of bacteria per gram of cecal content
require(ggplot2)
p <- ggplot(data = datalong2c, aes(x=groups, y=proportion)) 
p <- p + geom_boxplot(aes(fill = taxa))
p <- p + geom_point(aes(y=proportion, group=taxa), position = position_dodge(width=0.75))
p <- p + facet_wrap( ~ taxa, scales="fixed", ncol=2)
p <- p + xlab("Groups") + ylab("Proportion") 
p <- p + guides(fill=guide_legend(title="ASF Taxa"))
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

#Permanova analysis to decompose the variance and explain the ASF community structure
set.seed(4490)
bc2<-vegdist(df2, method="bray", binary=FALSE, na.rm=TRUE) 
#bc
set.seed(1041)
ada2=adonis(bc2~groups, data=data2a, permutations=999)
ada2
set.seed(1233)
ad2=adonis(bc2~asf*dss*pathobiont, data=data2a, permutations=999)
ad2

#Analysis of homogeneity of the Bray-curtis distances using the betadisper function 
#objective :analysis of multivariate homogeneity of group dispersions (variances). betadisper is a multivariate analogue of Levene's test for homogeneity of variances

mod2a <- with(data2a, betadisper(bc2, groups))
mod2a
plot(mod2a, main= "All Groups")
op <- par(mar = c(4,4,4,4) + 0.1)
boxplot(mod2a, cex.axis = 0.8, main= "All Groups", ylim= c(0,1), col=terrain.colors(2))
par(op)
anova(mod2a)
permutest(mod2a)
TukeyHSD(mod2a)# provide the p-values for group comparisons 

mod2b <- with(data2a, betadisper(bc2, asf))
mod2b
plot(mod2b, main= "Timing of ASF colonization")
boxplot(mod2b, main= "Timing of ASF colonization", ylim= c(0,1))
anova(mod2b)
permutest(mod2b)
TukeyHSD(mod2b)# provide the p-values for group comparisons 

mod2c <- with(data2a, betadisper(bc2, dss))
mod2c
plot(mod2c, main= "DSS status")
boxplot(mod2c, main= "DSS status", ylim= c(0,1))
anova(mod2c)
permutest(mod2c)
TukeyHSD(mod2c)# provide the p-values for group comparisons 

mod2d <- with(data2a, betadisper(bc2, pathobiont))
mod2d
plot(mod2d, main= "Timing of Pathobiont colonization")
boxplot(mod2d, main= "Timing of Pathobiont colonization", ylim= c(0,1))
anova(mod2d)
permutest(mod2d)
TukeyHSD(mod2d)# provide the p-values for group comparisons 

# ANOSIM 
library(vegan)
#dataframe from wide to long format for relative abundance plot
library(tidyr)

anosim11<- vegdist(df2,method="bray", binary=FALSE)
attach(datalong2a)
anosim12 <- anosim(anosim11, groups)
summary(anosim12)
plot(anosim12)
```

Does the parental generation (G1 vs G2) alter the ASF composition?

dataset comes from C3H/HeN mice colonized with H. bilis and ASF from birth
```{r, echo=TRUE}

data3<-read.csv("data3cont.csv")
sum(is.na(data3)) #checking for NA
str(data3)
head(data3)
tail(data3)
data3a<-na.omit(data3)
str(data3a)
sum(is.na(data3a))
data3a$generation=factor(data3a$generation)
str(data3a)
table(data3a$groups)
#PCA analysis
df3 <- data3a[c(6:13)]
str(df3)
pca3<-prcomp(df3, scale= TRUE)  #cor=TRUE
summary(pca3)
pca3
plot(pca3, type='l', main = "Scree plot")


p9<-autoplot(prcomp(df3, scale= TRUE), data = data3a, colour = 'generation', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (39.77%)", ylab="PC2 (21.4%)")
p10<-autoplot(prcomp(df3, scale = TRUE), data = data3a, colour = 'dss', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (39.77%)", ylab="PC2 (21.4%)")
p11<-autoplot(prcomp(df3, scale = TRUE), data = data3a, colour = 'groups', loadings= FALSE, loadings.label = FALSE, frame = TRUE, frame.type= 't', xlab="PC1 (39.77%)", ylab="PC2 (21.4%)",
label.size = 50, size = 3, label.n = 30) +
theme(axis.text=element_text(size=20, family="Times New Roman"),
        axis.title=element_text(size=20,face="bold", family="Times New Roman"),
        legend.text=element_text(size=16, family="Times New Roman"),
        legend.title=element_text(size=16,face="bold", family="Times New Roman"))



p11a<-p11+theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16,face="bold"))
p11a
#p11+theme_bw()
new('ggmultiplot', plots = list(p9, p10, p11))
new('ggmultiplot', plots = list(p8a, p11a))

#p9
#p10
#p11
#p11a<-p11 + theme(panel.background = element_rect(fill='white', colour='black'))
#p11a + theme(panel.grid.major = element_rect(color='black')
           #  , panel.grid.minor = element_blank())

#dataframe from wide to long format for relative abundance plot
library(tidyr)
datalong3a <- gather(data3a, taxa, abundance, ASF_356:ASF_519, factor_key=TRUE)
head(datalong3a)
str(datalong3a)

# calculate proportions
library(plyr)
require(plyr)
datalong3b <- ddply(datalong3a, .(groups, id, taxa))
head(datalong3b)

datalong3c<-ddply(datalong3b, .(groups, id), transform, total=sum(abundance))
#head(datalong3c)

datalong3c$proportion<-datalong3c$abundance/datalong3c$total
head(datalong3c, 20)
testnew3<-datalong3c[1:8,]
testnew3
sum(testnew3$proportion)
testnew3<-datalong3[9:16,]
testnew3
sum(testnew3$proportion)

#relative abundances for each ASF considering the total number of bacteria per gram of cecal content
require(ggplot2)
p <- ggplot(data = datalong3c, aes(x=groups, y=proportion)) 
p <- p + geom_boxplot(aes(fill = taxa))
p <- p + geom_point(aes(y=proportion, group=taxa), position = position_dodge(width=0.75))
p <- p + facet_wrap( ~ taxa, scales="fixed", ncol=2)
p <- p + xlab("Groups") + ylab("Proportion") 
p <- p + guides(fill=guide_legend(title="ASF Taxa"))
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

#Permanova analysis to decompose the variance and explain the ASF community structure
set.seed(445)
bc3<-vegdist(df3, method="bray", binary=FALSE, na.rm=TRUE) 
#bc
set.seed(1045)
ada4=adonis(bc3~groups, data=data3a, permutations=999)
ada4
set.seed(1443)
ad4=adonis(bc3~generation*dss, data=data3a, permutations=999)
ad4
set.seed(1442)
ad5=adonis(bc3~generation+dss, data=data3a, permutations=999)
ad5
anosim(df3, groups,permutations = 999, distance = "bray")
set.see(569)

library(vegan)
data(dune)
data(dune.env)
dune.dist <- vegdist(df3)
attach(datalong3a)
dune.ano <- anosim(dune.dist, groups)
summary(dune.ano)
par(mar = c(7,7,7,7))
plot(dune.ano, col=terrain.colors(5), horiz=FALSE, ylab="Dissimilarity Ranks", axes= FALSE)
box()
axis(2,
     las=2)
names1<-c("Between", "HbASFbirth1", "HbASFbirth1DSS","HbASFbirth2",
          "HbASFbirth2DSS")
axis(1,at=1:5, labels=names1,
     las=2, cex.axis=0.8) 
#Analysis of homogeneity of the Bray-curtis distances using the betadisper function 
#objective :analysis of multivariate homogeneity of group dispersions (variances). betadisper is a multivariate analogue of Levene's test for homogeneity of variances

mod3a <- with(data3a, betadisper(bc3, groups), type=c("centroid"))
mod3a
plot(mod3a, main= "All Groups")
op <- par(mar = c(4,4,4,4) + 0.1)
boxplot(mod3a, cex.axis = 0.8, main= "All Groups", ylim= c(0,1))
par(op)
anova(mod3a)
permutest(mod3a)
TukeyHSD(mod3a)# provide the p-values for group comparisons 

mod3b <- with(data3a, betadisper(bc3, generation))
mod3b
plot(mod3b, main= "Generation")
boxplot(mod3b, main= "Generation", ylim= c(0,1))
anova(mod3b)
permutest(mod3b)
TukeyHSD(mod3b)# provide the p-values for group comparisons 

mod3c <- with(data3a, betadisper(bc3, dss))
mod3c
plot(mod3c, main= "DSS status")
boxplot(mod3c, main= "DSS status", ylim= c(0,1))
anova(mod3c)
anova(mod3c)
permutest(mod3c)
TukeyHSD(mod3c)# provide the p-values for group comparisons 
```

ASF community structure analysis for C57BL/6 and Rag1-/- experiments 
ASF vs H. bilis-ASF groups treated or not with 2% DSS

```{r, echo=TRU}
datanew1<-read.csv("data1cont.csv")
sum(is.na(datanew1)) #checking for NA
str(datanew1)
head(datanew1)
tail(datanew1)
#C57BL/6

#1) need to subset the dataset 
datab6<-datanew1[c(45:101),]
datab6
str(datab6)
sub2 <- datab6[c(5:12)]
str(sub2)
pcab6<-prcomp(sub2, scale= TRUE) 
summary(pcab6)
pcab6
plot(pcab6, type='l', main="Scree plot")
plot1a<-autoplot(prcomp(sub2, scale=TRUE), data = datab6, colour = 'groups', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (35.79%)", ylab="PC2 (21.93%)")
plot2a<-autoplot(prcomp(sub2, scale=TRUE), data = datab6, colour = 'dss', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (35.79%)", ylab="PC2 (21.93%)")
plot3a<-autoplot(prcomp(sub2, scale=TRUE), data = datab6, colour = 'pathobiont', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (35.79%)", ylab="PC2 (21.93%)")
new('ggmultiplot', plots = list(plot1a, plot2a, plot3a))

#Permanova analysis to decompose the variance and explain the ASF community structure
set.seed(3495)
b6a<-vegdist(sub2, method="bray", binary=FALSE, na.rm=TRUE) 
#bc
set.seed(0948)
mod6a=adonis(b6a~groups, data=datab6, permutations=999)
mod6a
set.seed(5748)
mod6b=adonis(b6a~dss*pathobiont, data=datab6, permutations=999)
mod6b


# ANOSIM 
library(vegan)
#dataframe from wide to long format for relative abundance plot
library(tidyr)
datalongb6 <- gather(datab6, taxa, abundance, ASF_356:ASF_519, factor_key=TRUE)
head(datalongb6)
str(datalongb6)
anosimb6a <- vegdist(sub2,method="bray", binary=FALSE)
attach(datalongb6)
anosimb6b <- anosim(anosimb6a, groups)
summary(anosimb6b)

#BC variance
set.seed(4757)
bc.1a<-vegdist(sub2, method="bray", binary=FALSE, na.rm=TRUE)
mod.fit1a<- with(datab6, betadisper(bc.1a, groups))
mod.fit1a
plot(mod.fit1a, main= "groups")
boxplot(mod.fit1a, main= "groups", ylim= c(0,1))
anova(mod.fit1a)
permutest(mod.fit1a)
TukeyHSD(mod.fit1a)# provide the p-values for group comparisons 
##########################################################
#RAG1-/-
#1) need to subset the dataset 
datarag<-datanew1[c(102:142),]
datarag
rag1 <- datarag[c(5:12)]
str(rag1)
pcarag<-prcomp(rag1, scale= TRUE) 
summary(pcarag)
pcarag
plot(pcarag, type='l', main="Scree plot")
plot1b<-autoplot(prcomp(rag1, scale=TRUE), data = datarag, colour = 'groups', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (31.99%)", ylab="PC2 (21.24%)")
plot2b<-autoplot(prcomp(rag1, scale=TRUE), data =datarag, colour = 'dss', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (31.99%)", ylab="PC2 (21.24%)")
plot3b<-autoplot(prcomp(rag1, scale=TRUE), data = datarag, colour = 'pathobiont', loadings= TRUE, loadings.label = TRUE, frame = TRUE, frame.type= 't', xlab="PC1 (31.99%)", ylab="PC2 (21.24%)")
new('ggmultiplot', plots = list(plot1b, plot2b, plot3b))

#Permanova analysis to decompose the variance and explain the ASF community structure
set.seed(3567)
bc1a1<-vegdist(rag1, method="bray", binary=FALSE, na.rm=TRUE) 
#bc
set.seed(9056)
ada11=adonis(bc1a1~groups, data=datarag, permutations=999)
ada11
set.seed(5894)
ad12=adonis(bc1a1~dss*pathobiont, data=datarag, permutations=999)
ad12


# ANOSIM 
library(vegan)
#dataframe from wide to long format for relative abundance plot
library(tidyr)
datalongrag <- gather(datarag, taxa, abundance, ASF_356:ASF_519, factor_key=TRUE)
head(datalongrag)
str(datalongrag)
anosimrag <- vegdist(rag1,method="bray", binary=FALSE)
attach(datalongrag)
anosimrag1 <- anosim(anosimrag, groups)
summary(anosimrag1)

#BC variance
set.seed(4563)
bc.rag<-vegdist(rag1, method="bray", binary=FALSE, na.rm=TRUE)
mod.fit1bb <- with(datarag, betadisper(bc.rag, groups))
mod.fit1bb
plot(mod.fit1bb, main= "groups")
boxplot(mod.fit1bb, main= "groups", ylim= c(0,1))
anova(mod.fit1bb)
permutest(mod.fit1bb)
TukeyHSD(mod.fit1bb)# provide the p-values for group comparisons 
```

Correlation analysis of pathobiont and ASF members abundance in cecal luminal content and tissue-associated bacteria in control vs diseased C3H/HeN animals

```{r, echo=TRUE}
data4<-read.csv("data4cont.csv")
sum(is.na(data4)) #checking for NA
str(data4)
head(data4, 3)
tail(data4, 3)
data4a<-na.omit(data4)
str(data4a)
sum(is.na(data4a))
data4a=data4[,-c(1:5)]
str(data4a)
plot(data4a)

#subset C3H

data4b=data4[c(1:11),]
str(data4b)
data4b
sum(is.na(data4b))
data4b=na.omit(data4b)
sum(is.na(data4b))
datasub1a=data4b[,-c(1:5)]
str(datasub1a)
plot(datasub1a)

#Spearman correlation and heat map C3H

cormat1a<-round(cor(datasub1a, method="spearman", use = "complete.obs"),2)
heatmap(cormat1a)
head(cormat1a)
plot(pathobiont~ASF_356, data=datasub1a)
library(reshape2)
melted_cormat1a <- melt(cormat1a)
head(melted_cormat1a)

library(ggplot2)
ggplot(data = melted_cormat1a, aes(x=X1, y=X2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri1a<- get_upper_tri(cormat1a)
upper_tri1a


# Melt the correlation matrix
library(reshape2)
melted_cormat1b <- melt(upper_tri1a, na.rm = TRUE)
melted_cormat1b 
melted_cormat1b<-melted_cormat1b[order(melted_cormat1b$X1), ]
melted_cormat1b

#subseting from correlation matrix
head(melted_cormat1b, 10)
melted_cormat1c=melted_cormat1b[1:9,]
melted_cormat1c

# Heatmap
library(ggplot2)
ggheatmap1<-ggplot(data = melted_cormat1c, aes(X2, X1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman Correlation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 20, hjust = 0.7))+
  theme(axis.text.y = element_text(size=20))+
 coord_fixed()

#print(ggheatmap)

ggheatmap2<-ggheatmap1 + 
geom_text(aes(X2, X1, label = value), color = "black", size = 12) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(1, 1),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 12, barheight = 0.8,
                title.position = "top", title.hjust = 0.5,
                title.theme= element_text(size = 20, angle= 0)))

ggheatmap2 + ggtitle(expression(paste(italic('H. bilis'),"(pathobiont)-ASF-1.5% DSS")))

#subset HBASFnoDSS C3H/HeN
data5a=data4[c(42:51),]
str(data5a)
data5a
sum(is.na(data5a))
data4d=na.omit(data5a)
sum(is.na(data5a))
datasub5a=data5a[,-c(1:5)]
str(datasub5a)
plot(datasub5a)

#Pearson correlation and heat map C3H

cormat5a<-round(cor(datasub5a, method="spearman", use = "complete.obs"),2)
heatmap(cormat5a)
head(cormat5a)

library(reshape2)
melted_cormat5a <- melt(cormat5a)
head(melted_cormat5a)

library(ggplot2)
ggplot(data = melted_cormat5a, aes(x=X1, y=X2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri5a<- get_upper_tri(cormat5a)
upper_tri5a


# Melt the correlation matrix
library(reshape2)
melted_cormat5b <- melt(upper_tri5a, na.rm = TRUE)
melted_cormat5b 
melted_cormat5b<-melted_cormat5b[order(melted_cormat5b$X1), ]
melted_cormat5b

#subseting from correlation matrix
head(melted_cormat5b, 10)
melted_cormat5c=melted_cormat5b[1:9,]
melted_cormat5c

# Heatmap
library(ggplot2)
ggheatmap5a<-ggplot(data = melted_cormat5c, aes(X2, X1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman Correlation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 20, hjust = 0.7))+
  theme(axis.text.y = element_text(size=20))+
 coord_fixed()

#print(ggheatmap)

ggheatmap5b<-ggheatmap5a + 
geom_text(aes(X2, X1, label = value), color = "black", size = 12) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(1, 1),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 12, barheight = 0.8,
                title.position = "top", title.hjust = 0.5,   title.theme= element_text(size = 20, angle= 0)))


 

a<-ggheatmap5b + ggtitle(expression(paste(italic('H. bilis'),"(pathobiont)-ASF-no DSS")))+ 
      labs(x="C3H/HeN - Cecal content bacteria")+ 
      theme(plot.title = element_text(size = 22))+
      theme(axis.title.x = element_text(size = 22))


b<-ggheatmap2 + ggtitle(expression(paste(italic('H. bilis'),"(pathobiont)-ASF-1.5% DSS")))+ 
      labs(x="C3H/HeN - Cecal content bacteria")+ 
      theme(plot.title = element_text(size = 22))+
      theme(axis.title.x = element_text(size = 22))

#plist<-list(a,b)
#library(gridExtra)
#n <- length(plist)
#nCol <- floor(sqrt(n))
#do.call("grid.arrange", c(plist, ncol=nCol))
library(gridExtra)
grid.arrange(a,b, nrow=2, ncol=1, clip=TRUE)

#####################################################################
#cecal tissue associated bacteria
data5<-read.csv("data5cont.csv")
sum(is.na(data5)) #checking for NA
str(data5)


#subset HbASF no DSS
head(data5, 9)
datasub6a=data5[c(1:8),]
datasub6a
sum(is.na(datasub6a))
str(datasub6a)
#now subset the quantitative variables out of HbASF
datasub6b=datasub6a[,-c(1:5)]
str(datasub6b)
plot(datasub6b)

#Pearson correlation and heat map HbASF no DSS

cormat6a<-round(cor(datasub6b, method="spearman", use = "complete.obs"),2)
heatmap(cormat6a)
head(cormat6a)

library(reshape2)
melted_cormat6a <- melt(cormat6a)
head(melted_cormat6a)

library(ggplot2)
ggplot(data = melted_cormat6a, aes(x=X1, y=X2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri6a<- get_upper_tri(cormat6a)
upper_tri6a


# Melt the correlation matrix
library(reshape2)
melted_cormat6b <- melt(upper_tri6a, na.rm = TRUE)
melted_cormat6b 
melted_cormat6b<-melted_cormat6b[order(melted_cormat6b$X1), ]
melted_cormat6b

#subseting from correlation matrix
head(melted_cormat6b, 10)
melted_cormat6c=melted_cormat6b[1:8,]
melted_cormat6c

# Heatmap
library(ggplot2)
ggheatmap6a<-ggplot(data = melted_cormat6c, aes(X2, X1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman Correlation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 20, hjust = 0.7))+
  theme(axis.text.y = element_text(size=20))+
 coord_fixed()

#print(ggheatmap)

ggheatmap6b<-ggheatmap6a + 
geom_text(aes(X2, X1, label = value), color = "black", size = 12) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(1,1),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 12, barheight = 0.8,
                title.position = "top", title.hjust = 0.5,
                title.theme= element_text(size = 20, angle= 0)))

a1<-ggheatmap6b + ggtitle(expression(paste(italic('H. bilis'),"(pathobiont)-ASF-no DSS")))+ 
      labs(x="C3H/HeN - Cecal tissue associated bacteria")+ 
      theme(plot.title = element_text(size = 22))+
      theme(axis.title.x = element_text(size = 22))
a1


######################################################################
#subseting HbASF DSS for correlation matrix between pathobiont and ASF taxa in cecal tissue associated bacteria

str(data5)


#subset HbASF no DSS
data5
datasub7b=data5[c(9:17),]
datasub7b
sum(is.na(datasub7b))
str(datasub7b)
#now subset the quantitative variables out of HbASFdss
datasub7c=datasub7b[,-c(1:5)]
str(datasub7c)
plot(datasub7c)

#Pearson correlation and heat map HbASF no DSS

cormat7b<-round(cor(datasub7c, method="spearman", use = "complete.obs"),2)
heatmap(cormat7b)
head(cormat7b)

library(reshape2)
melted_cormat7b <- melt(cormat7b)
head(melted_cormat7b)

library(ggplot2)
ggplot(data = melted_cormat7b, aes(x=X1, y=X2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri7b<- get_upper_tri(cormat7b)
upper_tri7b


# Melt the correlation matrix
library(reshape2)
melted_cormat7c <- melt(upper_tri7b, na.rm = TRUE)
melted_cormat7c 
melted_cormat7c<-melted_cormat7c[order(melted_cormat7c$X1), ]
melted_cormat7c

#subseting from correlation matrix
head(melted_cormat7c, 10)
melted_cormat7c=melted_cormat7c[1:8,]
melted_cormat7c

# Heatmap
library(ggplot2)
ggheatmap7a<-ggplot(data = melted_cormat7c, aes(X2, X1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Spearman Correlation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 20, hjust = 0.7))+
  theme(axis.text.y = element_text(size=20))+
 coord_fixed()

#print(ggheatmap)

ggheatmap7b<-ggheatmap7a + 
geom_text(aes(X2, X1, label = value), color = "black", size = 12) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(1,1),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 12, barheight = 0.5,
                title.position = "top", title.hjust = 0.5,
                title.theme= element_text(size = 20, angle= 0)))

c1<-ggheatmap7b + ggtitle(expression(paste(italic('H. bilis'),"(pathobiont)-ASF-1.5% DSS")))+ 
      labs(x="C3H/HeN - Cecal tissue associated bacteria")+ 
      theme(plot.title = element_text(size = 22))+
      theme(axis.title.x = element_text(size = 22))
c1

library(gridExtra)
grid.arrange(a1,c1, nrow=2, clip=TRUE)
