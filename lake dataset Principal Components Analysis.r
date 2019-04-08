fish.dat=read.csv("~/Documents/GitHub/ACE/Fish summary.csv")

names(fish.dat)
dim(fish.dat)

library(vegan)

#Let's select only sectors with more than 10 fish
fish.rowsums=apply(fish.dat[,3:20],1,sum)
fish.dat=fish.dat[fish.rowsums>10,]
dim(fish.dat)
n=nrow(fish.dat)

####
# A block to set up for some graphics
# At least in RGui, You might want to run the following command so that plots appear side-by-side
par(mfrow=c(1,2))
# For ordinations below, the following commands set up some graphics for us.
# First we want to bring in a dataset of environmental variables for each lake.
# Then we match this with the sector-based fish dataset.
# Then set up the colors variable for plotting.
# Let's color the lake sector names in the ordination plots by lake type
lake.dat=read.csv("~/Documents/GitHub/ACE/MultiLake Data 1987 simplified.csv")

names(lake.dat)
match(fish.dat[,1], lake.dat[,1])
fish.lakes=cbind(fish.dat,lake.dat[match(fish.dat[,1], lake.dat[,1]),])
names(fish.lakes)
colors=ifelse(fish.lakes$Type=="LC",2,ifelse(fish.lakes$Type=="HC",3,ifelse(fish.lakes$Type=="LP",4,5)))
####


### PCA

#QUESTIONS FOR MONDAY
# Why do we take the means?
# What is centering, why do we do it? Why does it give us negative numbers, and why are the means 0?
#for the eigenvalues, do we caulcate a percentage of variance explained by an axes with those? When we go to plot arounf line 70, why do we multiple our eigen values by the raw data, why not just use the fish scores? 

# This block starts us off with Euclidean distance as a demo - simple, but incorrect for species data
# First I'm going to do it "by hand"
n=nrow(fish.dat)
p=ncol(fish.dat[,3:20])
fish.means=apply(fish.dat[,3:20],2,sum)/n
fish.means
means.matrix=t(matrix(nrow=p,ncol=n,fish.means))
fish.centered=fish.dat[,3:20]-means.matrix
# take a look at what this centering did
fish.centered[5,]
fish.centered.means=apply(fish.centered,2,sum)/n
fish.centered.means
# here is the covariance matrix then
#%*% is a matrix multiplication
#a correlation matrix would have diagonals all of 1- thats where the species compare to themselves, so they are perfectly correlated. 

(1/(n-1))*as.matrix(t(fish.centered))%*%as.matrix(fish.centered)

# Here is a function for the covariance matrix, shortcut
fish.cov=cov(fish.dat[,3:20])
fish.cov

# Eigen analysis with R base functions
fish.eigen=eigen(fish.cov)
fish.eigen
#^yields species scores coefficent for each pairwaise comparison and eigenvalues, one eigenvalue for each eigenvector. Higher eigenvalues tell you that more variation is explained in that axes. vectors are species scores, which are then used to calculate the site location? Values are the axes. 
fish.colsums=apply(fish.dat[,3:20],2,sum)
fish.propexpl=fish.eigen$values/sum(fish.eigen$values)
fish.propexpl
#^tells you the percent variation explained by each axes
fish.scores=as.matrix(fish.dat[,3:20])%*%as.matrix(fish.eigen$vectors)
head(fish.scores)
head(fish.dat)
head(fish.eigen$vectors)
plot(fish.scores[,1:2],col=colors)

## Ok here is the vegan way to do PCA
par(mfrow=c(1,1))
fish.pca=rda(fish.dat[,3:20])
class(fish.pca)
str(fish.pca)
summary(fish.pca)
eigenvals(fish.pca)
eigenvals(fish.pca)/sum(eigenvals(fish.pca))
plot(fish.pca)
# With PCA we go from this type of display
pairs(fish.dat[,3:20])
# ^what is this? 
# to this type of display
plot(fish.scores[,1:2],col=colors)
plot(fish.pca,type="none",xlab="PC1 (54%)",ylab="PC2 (34%)")
text(fish.pca,labels=fish.dat$LAKE,col=colors)
text(fish.pca,display="species",col=8)

# just to compare some scores and show the methods are giving the same relative positions among sites
# The "scores" command will extract the site and species scores
fish.pca.scores=scores(fish.pca,choices=1:18,display="sites")
fish.pca.spscores=scores(fish.pca,choices=1:18,display="species")
# To show that the function and the "by hand" method give the same results, although they are scaled differently the relative positions of points are identical
cor(cbind(fish.scores[,1:3],fish.pca.scores[,1:3]))
cor(cbind(fish.eigen$vectors[,1:3],fish.pca.spscores[,1:3]))
#^weren't we using covariance, not correlation?? What exactly are site scores? 

# What about comparing to PCoA?
fish.dist.euc=vegdist(fish.dat[,3:20],method="euclidean")
fish.euc.pcoa=wcmdscale(fish.dist.euc,eig=TRUE)
eigenvals(fish.euc.pcoa)/sum(eigenvals(fish.euc.pcoa))
fish.euc.pcoa.scores=scores(fish.euc.pcoa)
plot(fish.euc.pcoa,type="none",xlab="PCoA Axis 1 (54%)",ylab="PCoA Axis 2 (34%)")
text(fish.euc.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
cor(cbind(fish.euc.pcoa.scores[,1:3],fish.pca.scores[,1:3]))
#vegdist doesn't give species scores??

### Hellinger distance-based PCA
# Everything above was based on Euclidean distance - which we know is bad for species data
# Here is one better option that can be justified.
# Hellinger distance-based PCA
fish.hel=decostand(fish.dat[,3:20],method="hellinger")
fish.hel.pca=rda(fish.hel)
summary(fish.hel.pca)
fish.hel.pca.scores=scores(fish.hel.pca,choices=1:18,display="sites")
plot(fish.hel.pca,type="none",xlab="PC1 (35%)",ylab="PC2 (25%)")
text(fish.hel.pca,labels=fish.dat$LAKE,col=colors)
text(fish.hel.pca,display="species",col=8)
# playing around with adding environmental variable weighted averages to the plot
#^ why do you weight the averages?? 
#dominant species has been downweighted compared to euclidean, so less variation is explained by it.. the variation is parsed out more realistically
fish.hel.pca.envfit=envfit(fish.hel.pca.scores, fish.lakes[,24:28])
fish.hel.pca.envwa=wascores(fish.hel.pca.scores, fish.lakes[,24:28])
text(fish.hel.pca.envwa[,1:2],labels=rownames(fish.hel.pca.envwa),col=1)
ordihull(fish.hel.pca,fish.lakes[,23])
# The species are scaled quite differently, although the environmental variables appear to be similar, which is interesting and someone should explore it.  But its time to move on.

###Hellinger based PCoA and PCA yield the same results#
#how? they are calculated totally differently? 
# Is it the same ordination as found by Hellinger-based PCoA?
fish.dist.hel=vegdist(fish.hel,method="euclidean")
fish.hel.pcoa=wcmdscale(fish.dist.hel,eig=TRUE)
eigenvals(fish.hel.pcoa)/sum(eigenvals(fish.hel.pcoa))
fish.hel.pcoa.scores=scores(fish.hel.pcoa)
cor(cbind(fish.hel.pcoa.scores[,1:4],fish.hel.pca.scores[,1:4]))
#^what is this doing?? binding the same things?

########
# PCA based on a correlation matrix instead of covariance
# Lets take another look at the abiotic variables lake dataset. 
lake.dat=read.csv("~/Documents/GitHub/ACE/MultiLake Data 1987 simplified.csv")
names(lake.dat)
lake.pca=rda(lake.dat[,4:8])
summary(lake.pca)
# Does it make sense to perform this PCA using Euclidean distance directly on these variables?
# Nope, definitely not.
# We can use the "standardize" tranformation to instead base the PCA on the correlation matrix, which removes differences in scales among variables.
lake.pca.cor=rda(decostand(lake.dat[,4:8],method="standardize"))
summary(lake.pca.cor)
#how does correlation remove differences in scale between variables


### Correspondence Analysis

# Let's mess with the once popular Correspondence Analysis and Chi-square distance
# Chi square distance-based PCA
fish.chi=decostand(fish.dat[,3:20],method="chi.square")
fish.chi.pca=rda(fish.chi)
summary(fish.chi.pca)
plot(fish.chi.pca,type="none",xlab="PC1 (23%)",ylab="PC2 (17%)")
text(fish.chi.pca,labels=fish.dat$LAKE,col=colors)
text(fish.chi.pca,display="species",col=8)
fish.chi.pca.scores=scores(fish.chi.pca,choices=1:17,display="sites")

# Chi-square-based PCoA
fish.dist.chi=vegdist(fish.chi,method="euclidean")
fish.chi.pcoa=wcmdscale(fish.dist.chi,eig=TRUE)
eigenvals(fish.chi.pcoa)/sum(eigenvals(fish.chi.pcoa))
fish.chi.pcoa.scores=scores(fish.chi.pcoa)
plot(fish.chi.pcoa,type="none",xlab="Axis 1 (23%)",ylab="Axis 2 (17%)")
text(fish.chi.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)

# Classical method of Correspondence Analysis in Vegan
fish.ca=cca(fish.dat[,3:20])
summary(fish.ca)
plot(fish.ca,type="none",xlab="CA1 (20%)",ylab="CA2 (18%)")
text(fish.ca,labels=fish.dat$LAKE,col=colors)
text(fish.ca,display="species",col=8)
fish.ca.scores=scores(fish.ca,choices=1:17,display="sites")

# Chi-square-based PCoA, weighted by rowsums
fish.rowsums=apply(fish.dat[,3:20],1,sum)
fish.wchi.pcoa=wcmdscale(fish.dist.chi,eig=TRUE,w=as.numeric(fish.rowsums))
eigenvals(fish.wchi.pcoa)/sum(eigenvals(fish.wchi.pcoa))
fish.wchi.pcoa.scores=scores(fish.wchi.pcoa)
plot(fish.wchi.pcoa,type="none",xlab="Axis 1 (22%)",ylab="Axis 2 (14%)")
text(fish.wchi.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
#why does weighting change the ordination? How are you weighting? 

# Is it the same ordination as found by Chi-square-based PCA or the same as CA? 
cor(cbind(fish.chi.pcoa.scores[,1:4],fish.chi.pca.scores[,1:4],fish.ca.scores[,1:4]))
# Is it the same ordination as found by CA?
cor(cbind(fish.wchi.pcoa.scores[,1:4],fish.ca.scores[,1:4]))
#how do you interpret this correlation matrix? 





