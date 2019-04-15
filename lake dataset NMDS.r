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


#### Non-Metric MultiDimensional Scaling

# nMDS of Hellinger distance matrix
fish.hel=decostand(fish.dat[,3:20],method="hellinger")
fish.dist.hel=vegdist(fish.hel,method="euclidean")
fish.hel.nmds=metaMDS(fish.dist.hel, k=2, trymax=100, wascores=TRUE)
fish.hel.nmds$stress
stressplot(fish.hel.nmds)
goodness(fish.hel.nmds)
plot(fish.hel.nmds, type="none")
text(fish.hel.nmds, labels=fish.dat$LAKE,col=colors)
fish.hel.nmds.spwa=wascores(scores(fish.hel.nmds), fish.dat[,3:20])
text(fish.hel.nmds.spwa[,1:2],labels=rownames(fish.hel.nmds.spwa),col=8)
#ONE THING WE DON'T LIKE ABOUT NMDS IS THAT WE DON'T GET PERCENT VARIANCE EXPLAINED- THERE ARE NO EIGENVALUES TO CALCUATE FROM
# Compare to Hellinger-based PCA
fish.hel.pca=rda(fish.hel)
plot(fish.hel.pca,type="none",xlab="PC1 (35%)",ylab="PC2 (23%)")
text(fish.hel.pca,labels=fish.dat$LAKE,col=colors)
text(fish.hel.pca,display="species",col=8)
# Is nMDS doing its job as intended?
fish.hel.pca.dist=vegdist(scores(fish.hel.pca,display="sites",choices=1:2),method="euclidean")
fish.hel.nmds.dist=vegdist(scores(fish.hel.nmds),method="euclidean")
cor(cbind(fish.dist.hel,fish.hel.nmds.dist,fish.hel.pca.dist))
#HERE WE LOOK AT THE CORRELATION MATRIX FOR THE RAW, HELLINGER NMDS AND HELLINGER PCA. NMDS IS MORE SIMILAR TO RAW THAN PCA IS, WHICH IS WHAT WE EXPECT- FOR THE FIRST TWO AXES AT LEAST. NMDS DOES A BETTER JOB THAN THE PCA AT DEMONSTRATING THE DATA, BUT AGAIN, IT DOESN'T SHOW THE VARIANCE EXPLAINED.. VARIANCE EXPLAINED HELPS THE READER UNDERSATND THE ORDINATION'S VALUE BETTER. STRESS SEEMS LESS TELLING OF HOW WELL THE ORDINATION DOES ITS JOB

#YOU CAN FEED ANY DISTANCE MATRIX INTO NMDS, WHICH IS NICE
# I started doing some other nMDS types here.
# Bray-Curtis nMDS
fish.dist.bc=vegdist(fish.dat[,3:20], method="bray")
fish.bc.nmds=metaMDS(fish.dist.bc, k=2, trymax=100)
fish.bc.nmds$stress
plot(fish.bc.nmds, type="none")
text(fish.bc.nmds, labels=fish.dat$LAKE,col=colors)
text(matrix(c(0.3,0.45),nrow=1),labels="Stress=0.157")
#STRESS IS HIGHER, SO THIS ORDINATION ISN'T AS GOOD AS THE HELLINGER DISTANCE ORDINATION
#for fun, check this out
ordisurf(fish.bc.nmds, fish.lakes[,25])
#COOL, BUT WE CAN'T TELL HOW WELL THE TOPO LINES FIT... IT WOULD BE SO MUCH BETTER IF WE HAD AN R^2 OF SOME SORT...

# Compare to Bray-Curtis PCoA
plot(fish.bc.nmds, type="none")
text(fish.bc.nmds, labels=fish.dat$LAKE,col=colors)
text(matrix(c(0.3,0.45),nrow=1),labels="Stress=0.157")
fish.bc.pcoa=wcmdscale(fish.dist.bc,eig=TRUE)
eigenvals(fish.bc.pcoa)/sum(eigenvals(fish.bc.pcoa))
fish.bc.pcoa.scores=scores(fish.bc.pcoa)
plot(fish.bc.pcoa,type="none",xlab="Axis 1 (32%)",ylab="Axis 2 (25%)")
text(fish.bc.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)

# Is nMDS doing its job as intended?
fish.bc.nmds.dist=vegdist(scores(fish.bc.nmds),method="euclidean")
fish.bc.pcoa.dist=vegdist(scores(fish.bc.pcoa,display="sites",choices=1:2),method="euclidean")
cor(cbind(fish.dist.bc,fish.bc.nmds.dist,fish.bc.pcoa.dist))
#YES, NMDS IS DOING ITS JOB- GREATER CORRELATION WITH RAW

# Jaccard nMDS
fish.dist.jac=vegdist(fish.dat[,3:20],method="jaccard", binary=TRUE)
fish.jac.nmds.2=metaMDS(fish.dist.jac, k=2, trymax=100)
fish.jac.nmds.2$stress
plot(fish.jac.nmds.2, type="none")
text(fish.jac.nmds.2, labels=fish.dat$LAKE,col=colors)
text(matrix(c(0.3,0.4),nrow=1),labels="Stress=0.163")
#EVEN AFTER 100 STEPS, IT WASN'T ABLE TO CONVERGE. SAVES THE BEST. 
# Try varying the axes requested to see how it impacts the first two (the ones shown).
# This is a good demonstration of the problems with requesting more than two axes.
#SO, IF YOU INCREASE K TO 4, YOU INCREASE THE DIMENSIONS.. BUT DISPLAYING THE RESULTS IS IMPOSSIBLE..
fish.jac.nmds.4=metaMDS(fish.dist.jac, k=4, trymax=100)
fish.jac.nmds.4$stress
plot(fish.jac.nmds.4, type="none")
text(fish.jac.nmds.4, labels=fish.dat$LAKE,col=colors)
text(matrix(c(0,0.4),nrow=1),labels="I need 4 dimensions to show you this solution :-(")
# Compare to Jaccard PCoA
plot(fish.jac.nmds.2, type="none")
text(fish.jac.nmds.2, labels=fish.dat$LAKE,col=colors)
text(matrix(c(0.3,0.4),nrow=1),labels="Stress=0.163")
fish.jac.pcoa=wcmdscale(fish.dist.jac,eig=TRUE)
eigenvals(fish.jac.pcoa)/sum(eigenvals(fish.jac.pcoa))
fish.jac.pcoa.scores=scores(fish.jac.pcoa)
plot(fish.jac.pcoa,type="none",xlab="Axis 1 (30%)",ylab="Axis 2 (15%)")
text(fish.jac.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
# Is nMDS doing its job as intended?
fish.jac.pcoa.dist=vegdist(scores(fish.jac.pcoa,display="sites",choices=1:2),method="euclidean")
fish.jac.nmds.dist=vegdist(scores(fish.jac.nmds.2),method="euclidean")
cor(cbind(fish.dist.jac,fish.jac.nmds.dist,fish.jac.pcoa.dist))




