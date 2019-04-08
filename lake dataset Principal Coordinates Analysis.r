fish.dat=read.csv("~/Documents/GitHub/ACE/Fish summary.csv")

names(fish.dat)
dim(fish.dat)

library(vegan)

#Let's select only sectors with more than 10 fish
fish.rowsums=apply(fish.dat[,3:20],1,sum)
fish.dat=fish.dat[fish.rowsums>10,]
dim(fish.dat)
n=nrow(fish.dat)

# At least in RGui, You might want to run the following command so that plots appear side-by-side
par(mfrow=c(1,2))


### Distance metrics for Beta diversity analysis and ordination

# examples for illustration just on samples 1 through 5
# Euclidean distance
#vegdist calculated distances- feed it a cmmunity data matrix. default distance is Bray-curtis but you can specific others via the 'method'-you can base on P/A by setting binary to True, defualt is false. 
fish.dist.euc=vegdist(fish.dat[1:5,3:20],method="euclidean")
fish.dist.euc
class(fish.dist.euc)

# Hellinger distance using the square-root relative abundance transform implemented in decostand
#these distances are the A's in our formula
#decostand- transforms data into hellinger or chi-square- pre-transformation to use before you process like a euclidian???? You get hellinger distances by calculating euclidiean distance on hellinger transformed distances

fish.hel=decostand(fish.dat[,3:20],method="hellinger")
fish.dist.hel=vegdist(fish.hel[1:5,],method="euclidean")
fish.dist.hel

# Hellinger distance with the transform done with math functions
#this is wat decostand does...
fish.rowsums=apply(fish.dat[,3:20],1,sum)
fish.hel2=sqrt(fish.dat[,3:20]/fish.rowsums)
fish.dist.hel2=vegdist(fish.hel2[1:5,],method="euclidean")
fish.dist.hel2

# Bray-Curtis distance
fish.dist.bc=vegdist(fish.dat[1:5,3:20], method="bray")
fish.dist.bc

# Jaccard distance
# Note that binary is specified as TRUE or you get a non-binary version of Jaccard
fish.dist.jac=vegdist(fish.dat[1:5,3:20],method="jaccard", binary=TRUE)
fish.dist.jac

# Chao's Jaccard-type distance corrected for missing species and using abundance data
fish.dist.chao=vegdist(fish.dat[1:5,3:20],method="chao")
fish.dist.chao

# Let's compare them all by correlation (this is for illustration only - no testing can be done on this)
cor(cbind(fish.dist.bc,fish.dist.euc,fish.dist.jac,fish.dist.chao,fish.dist.hel,fish.dist.hel2))


# For ordinations below, the following commands set up some graphics for us.
# First we want to bring in a dataset of environmental variables for each lake.
# Then we match this with the sector-based fish dataset.
# Then set up the colors variable for plotting.
# Let's color the lake sector names in the ordination plots by lake type
lake.dat<- read.csv("~/Documents/GitHub/ACE/MultiLake Data 1987 simplified.csv")
names(lake.dat)
#match the environmental variables to the respective sites diversity data
match(fish.dat[,1], lake.dat[,1])
fish.lakes=cbind(fish.dat,lake.dat[match(fish.dat[,1], lake.dat[,1]),])
names(fish.lakes)
colors=ifelse(fish.lakes$Type=="LC",2,ifelse(fish.lakes$Type=="HC",3,ifelse(fish.lakes$Type=="LP",4,5)))


### Principal Coordinates Analysis (PCoA)
# also known as Metric Multidimensional Scaling (MDS)

# Euclidean-based PCoA
fish.dist.euc=vegdist(fish.dat[,3:20],method="euclidean")
fish.euc.pcoa=wcmdscale(fish.dist.euc,eig=TRUE)
#Let's look at some properties of the object created, which has various parts
class(fish.euc.pcoa)
str(fish.euc.pcoa)
#parts of this object can be extracted with $ syntax or with special functions
fish.euc.pcoa$eig
dim(fish.euc.pcoa$points)
# Calculate proportion variance explained by each axis
eigenvals(fish.euc.pcoa)/sum(eigenvals(fish.euc.pcoa))
# Plot the first two axes and you have an ordination
# The default vegan plotting is fairly useles, imho
plot(fish.euc.pcoa)
# We can also extract scores to make a more useful plot here, or in another program 
fish.euc.pcoa.scores=scores(fish.euc.pcoa)
plot(fish.euc.pcoa,type="none",xlab="Axis 1 (54%)",ylab="Axis 2 (34%)")
text(fish.euc.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
# But we know Euclidean distance is not appropriate for species data.
# Let's use some better distances for our purpose

# Hellinger-based PCoA
fish.dist.hel=vegdist(fish.hel,method="euclidean")
fish.hel.pcoa=wcmdscale(fish.dist.hel,eig=TRUE)
eigenvals(fish.hel.pcoa)/sum(eigenvals(fish.hel.pcoa))
fish.hel.pcoa.scores=scores(fish.hel.pcoa)
plot(fish.hel.pcoa,type="none",xlab="Axis 1 (35%)",ylab="Axis 2 (23%)")
text(fish.hel.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)

# Chi-square-based PCoA
#most heavily weighs the rare species
fish.chi=decostand(fish.dat[,3:20],method="chi.square")
fish.dist.chi=vegdist(fish.chi,method="euclidean")
fish.chi.pcoa=wcmdscale(fish.dist.chi,eig=TRUE)
eigenvals(fish.chi.pcoa)/sum(eigenvals(fish.chi.pcoa))
fish.chi.pcoa.scores=scores(fish.chi.pcoa)
plot(fish.chi.pcoa,type="none",xlab="Axis 1 (22%)",ylab="Axis 2 (14%)")
text(fish.chi.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)

# Chi-square-based PCoA, weighted by rowsums
#not working, skipped in class
fish.wchi.pcoa=wcmdscale(fish.dist.chi,eig=TRUE,w=fish.rowsums)
eigenvals(fish.wchi.pcoa)/sum(eigenvals(fish.wchi.pcoa))
fish.wchi.pcoa.scores=scores(fish.wchi.pcoa)
plot(fish.wchi.pcoa,type="none",xlab="Axis 1 (22%)",ylab="Axis 2 (14%)")
#text command adds points as text instead of just points. The line above sets up the rest of the graph, axis, etc.
text(fish.wchi.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)

# Bray-Curtis-based PCoA
#not based on euclidean distance format so no need for decostand fuction first
#gives negative eigenvalues because its not euclidean distance. Focus on the first few (positive with highest values)
fish.dist.bc=vegdist(fish.dat[,3:20],method="bray")
fish.bc.pcoa=wcmdscale(fish.dist.bc,eig=TRUE)
eigenvals(fish.bc.pcoa)/sum(eigenvals(fish.bc.pcoa))
fish.bc.pcoa.scores=scores(fish.bc.pcoa)
plot(fish.bc.pcoa,type="none",xlab="Axis 1 (32%)",ylab="Axis 2 (25%)")
text(fish.bc.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)

# Bray-Curtis-based PCoA, negative eigenvalues fixed
eigenvals(fish.bc.pcoa)
fish.dist.bc.noneg=as.matrix(vegdist(fish.dat[,3:20],method="bray",diag=TRUE,upper=TRUE))
n=nrow(fish.dist.bc.noneg)
for(i in 1:n){
	for(j in 1:n){
		fish.dist.bc.noneg[i,j]=ifelse(i==j,0,sqrt(fish.dist.bc.noneg[i,j]^2+2*0.488514))
}}
fish.bc.noneg.pcoa=wcmdscale(fish.dist.bc.noneg,eig=TRUE)
eigenvals(fish.bc.noneg.pcoa)/sum(eigenvals(fish.bc.noneg.pcoa))
fish.bc.noneg.pcoa.scores=scores(fish.bc.noneg.pcoa)
plot(fish.bc.noneg.pcoa,type="none",xlab="Axis 1 (32%)",ylab="Axis 2 (25%)")
text(fish.bc.noneg.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
#corrects for negative values so you can calculate percent variation explained. 



# Jaccard-based PCoA
fish.dist.jac=vegdist(fish.dat[,3:20],method="jaccard",binary=TRUE)
fish.jac.pcoa=wcmdscale(fish.dist.jac,eig=TRUE)
eigenvals(fish.jac.pcoa)/sum(eigenvals(fish.jac.pcoa))
fish.jac.pcoa.scores=scores(fish.jac.pcoa)
plot(fish.jac.pcoa,type="none",xlab="Axis 1 (30%)",ylab="Axis 2 (13%)")
text(fish.jac.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)

# Jaccard-based PCoA, negative eigenvalues fixed
eigenvals(fish.jac.pcoa)
fish.dist.jac.noneg=as.matrix(vegdist(fish.dat[,3:20],method="jaccard",binary=TRUE,diag=TRUE,upper=TRUE))
n=nrow(fish.dist.jac.noneg)
for(i in 1:n){
	for(j in 1:n){
		fish.dist.jac.noneg[i,j]=ifelse(i==j,0,sqrt(fish.dist.jac.noneg[i,j]^2+2*abs(min(eigenvals(fish.jac.pcoa)))))
}}
fish.jac.noneg.pcoa=wcmdscale(fish.dist.jac.noneg,eig=TRUE)
eigenvals(fish.jac.noneg.pcoa)/sum(eigenvals(fish.jac.noneg.pcoa))
fish.jac.noneg.pcoa.scores=scores(fish.jac.noneg.pcoa)
plot(fish.jac.noneg.pcoa,type="none",xlab="Axis 1 (30%)",ylab="Axis 2 (13%)")
text(fish.jac.noneg.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)


# Adding species to PCoA plots
# Remember, species scores on PCoA plots are post-hoc weighted averages based on sample scores 
# Using the corrected Bray-Curtis and Jaccard plots as examples

# Bray-Curtis PCoA, adding species
par(mfrow=c(1,2))
plot(fish.bc.noneg.pcoa,type="none",xlab="Axis 1 (32%)",ylab="Axis 2 (25%)")
text(fish.bc.noneg.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
fish.bc.noneg.pcoa.spwa=wascores(fish.bc.noneg.pcoa.scores, fish.dat[,3:20])
text(fish.bc.noneg.pcoa.spwa[,1:2],labels=rownames(fish.bc.noneg.pcoa.spwa),col=8)

# Jaccard PCoA, adding species
plot(fish.jac.noneg.pcoa,type="none",xlab="Axis 1 (30%)",ylab="Axis 2 (13%)")
text(fish.jac.noneg.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
fish.jac.noneg.pcoa.spwa=wascores(fish.jac.noneg.pcoa.scores, fish.dat[,3:20])
text(fish.jac.noneg.pcoa.spwa[,1:2],labels=rownames(fish.jac.noneg.pcoa.spwa),col=8)

# You could add anything that is measured on your samples this way. 
# Let's see what it does with Hellinger plots and species and environmental variables.
# Hellinger PCoA
plot(fish.hel.pcoa,type="none",xlab="Axis 1 (35%)",ylab="Axis 2 (23%)")
text(fish.hel.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
fish.hel.pcoa.spwa=wascores(fish.hel.pcoa.scores, fish.dat[,3:20])
text(fish.hel.pcoa.spwa[,1:2],labels=rownames(fish.hel.pcoa.spwa),col=8)
# Let's also try adding some environmental variables for kicks- post hoc weighted average scores
#tells you what weighted is or not or who?

fish.hel.pcoa.envwa=3*wascores(fish.hel.pcoa.scores, fish.lakes[,24:28])
text(fish.hel.pcoa.envwa[,1:2],labels=rownames(fish.hel.pcoa.envwa),col=1)
# Here is an example of the convex hull plots (they are straight lines connecting the outer most samples of one type)
plot(fish.hel.pcoa,type="none",xlab="Axis 1 (35%)",ylab="Axis 2 (23%)")
text(fish.hel.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
ordihull(fish.hel.pcoa,fish.lakes[,23])
# Here is an example of the environmental contour plots _fitted_ and then overlaid on the ordination.
plot(fish.hel.pcoa,type="none",xlab="Axis 1 (35%)",ylab="Axis 2 (23%)")
text(fish.hel.pcoa.scores[,1:2],labels=fish.dat$LAKE,col=colors)
ordisurf(fish.hel.pcoa.scores, fish.lakes[,24])

#watch out for hulls, they are convex and always angle inward
