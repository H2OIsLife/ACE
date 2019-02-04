# Read in the data and look at some basic features (make sure it read in correctly)
fish.dat=read.csv("C:\\Users\\cblackwo\\Documents\\Teaching\\2019_Spring_ACE\\lake survey dataset\\Fish summary.csv")
names(fish.dat)
dim(fish.dat)
head(fish.dat)
# Total caught per species
fish.colsums=apply(fish.dat[,3:20],2,sum)
# Total number of fish per sample
fish.rowsums=apply(fish.dat[,3:20],1,sum)

# Its so easy to calculate richness, right?  Well, maybe...
fish.rich=apply(fish.dat[,3:20]>0,1,sum)
# add total catch and raw richness to the dataset for future use, and look at the lowest number of fish caught in a sample
cbind(fish.dat[,1:2],fish.rich,fish.rowsums)
min(fish.rowsums)
as.data.frame(table(fish.rowsums))

library(vegan)

###Individual based rarefaction
# Using the rarefaction formula to find expected species richness
rarecurve(fish.dat[,3:20])
rarecurve(fish.dat[fish.rowsums>5,3:20])
# You can choose a minimum sample size to rarefy all samples to using the command below
rarefy(fish.dat[,3:20],sample=6)

# Shannon, Simpson, and 1/Simpson can be calculated by this command
diversity(fish.dat[,3:20],index="shannon")  
# Hill and Renyi numbers are calculated by this command
renyi(fish.dat[,3:20], hill=TRUE)
plot(renyi(fish.dat[,3:20], hill=TRUE))
plot(renyi(fish.dat[,3:20], scales=c(0,1,2,Inf), hill=TRUE))


# Individual-based rarefaction for other indices.  
# Writing our own simulation code to get average rarefied index values.
# Start by looking at one rarefaction "realization"
rrarefy(fish.dat[,3:20],sample=6)
rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6)
head(rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6))
head(rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6))

# Perform random sampling "nperm" times and then take averages
# First set some parameters and set up an empty matrix to store simulated communities
nperm=20
n=nrow(fish.dat[fish.rowsums>5,])
p=ncol(fish.dat[,3:20])
renyi.indperms=matrix(nrow=n,ncol=4*nperm)
# Simulate rarefied communities "nperm" times
for(i in 1:nperm){
	fish.rare6=rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6)
	cols=(1+(i-1)*4):(i*4)
	renyi.rare=renyi(fish.rare6, scales=c(0,1,2,Inf), hill=TRUE)
	renyi.indperms[,cols]=as.matrix(renyi.rare)
}
maxperm=nperm*4
# Below code is just for demo
apply(renyi.indperms[,seq(from=1,to=(maxperm-3),by=4)],1,mean)
apply(renyi.indperms[,seq(from=2,to=(maxperm-2),by=4)],1,mean)
apply(renyi.indperms[,seq(from=3,to=(maxperm-1),by=4)],1,mean)
apply(renyi.indperms[,seq(from=4,to=maxperm,by=4)],1,mean)
# Store the means across simulations
fish.hills=cbind(fish.dat[fish.rowsums>5,c(1,2)],
	apply(renyi.indperms[,seq(from=1,to=(maxperm-3),by=4)],1,mean),
	apply(renyi.indperms[,seq(from=2,to=(maxperm-2),by=4)],1,mean),
	apply(renyi.indperms[,seq(from=3,to=(maxperm-1),by=4)],1,mean),
	apply(renyi.indperms[,seq(from=4,to=maxperm,by=4)],1,mean))
colnames(fish.hills)=c("LAKE","SECTOR","H0","H1","H2","HInf")
fish.hills


# Sample-based rarefaction
# These functions work on an entire dataset. Further coding could be used to get results on subsets (like lakes).
sp.accum=specaccum(fish.dat[,3:20], method="exact", conditioned=T)
plot(sp.accum)
sp.accum=specaccum(fish.dat[,3:20], method="exact", conditioned=F)
plot(sp.accum)
sp.accum=specaccum(fish.dat[,3:20], method="rarefaction", conditioned=F)
plot(sp.accum, xvar="individuals")
# Sample-based rarefaction for Hill numbers.
# Plot shows each rarefaction realization.
renyiaccum(fish.dat[fish.rowsums>5,3:20], scales=c(0,1,2,Inf), permutations=100, raw=FALSE, hill=TRUE)
plot(renyiaccum(fish.dat[,3:20], scales=c(0,1,2,Inf), permutations=10, raw=TRUE, hill=TRUE))


# Sample-based rarefaction by sample subsets (here, lakes)
library(BiodiversityR)

# Method for rarefaction of richness
Accum.fish.uncond=accumcomp(fish.dat[,3:20], y=fish.dat[,1:2], factor='LAKE', method='exact', 
	conditioned =F, gamma = 'boot', permutations=100, legend=F, rainbow=T, ci=2, 
	ci.type='bar', cex=1, xlab='sites', scale='')
Accum.fish.uncond
Accum.fish.cond=accumcomp(fish.dat[,3:20], y=fish.dat[,1:2], factor='LAKE', method='exact', 
	conditioned =T, gamma = 'boot', permutations=100, legend=F, rainbow=T, ci=2, ci.type='bar', 
	cex=1, xlab='sites', scale='')
Accum.fish.cond
Accum.fish.rare=accumcomp(fish.dat[,3:20], y=cbind(fish.dat[,1:2],fish.rowsums), factor='LAKE', 
	method='rarefaction', conditioned =F, gamma = 'boot', permutations=100, legend=F, rainbow=T, 
	xvar="individuals", xlab='individuals', scale='fish.rowsums')
Accum.fish.rare

# Diversity index calculation on subsets
# NOTE: This does not perform rarefaction...so Feed this a rarefied dataset.
# Note difference below between pooled and mean methods
diversitycomp(fish.rare6, y=fish.dat[fish.rowsums>5,1:2], factor1='LAKE',
    index="inverseSimpson", method="pooled")
diversitycomp(fish.rare6, y=fish.dat[fish.rowsums>5,1:2], factor1='LAKE',
    index="inverseSimpson", method="mean")
# How would you transform this to Hill number equivalent?
diversitycomp(fish.rare6, y=fish.dat[fish.rowsums>5,1:2], factor1='LAKE',
    index="Shannon", method="mean")


# Using iNEXT
library(iNEXT)
# demo data
data(spider)
iNEXT(spider$Logged, q=0, datatype="abundance")
some.spider.data=iNEXT(spider, q=c(0,1,2), datatype="abundance")
ggiNEXT(some.spider.data, type=1, facet.var="site")

iNEXT(as.numeric(fish.dat[16,3:20]), q=c(0,1,2), datatype="abundance")
some.fish.data=iNEXT(as.numeric(fish.dat[16,3:20]), q=c(0,1,2), datatype="abundance")
ggiNEXT(some.fish.data, type=1, facet.var="site")
ggiNEXT(some.fish.data, type=3, facet.var="site")
some.fish.data=iNEXT(as.numeric(fish.dat[1,3:20]), q=c(0,1,2), datatype="abundance")
ggiNEXT(some.fish.data, type=1, facet.var="site")

fish.inext=apply(t(fish.dat[fish.rowsums>5,3:20]),2,as.abucount)
all.fish.data=iNEXT(fish.inext, q=0, datatype="abundance")
ggiNEXT(all.fish.data, type=1, facet.var="site")
estimateD(fish.inext,datatype="abundance",base="coverage")
estimateD(fish.inext,datatype="abundance",base="coverage",level=0.9)
estimateD(fish.inext,datatype="abundance",base="coverage",level=0.5)

estimateD(fish.inext,datatype="abundance",base="size")
estimateD(fish.inext,datatype="abundance",base="size",level=15)



#Load the lake dataset. 
lake.dat=read.csv("C:\\Users\\cblackwo\\Documents\\Teaching\\2015_Spring_AdvancedCommunityEcology\\lake survey dataset\\MultiLake Data 1987 simplified.csv")
lake.dat
names(lake.dat)
match(fish.hills[,1], lake.dat[,1])
fish.lakes=cbind(lake.dat[match(fish.hills[,1], lake.dat[,1]),],fish.hills[,2:6])
names(fish.lakes)

#Classical approach (not mixed)
H0.nomix=lm(fish.lakes[,11]~fish.lakes[,3])
anova(H0.nomix)
AIC(H0.nomix)

#Mixed model approach
#Must download the nlme package before using the nlme library
library(nlme)
H0.lake=lme(fixed=H0~Cobble, random=~1|factor(LAKE), data=fish.lakes)
anova(H0.lake)
AIC(H0.lake)
VarCorr(H0.lake)

