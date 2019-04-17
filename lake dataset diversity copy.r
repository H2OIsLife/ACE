# Read in the data and look at some basic features (make sure it read in correctly)
fish.dat<- read.csv("~/Desktop/ACE/Lab1/Fish summary.csv")
names(fish.dat)
dim(fish.dat)
head(fish.dat)
# Total caught per species
# the 2 means its applying this command by column, 1 is applying it by row 
fish.colsums=apply(fish.dat[,3:20],2,sum)
# Total number of fish per sample
fish.rowsums=apply(fish.dat[,3:20],1,sum)
#check out what we did
fish.colsums
fish.rowsums

# Its so easy to calculate richness, right?  Well, maybe...
fish.rich=apply(fish.dat[,3:20]>0,1,sum)
fish.rich
# add total catch and raw richness to the dataset for future use, and look at the lowest number of fish caught in a sample
#cbind stands for column bind, so we are combining the dataset info, richness and rowsums
cbind(fish.dat[,1:2],fish.rich,fish.rowsums)
min(fish.rowsums)
as.data.frame(table(fish.rowsums))
#start here on Thursday



library(vegan)

###Individual based rarefaction
# Using the rarefaction formula to find expected species richness
rarecurve(fish.dat[,3:20])
rarecurve(fish.dat[fish.rowsums>5,3:20])
# You can choose a minimum sample size to rarefy all samples to using the command below
#uses formulas from Gotelli et al
rarefy(fish.dat[,3:20],sample=6)

# Shannon, Simpson, and 1/Simpson can be calculated by this command
diversity(fish.dat[,3:20],index="shannon")  
# Hill and Renyi numbers are calculated by this command- exponentiates the shannon values
renyi(fish.dat[,3:20], hill=TRUE)
plot(renyi(fish.dat[,3:20], hill=TRUE))
plot(renyi(fish.dat[,3:20], scales=c(0,1,2,Inf), hill=TRUE))
#you would have to dump the NaN (dividing my 0 error) before you can do stats, you would also have to choose what hill number to run stats on and filter out the samples that are 'bad'
#

# Individual-based rarefaction for other indices.  
# Writing our own simulation code to get average rarefied index values.
# Start by looking at one rarefaction "realization"
###head commands are simulations- they will be different everytime- so you do it a bunch of times and then take an average. constrained by computing power.. so maybe 99 or 999 is enough- you can run at different nperm #'s and see how much your end values change
rrarefy(fish.dat[,3:20],sample=6)
rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6)
head(rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6))
head(rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6))
#
# Perform random sampling "nperm" times and then take averages
# First set some parameters and set up an empty matrix to store simulated communities
nperm=200
n=nrow(fish.dat[fish.rowsums>5,])
p=ncol(fish.dat[,3:20])
#matrix to store the results
renyi.indperms=matrix(nrow=n,ncol=4*nperm)
# Simulate rarefied communities "nperm" times, i is a variable that ends up being, 1,2,3... for each perm. 
#(1+(i-1)*4):(i*4) tells the loop where to store the data- we need four columns for every iteration
#rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6) simulates the community
#renyi(fish.rare6, scales=c(0,1,2,Inf), hill=TRUE) runs the calculation to get the hill numbers
#renyi.indperms[,cols]=as.matrix(renyi.rare)
for(i in 1:nperm){
	fish.rare6=rrarefy(fish.dat[fish.rowsums>5,3:20],sample=6)
	cols=(1+(i-1)*4):(i*4)
	renyi.rare=renyi(fish.rare6, scales=c(0,1,2,Inf), hill=TRUE)
	renyi.indperms[,cols]=as.matrix(renyi.rare)
}
maxperm=nperm*4

#look at what we did
head(renyi.indperms)
# seq(from=1,to=(maxperm-3),by=4) gives you 1,5,9.... to 797
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
#the product of this is what you would then run stats on... this would be averaged, rarified data. 
#how do you get the average for each hill number for each lake




# Sample-based rarefaction
#simulating a group of samples instead of a group of individuals
# These functions work on an entire dataset. Further coding could be used to get results on subsets (like lakes).
sp.accum=specaccum(fish.dat[,3:20], method="exact", conditioned=T)
plot(sp.accum)
#error bars are calculated with permutations- sampling simulations- lower variation^
sp.accum=specaccum(fish.dat[,3:20], method="exact", conditioned=F)
plot(sp.accum)
#more variation around the highest number, assumes this is a subset of the population^
sp.accum=specaccum(fish.dat[,3:20], method="rarefaction", conditioned=F)
plot(sp.accum, xvar="individuals")
# Sample-based rarefaction for Hill numbers.
# Plot shows each rarefaction realization.
renyiaccum(fish.dat[fish.rowsums>5,3:20], scales=c(0,1,2,Inf), permutations=100, raw=FALSE, hill=TRUE)
plot(renyiaccum(fish.dat[,3:20], scales=c(0,1,2,Inf), permutations=10, raw=TRUE, hill=TRUE))
#^ shows the 4 different hill number orders and each accumulation curve. The above graphs is a average of the accumulation curves. This would be more interesting if this was by lake... 
#inf plot is the Berger Parker hill numbers- its easier to estimate the evenness than it is to estimate richness, so Berger Parker reaches asymptope with few samples. High sample numbers are needed for richness analysis, but getting a good estimate for proportional abundance is not so dependent on rare taxa, which hardly make up any of the community anyway (contribute less to the Simposon(2) and Berger-Parker). Plots 1 (shannon) and 0(richness) are dependnet on high sampling efforts...
# Sample-based rarefaction by sample subsets (here, lakes)

#troubleshooting BiodiversityR package issues.....
install.packages("data.table", type = "source",
                 repos = "http://Rdatatable.github.io/data.table")
library('data.table')
install.packages('car')
library('car')
library(BiodiversityR)
library(tcltk)
library(vegan)

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


#### stopped here on Monday, 2/11
#Load the lake dataset. 
lake.dat <- read.csv("~/Desktop/ACE/Lab1/MultiLake Data 1987 simplified.csv")
#we pseudoreplicate on the response side (multiple sectors in the same lakes).. 
#we can pseudoreplicate on the explanatory side?
#we can take the average across the sectors for each lake and use that to avoid non-independence and get rid of pseudoreplication. 
lake.dat
names(lake.dat)
match(fish.hills[,1], lake.dat[,1])
fish.lakes=cbind(lake.dat[match(fish.hills[,1], lake.dat[,1]),],fish.hills[,2:6])
head(fish.hills)
head(fish.lakes)
names(fish.lakes)


#Classical approach (not mixed)
#dependent, independent. see if we can explain variation on the left (lake type) by variable on the right richness?
#ANOVA assumptions- independence, normality, equal variance- we are violating independence- lake sectors come from the same lakes...they are not independent- but averaging reduces the power of your test. so to get around this we use a mixed model with random effects. You basically tell statistics how you are violating independence and it accounts for that. corrects for the violation. 
H0.nomix=lm(fish.lakes[,11]~fish.lakes[,3])
anova(H0.nomix)
AIC(H0.nomix)

#Mixed model approach
#Must download the nlme package before using the nlme library
library(nlme)
H0.lake=lme(fixed=H0~Cobble, random=~1|factor(LAKE), data=fish.lakes)
H0.laketype=lme(fixed=H0~Type, random=~1|factor(LAKE), data=fish.lakes)
anova(H0.laketype)
AIC(H0.laketype)
VarCorr(H0.laketype)
#lower AIC means its a better model- choose the model with the lowest AIC
print("What's the pirate's coding language?")
print("R")
