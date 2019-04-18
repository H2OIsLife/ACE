library(vegan)

# load the data
# you will need to change the path for your computer, or import the data the way you like to
inv.dat=read.csv("C:\\Users\\cblackwo\\Dropbox\\1-files\\Teaching\\2015_Spring_AdvancedCommunityEcology\\manistee invert\\manistee invert factorial.csv")
dim(inv.dat)
names(inv.dat)

inv.colsums=apply(inv.dat[,11:2046],2,sum)
inv.rowsums=apply(inv.dat[,11:2046],1,sum)
min(inv.rowsums)
min=min(inv.rowsums)
n=nrow(inv.dat)

# Rarefaction and Hellinger distance transformation
inv.rare=cbind(inv.dat[,1:10],rrarefy(inv.dat[,11:2046],sample=min))
inv.rare.hel=cbind(inv.rare[,1:10],decostand(inv.rare[,11:2046],method="hellinger"))
inv.stands.hel.dist=vegdist(inv.rare.hel[,11:2046],method="euclidean")


# RDA - Redundancy Analysis

###
# All factors
# First, using Y,X,Z (with ,) format
rda.inv.rare.hel.all=rda(inv.rare.hel[,11:2046],inv.rare.hel[,6:10])
	summary(rda.inv.rare.hel.all)
	anova(rda.inv.rare.hel.all)
	RsquareAdj(rda.inv.rare.hel.all)
# Second, run this with formula interface (with ~) format
rda.inv.rare.hel.all=rda(inv.rare.hel[,11:2046]~inv.rare.hel$Ecosystem*inv.rare.hel$Habitat)
	summary(rda.inv.rare.hel.all)
	anova(rda.inv.rare.hel.all)
	RsquareAdj(rda.inv.rare.hel.all)

# Individual factors (there are some issues here)
# Can also easily get Type 1 anova output by terms using formula interface
	anova(rda.inv.rare.hel.all,by="terms")
# But note that for Type 1 anovas, ORDER MATTERS!
# This means that order of predictor variables changes results of tests for i.) unbalanced ANOVA designs or ii.) continuous predictors
rda.inv.rare.hel.all.unb=rda(inv.rare.hel[1:16,11:2046]~inv.rare.hel$Habitat[1:16]*inv.rare.hel$Ecosystem[1:16])
	anova(rda.inv.rare.hel.all.unb,by="terms")
rda.inv.rare.hel.all.unb=rda(inv.rare.hel[1:16,11:2046]~inv.rare.hel$Ecosystem[1:16]*inv.rare.hel$Habitat[1:16])
	anova(rda.inv.rare.hel.all.unb,by="terms")

###
# In contrast, order does not matter in Type 3 anova tests, also called marginal tests
# There are different ways to do this.
# The concept is that all other factors should be accounted for first, and variation is removed.
# Then the test factor is tested against remaining variation
# But this is where the RDA commands are acting weird :(

# Here is the one-line formula interface method.
	anova(rda.inv.rare.hel.all,by="margin")
# If interactions exist, this method only provides marginal tests for the interaction terms
# Let's try testing this using Y,X,Z format
rda.inv.rare.hel.int=rda(inv.rare.hel[,11:2046],inv.rare.hel[,c(9,10)],inv.rare.hel[6:8])
	anova(rda.inv.rare.hel.int)
# This way the same interactions are quite significant
# Another method of getting a marginal test provides similar results.
rda.inv.rare.hel.int=rda(inv.rare.hel[,11:2046] ~ inv.rare.hel$Ecosystem:inv.rare.hel$Habitat + Condition(inv.rare.hel$Ecosystem + inv.rare.hel$Habitat))
	anova(rda.inv.rare.hel.int)

# For completeness, here are tests of main effects using Y,X,Z format
# Obtain type 3 test for ecosystem type
rda.inv.rare.hel.eco=rda(inv.rare.hel[,11:2046],inv.rare.hel[,6:7],inv.rare.hel[8:10])
	anova(rda.inv.rare.hel.eco)
	RsquareAdj(rda.inv.rare.hel.eco)
# Obtain type 3 test for habitat
rda.inv.rare.hel.hab=rda(inv.rare.hel[,11:2046],inv.rare.hel[,8],inv.rare.hel[c(6,7,9,10)])
	anova(rda.inv.rare.hel.hab)
	RsquareAdj(rda.inv.rare.hel.hab)

###
# Plots and getting some more species info
# To get ordination plots of predicted values (group centroids), raw data on canonical axes, and species
plot(rda.inv.rare.hel.all)
# The above plot was not so great.
# Better plots can be made if you extract scores as shown below
rda.inv.rare.hel.all.spscores=scores(rda.inv.rare.hel.all,display="sp",choices=c(1:17))
head(rda.inv.rare.hel.all.spscores)
rda.inv.rare.hel.all.lcscores=scores(rda.inv.rare.hel.all,display="lc",choices=c(1:17))
head(rda.inv.rare.hel.all.lcscores)
rda.inv.rare.hel.all.wascores=scores(rda.inv.rare.hel.all,display="wa",choices=c(1:17))
head(rda.inv.rare.hel.all.wascores)

colors=ifelse(inv.rare.hel$Habitat=="Leaf",3,6)
symbols=ifelse(inv.rare.hel$Ecosys=="BOWO",1,ifelse(inv.rare.hel$Ecosys=="SMBW",2,0))
plot(rda.inv.rare.hel.all.wascores[,1:2],col=colors,pch=symbols,xlab="RDA1 (24%)",ylab="RDA2 (5%)")
points(rda.inv.rare.hel.all.lcscores[,1:2],col=colors,pch=symbols,cex=2,lwd=4)

# This command provides and R-squared for each taxon for each axis
goodness(rda.inv.rare.hel.all)
