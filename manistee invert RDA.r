#new vegan update installed
library(vegan)

# load the data
# you will need to change the path for your computer, or import the data the way you like to
inv.dat= read.csv("manistee invert factorial.csv")
#SM sugar maple varieties
#habitat describes leaf and soil data
#


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
#decostand does a hellinger transformation and preserved each species identity. We then create a distance matrix using the euclidean method (on the hellinger transformed data), so that's how we get hellinger distance

# RDA - Redundancy Analysis 
###
# All factors
#ANOVa design- interactive effects shows multiplicative effects between variables- so we feed it species data and environmental variables
# First, using Y,X,Z (with ,) format
#inv.rare.hel[,6:10] added in the environemtnal data as a factor string- two columns eco1 and eco2- coded as a legendre method where you can multiple eco1 and eco2 to get interactions 1 and 2- then those columns can be used to test for interactions- so its beyond 1s and 0s. this is numeric representations of 3 ecosystems and 1 habitat to look at the interactions between them- you need one less column than categories of factors.
#uses the xyz format- x=species data, y=column 6:10, z=partitioning of varience given at interia and proportion (constrained and unconstrined. So 60% explaned by the constrained(environmental variables), and 40% unconstrained.. great!) This comes form the fitted and unfitted values in the eigen analysis. 
#the eigenvalues block yield RDA1-5 and PC-, ignore the eigenvaue but the proportion explained is the R2 for that axes- the first axis is .3828. So, we get 5 constrained axes out because we loaded in 5 predictor variables- but there isn't a 1:1 correspondance to those variables. PC's are residual axes and are devived from the residuals of the eigen analysis. You can have 1 axes be RDA1 and one axis be PC1, use whatever is the highest. 
#you get the site scores for a given axes by multiplying the species scores by the raw data.-multiply all of the species by all of their species score and add them together for a site to get site scores.
#site constraints multiple species scores by fitted abundances  instead of raw abundances. fitted values are fitted to categories. they are the same for replicates in each category. Thats why the first 3 rows are identical. If you take the average of the first 3 site scores it should equal the replicate site constraint scores for those 3 rows. Blackwood like's plotting these site constraints. variance explained by RDA 1 is coming form those values, its the means coming from those groups. you can put error bars around those means. 
#biplot scores for constrained variables- consider the absolute values for thse- the HAB variable describes the most variation or plays the biggest role, but this is scale dependent..but these are coefficents just liek the specie scores- so 60% is explained in this data set by our factors 
#so, we use ANOVA to determine if the % explained is significant.. it is in this example at 0.01..uses random permutation of the environemntal variables and the species info. so this is saying that the environmental variables are significantly affecting the species... the independent variables are having a significant effect on the dependent variables. 
#so the r squares to the adjusted r squared drops substantially from 60 to 40 because there aren't very many variables- adjusted is the amount of variation we are explaining that we don't expect to explain with this number of variables- "its like the surpirse variation explained" 
rda.inv.rare.hel.all=rda(inv.rare.hel[,11:2046],inv.rare.hel[,6:10])
	summary(rda.inv.rare.hel.all)
	anova(rda.inv.rare.hel.all)
	RsquareAdj(rda.inv.rare.hel.all)
# Second, run this with formula interface (with ~) format
rda.inv.rare.hel.all=rda(inv.rare.hel[,11:2046]~inv.rare.hel$Ecosystem*inv.rare.hel$Habitat)
	summary(rda.inv.rare.hel.all)
	anova(rda.inv.rare.hel.all)
	RsquareAdj(rda.inv.rare.hel.all)
#^this output is the same, but the RDA line ises the ~ instead of the , so instead of using numeric code vegan does behind the scenes recoding...so no need to create the--- as long as all of the downstream things work with this formula...this boils down to the quality control of r being questionable because the 
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
	#the inv.rare.hel8:10 are the independent variables. you can parse out their variance first, and then look  for patterns in the residual variation in inv.rare.hel6:7
	#use varpart command to partition variance. you get an output that shows surprise variance, variance from intereaction can be positive or negative. varpart doesn't tst, only shows variation
	#one philosphy is to not include what's not significant in a model to assess what is more or less signnificant from what was distinguished as significant.
	
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
#^species scores
head(rda.inv.rare.hel.all.spscores)
rda.inv.rare.hel.all.lcscores=scores(rda.inv.rare.hel.all,display="lc",choices=c(1:17))
#^linear constrast scores
head(rda.inv.rare.hel.all.lcscores)
rda.inv.rare.hel.all.wascores=scores(rda.inv.rare.hel.all,display="wa",choices=c(1:17))
#^weighted average scores, closer to raw
head(rda.inv.rare.hel.all.wascores)

colors=ifelse(inv.rare.hel$Habitat=="Leaf",3,6)
symbols=ifelse(inv.rare.hel$Ecosys=="BOWO",1,ifelse(inv.rare.hel$Ecosys=="SMBW",2,0))
plot(rda.inv.rare.hel.all.wascores[,1:2],col=colors,pch=symbols,xlab="RDA1 (24%)",ylab="RDA2 (5%)")
points(rda.inv.rare.hel.all.lcscores[,1:2],col=colors,pch=symbols,cex=2,lwd=4)

# This command provides and R-squared for each taxon for each axis
goodness(rda.inv.rare.hel.all)
