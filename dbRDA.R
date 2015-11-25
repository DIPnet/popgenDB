# Example of dbRDA approach applied to Linkia using PhiST data based on ecoregions
library(vegan)

source("../config.R")

setwd("~/Documents/Cynthia-Atlas/Git_folders/DIPnet/popgenDB")

#read in geographical regionalizations from Treml
spatial<-read.table(spatial_path, header=T, sep="\t",stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")

#read in ipdb list of genetics distances (PhiST)
load("~/Desktop/DIPnet_structure_ecoregions_PhiST_071615.RData") 

#Subset for one species 
Linkia<-(diffstats$Linckia_laevigata_CO1)

#PRINCIPAL COORDINATES OF GENETIC DISTANCE MATRIX
PhiSt.pcoa<-cmdscale(Linkia, k=dim(as.matrix(Linkia))[1] - 1, eig=TRUE, add=TRUE)   #k = 23 for Linckia: 24 sites    #add=TRUE sets the Cailliez correction to avoid negative eigenvalues
PhiSt.scores<-PhiSt.pcoa$points  

#CREATE DISTANCE MATRIX FOR SPATIAL DATA AND FIND PRINCIPAL COORDINATES
#This is VERY HACKY. We obviously want much more sophisticated spatial predictors!! Here I am just yanking out longs and lats for representative pops in the ecoregion
#NOTE - important that the site locations in the geographic distance matrix follow EXACTLY the same order as the genetic distance matrix

sites<-row.names(as.matrix(Linkia))

ll<-unique(spatial[which(spatial$Genus_spec == "Linckia_laevigata_CO1"),c("ECOREGION", "decimalLat","decimalLon" )])
ll.sites<-data.frame(row.names=sites)
for (r in 1:nrow(ll.sites)) {
  ll.sites[r,"Latitude"]<-ll[ll$ECOREGION==row.names(ll.sites)[r],"decimalLat"][1]
  ll.sites[r,"Longitude"]<-ll[ll$ECOREGION==row.names(ll.sites)[r],"decimalLon"][1]
}

library(fossil)
geog.dists<-earth.dist(ll.sites)  #over earth distances..... not what we want to use!

#make the pcoa
#once we have geographic predictors we will have to examine to decide how many PCs to retain
geog.pcoa<-cmdscale(geog.dists, k=2, eig=TRUE, add=TRUE)     #add=TRUE sets the Cailliez correction to avoid negative eigenvalues; need to think about how many PCs are reasonable to retain
geog.scores<-geog.pcoa$points


#WE CAN ADD OTHER PREDICTIVE DISTANCE MATRICES 
#BARRIERS: use dummy variables: 0 for same side of barrier, 1 for different sides... set up as a distance matrix
#If we have many predictive variables, we will want to test for and deal with colinearity


##RDA
linckia.rda<-rda(PhiSt.scores,geog.scores)  #note that in vegan the Y matrix is the constraining matrix
#in reality we probably want geographic distance (overwater?) to be the conditioning matrix and to have additional predictive matrices
summary(linckia.rda)

#comments on this simple example: 
# The RDA eigenvectors reflect the amount of constrained variance explained, in this instance 48% (like R2 but bias correction follows later)
# Because PC1 is greater than RDA2, I would only take RDA1 seriously (41% of variance)
# "Species scores" - unclear in this context
# "Site scores" -  maps sites onto genetic PCoA axes (onto respose eigenvectors)
# "Site constraints"-  maps sites onto geographic PCoA axes (onto predictive eigenvectors)

#regresssion coefficients between explanatory variables and canonical axes
coef(linckia.rda)

#unadjusted R2 of model
R2<-RsquareAdj(linckia.rda)$r.squared

#adjusted R2
R2adj<-RsquareAdj(linckia.rda)$adj.r.squared

#ploting RDA results - 1 distance triplot
plot(linckia.rda, scaling =1)
#locations are mapped on by name, in red are the genetic PCoA axes, in blue are the geographic PCoA axes

#ploting RDA results - 2 corellation triplot
plot(linckia.rda, scaling =2)
#these two different ways of plotting the data confuse me and if I gain clarity on the subject I will let you know!

#Testing SIGNIFICANCE of RDA results
#global test
anova.cca(linckia.rda, step=1000)
#by canoncial axes
anova.cca(linckia.rda, by="axis", step=1000)  #why is this not working?



