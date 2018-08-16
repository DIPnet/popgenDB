# Cynthia Riginos - November 2016
# USING dbRDA TO INFER THE CONTRIBUTION OF BIOGEOGRAPHY AND IBD TO GENETIC STRUCTURE
# For now, using REALMS to define biogeographic hypotheses - Crandall to expand to other Ho's
# Also, focusing only on WIP vs CIP as a demonstration


library(dplyr)
library(vegan)

#set working directory to source location
source("config.R")

#read in geographical regionalizations from Eric Treml and reduce to unique fin100 id's
spatial<-read.table(spatial_path, header=T, sep="\t",stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")
fin100<-unique(spatial[,c("fn100id", "fn100_x","fn100_y", "REALM" )])
fin100<-fin100[order(fin100$fn100id),]  #order ascending to make sure it matches later on
fin100$fn100id<-as.factor(fin100$fn100id)  #so that dplyr join function does not have a fit later

#Make lists of WIP, CIP fin100 ID's
WIP<-droplevels(fin100[which(fin100$REALM == "Western Indo-Pacific"),])
CIP<-fin100[which(fin100$REALM == "Central Indo-Pacific"),]
#EIP<-fin100[which(fin100$REALM == "Eastern Indo-Pacific"),]

#read in ipdb list of genetics distances (PhiST) - grab from GoogleDrive 
load("~/google_drive/DIPnet_Gait_Lig_Bird/DIPnet_WG4_first_papers/statistics/By_Species/Pairwise_statistics/fn100id/DIPnet_structure_060715_PhiST_fn100id.Rdata")   #using fn100.id b/c they have x and y locations
##NOTE - WE SHOULD CHECK WITH EAT TO MAKE SURE THAT THESE PROJECTIONS ARE MOST SENSIBLE TO USE

species<-names(diffstats)

#To do - SET UP LOOPING STRUCTURE - LOOP OVER ALL SPECIES

#Example: Subset for one species 
s=180  #Linckia
single_species<-diffstats[["Linckia_laevigata_CO1"]] 

#Figure out which fin100 locations are in which biogeographic areas
#This initial mockup is focusing only on WIP-CIP 
#In order to subset the distance object, it is turned into a matrix, subsetted, and the turned back to a distance
single_species.matrix<-as.matrix(single_species)

yes.WIP<-which(row.names(single_species.matrix) %in% WIP$fn100id)
yes.CIP<-which(row.names(single_species.matrix) %in% CIP$fn100id)
WIP.CIP<-which(row.names(single_species.matrix) %in% WIP$fn100id | row.names(single_species.matrix) %in% CIP$fn100id)

#only proceed if minimum two pops per realm
#WIP vs CIP - Subset distance object, undertake PCOA, then dbRDA (following Legendre & Anderson 1999)
if (length(yes.WIP)>1 & length(yes.CIP)>1) {
  WIP_CIP<-single_species.matrix[WIP.CIP,WIP.CIP]
  WIP_CIP.dist<-as.dist(WIP_CIP)
  
  #PCOA
  PhiSt.pcoa<-cmdscale(WIP_CIP.dist, k=dim(as.matrix(single_species))[1] - 1, eig=TRUE, add=FALSE)   #k = 52 for Linckia: 53 sites   
  #ignore warnings for negative eigenvalues
  PhiSt.scores<-PhiSt.pcoa$points  #note that these are sorted by fn100id; these will be the response variables in the RDA
  
  #Set up for RDA
  sites<-as.data.frame(row.names(WIP_CIP))
  colnames(sites)<-"fn100id"
  sites<-left_join(sites, fin100)  #ignore warning about coercion - this means that some factors are being dropped from fin100
  RDA.res<-rda(PhiSt.scores~fn100_x+fn100_y+REALM, data=sites, scale=TRUE )
  # summary(RDA.res)
  # (R2adj <- RsquareAdj(RDA.res)$adj.r.squared)
  # anova.cca(RDA.res, step=1000)
  # anova(RDA.res, by="term", step=1000)
  
  #store results as list
  constrained.inertia<-summary(RDA.res)$constr.chi
  total.inertia<-summary(RDA.res)$tot.chi
  proportion.constrained.interia<-constrained.inertia/total.inertia
  adj.R2.total.model<-RsquareAdj(RDA.res)$adj.r.squared
  model.sig<-anova.cca(RDA.res, step=1000)
  terms.sig<-anova(RDA.res, by="term", step=1000)
  
  #need to create list of results with species name as key - EC, can you do this?
  
}

#THEN PROCEED WITH OTHER BIOGEOGRAPHIC DIVISIONS, etc





