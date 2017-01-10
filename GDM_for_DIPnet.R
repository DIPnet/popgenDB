##Generalized dissimilarity modeling - Cynthia Riginos
#Ferrier S, Manion G, Elith J, Richardson K, 2007. Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment. Diversity and Distributions 13: 252–264.


library(gdm)

#set working directory to source location
source("config.R")

#Much of "set up" is same as dbRDA.R script - modify both to suit

#read in geographical regionalizations from Eric Treml and reduce to unique fin100 id's
spatial<-read.table(spatial_path, header=T, sep="\t",stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")
fin100<-unique(spatial[,c("fn100id", "fn100_x","fn100_y", "REALM" )])
fin100<-fin100[order(fin100$fn100id),]  #order ascending to make sure it matches later on - might not have to do this
#fin100$fn100id<-as.factor(fin100$fn100id)  #so that dplyr join function does not have a fit later - might be able to skip this here

#To do - turn biogeographic assignment into series of dummy factors
#Hacky soln for now focusing on WIP vs CIP
#Make lists of WIP, CIP fin100 ID's
WIP<-fin100[which(fin100$REALM == "Western Indo-Pacific"),]
CIP<-fin100[which(fin100$REALM == "Central Indo-Pacific"),]


#read in ipdb list of genetics distances (PhiST) - grab from GoogleDrive and place where relevant
load("~/google_drive/DIPnet_Gait_Lig_Bird/DIPnet_WG4_first_papers/statistics/By_Species/Pairwise_statistics/fn100id/DIPnet_structure_060715_PhiST_fn100id.Rdata")
#using fn100.id b/c they have x and y locations
##NOTE - WE SHOULD CHECK WITH EAT TO MAKE SURE THAT THESE PROJECTIONS ARE MOST SENSIBLE TO USE
##ACTUALLY GDM CAN INCORPORATE ANY DISTANCE MEASURE SO MAYBE WE SHOULD USE OVER WATER AS A NULL....

species<-names(diffstats)

#To do - SET UP LOOPING STRUCTURE - LOOP OVER ALL SPECIES

#Example: Subset for one species 
s=180  #Linckia
single_species<-diffstats[["Linckia_laevigata_CO1"]]

##PREP FOR GDM - format predictive variables
##For now, just using spatial locations and biogeographic assignment, BUT distance measures and environmental rasters can be used
single_species.mat<-as.matrix(single_species)
single_species.df<-as.data.frame(single_species.mat)
single_species.df[single_species.df<0]<-0   #replace negative values with zero

#now subset within both sites and single_species.df to match the biogeographic regions
#yes.WIP<-which(row.names(single_species.df) %in% WIP$fn100id)
#yes.CIP<-which(row.names(single_species.df) %in% CIP$fn100id)

WIP_CIP<-which(row.names(single_species.df) %in% WIP$fn100id | row.names(single_species.df) %in% CIP$fn100id)
WIP_CIP.df<-single_species.df[WIP_CIP,WIP_CIP]
#sites$fn100id<-as.integer(as.character(sites$fn100id))

#subset fin100 to match locations for focal species and biogeographic regions
WIP_CIP.sites<-fin100[which(fin100$fn100id %in% row.names(WIP_CIP.df)),1:4] 
#set up dummy variables for biogeography
WIP_CIP.sites[WIP_CIP.sites$REALM == "Western Indo-Pacific","REALM"]<-0
WIP_CIP.sites[WIP_CIP.sites$REALM == "Central Indo-Pacific","REALM"]<-1
WIP_CIP.sites$REALM<-as.numeric(WIP_CIP.sites$REALM)

##PREP FOR GDM - format biological distances
#Using biological distance matrics - pairwise PhiSt - (bioData format #3: see GDM documentation)
fn100id<-as.integer(as.character(WIP_CIP.sites$fn100id))
speciesDissim<-cbind(fn100id, WIP_CIP.df)

#format site pairs and undertake GDM
#with biogeography
gdmTab.species <- formatsitepair(speciesDissim, bioFormat=3, XColumn="fn100_x", YColumn="fn100_y", predData=WIP_CIP.sites, siteColumn="fn100id")
gdm.species <- gdm(gdmTab.species, geo=T)
summary(gdm.species) #for the Linkia example WIPvsCIP biogeography not significant: all I spline coefficients = 0
plot(gdm.species)  #for Linkia, a biogeography plot is not produced b/c it is not a signif predictor
str(gdm.species)

#Stats to save: relative importance (=sum of coefficients for each predictor) and proportion of deviance explained by model
#See tables and Figures in Fitzpatrick MC, Sanders NJ, Normand S, Svenning JC, Ferrier S, Gove AD, Dunn RR, 2013. Environmental and historical imprints on beta diversity: insights from variation in rates of species turnover along gradients. Proc. Biol. Sci. 280: 20131201–20131201.
numb.geog.coeff<-gdm.species$splines[1]
numb.biog.coeff<-gdm.species$splines[2]

impt.geog<-sum(gdm.species$coefficients[1:numb.geog.coeff])   #SAVE THIS
impt.biog<-sum(gdm.species$coefficients[paste(numb.geog.coeff+1) : paste(numb.geog.coeff+numb.biog.coeff)])  #SAVE THIS

fullmodel.gdm.deviance<- gdm.species$gdmdeviance  #SAVE THIS
fulllmodel.explained.deviance<- gdm.species$explained #SAVE THIS
fullmodel.proportion.deviance.explained<- fullmodel.explained.deviance/ fullmodel.gdm.deviance   #probably save - see Fig 2 in Fitzpatrick et al. 2013 

#without biogeography
gdmTab.species.d <- formatsitepair(speciesDissim, bioFormat=3, XColumn="fn100_x", YColumn="fn100_y", predData=WIP_CIP.sites[,1:3], siteColumn="fn100id")
gdm.species.d <- gdm(gdmTab.species.d, geo=T)

numb.geog.coeff.d<-gdm.species.d$splines[1]

impt.geog.partialmodel<-sum(gdm.species.d$coefficients[1:numb.geog.coeff])   #SAVE THIS

partialmodel.gdm.deviance<- gdm.species.d$gdmdeviance  #SAVE THIS - should be same as full model
partialmodel.explained.deviance<- gdm.species.d$explained #SAVE THIS - for Linkia is same as full model but this will not usually be the case
partialmodel.proportion.deviance.explained<- partialmodel.explained.deviance/ partialmodel.gdm.deviance   #probably save - see Fig 2 in Fitzpatrick et al. 2013 



#Future - can use splines to interpolate genetic dissimilarities - could we add up rasters across species?


