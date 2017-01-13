# GDM Loop - Eric Crandall - 1/9/2017
# I am setting up an analytical loop using code that Cynthia developed
# in the file GDM_for_DIPnet.R. Here are the steps:
# 1. Import the IPDB and PhiST/FST table
# 2. Subsample for each species of interest, and filter based on Phi_ST table.
# 3. Calculate the overwater/great circle distance between all points
# 4. Create dummy distance matrices for each putative "barrier"
# 5. Run through gdm
# 6. Save coefficients
# 7. Loop back to 2.

# Initial config.
library(gdm)
library(plyr)
library(gdistance)
source("config.R")
source("DIPnet_Stats_Functions.R")

######################################################################
# Import the IPDB and Fst tables
ipdb<-read.table(ipdb_path,sep="\t",header=T,stringsAsFactors = F,quote="", na.strings=c("NA"," ","")) 


#read in geographical regionalizations from Treml
spatial<-read.table(spatial_path, header=T, sep="\t",stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")

#read in geographical regionalizations from Beger
spatial2<-read.table(spatial2_path, header=T,sep="\t", stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")

#read in ABGD groups
abgd<-read.table(abgd_path, header=T, sep="\t", stringsAsFactors = F)

#join spatial
ipdb<-join(ipdb,spatial, by = "IPDB_ID",type = "left")
ipdb<-join(ipdb,spatial2[,c(2,18:24)], by = "IPDB_ID", type = "left")

#join ABGD
ipdb<-join(ipdb,abgd[,c(1,3)], by = "IPDB_ID",type = "left")

# drop hybrids and divergent individuals
ipdb<-ipdb[ipdb$IPDB_ID %in% drops == FALSE, ] 

# read in the Fst/PhiSt table 
load("~/google_drive/DIPnet_Gait_Lig_Bird/DIPnet_WG4_first_papers/statistics/By_Species/Pairwise_statistics/sample/DIPnet_structure_060715_WC Theta_sample.Rdata")


######################################################################
# Subsample for each species of interest, and filter based on Phi_ST table.
esu_loci <- unique(ipdb$Genus_species_locus)
all.pops.table<-sapply(esu_loci, function(x) NULL)


for(gsl in esu_loci){ #gsl<-"Linckia_laevigata_CO1" 
  
  cat("\n","\n","\n","Now starting", gsl, "\n")
  
  # subsample IPDB, give it a name and order it alphabetically
  sp<-ipdb[which(ipdb$Genus_species_locus==gsl),]
  sp$sample<-paste(sp$locality,round(sp$decimalLatitude, digits=0),round(sp$decimalLongitude, digits=0),sep="_")  #sets up a variable that matches the name in Fst table
  sp<-sp[order(sp$sample),]
  # Not all localities are included in Veron's regionalization (e.g. Guam), so get their names and then zap NAs
   nonVeronpops<-unique(sp$sample[is.na(sp$VeronDivis)])
   sp<-sp[!is.na(sp$VeronDivis),]
  
  #subsample Fst and filter sp based on the localities that have Fst values
  gslFST<-diffstats[[gsl]]
  sp<-sp[sp$sample %in% labels(gslFST),]
  
  

  
  #create a locations data frame that has all the localities plus lats and longs and their Veron region.
  locs<-as.data.frame(unique(sp$sample))
  names(locs)<-"sample"
  #locs$Long<-sp$decimalLongitude[which(locs %in% sp$sample)]
  #can't do a unique on sample, lats and longs because some samples have non-unique lats and longs! So I do a join and take the first match.
  locs<-join(locs,sp[c("sample","decimalLongitude","decimalLatitude"
                       ,"VeronDivis")], by="sample", match="first")
  
  
  #make a matrix out of gslFST, convert negative values to zero
  gslFSTm<-as.matrix(gslFST)
  gslFSTm[which(gslFSTm<0)]<-0
  #zap the same na populations from the list of non existent pops from VeronDivis
  gslFSTm<-gslFSTm[-(rownames(gslFSTm) %in% nonVeronpops),-(colnames(gslFSTm) %in% nonVeronpops)]
  
  
  # convert to data frame with popsample names as first column
  gslFSTm<-cbind(sample=locs$sample,as.data.frame(gslFSTm))
  
######################################################################
# Calculate Great Circle Distance
  gcdist_km <- pointDistance(locs[,2:3],lonlat=T)/1000
  #cbind on the sample names
  gcdist_km <- cbind(sample=locs$sample,as.data.frame(gcdist_km))
  
####################################################################### Calculate Overwater Distances#
  #Save for later##
#######################################################################
#Create dummy distance matrices for each putative "barrier"
  
  #currently this creates matrices of zero for pairs of sites within a region and 1s for pairs of sites between 2 regions

  #Convert a logical test to 1s and 0s and then take the euclidean distance
  Sunda<-as.matrix(dist(as.numeric(locs$VeronDivis %in% c("Eastern Indian Ocean","Indian Ocean","Red Sea Plus","India"))))
  Sunda<-cbind(sample=locs$sample,as.data.frame(Sunda))
  
  EastCT<-as.matrix(dist(as.numeric(locs$VeronDivis %in% c("Eastern Indian Ocean","Indian Ocean","Red Sea Plus","India","Coral Triangle","Western Australia","Vietnam","South China Sea","Northern South China Sea"))))
  EastCT<-cbind(sample=locs$sample,as.data.frame(EastCT))

  locs$sample<-as.character(locs$sample)
  gslFSTm$sample<-as.character(gslFSTm$sample)
  Sunda$sample<-as.character(Sunda$sample)
  EastCT$sample<-as.character(EastCT$sample)
  gcdist_km$sample<-as.character(gcdist_km$sample)

  
  
#######################################################################
# Run through gdm

#format site pairs into column format
 gdm.format<-formatsitepair(bioData=gslFSTm, bioFormat=3, predData=locs[,1:3],XColumn = "decimalLongitude", YColumn = "decimalLatitude", siteColumn="sample", distPreds=list(gcdist_km,Sunda,EastCT))

 #run the model
gdm.species<-gdm(gdm.format)
 
summary(gdm.species)

#print the plots to a file
pdf(file=paste(gsl,"gdm_plot.pdf",sep=""))
plot(gdm.species)
dev.off()

#Stats to save: relative importance (=sum of coefficients for each predictor) and proportion of deviance explained by model
#See tables and Figures in Fitzpatrick MC, Sanders NJ, Normand S, Svenning JC, Ferrier S, Gove AD, Dunn RR, 2013. Environmental and historical imprints on beta diversity: insights from variation in rates of species turnover along gradients. Proc. Biol. Sci. 280: 20131201–20131201


numb.geog.coeff<-gdm.species$splines[1]
numb.Sunda.coeff<-gdm.species$splines[2]
numb.EastCT.coeff<-gdm.species$splines[3]

impt.geog<-sum(gdm.species$coefficients[1:numb.geog.coeff])
impt.Sunda<-sum(gdm.species$coefficients[(numb.geog.coeff+1) : (numb.geog.coeff+numb.Sunda.coeff)])
impt.EastCT<-sum(gdm.species$coefficients[(numb.geog.coeff+numb.Sunda.coeff+1):(numb.geog.coeff+numb.Sunda.coeff+numb.EastCT.coeff)])

fullmodel.gdm.deviance<- gdm.species$gdmdeviance  #SAVE THIS
fullmodel.explained.deviance<- gdm.species$explained #SAVE THIS
fullmodel.proportion.deviance.explained<- fullmodel.explained.deviance/ fullmodel.gdm.deviance   #probably save - see Fig 2 in Fitzpatrick et al. 2013 
}