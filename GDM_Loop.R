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
load("/Users/eric/google_drive/DIPnet_Gait_Lig_Bird/DIPnet_WG4_first_papers/statistics/By_Species/Pairwise_statistics/sample/DIPnet_structure_060715_WC Theta_sample.Rdata")


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
  
  #subsample Fst and filter sp based on the localities that have Fst values
  gslFST<-diffstats[[gsl]]
  sp<-sp[sp$sample %in% labels(gslFST),]
  
  
  # Not all localities are included in Veron's regionalization (e.g. Guam), so zap NAs
  sp<-sp[!is.na(sp$VeronDivis),]
  
  #create a locations data frame that has all the localities plus lats and longs and their Veron region.
  locs<-as.data.frame(unique(sp$sample))
  names(locs)<-"sample"
  locs$Long<-sp$decimalLongitude[which(locs %in% sp$sample)]
  #can't do a unique on sample, lats and longs because some samples have non-unique lats and longs! So I do a join and take the first match.
  locs<-join(locs,sp[c("sample","decimalLongitude","decimalLatitude"
                       ,"VeronDivis")], by="sample", match="first")

  #make a matrix out of gslFST with the first column sample names
  gslFSTm<-as.matrix(c(locs$sample,gslFST),nrow=length(locs$sample))
  gslFSTm<-cbind(locs$sample,gslFST)
  
######################################################################
# Calculate Great Circle Distance
  gcdist_km <- pointDistance(locs[,2:3],lonlat=T)/1000
  
####################################################################### Calculate Overwater Distances#
  #Save for later##
#######################################################################
#Create dummy distance matrices for each putative "barrier" 
  #currently this creates matrices of zero for pairs of sites within a region and 1s for pairs of sites between 2 regions
  #get the Veron regions for this species
  unique(locs$VeronDivis) 
  
  #Convert a logical test to 1s and 0s and then take the euclidean distance
  WA<-as.matrix(dist(as.numeric(locs$VeronDivis=="Western Australia")))
  CT<-as.matrix(dist(as.numeric(locs$VeronDivis=="Coral Triangle")))
  EIO<-as.matrix(dist(as.numeric(locs$VeronDivis=="Eastern Indian Ocean")))
  CP<-as.matrix(dist(as.numeric(locs$VeronDivis=="Central Pacific")))
  

  
  
#######################################################################
# Run through gdm
  gdmTab.species <- formatsitepair(speciesDissim, bioFormat=3, XColumn="fn100_x", YColumn="fn100_y", predData=WIP_CIP.sites, siteColumn="fn100id")
  
 gdm.format<-formatsitepair(bioData=gslFST, bioFormat=3, predData=locs,XColumn = "decimalLongitude", YColumn = "decimalLatitude", siteColumn="sample", distPreds=c(WA,CT,EIO,CP))
}