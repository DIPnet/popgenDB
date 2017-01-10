# GDM Loop - Eric Crandall - 1/9/2017
# I am setting up an analytical loop using code that Cynthia developed
# in the file GDM_for_DIPnet.R. Here are the steps:
# 1. Import the IPDB and PhiST/FST table
# 2. Subsample for each species of interest, and filter based on Phi_ST table.
# 3. Calculate the overwater distance between all points
# 4. Create dummy distance matrices for each putative "barrier"
# 5. Run through gdm
# 6. Save coefficients
# 7. Loop back to 2.

# Initial config.
library(gdm)
library(plyr)
source("config.R")

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
  locs<-unique(sp$sample)
  
  
  
  
  
}