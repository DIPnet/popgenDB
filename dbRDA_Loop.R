# GDM Loop2 - Eric Crandall - 4/20/2017
# I am revamping the original GDM loop that I wrote to follow an analysis 
# outlined by Michelle Gaither and myself.
# 1. Create a dummy column for each putative "barrier" (1s and 0s)
# 2. Import the IPDB and PhiST/FST table
# 3. Subsample for each species of interest, and filter based on Phi_ST table.
# 4. Calculate the overwater/great circle distance between all points
# 5. Calculate the principal coordinates for FST and geographic distance
# 6. Run through dbRDA and extract statistics
# 7. Loop back to 3.

##################################################################################
# Initial config.
library(gdm)
library(plyr)
library(dplyr)
library(gdistance)
library(vegan)

source("config.R")
source("DIPnet_Stats_Functions.R")

stats<-data.frame(Species_Locus=character(0),Barrier=character(0),WithBarrierDeviance=numeric(0),WithBarrierExplainedDeviance=numeric(0),ImportanceDistanceWithBarrier=numeric(0),ImportanceBarrierWithBarrier=numeric(0),NoBarrierDeviance=numeric(0),NoBarrierExplainedDeviance=numeric(0),ImportanceDistanceWithoutBarrier=numeric(0),DeltaDeviance=numeric(0),Pvalue=numeric(0),stringsAsFactors = F)

#######################################################################
# 1. Create a dummy column for each putative "barrier" (1s and 0s)
Indian_Ocean<-as.numeric(locs$VeronDivis %in% c("Indian Ocean","Red Sea Plus","India"))

Sunda<-as.numeric(locs$VeronDivis %in% c("Eastern Indian Ocean","Indian Ocean","Red Sea Plus","India"))

W_Australia<-as.numeric(locs$VeronDivis %in% c("Western Australia"))
N_Australia<-as.numeric(locs$VeronDivis %in% c("North Australia"))
SCS<-as.numeric(locs$VeronDivis %in% c("South China Sea"))
Vietnam<-as.numeric(locs$VeronDivis %in% c("Vietnam"))
Japan<-as.numeric(locs$VeronDivis %in% c("Northern Japan"))
NSCS<-as.numeric(locs$VeronDivis %in% c("Northern South China Sea","Northern Japan"))

East_CT<-as.numeric(locs$VeronDivis %in% c("Eastern Indian Ocean","Indian Ocean","Red Sea Plus","India","Coral Triangle","Western Australia","North Australia","Vietnam","South China Sea","Northern South China Sea","Northern Japan"))
Central_Pac<-as.numeric(locs$VeronDivis %in% c("Eastern Indian Ocean","Indian Ocean","Red Sea Plus","India","Coral Triangle","Western Australia","North Australia","Vietnam","South China Sea","Northern South China Sea","Northern Japan","Central Pacific",""))
Hawaii<-as.numeric(locs$VeronDivis %in% c("Hawaii"))


######################################################################
# 2. Import the IPDB and Fst tables
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


# Make an empty list to save gdm output for each species
esu_loci <- unique(ipdb$Genus_species_locus)
all.gsl.rda<-sapply(esu_loci, function(x) NULL)



###############################################################################
# 3. Subsample for each species of interest, and filter based on Phi_ST table.
for(gsl in esu_loci){ #gsl<-"Linckia_laevigata_CO1" "Tridacna_crocea_CO1" "Lutjanus_kasmira_CYB"
  
  cat("\n","\n","\n","Now starting", gsl, "\n")
  
  if(any(is.na(diffstats[[gsl]]))){cat("NAs in FST table, No gdm calculated"); next}
  
  if(diffstats[[gsl]]=="Fewer than 3 sampled populations after filtering. No stats calculated"){all.gsl.rda[[gsl]]<-"Fewer than 5 sampled populations after filtering."; cat("Fewer than 5 sampled populations after filtering.");next}
  sp<-ipdb[which(ipdb$Genus_species_locus==gsl),]
  #clean weird backslashes from names
  sp$locality<-gsub("\"","",sp$locality)
  
  sp$sample<-paste(sp$locality,round(sp$decimalLatitude, digits=0),round(sp$decimalLongitude, digits=0),sep="_")  #sets up a variable that matches the name in Fst table
  sp<-sp[order(sp$sample),]
  # Not all localities are included in Veron's regionalization (e.g. Guam), so get their names and then zap NAs
  nonVeronpops<-unique(sp$sample[is.na(sp$VeronDivis)])
  sp<-sp[!is.na(sp$VeronDivis),]
  
  #subsample Fst 
  gslFST<-diffstats[[gsl]]
  #make a matrix out of gslFST, convert negative values to zero
  gslFSTm<-as.matrix(gslFST)
  gslFSTm[which(gslFSTm<0)]<-0.0
  
  #zap weird slashes in the names
  rownames(gslFSTm)<-gsub("\"","",rownames(gslFSTm))
  colnames(gslFSTm)<-rownames(gslFSTm)
  
  #zap the same na populations from the list of non existent pops from VeronDivis
  if(any(rownames(gslFSTm) %in% nonVeronpops)){
    gslFSTm<-gslFSTm[-(which(rownames(gslFSTm) %in% nonVeronpops)),-(which(colnames(gslFSTm) %in% nonVeronpops))]
  }
  if(length(rownames(gslFSTm))<5){all.gsl.rda[[gsl]]<-"Fewer than 5 sampled populations";cat("Fewer than 5 sampled populations");next}
  
  #and filter sp based on the localities that have Fst values
  sp<-sp[sp$sample %in% rownames(gslFSTm),]
  
  #and vice versa
  
  gslFSTm<- gslFSTm[which(rownames(gslFSTm) %in% unique(sp$sample)),which(rownames(gslFSTm) %in% unique(sp$sample))]
  

  if(length(rownames(gslFSTm))<5){all.gsl.rda[[gsl]]<-"Fewer than 5 sampled populations";cat("Fewer than 5 sampled populations");next}
  
  #create a locations data frame that has all the localities plus lats and longs and their Veron region.
  locs<-as.data.frame(unique(sp$sample))
  names(locs)<-"sample"
  #locs$Long<-sp$decimalLongitude[which(locs %in% sp$sample)]
  #can't do a unique on sample, lats and longs because some samples have non-unique lats and longs! So I do a join and take the first match.
  locs<-join(locs,sp[c("sample","decimalLongitude","decimalLatitude"
                       ,"VeronDivis")], by="sample", match="first")
  
  #sort gslFSTm
  gslFSTm<-gslFSTm[order(rownames(gslFSTm)),order(colnames(gslFSTm))]
  # convert to data frame with popsample names as first column
  #gslFSTm<-cbind(sample=locs$sample,as.data.frame(gslFSTm))
  
  ######################################################################
  # 4. Calculate Great Circle Distance
  gcdist_km <- pointDistance(locs[,2:3],lonlat=T)/1000
  # and symmetricise it
  gcdist_km[upper.tri(gcdist_km)]<-0
  gcdist_km<-gcdist_km + t(gcdist_km)
   
  ####################################################################### Calculate Overwater Distances#
  #Save for later##


    

    
    ############################################################################
    # 5. Calculate the principal coordinates
    
    phiST.pcoa<-cmdscale(gslFSTm, k=dim(as.matrix(gslFSTm))[1] - 1, eig=TRUE, add=FALSE) 
    phiST.scores<-phiST.pcoa$points
    
    gcdist.pcoa<-cmdscale(gcdist_km, k=2, eig=TRUE, add=FALSE)
    gcdist.scores<-gcdist.pcoa$points
    gcdist.scores<-data.frame("pcx"=gcdist.pcoa$points[,1],"pcy"=gcdist.pcoa$points[,2])
    locs<-cbind(locs,gcdist.scores,Indian_Ocean,W_Australia,N_Australia,SCS,NSCS,Vietnam,Japan,Sunda,EastCT,Central_Pac,Hawaii)
    
    ###########################################################################
    # 6. Calculate the RDA and extract statistics
    RDA.res<-rda(phiST.scores~pcx+pcy+Indian_Ocean+Sunda+W_Australia+N_Australia+EastCT+SCS+Central_Pac+NSCS+Japan+Vietnam+Hawaii, data=locs, scale=TRUE )

    anova(RDA.res, by="term", step=1000)
    
    
    #Extract Statistics
    constrained.inertia<-summary(RDA.res)$constr.chi
    total.inertia<-summary(RDA.res)$tot.chi
    proportion.constrained.interia<-constrained.inertia/total.inertia
    
    adj.R2.total.model<-RsquareAdj(RDA.res)$adj.r.squared
    model.sig<-anova.cca(RDA.res, step=1000)
    terms.sig<-anova(RDA.res, by="term", step=1000)
    
    #save stats
     stats_model<-c(gsl,paste(barrier[1],barrier[2],sep="-"),gdm.barrier.deviance,gdm.barrier.explained, impt.dist.gdm.barrier, impt.barrier.gdm.barrier, gdm.no.barrier.deviance, gdm.no.barrier.explained,impt.dist.gdm.no.barrier, (gdm.barrier.deviance-gdm.no.barrier.deviance),pvalue)
    
    stats[nrow(stats)+1,]<-stats_model
  }
  

  
  
}