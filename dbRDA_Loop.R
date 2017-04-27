# GDM Loop2 - Eric Crandall - 4/20/2017
# I am revamping the original GDM loop that I wrote to follow an analysis 
# outlined by Michelle Gaither and myself.
# 1. Import the IPDB and PhiST/FST table
# 2. Subsample for each species of interest, and filter based on Phi_ST table.
# 3. Calculate the overwater/great circle distance between all points
# 4. Calculate the principal coordinates for FST and geographic distance
# 5. Run through dbRDA and extract statistics
# 6. Loop back to 3.

##################################################################################
# Initial config.
library(gdm)
library(plyr)
library(dplyr)
library(gdistance)
library(vegan)

source("config.R")
source("DIPnet_Stats_Functions.R")

barriers1<-c("Red Sea Plus","Indian Ocean","Eastern Indian Ocean","Coral Triangle","Coral Triangle","Western Australia","North Australia","Coral Triangle","Coral Triangle","Northern South China Sea","Coral Triangle","Northern South China Sea","Central Pacific","Central Pacific","French Polynesia","Northern Japan")
barriers2<-c("Indian Ocean","Eastern Indian Ocean","Coral Triangle","North Australia","Western Australia","North Australia","Central Pacific","South China Sea","Northern South China Sea","Northern Japan","Central Pacific","Central Pacific","French Polynesia","Hawaii","Hawaii","Hawaii")

barriers<-cbind(barriers1,barriers2)

stats<-data.frame(Species_Locus=character(0),Barrier=character(0),constrained.inertia=numeric(0),totalInertia=numeric(0),ProportionConstrained=numeric(0),adj.R2.total=numeric(0),modelF=numeric(0),modelPvalue=numeric(0),pcx_Var=numeric(0),pcx_p=numeric(0),pcy_Var=numeric(0),pcy_p=numeric(0),barrier_Var=numeric(0),barrier_p=numeric(0),stringsAsFactors = F)
stats$Species_Locus<-as.character(stats$Species_Locus)
stats$Barrier<-as.character(stats$Barrier)


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
# 2. Subsample for each species of interest, and filter based on Phi_ST table.
for(gsl in esu_loci){ #gsl<-"Linckia_laevigata_CO1" "Tridacna_crocea_CO1" "Lutjanus_kasmira_CYB"
  
  cat("\n","\n","\n","Now starting", gsl, "\n")
  
  if(any(is.na(diffstats[[gsl]]))){cat("NAs in FST table, No dbRDA calculated"); next}
  
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

  
  ######################################################################
  # 4. Calculate Great Circle Distance
  gcdist_km <- pointDistance(locs[,2:3],lonlat=T)/1000
  # and symmetricise it
  gcdist_km[upper.tri(gcdist_km)]<-0
  gcdist_km<-gcdist_km + t(gcdist_km)
   
  ####################################################################### Calculate Overwater Distances#
  #Save for later##

  ##############################################################################
  #5. Create a subset of the distance matrices including only the localities from
  #    two neighboring Veron regions.

  for(j in 1:16){
    barrier<-c(barriers[j,1],barriers[j,2])
    subset_locs<-which(locs$VeronDivis==barrier[1] | locs$VeronDivis==barrier[2])
    locs2<-locs[subset_locs,]
    
    cat("Now Starting",barrier,"\n")
    
    gcdist_km2<-gcdist_km[subset_locs,subset_locs]
    gslFSTm2<-gslFSTm[subset_locs,subset_locs]
    
    #######################################################################
    # 6. Create a dummy vector for each putative "barrier" 
    #     between the two regions (1s on one side and 0s on the other)
    
    barrier_test<-length(unique(locs2$VeronDivis))>1
    
    #if there aren't samples on either side of this barrier, then go to next barrier
    if(!barrier_test){cat("Not testable \n");next}
    
    
    ############################################################################
    # 7. Calculate the principal coordinates
    
    FST.pcoa<-cmdscale(gslFSTm2, k=dim(as.matrix(gslFSTm2))[1] - 1, eig=TRUE, add=FALSE) 
    FST.scores<-FST.pcoa$points
    
    gcdist.pcoa<-cmdscale(gcdist_km2, k=2, eig=TRUE, add=FALSE)
    gcdist.scores<-gcdist.pcoa$points
    gcdist.scores<-data.frame("pcx"=gcdist.pcoa$points[,1],"pcy"=gcdist.pcoa$points[,2])
    locs2<-cbind(locs2,gcdist.scores)
    
    ###########################################################################
    # 8. Calculate the RDA and extract statistics
    RDA.res<-rda(FST.scores~pcx+pcy+VeronDivis, data=locs2, scale=TRUE )
    
    #Extract Statistics
    constrained.inertia<-summary(RDA.res)$constr.chi
    total.inertia<-summary(RDA.res)$tot.chi
    proportion.constrained.inertia<-constrained.inertia/total.inertia
    
    adj.R2.total.model<-RsquareAdj(RDA.res)$adj.r.squared
    model.sig<-anova.cca(RDA.res, step=1000)
    modelF<-model.sig$F[1]
    modelPvalue<-model.sig$`Pr(>F)`[1]
    
    terms.sig<-anova(RDA.res, by="term", step=1000)
    pcx_Var<-terms.sig$Variance[1]
    pcx_p<-terms.sig$`Pr(>F)`[1]
    pcy_Var<-terms.sig$Variance[2]
    pcy_p<-terms.sig$`Pr(>F)`[2]
    barrier_Var<-terms.sig$Variance[3]
    barrier_p<-terms.sig$`Pr(>F)`[3]
    
    
    
    #save stats
    
    stats_model<-c(gsl,paste(barrier[1],barrier[2],sep="-"),constrained.inertia,total.inertia,proportion.constrained.inertia,adj.R2.total.model,modelF,modelPvalue,pcx_Var,pcx_p,pcy_Var,pcy_p,barrier_Var,barrier_p)
    
    stats[nrow(stats)+1,]<-stats_model
    all.gsl.rda[[gsl]]<-RDA.res
    
    
 
  }

    


  }
  

  
  
}