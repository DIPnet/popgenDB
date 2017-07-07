# GDM Loop2 - Eric Crandall - 4/20/2017
# I am revamping the original GDM loop that I wrote to follow an analysis 
# outlined by Michelle Gaither and myself.

# 1. Import the IPDB and PhiST/FST table
# 2. Subsample for each species of interest, and filter based on Phi_ST table.
# 3. Calculate the overwater/great circle distance between all points
# 4. Create a subset of the distance matrices including only the localities from
#    two neighboring Veron regions.
# 5. Create a dummy distance matrix for each putative "barrier" 
#     between the two regions (1s and 0s)
# 6. Run through gdm with the barrier and without. Save the deviance values.
# 7. Perform Monte-Carlo permutations to develop a null distribution 
#    of deviance values and determine significance (pvalue)
# 8. Loop back to 4, over all putative barriers.
# 9. Loop back to 2.

# Initial config.
library(gdm)
library(plyr)
library(gdistance)
source("config.R")
source("DIPnet_Stats_Functions.R")

barriers<-read.csv("VeronBarriers.csv",header=F,stringsAsFactors = F)

stats<-data.frame(Species_Locus=character(0),Barrier=character(0),WithBarrierDeviance=numeric(0),WithBarrierExplainedDeviance=numeric(0),ImportanceDistanceWithBarrier=numeric(0),ImportanceBarrierWithBarrier=numeric(0),NoBarrierDeviance=numeric(0),NoBarrierExplainedDeviance=numeric(0),ImportanceDistanceWithoutBarrier=numeric(0),DeltaDeviance=numeric(0),Pvalue=numeric(0),stringsAsFactors = F)




######################################################################
# 1. Import the IPDB and Fst tables
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
load("~/google_drive/DIPnet_Gait_Lig_Bird/DIPnet_WG4_first_papers/statistics/By_Species/Pairwise_statistics/sample/DIPnet_structure_sample_PhiST_042817.Rdata")


# Make an empty list to save gdm output for each species
esu_loci <- unique(ipdb$Genus_species_locus)
all.pops.table<-sapply(esu_loci, function(x) NULL)
solution<-list()
nosolution<-list()

###############################################################################
# 2. Subsample for each species of interest, and filter based on Phi_ST table.
for(gsl in esu_loci){ #gsl<-"Linckia_laevigata_CO1" "Tridacna_crocea_CO1" "Lutjanus_kasmira_CYB" "Acanthaster_planci_CO1"
  
  cat("\n","\n","\n","Now starting", gsl, "\n")
  
  if(any(is.na(diffstats[[gsl]]))){cat("NAs in FST table, No gdm calculated"); next}
  
  if(diffstats[[gsl]]=="Fewer than 3 sampled populations after filtering. No stats calculated"){all.pops.table[[gsl]]<-"Fewer than 5 sampled populations after filtering."; cat("Fewer than 5 sampled populations after filtering.");next}
  
  #pull out the data for this genus-species-locus (gsl)
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
  gslFSTm[which(gslFSTm<=0)]<-0.00001
  
  #zap weird slashes in the names
  rownames(gslFSTm)<-gsub("\"","",rownames(gslFSTm))
  colnames(gslFSTm)<-rownames(gslFSTm)
  
  #zap the same na populations from the list of non existent pops from VeronDivis
  if(any(rownames(gslFSTm) %in% nonVeronpops)){
    gslFSTm<-gslFSTm[-(which(rownames(gslFSTm) %in% nonVeronpops)),-(which(colnames(gslFSTm) %in% nonVeronpops))]
  }
  if(length(rownames(gslFSTm))<5){all.pops.table[[gsl]]<-"Fewer than 5 sampled populations";cat("Fewer than 5 sampled populations");next}
  
  #and filter sp based on the localities that have Fst values
  sp<-sp[sp$sample %in% rownames(gslFSTm),]
  
  #and vice versa
  
  gslFSTm<- gslFSTm[which(rownames(gslFSTm) %in% unique(sp$sample)),which(rownames(gslFSTm) %in% unique(sp$sample))]
  

  if(length(rownames(gslFSTm))<5){all.pops.table[[gsl]]<-"Fewer than 5 sampled populations";cat("Fewer than 5 sampled populations");next}
  
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
  gslFSTm<-cbind(sample=locs$sample,as.data.frame(gslFSTm))
  
  ######################################################################
  # 3. Calculate Great Circle Distance
  gcdist_km <- pointDistance(locs[,2:3],lonlat=T)/1000
  #symmetricize the matrix
  gcdist_km[upper.tri(gcdist_km)]<-t(gcdist_km)[upper.tri(gcdist_km)]
  
  #cbind on the sample names
  gcdist_km <- cbind(sample=locs$sample,as.data.frame(gcdist_km))
  
  ####################################################################### Calculate Overwater Distances#
  #Save for later##
  #######################################################################
  # 4. Create a subset of the distance matrices including only the localities from
  #    two neighboring Veron regions.
  for(j in 1:16){
    barrier<-c(barriers[j,1],barriers[j,2])
    subset_locs<-which(locs$VeronDivis==barrier[1] | locs$VeronDivis==barrier[2])
    locs2<-locs[subset_locs,]
    
    cat("Now Starting",barrier,"\n")
    
    gcdist_km2<-gcdist_km[subset_locs,c(1,subset_locs+1)]
    gslFSTm2<-gslFSTm[subset_locs,c(1,subset_locs+1)]
    
    #######################################################################
    # 5. Create a dummy distance matrix for each putative "barrier" 
    #     between the two regions (1s and 0s)
    
    barrier_m2<-as.matrix(dist(as.numeric(locs2$VeronDivis %in% barrier[1])))
    barrier_test<-sum(barrier_m2)>0
    
    barrier_m2 <- cbind(sample=locs2$sample,as.data.frame(barrier_m2))
    
    #if there aren't samples on either side of this barrier, then go to next barrier
    if(!barrier_test){cat("Not testable \n");next}
    
    ############################################################################
    # 6. Run through gdm with the barrier and without. Save the deviance values.
    
    locs2$sample<-as.character(locs2$sample)
    gslFSTm2$sample<-as.character(gslFSTm2$sample)
    gcdist_km2$sample<-as.character(gcdist_km2$sample)
    
    gdm.format<-formatsitepair(bioData=gslFSTm2, bioFormat=3, predData=locs2[,1:3],XColumn = "decimalLongitude", YColumn = "decimalLatitude", siteColumn="sample", distPreds=list(gcdist_km2,barrier_m2))
    
    #run gdm with and without the barrier
    gdm.barrier<-gdm(gdm.format)
    gdm.no.barrier<-gdm(gdm.format[,-grep("matrix_2",names(gdm.format))])
    
    #TROUBLESHOOTING: save fst matrices from gdm models that obtain no solution
    if(is.null(gdm.barrier) | is.null(gdm.no.barrier))
      {cat("No Solution Obtained \n");
      nosolution[[paste(gsl,barrier[1],barrier[2],sep=",")]]<-list(locs2,gdm.format);
      next}
    
    #difference in deviance is the more complex model - less complex model
    deltadev<-gdm.no.barrier$gdmdeviance-gdm.barrier$gdmdeviance
    
    ##############################################################################
    # 7. Perform Monte-Carlo permutations to develop a null distribution 
    #    of deviance values and determine significance (pvalue)
    gdm.format.rand<-gdm.format
    rand.deltas<-vector() 
    
    while(length(rand.deltas) < 100) {
      gdm.format.rand$distance<-sample(gdm.format.rand$distance,size=length(gdm.format.rand$distance))
      gdm.barrier.rand<-gdm(gdm.format.rand)
      gdm.no.barrier.rand<-gdm(gdm.format.rand[,-grep("matrix_2",
                                                      names(gdm.format))])
      if(is.null(gdm.barrier.rand) | is.null(gdm.no.barrier.rand)){next}
      deltadev.rand<-gdm.no.barrier.rand$gdmdeviance-gdm.barrier.rand$gdmdeviance
      rand.deltas<-c(rand.deltas,deltadev.rand)
    }
    pvalue<-length(which(abs(deltadev) < abs(rand.deltas)))/length(rand.deltas)
    
    #TROUBLESHOOTING: code to keep track of which ones worked
    cat("Good Solution \n");
    solution<-rbind(solution,c(gsl,barrier))
    
    #save stats
    gdm.barrier.deviance<-gdm.barrier$gdmdeviance
    gdm.barrier.explained<-gdm.barrier$explained

    gdm.no.barrier.deviance<-gdm.no.barrier$gdmdeviance
    gdm.no.barrier.explained<-gdm.no.barrier$explained

    impt.dist.gdm.barrier<-sum(gdm.barrier$coefficients[1:gdm.barrier$splines[1]])
    impt.barrier.gdm.barrier<-sum(gdm.barrier$coefficients[gdm.barrier$splines[1]+1:gdm.barrier$splines[1]])
    impt.dist.gdm.no.barrier<-sum(gdm.no.barrier$coefficients[1:gdm.no.barrier$splines[1]])
   
     stats_model<-c(gsl,paste(barrier[1],barrier[2],sep="-"),gdm.barrier.deviance,gdm.barrier.explained, impt.dist.gdm.barrier, impt.barrier.gdm.barrier, gdm.no.barrier.deviance, gdm.no.barrier.explained,impt.dist.gdm.no.barrier, deltadev,pvalue)
    
    stats[nrow(stats)+1,]<-stats_model
  }
  

  
  
}


# do some figuring with the results - see how often each barrier is significant at alpha of 0.05 - will get more sophisticated later.
library(dplyr)

goodbarriers<-stats %>% filter(Pvalue < 0.05) %>% group_by(Barrier) %>% summarize(goodbarriers = n())

allbarriers<-stats %>% group_by(Barrier) %>% summarize(barrier_tests = n())

barrier_ratios<-left_join(allbarriers,goodbarriers,by="Barrier")

barrier_ratios$goodbarriers[which(is.na(barrier_ratios$goodbarriers))]<-0

barrier_ratios<-mutate(barrier_ratios, goodbarriers/barrier_tests)
