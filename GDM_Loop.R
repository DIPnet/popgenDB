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


# Make an empty list to save gdm output for each species
esu_loci <- unique(ipdb$Genus_species_locus)
all.pops.table<-sapply(esu_loci, function(x) NULL)


for(gsl in esu_loci){ #gsl<-"Linckia_laevigata_CO1" 
  
  cat("\n","\n","\n","Now starting", gsl, "\n")
  
  
  if(any(is.na(diffstats[[gsl]]))){"NAs in FST table, No gdm calculated"; next}

  if(diffstats[[gsl]]=="Fewer than 3 sampled populations after filtering. No stats calculated"){all.pops.table[[gsl]]<-"Fewer than 5 sampled populations after filtering."; next}
  
  
  ######################################################################
  #Subsample for each species of interest, and filter based on Phi_ST table.
  # subsample IPDB, give it a name and order it alphabetically
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
  gslFSTm[which(gslFSTm<0)]<-0
  
  #zap weird slashes in the names
  rownames(gslFSTm)<-gsub("\"","",rownames(gslFSTm))
  colnames(gslFSTm)<-rownames(gslFSTm)
  
  #zap the same na populations from the list of non existent pops from VeronDivis
  if(any(rownames(gslFSTm) %in% nonVeronpops)){
    gslFSTm<-gslFSTm[-(which(rownames(gslFSTm) %in% nonVeronpops)),-(which(colnames(gslFSTm) %in% nonVeronpops))]
  }
  if(length(rownames(gslFSTm))<5){all.pops.table[[gsl]]<-"Fewer than 5 sampled populations";next}
  
  #and filter sp based on the localities that have Fst values
  sp<-sp[sp$sample %in% rownames(gslFSTm),]
  
  #and vice versa
  
  gslFSTm<- gslFSTm[which(rownames(gslFSTm) %in% unique(sp$sample)),which(rownames(gslFSTm) %in% unique(sp$sample))]
  
  #
  if(length(rownames(gslFSTm))<5){all.pops.table[[gsl]]<-"Fewer than 5 sampled populations";next}
  
  #create a locations data frame that has all the localities plus lats and longs and their Veron region.
  locs<-as.data.frame(unique(sp$sample))
  names(locs)<-"sample"
  #locs$Long<-sp$decimalLongitude[which(locs %in% sp$sample)]
  #can't do a unique on sample, lats and longs because some samples have non-unique lats and longs! So I do a join and take the first match.
  locs<-join(locs,sp[c("sample","decimalLongitude","decimalLatitude"
                       ,"VeronDivis")], by="sample", match="first")
  
  
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
  
  #run all the models
  gdm1.all<-gdm(gdm.format)
  gdm2.NO.Geog<-gdm(gdm.format[-grep("matrix_1",names(gdm.format)),])
  gdm3.NO.Sunda<-gdm(gdm.format[-grep("matrix_2",names(gdm.format)),])
  gdm4.NO.EastCT<-gdm(gdm.format[-grep("matrix_3",names(gdm.format)),])
  gdm5.ONLY.Geog<-gdm(gdm.format[-grep("matrix_[23]",names(gdm.format)),])
  gdm6.ONLY.Sunda<-gdm(gdm.format[-grep("matrix_[13]",names(gdm.format)),])
  gdm7.ONLY.EastCT<-gdm(gdm.format[-grep("matrix_[12]",names(gdm.format)),])
  
  if(any(is.null(gdm1.all),is.null(gdm2.NO.Geog),is.null(gdm3.NO.Sunda),is.null(gdm4.NO.EastCT),is.null(gdm5.ONLY.Geog),is.null(gdm6.ONLY.Sunda),is.null(gdm7.ONLY.EastCT))) {paste("some null models for",gsl, "skipping"); all.pops.table[[gsl]]<-"Problem with model fitting, no solution obtained. Sum of spline coefficients = 0. Deviance explained = NULL."; next}
  
  #summary(gdm1.all)
  
  #print the plots for gdm_all to a file
  #pdf(file=paste(gsl,"gdm_plot.pdf",sep=""))
  #plot(gdm1.all)
  #dev.off()
  
  ###############################################################################################
  #Stats to save: relative importance (=sum of coefficients for each predictor) and proportion of deviance explained by model
  #See tables and Figures in Fitzpatrick MC, Sanders NJ, Normand S, Svenning JC, Ferrier S, Gove AD, Dunn RR, 2013. Environmental and historical imprints on beta diversity: insights from variation in rates of species turnover along gradients. Proc. Biol. Sci. 280: 20131201–20131201
  ####
  
  
  # GDM 1 - Full Model
  #importance stats
  numb.1.coeff<-gdm1.all$splines[1]
  numb.2.coeff<-gdm1.all$splines[2]
  numb.3.coeff<-gdm1.all$splines[3]
  
  impt.Geog.gdm1.all<-sum(gdm1.all$coefficients[1:numb.1.coeff])
  impt.Sunda.gdm1.all<-sum(gdm1.all$coefficients[(numb.1.coeff+1) : (numb.1.coeff+numb.2.coeff)])
  impt.EastCT.gdm1.all<-sum(gdm1.all$coefficients[(numb.1.coeff+numb.2.coeff+1):(numb.1.coeff+numb.2.coeff+numb.3.coeff)])
  
  #deviance stats
  gdm1.all.gdm.deviance<- gdm1.all$gdmdeviance  
  gdm1.all.explained.deviance<- gdm1.all$explained 
  gdm1.all.proportion.deviance.explained<- gdm1.all.explained.deviance/ gdm1.all.gdm.deviance   #see Fig 2 in Fitzpatrick et al. 2013 
  
  
  # GDM 2 - No geog distance
  #importance stats
  numb.1.coeff<-gdm2.NO.Geog$splines[1]
  numb.2.coeff<-gdm2.NO.Geog$splines[2]
  
  impt.Sunda.gdm2.NO.Geog<-sum(gdm2.NO.Geog$coefficients[1 : (numb.1.coeff)])
  impt.EastCT.gdm2.NO.Geog<-sum(gdm2.NO.Geog$coefficients[(numb.1.coeff+1):(numb.1.coeff+numb.2.coeff)])
  
  #deviance stats
  gdm2.NO.Geog.gdm.deviance<- gdm2.NO.Geog$gdmdeviance  
  gdm2.NO.Geog.explained.deviance<- gdm2.NO.Geog$explained 
  gdm2.NO.Geog.proportion.deviance.explained<- gdm2.NO.Geog.explained.deviance/ gdm2.NO.Geog.gdm.deviance   #see Fig 2 in Fitzpatrick et al. 2013 
  
  
  
  # GDM 3 - No Sunda barrier
  #importance stats
  numb.1.coeff<-gdm3.NO.Sunda$splines[1]
  numb.2.coeff<-gdm3.NO.Sunda$splines[2]
  
  impt.Geog.gdm3.NO.Sunda<-sum(gdm3.NO.Sunda$coefficients[1 : (numb.1.coeff)])
  impt.EastCT.gdm3.NO.Sunda<-sum(gdm3.NO.Sunda$coefficients[(numb.1.coeff+1):(numb.1.coeff+numb.2.coeff)])
  
  #deviance stats
  gdm3.NO.Sunda.gdm.deviance<- gdm3.NO.Sunda$gdmdeviance  
  gdm3.NO.Sunda.explained.deviance<- gdm3.NO.Sunda$explained 
  gdm3.NO.Sunda.proportion.deviance.explained<- gdm3.NO.Sunda.explained.deviance/ gdm3.NO.Sunda.gdm.deviance   #see Fig 2 in Fitzpatrick et al. 2013 
  
  
  # GDM 4 - No EastCT barrier
  #importance stats
  numb.1.coeff<-gdm4.NO.EastCT$splines[1]
  numb.2.coeff<-gdm4.NO.EastCT$splines[2]
  
  impt.Geog.gdm4.NO.EastCT<-sum(gdm4.NO.EastCT$coefficients[1 : (numb.1.coeff)])
  impt.Sunda.gdm4.NO.EastCT<-sum(gdm4.NO.EastCT$coefficients[(numb.1.coeff+1):(numb.1.coeff+numb.2.coeff)])
  
  #deviance stats
  gdm4.NO.EastCT.gdm.deviance<- gdm4.NO.EastCT$gdmdeviance  
  gdm4.NO.EastCT.explained.deviance<- gdm4.NO.EastCT$explained 
  gdm4.NO.EastCT.proportion.deviance.explained<- gdm4.NO.EastCT.explained.deviance/ gdm4.NO.EastCT.gdm.deviance   #see Fig 2 in Fitzpatrick et al. 2013 
  
  # GDM 5 - ONLY Geog
  #importance stats
  numb.1.coeff<-gdm5.ONLY.Geog$splines[1]
  
  impt.Geog.gdm5.ONLY.Geog<-sum(gdm5.ONLY.Geog$coefficients[1 : (numb.1.coeff)])
  
  #deviance stats
  gdm5.ONLY.Geog.gdm.deviance<- gdm5.ONLY.Geog$gdmdeviance  
  gdm5.ONLY.Geog.explained.deviance<- gdm5.ONLY.Geog$explained 
  gdm5.ONLY.Geog.proportion.deviance.explained<- gdm5.ONLY.Geog.explained.deviance/ gdm5.ONLY.Geog.gdm.deviance   #see Fig 2 in Fitzpatrick et al. 2013 
  
  
  # GDM 6 - ONLY Sunda
  #importance stats
  numb.1.coeff<-gdm6.ONLY.Sunda$splines[1]
  
  impt.Sunda.gdm6.ONLY.Sunda<-sum(gdm6.ONLY.Sunda$coefficients[1 : (numb.1.coeff)])
  
  #deviance stats
  gdm6.ONLY.Sunda.gdm.deviance<- gdm6.ONLY.Sunda$gdmdeviance  
  gdm6.ONLY.Sunda.explained.deviance<- gdm6.ONLY.Sunda$explained 
  gdm6.ONLY.Sunda.proportion.deviance.explained<- gdm6.ONLY.Sunda.explained.deviance/ gdm6.ONLY.Sunda.gdm.deviance   #see Fig 2 in Fitzpatrick et al. 2013 
  
  
  
  # GDM 7 - ONLY EastCT
  #importance stats
  numb.1.coeff<-gdm7.ONLY.EastCT$splines[1]
  
  impt.EastCT.gdm7.ONLY.EastCT<-sum(gdm7.ONLY.EastCT$coefficients[1 : (numb.1.coeff)])
  
  #deviance stats
  gdm7.ONLY.EastCT.gdm.deviance<- gdm7.ONLY.EastCT$gdmdeviance  
  gdm7.ONLY.EastCT.explained.deviance<- gdm7.ONLY.EastCT$explained 
  gdm7.ONLY.EastCT.proportion.deviance.explained<- gdm7.ONLY.EastCT.explained.deviance/ gdm7.ONLY.EastCT.gdm.deviance   #see Fig 2 in Fitzpatrick et al. 2013 
  
  #Cat all the stats together
  catstats<-c(impt.Geog.gdm1.all,impt.Sunda.gdm1.all,impt.EastCT.gdm1.all,gdm1.all.gdm.deviance,gdm1.all.explained.deviance,gdm1.all.proportion.deviance.explained,
              impt.Sunda.gdm2.NO.Geog,impt.EastCT.gdm2.NO.Geog,gdm2.NO.Geog.gdm.deviance,gdm2.NO.Geog.explained.deviance,gdm2.NO.Geog.proportion.deviance.explained,
              impt.Geog.gdm3.NO.Sunda,impt.EastCT.gdm3.NO.Sunda,gdm3.NO.Sunda.gdm.deviance,gdm3.NO.Sunda.explained.deviance,gdm3.NO.Sunda.proportion.deviance.explained,
              impt.Geog.gdm4.NO.EastCT,impt.Sunda.gdm4.NO.EastCT,gdm4.NO.EastCT.gdm.deviance,gdm4.NO.EastCT.explained.deviance,gdm4.NO.EastCT.proportion.deviance.explained,
              impt.Geog.gdm5.ONLY.Geog,gdm5.ONLY.Geog.gdm.deviance,gdm5.ONLY.Geog.explained.deviance,gdm5.ONLY.Geog.proportion.deviance.explained,
              impt.Sunda.gdm6.ONLY.Sunda,gdm6.ONLY.Sunda.gdm.deviance,gdm6.ONLY.Sunda.explained.deviance,gdm6.ONLY.Sunda.proportion.deviance.explained,
              impt.EastCT.gdm7.ONLY.EastCT,gdm7.ONLY.EastCT.gdm.deviance,gdm7.ONLY.EastCT.explained.deviance,gdm7.ONLY.EastCT.proportion.deviance.explained)
  
  #add names
  names(catstats)<-c("impt.Geog.gdm1.all","impt.Sunda.gdm1.all","impt.EastCT.gdm1.all","gdm1.all.gdm.deviance","gdm1.all.explained.deviance","gdm1.all.proportion.deviance.explained","impt.Sunda.gdm2.NO.Geog","impt.EastCT.gdm2.NO.Geog","gdm2.NO.Geog.gdm.deviance","gdm2.NO.Geog.explained.deviance","gdm2.NO.Geog.proportion.deviance.explained","impt.Geog.gdm3.NO.Sunda","impt.EastCT.gdm3.NO.Sunda","gdm3.NO.Sunda.gdm.deviance","gdm3.NO.Sunda.explained.deviance","gdm3.NO.Sunda.proportion.deviance.explained","impt.Geog.gdm4.NO.EastCT","impt.Sunda.gdm4.NO.EastCT","gdm4.NO.EastCT.gdm.deviance","gdm4.NO.EastCT.explained.deviance","gdm4.NO.EastCT.proportion.deviance.explained","impt.Geog.gdm5.ONLY.Geog","gdm5.ONLY.Geog.gdm.deviance","gdm5.ONLY.Geog.explained.deviance","gdm5.ONLY.Geog.proportion.deviance.explained","impt.Sunda.gdm6.ONLY.Sunda","gdm6.ONLY.Sunda.gdm.deviance","gdm6.ONLY.Sunda.explained.deviance","gdm6.ONLY.Sunda.proportion.deviance.explained","impt.EastCT.gdm7.ONLY.EastCT","gdm7.ONLY.EastCT.gdm.deviance","gdm7.ONLY.EastCT.explained.deviance","gdm7.ONLY.EastCT.proportion.deviance.explained")
  
  all.pops.table[[gsl]]<-catstats
}

