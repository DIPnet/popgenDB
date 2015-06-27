#Modified the stats functions to push out FABOX formatted FASTA files for conversion to Arlequin format. These files will be used to test our statistics against results from Arlequin
#Eric Crandall June 3, 2015

esu_loci<-c("Acanthurus_nigrofuscus_CO1","Zebrasoma_flavescens_CYB","Sphyrna_lewini_CR","Halichoeres_hortulanus_CR","Diadema_paucispinum_A68","Stenopus_hispidus_CO1","Stenella_longirostris_CYB","Diadema_savignyi_A68")
regionalization<-"sample"
minseqs = 5
minsamps = 3
mintotalseqs = 0

for(gsl in esu_loci){ #gsl<-"Zebrasoma_flavescens_CYB" 
  
  cat("\n","\n","\n","Now starting", gsl, "\n")
  
  #POPULATION SAMPLES
  
  #modify/replace the following once we have the script that assigns populations
  #subset the genus-species-locus and order it alphabetically by locality name
  sp<-subset(ipdb, Genus_species_locus == gsl )
  sp$sample<-paste(sp$locality,round(sp$decimalLatitude, digits=0),round(sp$decimalLongitude, digits=0),sep="_")  #sets up a variable that matches assignsamp function outcome
  sp<-sp[order(sp$sample),]
  
  #this code is not currently necessary given the sp$sample field above may come in handy when we start using ecoregions etc.
  #create a set of unique samples, as denoted by lat+long+pi, and sort it alphabetically. 
  #This is because the genind object will return things in alphaetical order, so we want to be on the same page with it
  #samps<-sort(unique(paste(sp$locality,sp$decimalLatitude,sp$decimalLongitude,sep="_")))
  ##assign each individual to its population sample##
  #sp$sample<-sapply(1:length(sp[,1]), assignsamp) 
  
  #FILTER
  cat("filtering out population samples with n <", minseqs,"and species with fewer than", minsamps,"total populations \n")
  #filter out populations with low numbers (n<5 for now). Skip the species if this procedure deletes all samples or if there are less than 3 populations total.
  sampN<-table(sp[[regionalization]])
  lowsamps<-names(sampN[sampN < minseqs])
  if(length(lowsamps)>0){sp<-sp[-which(sp[[regionalization]] %in% lowsamps),]}
  if(length(sampN) - length(lowsamps) < minsamps){
                                                  paste("Fewer than", minsamps, "sampled populations after filtering. No stats calculated")
                                                  next}
  if(length(sp[,1])<mintotalseqs){
                                  cat("fewer than",mintotalseqs,"samples left after filtering. No stats calculated")
                                  next}
  cat("Removed the following population samples:", lowsamps, "\n")
  
  #FORMAT CONVERSION
  cat("converting data to various R DNA formats", "\n")
  #convert to seqinr alignment (need to specify that we are using the as.alignment() function from seqinr rather than ape)
  spseqs<-seqinr::as.alignment(nb=length(sp[,1]),nam=sp$materialSampleID, seq=sp$sequence)
  
  faboxnames<-paste(spseqs$nam,"F",sp[[regionalization]],sep="|")
  #write to fasta format for input to FaBox
  seqinr::write.fasta(as.list(spseqs$seq),faboxnames,file.out=paste(gsl,regionalization,"fabox.fasta",sep="_"))
  
}



load("/Users/eric/Google Drive/!DIPnet_Gait_Lig_Bird/DIPnet_WG4/statistics/By_Species/Diversity_statistics/sample/DIPnet_stats_061015sample.Rdata")

for(gsl in esu_loci){#gsl<-"Zebrasoma_flavescens_CYB" 
  
    level<-"sample"
    cat(gsl,level,"\n \n")
    arlstats<-read.csv(paste("/Users/eric/Google Drive/!DIPnet_Gait_Lig_Bird/DIPnet_WG4/statistics/Arlequin_crosscheck/",gsl,"_",level,"_stats.csv", sep=""),fill = T)
    
    dipnetstats<-divstats[[gsl]]
    
    tol<-4
    
    for(stat in c("sampleN","UniqHapNum","HaploDiv","NucDivSite","ThetaS","TajD")) { #stat<-"ThetaS"
      cat("Now checking", stat,"\n")
      diffs<-setdiff(round(arlstats[[stat]],tol), round(dipnetstats[[stat]],tol))
      if(length(diffs > 0)){ cat("different values for", stat,
                                  ": Arlequin =", round(arlstats[[stat]],tol)[which(round(arlstats[[stat]],tol) %in% diffs)], 
                                    "and DIPnet =", round(dipnetstats[[stat]],tol)[which(round(arlstats[[stat]],tol) %in% diffs)],
                                  "for",rownames(dipnetstats)[which(round(arlstats[[stat]],tol) %in% diffs)],"\n","\n")
                              
                              }
    }
}
  