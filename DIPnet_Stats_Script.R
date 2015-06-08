# Eric Crandall March 2015 - These commands should be helpful in running the scripts in DIPnet_Stats_Functions.R
# Relevant Links:
# Interesting discussion about calculating Fst in R: http://www.molecularecologist.com/2012/05/calculating-pair-wise-unbiased-fst-with-r/

setwd("~/github/popgenDB/output")

#Load libraries and source the functions
library(seqinr)  #might have to download .tgz file directly from CRAN site and install locally, not directly from CRAN repository
library(ape)  
library(pegas)
library(mmod)
library(adegenet)
library(plyr)
library(strataG)
##Do the first time to install iNEXT##
# install.packages('devtools')
# library(devtools)
# install_github('JohnsonHsieh/iNEXT')
library(iNEXT)

#set the working directory and source stuff

source("../DIPnet_Stats_Functions.R")
source("../config.R")

setwd(working_directory)



##READING IN THE DATA. ##
#need to turn off the quoting with quote="" for it to read correctly. 

ipdb<-read.table(ipdb_path,sep="\t",header=T,stringsAsFactors = F,quote="", na.strings=c("NA"," ","")) 


#read in geographical regionalizations from Treml
spatial<-read.table(spatial_path, header=T, sep="\t",stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")

#read in ABGD groups
abgd<-read.table(abgd_path, header=T, sep="\t", stringsAsFactors = F)

#join spatial
ipdb<-join(ipdb,spatial, by = "IPDB_ID",type = "left")

#join ABGD
ipdb<-join(ipdb,abgd[,c(1,3)], by = "IPDB_ID",type = "left")

#CHECK FOR DUPLICATES
#dups<-ipdb[duplicated(ipdb),]
#ipdb<-ipdb[!duplicated(ipdb),]  #remove duplicates if you need to

#the following 13 IDs should be removed if not using the ABGD dataset per Libby Liggins. They are highly divergent individuals, likely to be due to hybridization, cryptic species etc.
drops <- c("Acanthurus_sohal_CO1_A8","Acanthurus_sohal_CO1_RS3407","Acanthurus_sohal_CO1_RS3406","Linckia_laevigata_CO1_2130.14","Linckia_laevigata_CO1_1195.01","Linckia_laevigata_CO1_1195.02","Labroides_dimidiatus_CR_LB10","Labroides_dimidiatus_CR_LB01","Labroides_dimidiatus_CR_LB02","Labroides_dimidiatus_CR_LB09","Labroides_dimidiatus_CR_LP03","Labroides_dimidiatus_CR_LP01","Labroides_dimidiatus_CR_LP04","Labroides_dimidiatus_CR_LP02","Katsuwonus_pelamis_CR_3637_01","Katsuwonus_pelamis_CR_3637_1","Periclimenes_soror_CO1_2227.07","Periclimenes_soror_CO1_2227.01","Periclimenes_soror_CO1_2227.02","Periclimenes_soror_CO1_2227.11","Periclimenes_soror_CO1_2227.12","Periclimenes_soror_CO1_2227.1","Periclimenes_soror_CO1_2227.16","Lutjanus_kasmira_CYB_Lka_OK4", "Lutjanus_kasmira_CYB_Lka_OK3")
ipdb<-ipdb[ipdb$IPDB_ID %in% drops == FALSE, ] 

#CHECK FOR MISSING LOCALITIES
#table(ipdb$principalInvestigator[which(is.na(ipdb$locality))])
#ipdb$locality[which(is.na(ipdb$locality))]<-"no_name"  #give "no_name" to NA localities

#subset by species (or whatever else) if necessary
#zebfla<-subset(ipdb, Genus_species_locus == "Zebrasoma_flavescens_CYB" )

####Examples###
###Diversity Stats Function###
#Computes diversity stats by species and population for a flatfile of mtDNA sequences and metadata (with required fields $Genus_species_locus and $loc_lat_long)
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
# To be added: rarefaction, Fus Fs, Fu and Li's D
divstats<-genetic.diversity.mtDNA.db(ipdb=ipdb, minseqs = 5, minsamps= 3, mintotalseqs= 0, regionalization = "sample")


###Pairwise Genetic Structure Function###
#Computes genetic differentiation statistics by species and population for a flatfile of mtDNA sequences and metadata (with required fields $Genus_species_locus and $loc_lat_long)
# gdist = You must choose one genetic distance to calculate: choose gdist from: "Nei FST","Nei GST", "Hedrick G'ST", "Jost D", "WC Theta", "PhiST", "Chi2", or "NL dA"
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
# nreps = number of resampling permutations for WC Theta, PhiST, Chi2, and NL dA (strataG package). This is working, but you will not currently see p-values in the output, so not much use for now
# num.cores = number of computer cores to devote to computations for WC Theta, PhiST, Chi2, and NL dA (strataG package)
# To be added: option to output square matrices with p-values.
diffstats<-pairwise.structure.mtDNA.db(ipdb=ipdb,gdist = "PhiST",minseqs = 5,minsamps=3, mintotalseqs=0, num.cores=2, regionalization = "ECOREGION")


###Hierarchical Genetic Structure Function###
#Computes heigenetic differentiation statistics by species and population for a flatfile of mtDNA sequences and metadata (with required fields $Genus_species_locus and $loc_lat_long)
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
# nperm = number of AMOVA permutations
# model= model of molecular evolution to be passed to dna.dist() = c("raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock")
# To be added: capability to do 3-level F-statistics and option for Fst (distance matrix =1)
hierstats<-hierarchical.structure.mtDNA.db(ipdb = ipdb,level1 = "sample",level2="ECOREGION",model="raw",nperm=100)


#Save diversity stats to file(s)
save(divstats,file="DIPnet_stats_samples_060815.Rdata")
write.stats(divstats,filename="DIPnet_stats_samples_060815.csv",structure=F) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
#write.stats(divstats,"DIPnet_stats_ecoregions_032415.csv",structure=F)

#Save differentiation stats to files
save(diffstats,file="DIPnet_structure_ecoregions_PhiST_ABGD_032415.Rdata") # for an R object
write.stats(diffstats,filename="DIPnet_structure_ecoregions_PhiST_ABGD_032415.csv",structure=T) # for an excel-readable csv. Ignore warnings. structure=T for triangular matrices. Note this function will not overwrite, it will append to existing files





# Loop through all regionalizations and calculate the statistics
for(r in c("sample","fn100id", "fn500id", "ECOREGION", "PROVINCE", "REALM", "EEZ")){
  divstats<-genetic.diversity.mtDNA.db(ipdb=ipdb, minseqs = 5, minsamps= 3, mintotalseqs= 0, regionalization = r, ABGD=T)
  dir.create(file.path("./",r))
  save(divstats,file=file.path("./",r,paste("DIPnet_stats_060315_ABGD_no_coverage_correction",r,".Rdata",sep="")))
  write.stats(divstats,filename=file.path("./",r,paste("DIPnet_stats_060315_ABGD_no_coverage_correction",r,".csv",sep="")),structure=F) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
  
}

for(r in c("sample","fn100id", "fn500id", "ECOREGION", "PROVINCE", "REALM", "EEZ")){
  for(g in c("WC Theta","PhiST", "Jost D")){
  diffstats<-pairwise.structure.mtDNA.db(ipdb=ipdb, gdist = g, regionalization = r, minseqs = 5, minsamps= 3, mintotalseqs= 0, num.cores = 2)
  #dir.create(file.path("./",r))
  save(diffstats,file=file.path("./",r,paste("DIPnet_structure_060315_",g,"_",r,".Rdata",sep="")))
  write.stats(diffstats,filename=file.path("./",r,paste("DIPnet_structure_060315_",g,"_",r,".csv",sep="")),structure=T) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
  }
}

