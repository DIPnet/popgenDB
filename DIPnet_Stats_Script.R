# Eric Crandall March 2015 - These commands should be helpful in running the scripts in DIPnet_Stats_Functions.R
# Relevant Links:
# Interesting discussion about calculating Fst in R: http://www.molecularecologist.com/2012/05/calculating-pair-wise-unbiased-fst-with-r/


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

source("DIPnet_Stats_Functions.R")
source("config.R")

setwd(working_directory)



##READING IN THE DATA. ##
#need to turn off the quoting with quote="" for it to read correctly. 
#Point this to keepsdata_ABGD files if you want to calculate stats by ESUs rather than by species
#Add "_ABGD" to file name to load the ABGD assignments
ipdb<-read.table(ipdb_path,sep="\t",header=T,stringsAsFactors = F,quote="") 

#windows version of reading in a file
#ipdb<-read.table("C:/Users/Chris/Google Drive/!DIPnet_DataQC/Reunite_metadata_and_alignments/2015-03-24/2015-03-24_keepsdata.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
 

#read in geographical regionalizations from Treml
spatial<-read.csv(spatial_path, stringsAsFactors = F, na.strings=c("NA"," ",""))
#spatial<-read.csv("~/Google Drive/!DIPnet_DataQC/Reunite_metadata_and_alignments/2015-03-24/KD2_24Mar15_AllJoin_cln.csv", stringsAsFactors = F, na.strings=c("NA"," ",""))
#windows version
#spatial<-read.csv("C:/Users/Chris/Google Drive/!DIPnet_DataQC/Reunite_metadata_and_alignments/2015-03-24/KD2_24Mar15_AllJoin_cln.csv", stringsAsFactors = F, na.strings=c("NA"," ",""))


#Remove commas from fn100id and fn500id (also may need to fix reunion accent issue manually)
spatial$fn100id<-gsub(",","", spatial$fn100id)
spatial$fn500id<-gsub(",","", spatial$fn500id)
ipdb<-join(ipdb,spatial, by = "IPDB_ID")
#CHECK FOR DUPLICATES
dups<-ipdb[duplicated(ipdb),]
#ipdb<-ipdb[!duplicated(ipdb),]  #remove duplicates if you need to

#the following IDs should be removed if not using the ABGD dataset per Libby Liggins
drops <- c("Acanthurus_sohal_CO1_A8","Acanthurus_sohal_CO1_RS3407","Acanthurus_sohal_CO1_RS3406","Linckia_laevigata_CO1_2130.14","Linckia_laevigata_CO1_1195.01","Linckia_laevigata_CO1_1195.02","Labroides_dimidiatus_CR_LB10","Labroides_dimidiatus_CR_LB01","Labroides_dimidiatus_CR_LB02","Labroides_dimidiatus_CR_LB09","Labroides_dimidiatus_CR_LP03","Labroides_dimidiatus_CR_LP01","Labroides_dimidiatus_CR_LP04","Labroides_dimidiatus_CR_LP02","Katsuwonus_pelamis_CR_3637_01","Katsuwonus_pelamis_CR_3637_1","Periclimenes_soror_CO1_2227.07","Periclimenes_soror_CO1_2227.01","Periclimenes_soror_CO1_2227.02","Periclimenes_soror_CO1_2227.11","Periclimenes_soror_CO1_2227.12","Periclimenes_soror_CO1_2227.1","Periclimenes_soror_CO1_2227.16","Lutjanus_kasmira_CYB_Lka_OK4", "Lutjanus_kasmira_CYB_Lka_OK3")
ipdb<-ipdb[ipdb$IPDB_ID %in% drops == FALSE, ] 

#Also found some duplicate MSIDs in the COTS data removed thusly (apparently Toonen and Worheide contributed the same data)
#cots_dups<-rownames(sp[which(duplicated(sp$materialSampleID,fromLast=F)),1:7])
cots_dups<-c("1536","1539","1541","1542","1543","1546","1548","1550","1551","490","477","478","479","480","481","482","483","484","485")
ipdb<-ipdb[rownames(ipdb) %in% cots_dups==FALSE,]

#CHECK FOR MISSING LOCALITIES
#note there are 1162 individual sequences with no locality! Distributed thusly:
table(ipdb$principalInvestigator[which(is.na(ipdb$locality))])
ipdb$locality[which(is.na(ipdb$locality))]<-"no_name"  #give "no_name" to NA localities

#Remove this later on - lat and long are replaced with rounded values to aid with population reduction
#ipdb$decimalLatitude<-round(ipdb$decimalLatitude, digits=0)      
#ipdb$decimalLongitude<-round(ipdb$decimalLongitude, digits=0)

#subset by species (or whatever else) if necessary
zebfla<-subset(ipdb, Genus_species_locus == "Zebrasoma_flavescens_CYB" )

####RUN IT###
###Diversity Stats Function###
#Computes diversity stats by species and population for a flatfile of mtDNA sequences and metadata (with required fields $Genus_species_locus and $loc_lat_long)
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
# To be added: rarefaction, Fus Fs, Fu and Li's D
divstats<-genetic.diversity.mtDNA.db(ipdb=ipdb, minseqs = 5, minsamps= 3, mintotalseqs= 0, regionalization = "meow_ecoregion")


###Pairwise Genetic Structure Function###
#Computes genetic differentiation statistics by species and population for a flatfile of mtDNA sequences and metadata (with required fields $Genus_species_locus and $loc_lat_long)
# gdist = You must choose one genetic distance to calculate: choose gdist from: "Nei FST","Nei GST", "Hedrick G'ST", "Jost D", "WC Theta", "PhiST", "Chi2", or "NL dA"
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
# nreps = number of resampling permutations for WC Theta, PhiST, Chi2, and NL dA (strataG package). This is working, but you will not currently see p-values in the output, so not much use for now
# num.cores = number of computer cores to devote to computations for WC Theta, PhiST, Chi2, and NL dA (strataG package)
# To be added: option to output square matrices with p-values.
diffstats<-pairwise.structure.mtDNA.db(ipdb=ipdb,gdist = "PhiST",minseqs = 5,minsamps=3, mintotalseqs=0, num.cores=2, regionalization = "meow_ecoregion")


###Hierarchical Genetic Structure Function###
#Computes heigenetic differentiation statistics by species and population for a flatfile of mtDNA sequences and metadata (with required fields $Genus_species_locus and $loc_lat_long)
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
# nperm = number of AMOVA permutations
# model= model of molecular evolution to be passed to dna.dist() = c("raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock")
# To be added: capability to do 3-level F-statistics and option for Fst (distance matrix =1)
hierstats<-hierarchical.structure.mtDNA.db(ipdb = ipdb,level1 = "sample",level2="meow_ecoregion",model="raw",nperm=100)


#Save diversity stats to file(s)
save(divstats,file="../DIPnet_stats_ecogregions_032415.Rdata")
write.stats(divstats,filename=".../DIPnet_stats_ecoregions_032415.csv",structure=F) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
#write.stats(divstats,"DIPnet_stats_ecoregions_032415.csv",structure=F)

#Save differentiation stats to files
save(diffstats,file="./For_NESCent/DIPnet_structure_ecoregions_PhiST_ABGD_032415.Rdata") # for an R object
write.stats(diffstats,filename="./For_NESCent/DIPnet_structure_ecoregions_PhiST_ABGD_032415.csv",structure=T) # for an excel-readable csv. Ignore warnings. structure=T for triangular matrices. Note this function will not overwrite, it will append to existing files





# Loop through all regionalizations and calculate the statistics
for(r in c("sample","meow_ecoregion", "meow_prov_name", "meow_rlm_name","fn100id", "fn500id","EEZ_country")){
  divstats<-genetic.diversity.mtDNA.db(ipdb=ipdb, minseqs = 5, minsamps= 3, mintotalseqs= 0, regionalization = r)
  dir.create(file.path("./For_NESCent",r))
  save(divstats,file=file.path("./For_NESCent",r,paste("DIPnet_stats_032415_",r,".Rdata",sep="")))
  write.stats(divstats,filename=file.path("./For_NESCent",r,paste("DIPnet_stats_032415_",r,".csv",sep="")),structure=F) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
  
}

for(r in c("sample","meow_ecoregion", "meow_prov_name", "meow_rlm_name","fn100id", "fn500id","EEZ_country")){
  for(g in c("WC Theta","PhiST", "Jost D")){
  diffstats<-pairwise.structure.mtDNA.db(ipdb=ipdb, gdist = g, regionalization = r, minseqs = 5, minsamps= 3, mintotalseqs= 0, num.cores = 2)
  #dir.create(file.path("./For_NESCent",r))
  save(diffstats,file=file.path("./For_NESCent",r,paste("DIPnet_structure_032415_ABGD_",g,"_",r,".Rdata",sep="")))
  write.stats(diffstats,filename=file.path("./For_NESCent",r,paste("DIPnet_structure_032415_",g,"_",r,".csv",sep="")),structure=T) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
  }
}

