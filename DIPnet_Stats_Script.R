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
# source("DIPnet_Stats_Functions.R")
# source("config.R")

setwd(working_directory)



##READING IN THE DATA. ##
#need to turn off the quoting with quote="" for it to read correctly. 

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
divstats<-genetic.diversity.mtDNA.db(ipdb=ipdb, basic_diversity = T, sequence_diversity = T, coverage_calc = F, coverage_correction = F, minseqs = 6, minsamps = 3, mintotalseqs = 0, ABGD=F,regionalization = "sample", keep_all_gsls=F, mincoverage = 0.4, hill.number = 0)

###Pairwise Genetic Structure Function###
#Computes genetic differentiation statistics by species and population for a flatfile of mtDNA sequences and metadata (with required fields $Genus_species_locus and $loc_lat_long)
# gdist = You must choose one genetic distance to calculate: choose gdist from:"Nei GST", "Hedrick G'ST", "Jost D", "WC Theta", "PhiST", "Chi2", or "NL dA"
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
# nreps = number of resampling permutations for WC Theta, PhiST, Chi2, and NL dA (strataG package). This is working, but you will not currently see p-values in the output, so not much use for now
# num.cores = number of computer cores to devote to computations for WC Theta, PhiST, Chi2, and NL dA (strataG package)
# To be added: option to output square matrices with p-values.
diffstats<-pairwise.structure.mtDNA.db(ipdb=ipdb, gdist = "PhiST", minseqs = 5, minsamps = 3, mintotalseqs = 0, nrep = 0, num.cores = 1, ABGD = F, regionalization = "sample")


###Hierarchical Genetic Structure Function###
#Computes heigenetic differentiation statistics by species and population for a flatfile of mtDNA sequences and metadata (with required fields $Genus_species_locus and $loc_lat_long)
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
# nperm = number of AMOVA permutations
# model= model of molecular evolution to be passed to dna.dist() = c("raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock")
# model defaults to "N" which is the raw count of differences ("raw" is the proportion- same thing). If you use model = "none" you will get all distances between haplotypes = 1, which is the same as "regular" FST
# levels can be one of c("sample","fn100id", "fn500id", "ECOREGION", "PROVINCE", "REALM", "EEZ") or new regionalizations as they are added.
hierstats<-hierarchical.structure.mtDNA.db(ipdb = ipdb,level1 = "sample",level2="ECOREGION",model="raw",nperm=1)


#Save diversity stats to file(s)
save(divstats,file="DIPnet_stats_samples_060815.Rdata")
write.stats(divstats,filename="DIPnet_stats_samples_060815.csv",structure=F) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
#write.stats(divstats,"DIPnet_stats_ecoregions_032415.csv",structure=F)

#Save differentiation stats to files
save(diffstats,file="DIPnet_structure_ecoregions_PhiST_070215.Rdata") # for an R object
write.stats(diffstats,filename="DIPnet_structure_ecoregions_PhiST_070215.csv",structure=T) # for an excel-readable csv. Ignore warnings. structure=T for triangular matrices. Note this function will not overwrite, it will append to existing files





# Loop through all regionalizations and calculate the statistics
for(r in c("sample","ECOREGION", "PROVINCE", "REALM", "EEZ", "fn100id", "fn500id")){
  divstats<-genetic.diversity.mtDNA.db(ipdb=ipdb, basic_diversity = T, sequence_diversity = T, coverage_calc = T, coverage_correction = T, minseqs = 6, minsamps = 3, mintotalseqs = 0, ABGD=F,regionalization = r, keep_all_gsls=F, mincoverage = 0.4, hill.number = 0)
  dir.create(file.path("./",r))
  save(divstats,file=file.path("./",r,paste("DIPnet_stats_Hill0_072115_",r,".Rdata",sep="")))
  write.stats(divstats,filename=file.path("./",r,paste("DIPnet_stats_Hill0_072115_",r,".csv",sep="")),structure=F) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
  
}

for(r in c("sample","ECOREGION", "PROVINCE", "REALM", "EEZ", "fn100id", "fn500id")){
  for(g in c("WC Theta","PhiST", "Jost D")){
  diffstats<-pairwise.structure.mtDNA.db(ipdb=ipdb, gdist = g, regionalization = r, minseqs = 6, minsamps= 3, mintotalseqs= 0, num.cores = 2)
  #dir.create(file.path("./",r))
  save(diffstats,file=file.path("./",r,paste("DIPnet_structure_060315_",g,"_",r,".Rdata",sep="")))
  write.stats(diffstats,filename=file.path("./",r,paste("DIPnet_structure_060315_",g,"_",r,".csv",sep="")),structure=T) # for an excel-readable csv. Ignore warnings. Note this function will not overwrite, it will append to existing files
  }
}


#remove anything not included in the ecoregions scheme (some dolphins, some COTS from Kingman and Madagascar(?), some A. nigros from Kiribati, som C. auriga from Fakareva, hammerheads from Western Australia, and West Africa, and some dolphins from the middle of the eastern tropical pacific

ipdb_ecoregions<-ipdb[-which(is.na(ipdb$ECOREGION)),]

#remove anything that doesn't occur in the 3 Indo-Pacific realms
ipdb_ip<-ipdb_ecoregions[which(ipdb_ecoregions$REALM %in% c("Central Indo-Pacific","Western Indo-Pacific","Eastern Indo-Pacific")),]


# Loop through hypotheses, calculating AMOVA
amova_list<-list()
for(h in c("Lat_Zone","VeronDivis","Kulbicki_b","Kulbicki_r","Bowen","Keith","ECOREGION", "PROVINCE","REALM")){
  ipdb_trim<-ipdb[-which(is.na(ipdb[[h]])),]

  hierstats<-hierarchical.structure.mtDNA.db(ipdb = ipdb_trim,level1 = "sample",level2=h,model="raw",nperm=1)
  amova_list[[h]]<-hierstats
}
  

#Summarize AMOVA results
stat.list<-list()
for(h in c("Lat_Zone","VeronDivis","Kulbicki_b","Kulbicki_r","Bowen","Keith","ECOREGION", "PROVINCE","REALM")){

  #Create an empty table the length of all the gsls in the dataset
    stat.table<-data.frame(row.names=names(amova_list[[h]]),level1_k=integer(length(names(amova_list[[h]]))),level2_k=integer(length(names(amova_list[[h]]))),FCT=numeric(length(names(amova_list[[h]]))),FSC=numeric(length(names(amova_list[[h]]))),FST=numeric(length(names(amova_list[[h]]))))
  
  for(gsl in names(amova_list[[1]])){
    # loop through all the gsls pulling out various stats from each one into a data table
    
    amova<-amova_list[[h]][[gsl]]
  
    #some tests to skip over non-existant of single level AMOVAs that don't have FCT  
    if(is.null(amova)){next} # skip to the next gsl if this one has null results
    if(grepl(pattern="fewer",x = amova[1],ignore.case = T)){next} # skip to the next gsl if this one has no results
    if(length(amova$level2_names)<=1){next} # skip to the next gsl if level 2 of the AMOVA has 0 or 1 levels
 
    #pull out the stats 
    level1_k<-length(amova$level1_names) 
    level2_k<-length(amova$level2_names)
    FCT<-amova$FCT
    FSC<-amova$FSC
    FST<-amova$FST
    
    #tie it all up - place it in the appropriate line of the data frame
    stats<-c(level1_k,level2_k,FCT,FSC,FST)
    stat.table[gsl,]<-stats
    
    #all skipped gsls replaced with NA
    stat.table[which(stat.table$level1_k==0),]<-NA

  }
    
  
stat.list[[h]]<-stat.table
}
melted<-melt(stat.list)
FCT<-melted[which(melted$variable=="FCT"),]

plot<-ggplot(data=FCT, aes(x=L1,y=value))
plot<-plot + geom_boxplot() + ylim(c(-0.1,1))
  
  



