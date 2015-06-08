#Copied from DIPnet_Stats_Script.R
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

source("coverage_function_40percentSCtrigger.R")
source("config.R")

setwd(working_directory)



##READING IN THE DATA. ##
#need to turn off the quoting with quote="" for it to read correctly. 
#Point this to keepsdata_ABGD files if you want to calculate stats by ESUs rather than by species
#Add "_ABGD" to file name to load the ABGD assignments
ipdb<-read.table(ipdb_path,sep="\t",header=T,stringsAsFactors = F,quote="") 

#read in ABGD groups
abgd<-read.table(abgd_path, header=T, sep="\t", stringsAsFactors = F)

#read in geographical regionalizations from Treml
spatial<-read.table(spatial_path, header=T, sep="\t",stringsAsFactors = F, na.strings=c("NA"," ",""), quote="")

#join spatial
ipdb<-join(ipdb,spatial, by = "IPDB_ID",type = "left")

#join ABGD
ipdb<-join(ipdb,abgd[,c(1,3)], by = "IPDB_ID",type = "left")

#the following 13 IDs should be removed if not using the ABGD dataset per Libby Liggins. They are highly divergent individuals, likely to be due to hybridization, cryptic species etc.
drops <- c("Acanthurus_sohal_CO1_A8","Acanthurus_sohal_CO1_RS3407","Acanthurus_sohal_CO1_RS3406","Linckia_laevigata_CO1_2130.14","Linckia_laevigata_CO1_1195.01","Linckia_laevigata_CO1_1195.02","Labroides_dimidiatus_CR_LB10","Labroides_dimidiatus_CR_LB01","Labroides_dimidiatus_CR_LB02","Labroides_dimidiatus_CR_LB09","Labroides_dimidiatus_CR_LP03","Labroides_dimidiatus_CR_LP01","Labroides_dimidiatus_CR_LP04","Labroides_dimidiatus_CR_LP02","Katsuwonus_pelamis_CR_3637_01","Katsuwonus_pelamis_CR_3637_1","Periclimenes_soror_CO1_2227.07","Periclimenes_soror_CO1_2227.01","Periclimenes_soror_CO1_2227.02","Periclimenes_soror_CO1_2227.11","Periclimenes_soror_CO1_2227.12","Periclimenes_soror_CO1_2227.1","Periclimenes_soror_CO1_2227.16","Lutjanus_kasmira_CYB_Lka_OK4", "Lutjanus_kasmira_CYB_Lka_OK3")
ipdb<-ipdb[ipdb$IPDB_ID %in% drops == FALSE, ] 

# #Remove commas from fn100id and fn500id (also may need to fix reunion accent issue manually)
# spatial$fn100id<-gsub(",","", spatial$fn100id)
# spatial$fn500id<-gsub(",","", spatial$fn500id)
# ipdb<-join(ipdb,spatial, by = "IPDB_ID")
# #CHECK FOR DUPLICATES
# dups<-ipdb[duplicated(ipdb),]
# #ipdb<-ipdb[!duplicated(ipdb),]  #remove duplicates if you need to
# 
# #the following IDs should be removed if not using the ABGD dataset per Libby Liggins
# drops <- c("Acanthurus_sohal_CO1_A8","Acanthurus_sohal_CO1_RS3407","Acanthurus_sohal_CO1_RS3406","Linckia_laevigata_CO1_2130.14","Linckia_laevigata_CO1_1195.01","Linckia_laevigata_CO1_1195.02","Labroides_dimidiatus_CR_LB10","Labroides_dimidiatus_CR_LB01","Labroides_dimidiatus_CR_LB02","Labroides_dimidiatus_CR_LB09","Labroides_dimidiatus_CR_LP03","Labroides_dimidiatus_CR_LP01","Labroides_dimidiatus_CR_LP04","Labroides_dimidiatus_CR_LP02","Katsuwonus_pelamis_CR_3637_01","Katsuwonus_pelamis_CR_3637_1","Periclimenes_soror_CO1_2227.07","Periclimenes_soror_CO1_2227.01","Periclimenes_soror_CO1_2227.02","Periclimenes_soror_CO1_2227.11","Periclimenes_soror_CO1_2227.12","Periclimenes_soror_CO1_2227.1","Periclimenes_soror_CO1_2227.16","Lutjanus_kasmira_CYB_Lka_OK4", "Lutjanus_kasmira_CYB_Lka_OK3")
# ipdb<-ipdb[ipdb$IPDB_ID %in% drops == FALSE, ] 
# 
# #Also found some duplicate MSIDs in the COTS data removed thusly (apparently Toonen and Worheide contributed the same data)
# #cots_dups<-rownames(sp[which(duplicated(sp$materialSampleID,fromLast=F)),1:7])
# cots_dups<-c("1536","1539","1541","1542","1543","1546","1548","1550","1551","490","477","478","479","480","481","482","483","484","485")
# ipdb<-ipdb[rownames(ipdb) %in% cots_dups==FALSE,]
# 
# #CHECK FOR MISSING LOCALITIES
# #note there are 1162 individual sequences with no locality! Distributed thusly:
# table(ipdb$principalInvestigator[which(is.na(ipdb$locality))])
# ipdb$locality[which(is.na(ipdb$locality))]<-"no_name"  #give "no_name" to NA localities

####RUN IT###
#source("coverage_function_UniqueHaplotypesTrigger.R")  #any gsl with pops of only unique haplotypes gets resampled by tranvservsion defined haplotypes
source("coverage_function_40percentSCtrigger.R")  #any gsl with pops with Standardized Coverage < 0.4 gets resampled by tranvservsion defined haplotypes


## minsamps = 6 seems reasonable based on coverage sensitivity analyses below ##
coverage<-coverage.db(ipdb=ipdb, minseqs = 6, minsamps= 3, mintotalseqs= 0, regionalization = "ECOREGION")


coverage_stats<-rbind(coverage_stats, coverage)







## CHECKING EFFECTS OF COVERAGE SENSITIVITY ##
###Coverage Sensitivity###
#Computes coverage and sample size depending on minsamp cut off
# minseqs = minimum sequences per sampled population, 
# minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
# mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)


coverage_stats<-data.frame()
min_minseqs<-5
max_minseqs<-15
for (s in min_minseqs:max_minseqs) {
coverage<-coverage.db(ipdb=ipdb, minseqs = s, minsamps= 3, mintotalseqs= 0, regionalization = "meow_ecoregion")
coverage_stats<-rbind(coverage_stats, coverage)
}

##Summarizing results
samplesize_sum<-data.frame()
# min_minseqs<-5
# max_minseqs<-15
#par(mfrow=c(1,1))

def.par <- par(no.readonly = TRUE) # save default, for resetting...
#pdf("CoverageSamplingEffort/ActualCoverage_minseqs_standardhaplos.pdf")
#pdf("CoverageSamplingEffort/ActualCoverage_minseqs_TVhaplos_UniqHapTrigger.pdf")
pdf("CoverageSamplingEffort/ActualCoverage_minseqs_TVhaplos_40percentSCtrigger.pdf")
par(mfrow=c(4,4))

for (s in min_minseqs:max_minseqs) {
  coverage_samp<-subset(coverage_stats, minseqs == s)
  samplesize_sum[s-min_minseqs+1,"minseqs"]<-s
  samplesize_sum[s-min_minseqs+1,"retained_gsl"]<-length(unique(coverage_samp$gsl))
  samplesize_sum[s-min_minseqs+1,"retained_pops"]<-nrow(coverage_samp)
  samplesize_sum[s-min_minseqs+1,"mean_actual_coverage"]<-mean(coverage_samp$CoverageforActualSampleSize, na.rm=TRUE)
  samplesize_sum[s-min_minseqs+1,"mean_standardized_coverage"]<-mean(coverage_samp$SC, na.rm=TRUE)
  hist(coverage_samp$CoverageforActualSampleSize, main = paste("minseqs = ", s), xlab = "Actual coverage")
  }
dev.off()

def.par <- par(no.readonly = TRUE) # save default, for resetting...
#pdf("CoverageSamplingEffort/StandardizedCoverage_minseqs_standardhaplos.pdf")
#pdf("CoverageSamplingEffort/StandardizedCoverage_minseqs_TVhaplos_UniqHapTrigger.pdf")
pdf("CoverageSamplingEffort/StandardizedCoverage_minseqs_TVhaplos_40percentSCtrigger.pdf")

par(mfrow=c(4,4))
for (s in min_minseqs:max_minseqs) {
  hist(coverage_samp$SC, main = paste("minseqs = ", s), xlab = "Standardized coverage")
}
dev.off()

#Plots of sample size effect
def.par <- par(no.readonly = TRUE) # save default, for resetting...
#pdf("CoverageSamplingEffort/Summary_effect_minseqs_standardhaplos.pdf")
#pdf("CoverageSamplingEffort/Summary_effect_minseqs_TVhaplos_UniqHapTrigger.pdf")
pdf("CoverageSamplingEffort/Summary_effect_minseqs_TVhaplos_40percentSCtrigger.pdf")

par(mfrow=c(2,2))
with(samplesize_sum, plot(minseqs, retained_gsl, main = "Genus Species Locus"))
with(samplesize_sum, plot(minseqs, retained_pops, main = "Total Populations"))
with(samplesize_sum, plot(minseqs, mean_actual_coverage, main = "Mean Actual Coverage"))
with(samplesize_sum, plot(minseqs, mean_standardized_coverage, main = "Mean Standardized Coverage"))
dev.off()

#write.csv(coverage,"CoverageSamplingEffort/CoverageResults_StandardHaplotype.csv")
#write.csv(coverage,"CoverageSamplingEffort/CoverageResults_TVhaplos_UniqHapTrigger.csv")
write.csv(coverage,"CoverageSamplingEffort/CoverageResults_TVhaplos_40percentSCtrigger.csv")
