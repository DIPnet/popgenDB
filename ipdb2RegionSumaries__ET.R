# To summarise ipdb by region for sampled species
#  Create a site (row) by species (column) matrix 
#  of samples per species per region (ECO_CODE first)
#
#  etreml 1Dec2015 
library(labdsv)
library(vegan)   #for sites (rows) by species (columns) community analysis
library(bipartite)
library(ggplot2)
minSamp<-3

#load data from latest spatial ipdb (35,897 samples, dated 5June2015)
#DP.datfile<-'C:/etreml/DIPnet/IPDB/ipdb_sp.tsv'
#dat<-read.table(DP.datfile, header = TRUE, sep = "\t", stringsAsFactors=F,na.strings=c("","NA","#N/A"),quote="",comment.char="",strip.white=T,fill=T)


#using the ipdb object from DIPnet_Stats_Script.R
summary(ipdb)
reg<-'VeronDivis2'                                     # set region of interest
#length(unique(ipdb$fn500id))                         # 238 
#length(unique(ipdb$fn100id))                         # 571 
#length(unique(ipdb$VeronDivis))                      # 17
#reg.num<-length(unique(ipdb[,reg]))                  # 114 ECO_CODE (some 0s)
sp.num<-length(unique(ipdb$Genus_species_locus))              # 281 

#preallocate matrix and build RegionBySpecies mat
regIDs<-sort(unique(ipdb[,reg]))                  # unique reg IDs
#regIDs<-sort(as.numeric(unique(ipdb[,reg])))      # use this if you are using numeric identifiers for the regions
spNames<-sort(unique(ipdb$Genus_spec))
regIDs<-regIDs[!regIDs %in% 0]                      # remove any zero IDs
reg.num<-length(regIDs)
regSp.mat<-matrix(data=0,nrow=reg.num,ncol=sp.num)
colnames(regSp.mat)<-spNames
rownames(regSp.mat)<-regIDs

i<-1;j<-1;                                          # row/col index
for (r in regIDs) {                                 # for each region
  sub<-ipdb[ipdb[,reg]==r,]
  for (sp in spNames){                              # for each sp
    #sp<-'Sphyrna_lewini_CR'
    cnt<-(sub$Genus_spec==sp)
    regSp.mat[i,j]<-table(cnt)["TRUE"]
    j<-j+1
  }
  i<-i+1
  j<-1
}

#QA/QC with Regionalizations_pivot.xlsx
# regSp.mat[20:21,49:50]
# regIDs[20:21]
# spNames[49:50]
# Qsp<-ipdb[(ipdb$Genus_spec=='Cephalopholis_argus_CO1' & ipdb$ECO_CODE==20088),]
# length(Qsp[,1])

#sort regSp.mat into upper left
regSp.mat[is.na(regSp.mat)]<-0
labdsv::abuocc(regSp.mat,minabu=minSamp,panel=2)    # species/reg

#play with vegan
reg.dis<-vegdist(regSp.mat)
reg.mds0<-isoMDS(reg.dis)
stressplot(reg.mds0,reg.dis)
nestedchecker(regSp.mat)
reg.temp<-nestedtemp(regSp.mat)
plot(reg.temp,kind='incidence', names=TRUE)

#test nestedness (from LL)
#reg.nested<-nestednodf(regSp.mat, order=TRUE, weight=TRUE)
#plot(reg.nested) #weird
#reg.nested$statistic[3]

reg.nr<-bipartite::nestedrank(regSp.mat, method="NODF", return.matrix=TRUE)
dim(reg.nr$nested.matrix)

#plot the nested rank region by species matrix (reg.nr.mat)
reg.nr.mat<-reg.nr$nested.matrix
reg.nr.mat.na<-reg.nr.mat
reg.nr.mat.na[reg.nr.mat.na==0]<-NA
heatmap(reg.nr.mat.na, Rowv=NA, Colv=NA, col=heat.colors(256))


heatmap(reg.nr.mat.na[,1:50], Rowv=NA, Colv=NA, col=heat.colors(256))

#sp.abun<-colSums(reg.nr.mat.na, na.rm=T)

#test if same as vegan. 
nestednodf(reg.nr.mat, order=TRUE, weight=TRUE)


#ToDo
#  1. decide if we need to filter/threshold ecoregions based on
#      minimum number of species and/or minimum number of samples
#      or both.
#  2. compare matrix nestedness between regionaisations - is this appropriate
#     with different numbers of regions per regionalisation?
#  3. extract data data based on thresholds and export to format ready for EC
#     what format is this?

          
