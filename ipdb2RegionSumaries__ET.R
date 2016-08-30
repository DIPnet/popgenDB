# To summarise ipdb by region for sampled species
#  Create a site (row) by species (column) matrix 
#  of samples per species per region (ECO_CODE first)
#
#  etreml 1Dec2015 
library(labdsv)
library(vegan)   #for sites (rows) by species (columns) community analysis
library(bipartite)
library(ggplot2)
library(reshape2)
minSamp<-3


regionalizations<-c("ECOREGION", "PROVINCE", "REALM", "EEZ","Bowen","Keith","Kulbicki_r","Kulbicki_b", "VeronDivis")

reg<-"ECOREGION"

# Make a new regionalization out of VeronDivis that lumps India and Red Sea into Indian Ocean
ipdb$VeronDivis2<-ipdb$VeronDivis
ipdb$VeronDivis2[which(ipdb$VeronDivis=="India")]<-"Indian Ocean"
ipdb$VeronDivis2[which(ipdb$VeronDivis=="Red Sea Plus")]<-"Indian Ocean"

#using the ipdb object from DIPnet_Stats_Script.R
#summary(ipdb)
#reg<-'VeronDivis'                                     # set region of interest
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

#remove NAs
regSp.mat[is.na(regSp.mat)]<-0
    # species/reg

#play with vegan
#reg.dis<-vegdist(regSp.mat)
#reg.mds0<-isoMDS(reg.dis)
#stressplot(reg.mds0,reg.dis)
#nestedchecker(regSp.mat)
#reg.temp<-nestedtemp(regSp.mat)
#plot(reg.temp,kind='incidence', names=TRUE)

#test nestedness (from LL)
#reg.nested<-nestednodf(regSp.mat, order=TRUE, weight=TRUE)
#plot(reg.nested) #weird
#reg.nested$statistic[3]

#sort regSp.mat into upper left
reg.nr<-bipartite::nestedrank(regSp.mat, method="NODF", return.matrix=TRUE)
dim(reg.nr$nested.matrix)

#plot the nested rank region by species matrix (reg.nr.mat)
reg.nr.mat<-reg.nr$nested.matrix
reg.nr.mat.na<-reg.nr.mat
reg.nr.mat.na[reg.nr.mat.na==0]<-NA

#heatmap with GGplot

reg.nr.mat.melt<-melt(reg.nr.mat.na)
colnames(reg.nr.mat.melt)<-c("Region","Species","n")
reg.nr.mat.melt$n[which(reg.nr.mat.melt$n>200)]<-200

pdf(file = paste(reg,".pdf",sep=""),paper = "letter")

labdsv::abuocc(regSp.mat,minabu=minSamp,panel=2)

region_n<-ggplot(reg.nr.mat.melt,aes(y=Region,x=Species,fill=n))
region_n<-region_n+geom_tile()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylim(levels(reg.nr.mat.melt$Region))+xlim(levels(reg.nr.mat.melt$Species)[1:50])+
  scale_fill_gradient(low = "white", high = "red", space = "Lab",
                      na.value = "grey50", guide = "colourbar")
region_n
dev.off()


#Pick out species with complete samples across a given set of regions from Veron (but with Red Sea and India lumped into "Indian Ocean")
IP<-c("Indian Ocean","Coral Triangle","Central Pacific","French Polynesia","Northern South China Sea", "Eastern Indian Ocean")  #for VeronDivis2



reg.IP<-reg.nr.mat.na[which(rownames(reg.nr.mat.na) %in% IP),]
reg.IP<-t(reg.IP)
reg.IP<-reg.IP[which(complete.cases(reg.IP)),]
# leaves only 3 species, two of which (N. albicilla and A. plancii) likely consist of two cryptic species)

# Indian Ocean Coral Triangle Central Pacific French Polynesia Eastern Indian Ocean
# Acanthaster_planci_CO1          104             10              26               22                   21
# Nerita_plicata_CO1              125            104             242              118                   36
# Nerita_albicilla_CO1            127             99             151               30                   59

#take Eastern Indian Ocean out - we get 11 datasets (two of which are Stenella)
IP<-c("Indian Ocean","Coral Triangle","Central Pacific","French Polynesia","Northern South China Sea")  
reg.IP<-reg.nr.mat.na[which(rownames(reg.nr.mat.na) %in% IP),]
reg.IP<-t(reg.IP)
reg.IP<-reg.IP[which(complete.cases(reg.IP)),]

# Indian Ocean Coral Triangle Central Pacific French Polynesia Northern South China Sea
# Acanthaster_planci_CO1              104             10              26               22                        6
# Panulirus_penicillatus_CO1           26              4             149              120                       13
# Stenella_longirostris_CR             15             12              16               50                        4
# Stenella_longirostris_CYB            15             13              15               47                        4
# Nerita_plicata_CO1                  125            104             242              118                       30
# Nerita_albicilla_CO1                127             99             151               30                       51
# Tripneustes_gratilla_CO1             24            346              31               10                       10
# Naso_unicornis_CR                    21              3              62               18                        1
# Scarus_psittacus_CR                  20              2              21               16                       12
# Diadema_savignyi_A68                 14             10               9               23                       12
# Dascyllus_aruanus_CYB                52             64              38               24                       52

# remove N South China Sea
IP<-c("Indian Ocean","Coral Triangle","Central Pacific","French Polynesia")  
reg.IP<-reg.nr.mat.na[which(rownames(reg.nr.mat.na) %in% IP),]
reg.IP<-t(reg.IP)
reg.IP<-reg.IP[which(complete.cases(reg.IP)),]

# Indian Ocean Coral Triangle Central Pacific French Polynesia
# Acanthaster_planci_CO1              104             10              26               22
# Panulirus_penicillatus_CO1           26              4             149              120
# Stenella_longirostris_CR             15             12              16               50
# Stenella_longirostris_CYB            15             13              15               47
# Nerita_plicata_CO1                  125            104             242              118
# Nerita_albicilla_CO1                127             99             151               30
# Tripneustes_gratilla_CO1             24            346              31               10
# Naso_unicornis_CR                    21              3              62               18
# Scarus_psittacus_CR                  20              2              21               16
# Diadema_savignyi_A68                 14             10               9               23
# Cephalopholis_argus_CYB             104             23              79              142
# Dascyllus_aruanus_CYB                52             64              38               24
# Dascyllus_trimaculatus_CR            18             65              41                9
# Naso_brevirostris_CR                 19              8              30                7


# add Hawaii 
IP<-c("Indian Ocean","Coral Triangle","Central Pacific","French Polynesia", "Hawaii")  
reg.IP<-reg.nr.mat.na[which(rownames(reg.nr.mat.na) %in% IP),]
reg.IP<-t(reg.IP)
reg.IP<-reg.IP[which(complete.cases(reg.IP)),]

# Indian Ocean Coral Triangle Central Pacific French Polynesia Hawaii
# Acanthaster_planci_CO1              104             10              26               22     13
# Panulirus_penicillatus_CO1           26              4             149              120    331
# Stenella_longirostris_CR             15             12              16               50     37
# Stenella_longirostris_CYB            15             13              15               47     36
# Tripneustes_gratilla_CO1             24            346              31               10     10
# Naso_unicornis_CR                    21              3              62               18     24
# Scarus_psittacus_CR                  20              2              21               16     20
# Cephalopholis_argus_CYB             104             23              79              142    236



#Heatmaps from package stats
#heatmap(reg.nr.mat.na, Rowv=NA, Colv=NA, col=heat.colors(256))
#heatmap(reg.nr.mat.na[,1:50], Rowv=NA, Colv=NA, col=heat.colors(256))

#sp.abun<-colSums(reg.nr.mat.na, na.rm=T)

#test if same as vegan. 
nestednodf(reg.nr.mat, order=TRUE, weight=TRUE)

CT<-c("Western Coral Triangle", "Eastern Coral Triangle", "Sunda Shelf", "Andaman", "Java Transitional", "Sahul Shelf", "Northwest Australian Shelf") - no cosampled species # for PROVINCE


CT<-c("Western Coral Triangle", "Eastern Coral Triangle", "Sunda Shelf", "Andaman", "Java Transitional", "Sahul Shelf")
# 
# Western Coral Triangle Java Transitional Sunda Shelf Eastern Coral Triangle Andaman
# Tridacna_maxima_CO1                    379               120          28                      5      48
# Tridacna_crocea_CO1                    404                78          21                      4      61
# Sahul Shelf
# Tridacna_maxima_CO1           4
# Tridacna_crocea_CO1           4

CT<-c("Western Coral Triangle", "Eastern Coral Triangle", "Sunda Shelf", "Andaman", "Java Transitional")

# Western Coral Triangle Java Transitional Sunda Shelf Eastern Coral Triangle
# Nerita_albicilla_CO1                      80                38          19                     10
# Holothuria_atra_CO1                      117                39          18                      5
# Linckia_laevigata_CO1                    610               187          10                     59
# Tridacna_maxima_CO1                      379               120          28                      5
# Tridacna_crocea_CO1                      404                78          21                      4
# Andaman
# Nerita_albicilla_CO1       11
# Holothuria_atra_CO1        27
# Linckia_laevigata_CO1      46
# Tridacna_maxima_CO1        48
# Tridacna_crocea_CO1        61

CT<-c("Western Coral Triangle", "Eastern Coral Triangle", "Sunda Shelf", "Java Transitional")

# Western Coral Triangle Java Transitional Sunda Shelf
# Acanthaster_planciPac_CR                         264                32          48
# Nerita_albicilla_CO1                              80                38          19
# Holothuria_atra_CO1                              117                39          18
# Linckia_laevigata_CO1                            610               187          10
# Tridacna_maxima_CO1                              379               120          28
# Tridacna_crocea_CO1                              404                78          21
# Tridacna_squamosa_CO1                            242                61           6
# Caesio_cuning_CR                                 238                33          40
# Haptosquilla_glyptocercus_CO1                    194                34          10
# Haptosquilla_pulchella_CO1                       289                27          45
# Gonodactylellus_viridis_CO1                      159                20          39

# Eastern Coral Triangle
# Acanthaster_planciPac_CR                          14
# Nerita_albicilla_CO1                              10
# Holothuria_atra_CO1                                5
# Linckia_laevigata_CO1                             59
# Tridacna_maxima_CO1                                5
# Tridacna_crocea_CO1                                4
# Tridacna_squamosa_CO1                              2
# Caesio_cuning_CR                                   2
# Haptosquilla_glyptocercus_CO1                     11
# Haptosquilla_pulchella_CO1                        20
# Gonodactylellus_viridis_CO1                       10



reg.CT<-reg.nr.mat.na[which(rownames(reg.nr.mat.na) %in% CT),]
reg.CT<-t(reg.CT)
reg.CT<-reg.CT[which(complete.cases(reg.CT)),]




#ToDo
#  1. decide if we need to filter/threshold ecoregions based on
#      minimum number of species and/or minimum number of samples
#      or both.
#  2. compare matrix nestedness between regionaisations - is this appropriate
#     with different numbers of regions per regionalisation?
#  3. extract data data based on thresholds and export to format ready for EC
#     what format is this?

          
