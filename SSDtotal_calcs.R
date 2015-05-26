#Eric Crandall 5/20 - just playing with some ways to reduce diversity by testing each haplotype for its effect on SSD total when removed from the dataset. Based on an email from Peter Smouse.
#This may be an interesting way to reduce a dataset by removing the least phylogeographically informative haplotypes, or at least sorting them based on their information to an AMOVA.
#Putting it aside for now.

dists_full<<-dist.dna(spseqsbin, model="N")
level1factor_full<<-as.factor(sp[["meow_ecoregion"]]) 
amova_out_full<-pegas::amova(dists_full~level1factor_full, nperm=0)
FST_full<-amova_out_full$varcomp[1]/sum(amova_out_full$varcomp)
SSDtotal_full<-amova_out_full$tab[3,1]


hspseqsbin<-haplotype(spseqsbin)
#hspseqsbin<-sort(haplotype(spseqsbin))

haptable<-data.frame()

for(hap in 1:length(attr(hspseqsbin,"index"))){
  hapseqnums<-attr(hspseqsbin,"index")[[hap]]
  spseqsbin_1out<-spseqsbin[-hapseqnums,]
  sp_1out<-sp[-hapseqnums,]
  dists<<-dist.dna(spseqsbin_1out, model="N")
  
  level1factor<<-as.factor(sp_1out[["sample"]]) 
  amova_out<-pegas::amova(dists~level1factor, nperm=0)
  haptable[hap,"hapnumber"]<-hap
  haptable[hap,"label"]<-attr(hspseqsbin,"dimnames")[[1]][hap]
  haptable[hap,"abundance"]<-length(hapseqnums)
  haptable[hap,"FST"]<-amova_out$varcomp[1]/sum(amova_out$varcomp)
  haptable[hap,"SSDtotal"]<-amova_out$tab[3,1]
  haptable[hap,"deltaSSD"]<-amova_out$tab[3,1]-SSDtotal_full
}

haptable$zscore<-scale(haptable$deltaSSD, center=T, scale=T)
haptable[order(haptable$deltaSSD),]

#expendables<-haptable2[which(abs(haptable2$deltaSSD)<2),"hapnumber"]
#expendables<-haptable2[which(haptable2$deltaSSD>0),"hapnumber"]
#expendables<-haptable2[which(haptable2$deltaSSD<0),"hapnumber"]
#expendables<-haptable[which(abs(haptable$zscore)<1),"hapnumber"]
expendables<-haptable[which(abs(haptable$deltaSSD)<10),"hapnumber"]
expendables<-haptable[which(haptable$deltaSSD>-100 & haptable$deltaSSD<0),"hapnumber"]


hapseqnums<-attr(hspseqsbin,"index")[expendables]
hapseqnums<-unlist(hapseqnums)
spseqsbin_expendables<-spseqsbin[-hapseqnums,]
sp_expendables<-sp[-hapseqnums,]

#spseqsbin<-spseqsbin_expendables
#sp<-sp_expendables
#write.csv(pop.data,"../Nplicata_ecoregion_deltaSSDgreaterthanneg100.csv")
