#####Eric Crandall and Cynthia Riginos
#####Started: March 2015

genetic.diversity.mtDNA.db<-function(ipdb=ipdb, basic_diversity = T, sequence_diversity = T, coverage_calc = F, coverage_correction = F, minseqs = 6, minsamps = 3, mintotalseqs = 0, ABGD=F,regionalization = c("sample","fn100id", "fn500id", "ECOREGION", "PROVINCE", "REALM", "EEZ"), keep_all_gsls=F, mincoverage = 0.4, hill.number = 0){
  ###Diversity Stats Function###
  #Computes diversity stats by species and population for a flatfile of mtDNA sequences and metadata (with required field $Genus_species_locus)
  # minseqs = minimum sequences per sampled population, 
  # minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
  # keep_all_gsls = set to T if you want to keep the name of all gsls in the output, even if they don't meet the above thresholds
  # mincoverage = minimum coverage for a single population at which the dataset will be reduced to transversions for haplotype diversity computation
  #               set mincoverage to 1 to apply this correction to all species. Set to 0 to apply it to no species.
  # Troubleshooting settings:
  # basic_diversity = T; sequence_diversity = T; coverage_calc = T; coverage_correction = T; minseqs = 6; minsamps = 3; mintotalseqs = 0; ABGD=F;regionalization = "sample"; keep_all_gsls=F; mincoverage = 0.4; hill.number = 0
  
  
  require(seqinr)
  require(ape)
  require(adegenet)
  require(pegas)
  require(mmod)
  require(hierfstat)
  require(iNEXT)
  
  #If ABGD=T and there is an ABGD field in the data frame, then replace Genus_species_locus with this ABGD version.
  if(ABGD==T & "ABGD_Genus_species_locus" %in% colnames(ipdb)){ipdb$Genus_species_locus <- ipdb$ABGD_Genus_species_locus}
  if(ABGD==T & !"ABGD_Genus_species_locus" %in% colnames(ipdb)){cat("No ABGD groupings provided. Using standard taxonomy")}
  
  #create an empty list with headers named for each species/locus combo (or ESU/locus combo)
  esu_loci <- unique(ipdb$Genus_species_locus)
  
  all.pops.table<-sapply(esu_loci, function(x) NULL) 
  
  # LOOP through all gsl combos
  for(gsl in esu_loci){ #gsl<-"Chaetodon_auriga_CYB" 
    
    cat("\n","\n","\n","Now starting", gsl, "\n")
    
    #POPULATION SAMPLES
    
    #modify/replace the following once we have the script that assigns populations
    #subset the genus-species-locus and order it alphabetically by locality name
    sp<-ipdb[which(ipdb$Genus_species_locus==gsl),]
    sp$sample<-paste(sp$locality,round(sp$decimalLatitude, digits=0),round(sp$decimalLongitude, digits=0),sep="_")  #sets up a variable that matches assignsamp function outcome
    sp<-sp[order(sp$sample),]
    
    #this code is not currently necessary given the sp$sample field above may come in handy when we start using ecoregions etc.
    #create a set of unique samples, as denoted by lat+long+pi, and sort it alphabetically. 
    #This is because the genind object will return things in alphaetical order, so we want to be on the same page with it
    #samps<-sort(unique(paste(sp$locality,sp$decimalLatitude,sp$decimalLongitude,sep="_")))
    ##assign each individual to its population sample##
    #sp$sample<-sapply(1:length(sp[,1]), assignsamp) 
    #To be added: add option for rarefaction, fu's Fs and Fu and Li's D
    
    #FILTER
    cat("filtering out population samples with n <", minseqs,"and species with fewer than", minsamps,"total populations \n")
    sampN<-table(sp[[regionalization]])
    lowsamps<-names(sampN[sampN < minseqs])
    if(length(lowsamps)>0){sp<-sp[-which(sp[[regionalization]] %in% lowsamps),]}
    if(length(sampN) - length(lowsamps) < minsamps){cat("Fewer than", minsamps, "sampled populations after filtering. No stats calculated")
      if(keep_all_gsls==T){all.pops.table[[gsl]]<-paste("Fewer than", minsamps, "sampled populations after filtering. No stats calculated")} 
      next}
    if(length(sp[,1])<mintotalseqs){cat("fewer than",mintotalseqs,"samples left after filtering. No stats calculated")
      if(keep_all_gsls==T){all.pops.table[[gsl]]<-paste("fewer than",mintotalseqs,"samples left after filtering. No stats calculated")}
      next}
    cat("Removed the following population samples:", lowsamps, "\n")
    
    #FORMAT CONVERSION
    cat("converting data to various R DNA formats", "\n")
    #convert to seqinr alignment (need to specify that we are using the as.alignment() function from seqinr rather than ape)
    spseqs<-seqinr::as.alignment(nb=length(sp[,1]),nam=sp$materialSampleID, seq=sp$sequence)
    
    #convert to ape DNAbin format (ape)
    spseqsbin<-as.DNAbin(spseqs)
    
    #convert to genind format (adegenet and mmod) with various pieces of information accompanying it (may not be important though)
    spseqs.genind<-as.genind.DNAbin(x = spseqsbin, pops = sp[[regionalization]])
    spseqs.genind@pop <- factor(as.character((sp[[regionalization]]))) #convert to character then to factor
    spseqs.genind@ploidy<-rep(1L,length.out=length(sp[,1]))  #L needed to indicate an integer
    spseqs.genind@other <- sp[,c("decimalLatitude","decimalLongitude")] #and the lat longs, why not
    
    #convert to loci format (pegas) and a data-frame format that works for hierfstat
    spseqs.loci<-as.loci(spseqs.genind)
    spseqs_wc<-as.data.frame(spseqs.loci) #convert locus format to format for hierfstat by making  
    spseqs_wc[,1]<-as.integer(spseqs_wc[,1]) #a vector of integers for population
    spseqs_wc[,2]<-as.numeric(as.character(spseqs_wc[,2])) #and a vector of numeric for haplotype
    spseqs_wc$dummy<-1 #OMFG! You need to have two loci for betai to work, this is a dummy matrix of 1s
    
    
    #Set up the pop.data data frame
    spsummary<-summary(spseqs.genind) # the summary of the genind object does some legwork for us
    
    pop.data<-data.frame(popname=sort(unique(sp[[regionalization]])),sampleN=spsummary$pop.eff) 
    
    populations<-levels(pop(spseqs.genind))  #Returns in same order as used to create pop.data
    
    #BASIC DIVERSITY STATS CALCULATION - Stats that can be calculated for all populations at the same time
    #create a data frame alphabetically sorted by locality to populate with popgen statistics
    #start with sampleN and Unique Haps that are already calculated in the genind object
    if(basic_diversity == T){
      cat("Calculating Basic Diversity Statistics: Haplotype diversity, Shannon-Wiener diversity, Effective Number of Haps, Local FST \n")
      #unique haplotypes
      pop.data$UniqHapNum<-spsummary[[4]]
      #haplotype diversity, calculated by Hs in adegenet corrected for sample size
      pop.data$HaploDiv<-(pop.data$sampleN/(pop.data$sampleN-1))*Hs(spseqs.genind, truenames=TRUE)   #with sample size correction
      #Shannon-Wiener Diversity based on the modified Hs function above
      pop.data$SWdiversity<-shannon.wiener.d(spseqs.genind, truenames=TRUE)
      #Effective number of haplotpyes
      pop.data$EffNumHaplos<-1/(1-Hs(spseqs.genind, truenames=TRUE))   #No sample size correction - based on Crow & Kimura 1964, eq 21. See also Jost 2008 eq 5
      #local Fst (Beta of Weir and Hill 2002 NOT of Foll and Gaggiotti 2006)
      betaWH<-betai_haploid(spseqs_wc)
      pop.data$localFST<-betaWH$betaiov 
    }
    
    #SEQUENCE-BASED DIVERSITY STATS CALCULATION - Stats that need to be calculated one population at a time.
    #now loop through the populations (as delimited by the pop.names in the genind object) 
    #to calculate additional stats that can't be done on a per-population basis as above 
    if(sequence_diversity == T){
      cat("Calculating Sequence-Based Diversity Statistics: Nucleotide diversity, ThetaS, Tajima's D \n")  
      
      for (p in 1:length(populations)) {
        singlepop<-spseqsbin[(spseqs.genind@pop == populations[p] & !is.na(spseqs.genind@pop == populations[p])),]  #DNAbin object containing only the sequences from population p
        #nucleotide diversity, pi (percent)  - based on Nei 1987
        pop.data[p, "NucDivSite"] <- nuc.div( singlepop, variance = FALSE, pairwise.deletion = FALSE)[1]
        pop.data[p, "NucDivSiteVAR"] <- nuc.div( singlepop, variance = TRUE, pairwise.deletion = FALSE)[2]
        pop.data[p, "NucDivLocus"] <- nuc.div( singlepop, variance = FALSE, pairwise.deletion = FALSE)[1]* nchar(sp$sequence[1])
        pop.data[p, "NucDivLocusVAR"] <- nuc.div( singlepop, variance = TRUE, pairwise.deletion = FALSE)[2]* nchar(sp$sequence[1])
        #thetaS - based on Watterson 1975
        pop.data[p, "ThetaS"] <- theta.s(s=length(seg.sites(singlepop)),n=pop.data[p,"sampleN"])
        pop.data[p, "ThetaSVAR"] <- theta.s(s=length(seg.sites(singlepop)),n=pop.data[p,"sampleN"], variance=TRUE)[2]
        tajD<-tajima.test(singlepop)
        pop.data[p, "TajD"] <- tajD[[1]]
        pop.data[p, "TajDp"] <- tajD[[3]] #Pval.beta - p-value assuming that D follows a beta distribution (Tajima, 1989)
      }
    }
    
    ##COVERAGE CALCULATION##
    
    cat("Calculating Coverage \n")
    #list to hold frequency distribution of haplotypes for coverage adjustments
    hap_freq_dist<-list()  
    for (p in 1:length(populations)) {
      hapfreq<-as.data.frame(table(spseqs.loci[(spseqs.genind@pop == populations[p] & !is.na(spseqs.genind@pop == populations[p])),2]))
      hap_freq_dist[p][[1]]<-subset(hapfreq[,2], hapfreq[,2]>0)  #non zero haplotype occurances added to item p in hap_freq_dist list
      f1<-length(which(hap_freq_dist[[p]]==1))
      f2<-length(which(hap_freq_dist[[p]]==2))
      n<-sum(hap_freq_dist[[p]])
      coverage<-1-(f1/n)*(((n-1)*f1)/(((n-1)*f1)+2*f2)) ##Chao & Jost 2012
      if(coverage_calc == T){
        pop.data[p, "CoverageforActualSampleSize"] <- coverage
      }
    }
    
    
    ##COVERAGE STANDARDIZED DIVERSITY##
    if(coverage_correction == T){
      cat("Calculating Coverage Standardized Diversity based on original data \n")
      #find entries in hap_freq_dist with length=1
      singlehap.pops<-which(sapply(hap_freq_dist,length)==1)
      #replace these entries with c(1,1) just as a placeholder so iNEXT will work
      hap_freq_dist[singlehap.pops]<-lapply(1:length(singlehap.pops),function(x) c(1,1))
      
      #calculate coverage and "species" (haplotype) richness
      coverage <- iNEXT(hap_freq_dist, q=c(hill.number))
      #create vector to hold max SC (species coverage) values and loop through list of dataframes
      max_coverage<-vector()
      for (p in 1:length(populations)) {
        #coverage[3][[1]][[1]][7]   #[3] list of results; [[1]][p] population p; [7] is SC
        max_coverage<-(c(max_coverage, max(coverage[3][[1]][[p]][7])))
      }
      min_SC<-min(max_coverage) #this is the smallest value of maximal coverage achievable across the sampled popualations
      rm(max_coverage)
      #Loop through the dataframes from the coverage output (again), extract the row with the max SC < min_SC and add to temporary dataframe
      SC_standardized_res<-data.frame()
      for (p in 1:length(populations)) {
        popdataframe<-coverage[3][[1]][[p]]
        popdataframe<-subset(popdataframe, SC<=min_SC)
        ifelse(nrow(popdataframe)==0,
               SC_standardized_res<-rbind(SC_standardized_res, coverage[3][[1]][[p]][1,]),
               SC_standardized_res<-rbind(SC_standardized_res, popdataframe[nrow(popdataframe),])  )
        #m is interpolated sample size, qD is diversity, SC is coverage - see iNEXT documentation
        #pop.data[p, "HaploSimplification"] <- "Orig_haplos"
      }
      
      SC_standardized_res[singlehap.pops,]<-data.frame(1,"singleHaplotype",hill.number,1,0,1,1,0,1,stringsAsFactors = F)  #replace the pops with single haplotypes with appropriate values. m=1,method="singleHaplotype,q=1,SC=1)
      pop.data<-cbind(pop.data, SC_standardized_res)  
      #TV = F
      SC<-SC_standardized_res$SC
      
      ##REDUCE TO TRANSVERSIONS ONLY GSL'S CONTAINING AT LEAST ONE POP WITH STANDARDIZED COVERAGE < mincoverage ##  code from LL
      if(any((sapply(SC,max,na.rm=T) < mincoverage) == TRUE)) {
        
        #TV = T
        
        #Loop by population - just need to reduce haplotypes within each population
        for (p in 1:length(populations)) {
          spseqsbin_pop<-spseqsbin[(spseqs.genind@pop == populations[p] & !is.na(spseqs.genind@pop == populations[p])),]
          h<-haplotype(spseqsbin_pop)
          #dist.TV will give you a lower triangular matrix of transversions between every sequence in the dataset
          dist.TV <- dist.dna(spseqsbin_pop, model = "TV") #calculate the distance among sequences in the datasets based on transversions only
          #put all the pairwise comparisons among sequences into a dataframe
          m <- nrow(spseqsbin_pop)  #m is the number of sequences in the dataset
          prs <- cbind(rep(1:m, each = m), 1:m)  #prs is all of the possible pairwise comparisons shifted into columns (as above, just different format)
          prs.comp <- prs[prs[, 1] < prs[, 2], , drop = FALSE]  #this drops the diagonal comparisons, i.e. sequence 1 with sequence 1
          #rename the elements in the column so that they are the actual names of the sequences, rather than an integer telling us their position. Append the pairwise distance among sequences based on transversions only
          seq.df <- data.frame(t(apply(prs.comp, 1, function(x) labels(dist.TV)[x])), TV.dist = as.vector(dist.TV))
          h.lib <- data.frame(labels(dist.TV), 0) #allocating sequences with 0 distance between them to haplotypes
          
          for (i in 1:m){
            if(h.lib[i,2] == 0){
              index <- labels(dist.TV)[i]
              tmp <- seq.df[seq.df[,1] == index, ]    
              tmp.2 <- c(index, levels(droplevels(tmp[tmp[,3] == 0, 2])))   #identifies haplotype labels with zero distance
              h.lib[h.lib[,1] %in% tmp.2, 2] <- i                         #identical haplotypes are labeled as "i"
            }
          }
          
          TV_defined_hapfreqs<-as.data.frame(with(h.lib,table(X0)))
          hap_freq_dist[p][[1]]<-TV_defined_hapfreqs[,2]       
        }
        
        ##REDO COVERAGE AND STANDARDIZED COVERAGE CALCS FOR TV HAPLOTYPES
        for (p in 1:length(populations)) {
          f1<-length(which(hap_freq_dist[[p]]==1))
          f2<-length(which(hap_freq_dist[[p]]==2))
          n<-sum(hap_freq_dist[[p]])
          coverage<-1-(f1/n)*(((n-1)*f1)/(((n-1)*f1)+2*f2)) ##Chao & Jost 2012
          pop.data[p, "CoverageforActualSampleSize_TV"] <- coverage
        }
        
        ##COVERAGE STANDARDIZED DIVERSITY##
        cat("Calculating Coverage Standardized Diversity based on transversion defined haplotypes\n")
        #find entries in hap_freq_dist with length=1
        singlehap.pops<-which(sapply(hap_freq_dist,length)==1)
        #replace these entries with c(1,1) just as a placeholder so iNEXT will work
        hap_freq_dist[singlehap.pops]<-lapply(1:length(singlehap.pops),function(x) c(1,1))
        #calculate coverage and "species" (haplotype) richness
        TVcoverage <- iNEXT(hap_freq_dist, q=hill.number) 
        #ERIC - maybe there is a more efficient way to do this loop with an apply function?
        #create vector to hold max SC (species coverage) values and loop through list of dataframes
        max_coverage<-vector()
        for (p in 1:length(populations)) {
          #coverage[3][[1]][[1]][7]   #[3] list of results; [[1]][p] population p; [7] is SC
          max_coverage<-(c(max_coverage, max(TVcoverage[3][[1]][[p]][7])))
        }
        
        
        min_SC<-min(max_coverage) #this is the smallest value of maximal coverage achievable across the sampled popualations
        rm(max_coverage)
        #Loop through the dataframes from the coverage output (again), extract the row with the max SC < min_SC and add to temporary dataframe
        SC_TV_standardized_res<-data.frame()
        for (p in 1:length(populations)) {
          popdataframe<-TVcoverage[3][[1]][[p]]
          popdataframe<-subset(popdataframe, SC<=min_SC)
          ifelse(nrow(popdataframe)==0,
                 SC_TV_standardized_res<-rbind(SC_TV_standardized_res, TVcoverage[3][[1]][[p]][1,]),
                 SC_TV_standardized_res<-rbind(SC_TV_standardized_res, popdataframe[nrow(popdataframe),])  )
          #m is interpolated sample size, qD is diversity, SC is coverage - see iNEXT documentation
          
        }  
        names(SC_TV_standardized_res)<-paste(names(SC_TV_standardized_res),"TV",sep="_")
        
        SC_TV_standardized_res[singlehap.pops,]<-data.frame(1,"singleHaplotype",hill.number,1,1,1,1,0,1,stringsAsFactors = F)  #replace the pops with single haplotypes with appropriate values. m=1,method="singleHaplotype,q=1,SC=1)
        pop.data<-cbind(pop.data, SC_TV_standardized_res)
        
      } #ends (any((sapply(SC_standardized_res$SC,max) < mincoverage) == TRUE)
      
      
      
      
      
      
    } #end coverage_correction 
    
    all.pops.table[[gsl]]<-pop.data
    
  }  #end gsl esu_loci
  if(keep_all_gsls==F) {all.pops.table<-all.pops.table[!sapply(all.pops.table, is.null)]} # remove the NULL gsls if they were not requested
  return(all.pops.table)
}  #end genetic.diversity.mtDNA.db




pairwise.structure.mtDNA.db<-function(ipdb=ipdb, gdist = c("Nei GST", "Hedrick G'ST", "Jost D", "WC Theta", "PhiST", "Chi2", "NL dA"), minseqs = 5, minsamps = 3, mintotalseqs = 0, nrep = 0, num.cores = 1, ABGD = F, regionalization = c("sample","fn100id", "fn500id", "ECOREGION", "PROVINCE", "REALM", "EEZ")){
  ###Genetic Structure Function###
  #Computes genetic differentiation statistics by species and population for a flatfile of mtDNA sequences and metadata (with required field $Genus_species_locus)
  # gdist = You must choose one genetic distance to calculate
  # minseqs = minimum sequences per sampled population, 
  # minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
  # mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
  # nreps = number of resampling permutations for WC Theta, PhiST, Chi2, and NL dA (strataG package). This is working, but you will not currently see p-values in the output, so not much use for now
  # num.cores = number of computer cores to devote to computations for WC Theta, PhiST, Chi2, and NL dA (strataG package)
  # To be added: option to output square matrices with p-values.
  
  require(seqinr)
  require(ape)
  require(adegenet)
  require(pegas)
  require(mmod)
  require(hierfstat)
  require(strataG)
  
  ###LOOP THROUGH THE SPECIES###
  gdistlist<-c("Nei GST", "Hedrick G'ST", "Jost D", "WC Theta", "PhiST", "Chi2", "NL dA")
  if(!gdist %in% gdistlist){
    stop("Please select a genetic distance from the following", gdistlist)
  }
  
  #If ABGD=T and there is an ABGD field in the data frame, then replace Genus_species_locus with this ABGD version.
  if(ABGD==T & "ABGD_Genus_species_locus" %in% colnames(ipdb)){ipdb$Genus_species_locus <- ipdb$ABGD_Genus_species_locus}
  if(ABGD==T & !"ABGD_Genus_species_locus" %in% colnames(ipdb)){cat("No ABGD groupings provided. Using standard taxonomy")}
  
  #create an empty list with headers named for each species/locus combo (or ESU/locus combo)
  esu_loci <- unique(ipdb$Genus_species_locus)
  
  all.pops.table<-sapply(esu_loci, function(x) NULL) 
  
  # LOOP through all gsl combos  
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
    if(length(sampN) - length(lowsamps) < minsamps){all.pops.table[[gsl]]<-paste("Fewer than", minsamps, "sampled populations after filtering. No stats calculated")
                                                    paste("Fewer than", minsamps, "sampled populations after filtering. No stats calculated")
                                                    next}
    if(length(sp[,1])<mintotalseqs){all.pops.table[[gsl]]<-paste("fewer than",mintotalseqs,"samples left after filtering. No stats calculated")
                                   cat("fewer than",mintotalseqs,"samples left after filtering. No stats calculated")
                                   next}
    cat("Removed the following population samples:", lowsamps, "\n")
    
    #FORMAT CONVERSION
    cat("converting data to various R DNA formats", "\n")
    #convert to seqinr alignment (need to specify that we are using the as.alignment() function from seqinr rather than ape)
    spseqs<-seqinr::as.alignment(nb=length(sp[,1]),nam=sp$materialSampleID, seq=sp$sequence)
    
    #convert to ape DNAbin format (ape)
    spseqsbin<-as.DNAbin(spseqs)
    
    #convert to genind format (adegenet and mmod) with various pieces of information accompanying it (may not be important though)
    spseqs.genind<-as.genind.DNAbin(x = spseqsbin, pops = sp[[regionalization]])
    spseqs.genind@pop <- factor(as.character((sp[[regionalization]]))) #convert to character then to factor
    spseqs.genind@ploidy<-rep(1L,length.out=length(sp[,1]))
    spseqs.genind@other <- sp[,c("decimalLatitude","decimalLongitude")] #and the lat longs, why not
    
    #convert to gtypes format (strataG)
    seqs<-sp$sequence
    names(seqs)<-sp$materialSampleID
    gseqs<-as.dna.seq(seqs)
    gseqhaps<-label.haplotypes(gseqs,prefix="H")
    seq.gtype<-gtypes(gen.data=data.frame(sp$materialSampleID,sp[[regionalization]],gseqhaps$haps),id.col=1,strata.col=2,locus.col=3,dna.seq=gseqhaps$hap.seqs)
    
    #DIFFERENTIATION STATS CALCULATION
    cat("Calculating", gdist)

      
    #Nei's Gst (Nei 1973, Nei and Chesser 1983) - mmod package
    if(gdist=="Nei GST"){
        diffs<-pairwise_Gst_Nei(spseqs.genind)
      }
        
    #Hedrick's G'st (Hedrick 2005) mmod package
    if(gdist=="Hedrick G'ST"){
        diffs<-pairwise_Gst_Hedrick(spseqs.genind)
      }
    
    #Jost's D (Jost 2008) - mmod package
    if(gdist=="Jost D"){
      diffs<-pairwise_D(spseqs.genind)
    }
    
    #Weir-Cockerhams (1984) theta - StrataG package (also commented code for hierfstat package - results are the same)
    if(gdist=="WC Theta"){
      pairwise<-pairwise.test(seq.gtype,stats="fst",nrep=nrep,num.cores=num.cores,quietly=T)
      diffs<-as.dist(t(pairwise$pair.mat$Fst))
    }
    
    
    #PhiST - Excoffier et al. 1992 - StrataG package
    if(gdist=="PhiST"){
      pairwise<-pairwise.test(seq.gtype,stats="phist",nrep=nrep,num.cores=num.cores,quietly=T)
      diffs<-as.dist(t(pairwise$pair.mat$PHIst))
    }
    
    #ChiSq - Raymond and Rousset 1995? I think? - StrataG package
    if(gdist=="Chi2"){
      pairwise<-pairwise.test(seq.gtype,stats="Chi2",nrep=nrep,num.cores=num.cores,quietly=T)
      diffs<-as.dist(t(pairwise$pair.mat$Chi2))
    }
    
    #neis dA Nei and Li 1979 - net divergence (within pop divergence removed) - StrataG package
    if(gdist=="NL dA") {
      nucdiv<-nucleotide.divergence(seq.gtype)
      diffs<-matrix(nrow=19,ncol=19)
      diffs[lower.tri(diffs)]<-nucdiv$between$dA
      diffs<-as.dist(diffs)
      attr(diffs, "Labels") <- names(sampN)
    }
    
    all.pops.table[[gsl]]<-diffs
    
  }
return(all.pops.table)
}



hierarchical.structure.mtDNA.db<-function(ipdb=ipdb, level1=NULL, level2=NULL, level3=NULL, minseqs = 6, minsamps = 3, mintotalseqs = 0, ABGD = F, nperm=10, model="N"){
  #add the ability to filter based on sample size at levels 2 and 3?
  
  ###Hierarchical Genetic Structure Function###
  #Computes hierarchical genetic differentiation statistics by species and population for a flatfile of mtDNA sequences and metadata (with required field $Genus_species_locus)
  # minseqs = minimum sequences per sampled population, 
  # minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
  # mintotalseqs = minimum sampled sequences per species (after pops with n < minseqs have been removed)
  # nperm = number of AMOVA permutations
  # model= model of molecular evolution to be passed to dna.dist() = c("none","N", "raw", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock")
  # model defaults to "N" which is the raw count of differences ("raw" is the proportion- same thing). If you use model = "none" you will get all distances between haplotypes = 1, which is the same as "regular" FST
  # levels can be one of c("sample","fn100id", "fn500id", "ECOREGION", "PROVINCE", "REALM", "EEZ") or new regionalizations as they are added.
  
  require(seqinr)
  require(ape)
  require(adegenet)
  require(pegas)
  

  ##Code to figure out the number of levels in the analysis
  if(is.null(level1)){stop("You must enter at least one level - try using 'sample'")}
  if(is.null(level2) && is.null(level3)){level<-1}
  if(!is.null(level2) && is.null(level3)){level<-2}
  if(!is.null(level2) && !is.null(level3)){level<-3}
  
  #If ABGD=T and there is an ABGD field in the data frame, then replace Genus_species_locus with this ABGD version.
  if(ABGD==T & "ABGD_Genus_species_locus" %in% colnames(ipdb)){ipdb$Genus_species_locus <- ipdb$ABGD_Genus_species_locus}
  if(ABGD==T & !"ABGD_Genus_species_locus" %in% colnames(ipdb)){cat("No ABGD groupings provided. Using standard taxonomy")}
  
  #create an empty list with headers named for each species/locus combo (or ESU/locus combo)
  esu_loci <- unique(ipdb$Genus_species_locus)
  
  all.pops.table<-sapply(esu_loci, function(x) NULL) 
  
  #LOOP through all gsl combos
  for(gsl in esu_loci){ #gsl<-"Zebrasoma_flavescens_CYB" 
    
    cat("\n","\n","\n","Now starting", gsl, "\n")
    
    # Pass in the intended number of levels
    l<-level
    l1<-level1
    l2<-level2
    l3<-level3
    
    #POPULATION SAMPLES
    
    #modify/replace the following once we have the script that assigns populations
    #subset the genus-species-locus and order it alphabetically by locality name
    sp<-ipdb[which(ipdb$Genus_species_locus==gsl),]
    sp$sample<-paste(sp$locality,round(sp$decimalLatitude, digits=0),round(sp$decimalLongitude, digits=0),sep="_")#sets up a variable that matches assignsamp function outcome
    sp<-sp[order(sp$sample),]
    
    
    #FILTER
    cat("filtering out population samples with n <", minseqs,"and species with fewer than", minsamps,"total populations \n")
    #filter out populations with low numbers (n<5 for now). Skip the species if this procedure deletes all samples or if there are less than 3 populations total.
    sampN<-table(sp[[l1]])
    lowsamps<-names(sampN[sampN < minseqs])
    if(length(lowsamps)>0){sp<-sp[-which(sp[[l1]] %in% lowsamps),]}
    if(length(sampN) - length(lowsamps) < minsamps){all.pops.table[[gsl]]<-paste("Fewer than", minsamps, "sampled populations after filtering. No stats calculated")
                                                    paste("Fewer than", minsamps, "sampled populations after filtering. No stats calculated")
                                                    next}
    if(length(sp[,1])<mintotalseqs){all.pops.table[[gsl]]<-paste("fewer than",mintotalseqs,"samples left after filtering. No stats calculated")
                                    cat("fewer than",mintotalseqs,"samples left after filtering. No stats calculated")
                                    next}
    cat("Removed the following population samples:", lowsamps, "\n")
    

    #Adjust l for the number of groups at each hierarchical level. If only one group at a certain level, reduce l by 1. If the same number at two levels, then reduce l by 1
    if(l==3 && length(unique(sp[[l3]]))==1){l<-2}
    if(l==2 && length(unique(sp[[l2]]))==1){l<-1}
    if(l==3 && length(unique(sp[[l3]]))>=length(unique(sp[[l2]]))){l<-2}  #if there are the same number groups at level3 as level2, then just treat it as a two level AMOVA
    if(l==2 && length(unique(sp[[l2]]))>=length(unique(sp[[l1]]))){l<-1} #if there are the same number groups at level2 as level1, then just treat it as a one level AMOVA
    if(l==3 && length(unique(sp[[l2]]))>=length(unique(sp[[l1]]))){l<-2; l1<-level2; l2<-level3}
    
    #remove NAs
    if(l==2){sp<-sp[!(is.na(sp[[l2]])),]} # remove any individuals that have NAs for level2
    
    if(l==3){sp<-sp[!(is.na(sp[[l3]])),]} # remove any individuals that have NAs for level3
    # or possibly keep them in "other?"

    #FORMAT CONVERSION
    cat("converting data to various R DNA formats", "\n")
    #convert to seqinr alignment (need to specify that we are using the as.alignment() function from seqinr rather than ape)
    spseqs<-seqinr::as.alignment(nb=length(sp[,1]),nam=sp$materialSampleID, seq=sp$sequence)
    
    #convert to ape DNAbin format (ape)
    spseqsbin<-as.DNAbin(spseqs)

    
    
    #HEIRARCHICAL DIFFERENTIATION STATS CALCULATION
    #pie charts? stacked bars?
    #for eventual stacked bar function barplot(height = amova_out$varcomp[,1]/sum(amova_out$varcomp[,1]),beside = T)
    
    #Calculate distances among individuals for PhiST. Need to write it to the global environment, because that is where pegas::amova() will look for it
    #change dists to all 1's if Fst is requested (model = "none")
    if(model=="none"){dists<<-dist.dna(spseqsbin, model="N"); .GlobalEnv$dists[.GlobalEnv$dists>0]<-1} else{dists<<-dist.dna(spseqsbin, model=model)}
   
    # add evolutionary distance along a tree?
    
    #1 level AMOVA
    if(l==1){
      level1factor<<-as.factor(sp[[l1]]) #Need to write it to the global environment, because that is where pegas::amova() will look for it
      #amova_df<-data.frame(level1factor)
      amova_out<-pegas::amova(dists~level1factor, nperm=nperm, is.squared=T)
      FST<-amova_out$varcomp[1,1]/sum(amova_out$varcomp[,1])
      diffs<-list("raw_amova_output"=amova_out,"level1_names"=levels(level1factor),"FST"=FST)
      
      rm(dists, level1factor, envir=.GlobalEnv)
   }
    
    #2 level AMOVA (per Arlequin 3.5 manual)
    if(l==2){
      
      level1factor<<-as.factor(sp[[l1]])
      level2factor<<-as.factor(sp[[l2]])
      amova_out<-pegas::amova(dists~level2factor/level1factor, nperm=nperm, is.squared=T)
      
      FCT<-amova_out$varcomp[1,1]/sum(amova_out$varcomp[,1])
      FCTp<-amova_out$varcomp[1,2]
     
      FSC<-amova_out$varcomp[2,1]/(amova_out$varcomp[2,1]+amova_out$varcomp[3,1])
      FSCp<-amova_out$varcomp[2,2]
      
      FST<-(amova_out$varcomp[1,1]+amova_out$varcomp[2,1])/sum(amova_out$varcomp[,1])
      FSTp<-NA
      
      diffs<-list("raw_amova_output"=amova_out,"level2_names"=levels(level2factor), "level1_names"=levels(level1factor), "FCT"=FCT, "FSC"=FSC,"FST"=FST)
     # if(round((1-FST),digits=4) != round((1-FSC)*(1-FCT),digits = 4)){warning("Variance components don't meet expectation of (1-FST)==(1-FSC)*(1-FCT) for ", gsl)}
    
      rm(dists, level1factor, level2factor, envir=.GlobalEnv)
    }
    
  
    #3 level AMOVA
    if(l==3){
    
      level1factor<<-as.factor(sp[[l1]])
      level2factor<<-as.factor(sp[[l2]])
      level3factor<<-as.factor(sp[[l3]])
      amova_out<-pegas::amova(dists~level3factor/level2factor/level1factor, nperm=nperm, is.squared=T)
      
      FRT<-amova_out$varcomp[1,1]/sum(amova_out$varcomp[,1])
      FRTp<-amova_out$varcomp[1,2]
      
      FCR<-amova_out$varcomp[2,1]/(amova_out$varcomp[2,1]+amova_out$varcomp[3,1]+amova_out$varcomp[4,1])
      FCRp<-amova_out$varcomp[2,2]
      
      FSC<-amova_out$varcomp[3,1]/(amova_out$varcomp[3,1]+amova_out$varcomp[4,1])
      FSCp<-amova_out$varcomp[3,2]
      
      FST<-(amova_out$varcomp[1,1]+amova_out$varcomp[2,1]+amova_out$varcomp[3,1])/sum(amova_out$varcomp[,1])
      FSTp<-NA
      
      
      diffs<-list("raw_amova_output"=amova_out,"level3_names"=levels(level3factor),"level2_names"=levels(level2factor), "level1_names"=levels(level1factor), "FRT"=FRT, "FCR"=FCR,"FSC"=FSC, "FST"=FST)
      #if(round((1-FST),digits=4) != round((1-FSC)*(1-FCR)*(1-FRT),digits = 4)){warning("Variance components don't meet expectation of (1-FST)==(1-FSC)*(1-FCT) for ", gsl)}
   
      rm(dists, level1factor, level2factor, level3factor, envir=.GlobalEnv) 
    }
    
    #change this
    all.pops.table[[gsl]]<-diffs
  
    
  }
  return(all.pops.table)
}



# this function, assignsamp, will lookup the sample index number in samps and sapply it across all individuals in sp
assignsamp<-function(x) {which(samps==paste(sp$locality[x],sp$decimalLatitude[x],sp$decimalLongitude[x],sep="_"))
} 

#Adapt the Hs function to make Shannon-Wiener diversity         
shannon.wiener.d<-function (x, truenames = TRUE) 
{
  if (is.genind(x)) {
    x <- genind2genpop(x, quiet = TRUE)
  }
  if (!is.genpop(x)) 
    stop("x is not a valid genpop object")
  if (x@type == "PA") 
    stop("not implemented for presence/absence markers")
  x.byloc <- seploc(x)
  lX <- lapply(x.byloc, function(e) makefreq(e, quiet = TRUE))
  lres<-lapply(lX, function(X) apply((-X+0.0000000001)*log(X+0.0000000001), 1, sum))
  res <- apply(as.matrix(data.frame(lres)), 1, mean)
  return(res)
}

#adapt the betai function from heirfstat for haploid data (add diploid=F to pop.freq() function call in this function)
betai_haploid<-function (gendata) 
{
  nloc <- dim(gendata)[2] - 1
  npop <- dim(table(gendata[, 1]))
  ninds <- 2 * ind.count(gendata)
  p <- pop.freq(gendata, diploid=F)
  pb <- pop.freq(cbind(rep(1, dim(gendata)[1]), gendata[, -1]), diploid=F)
  n2 <- ninds^2
  n2 <- sweep(n2, 2, apply(ninds, 2, sum), FUN = "/")
  nic <- ninds - n2
  snic <- apply(nic, 2, sum)
  betas <- array(dim = c(npop, npop, nloc))
  nums <- array(dim = c(npop, npop, nloc))
  lden <- numeric(nloc)
  for (il in 1:nloc) {
    dum1 <- sweep(sweep(p[[il]], 1, pb[[il]], FUN = "-")^2, 
                  2, ninds[, il], FUN = "*")
    dum2 <- sweep(p[[il]] * (1 - p[[il]]), 2, nic[, il], 
                  FUN = "*")
    sden <- sum(dum1 + dum2)
    lden[il] <- sden
    for (ip1 in 1:npop) {
      for (ip2 in ip1:npop) {
        if (ip1 == ip2) {
          dum1 <- p[[il]][, ip1] * (1 - p[[il]][, ip1])
          nums[ip1, ip1, il] <- snic[il] * ninds[ip1, 
                                                 il]/(ninds[ip1, il] - 1) * sum(dum1)
          betas[ip1, ip1, il] <- 1 - nums[ip1, ip1, il]/sden
        }
        else {
          dum1 <- p[[il]][, ip1] * (1 - p[[il]][, ip2])
          dum2 <- p[[il]][, ip2] * (1 - p[[il]][, ip1])
          nums[ip1, ip2, il] <- snic[il] * sum(dum1 + 
                                                 dum2)
          betas[ip1, ip2, il] <- 1 - nums[ip1, ip2, il]/2/sden
        }
      }
    }
  }
  betai <- t(apply(betas, 3, diag))
  betaio <- 1 - apply(apply(nums, 3, diag), 1, sum)/sum(lden)
  betaw <- 1 - apply(apply(nums, 3, diag) * ninds, 2, sum)/(apply(ninds, 
                                                                  2, sum) * lden)
  return(list(betaiip = betas, nic = snic, betai = betai, betaiov = betaio, 
              betaw = betaw))
}

#Write the output of genetic.diversity.mtDNA.db and genetic.structure.mtDNA.db to a csv file readable by Excel
write.stats<-function(x=divstats,filename=NULL,structure=F){
  outfile<-file(filename)
  for(gsl in names(x)){
    
    if(structure){
      dm<-as.matrix(x[[gsl]])
      dm[upper.tri(dm)]<-NA
    }
    else{dm<-x[[gsl]]}
    
    dm<-as.data.frame(dm)
    open(outfile,open="a")
    writeLines(paste("\n",gsl,sep=""),outfile,sep="\n")
    if(structure){write.table(dm, file=outfile, sep=",",append=T, quote=F, col.names=NA, row.names=T, na="")} else
      {write.table(dm, file=outfile, sep=",",append=T, quote=F, col.names=T, row.names=F, na="")}
  }
  close(outfile)
}

# another function that uses sink()?


