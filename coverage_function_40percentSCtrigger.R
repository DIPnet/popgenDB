coverage.db<-function(ipdb=ipdb, minseqs = 5, minsamps = 3, mintotalseqs = 0, regionalization = c("sample","fn100id", "fn500id", "meow_ecoregion", "meow_prov_name", "meow_rlm_name", "EEZ_country")){
  
  ###Based on Diversity Stats Function###
  #Computes coverage stats by species and population for a flatfile of mtDNA sequences and metadata (with required field $Genus_species_locus)
  # minseqs = minimum sequences per sampled population, 
  # minsamps = minimum sampled populations per species (after pops with n < minseqs have been removed)
   
  require(seqinr)
  require(ape)
  require(adegenet)
  require(pegas)
  require(mmod)
  require(hierfstat)
  require(iNEXT)
  #require(plyr)
  
  #create an empty list with headers named for each species/locus combo (or ESU/locus combo)
  esu_loci <- unique(ipdb$Genus_species_locus) 
  coverage_results<-data.frame()
   
  for(gsl in esu_loci){ #gsl<-"Chaetodon_auriga_CYB" 
    
    cat("\n","\n","\n","Now starting", gsl, "\n")
    
    #POPULATION SAMPLES
    
    #modify/replace the following once we have the script that assigns populations
    #subset the genus-species-locus and order it alphabetically by locality name
    sp<-subset(ipdb, Genus_species_locus == gsl )
    sp$sample<-paste(sp$locality,round(sp$decimalLatitude, digits=0),round(sp$decimalLongitude, digits=0),sep="_")  #sets up a variable that matches assignsamp function outcome
    sp<-sp[order(sp$sample),]
    #FILTER
    cat("filtering out population samples with n <", minseqs,"and species with fewer than", minsamps,"total populations \n")
    sampN<-table(sp[[regionalization]])
    lowsamps<-names(sampN[sampN < minseqs])
    if(length(lowsamps)>0){sp<-sp[-which(sp[[regionalization]] %in% lowsamps),]}
    if(length(sampN) - length(lowsamps) < minsamps){ 
                                                    next}
    if(length(sp[,1])<mintotalseqs){
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
    spseqs.genind@ploidy<-1L  #L needed to indicate an integer
    spseqs.genind@loc.names<-strsplit(sp$Genus_species_locus[1],"_")[[1]][3] #add in the locus name, why not
    spseqs.genind@pop.names <- sort((as.character(unique(sp[[regionalization]])))) # I will add the sample name into the pop.names slot here anyway.
    spseqs.genind@other <- sp[,c("decimalLatitude","decimalLongitude")] #and the lat longs, why not
    
    #convert to genpop (adegenet)
    #spseqs.pop<-genind2genpop(spseqs.genind) 
    
    #convert to loci format (pegas) and a data-frame format that works for hierfstat
    spseqs.loci<-as.loci(spseqs.genind)
    spseqs_wc<-as.data.frame(spseqs.loci) #convert locus format to format for hierfstat by making  
    spseqs_wc[,1]<-as.integer(spseqs_wc[,1]) #a vector of integers for population
    spseqs_wc[,2]<-as.numeric(as.character(spseqs_wc[,2])) #and a vector of numeric for haplotype
    spseqs_wc$dummy<-1 #OMFG! You need to have two loci for betai to work, this is a dummy matrix of 1s
    
    
    #DIVERSITY STATS CALCULATION - Stats that can be calculated for all populations at the same time
    #create a data frame alphabetically sorted by locality to populate with popgen statistics
    #start with sampleN and Unique Haps that are already calculated in the genind object
    spsummary<-summary(spseqs.genind) # the summary of the genind object does some legwork for us
    cat("Creating Data Frame")
    pop.data<-data.frame(row.names=NULL, popname=sort(unique(sp[[regionalization]])),sampleN=spsummary$pop.eff) 
    pop.data$UniqHapNum<-spsummary[[4]]
    
    ##HAPLOTYPES PER POPULATION##
    cat("Calculating Haplotypes per Population \n")  
    #list to hold frequency distribution of haplotypes for coverage adjustments
    hap_freq_dist<-list()  
    populations<-spseqs.genind$pop.names  #Returns in same order as used to create pop.data
    for (p in 1:length(populations)) {
      hapfreq<-as.data.frame(table(spseqs.loci[(spseqs.genind@pop == populations[p] & !is.na(spseqs.genind@pop == populations[p])),2]))
      hap_freq_dist[p][[1]]<-subset(hapfreq[,2], hapfreq[,2]>0)  #non zero haplotype occurances added to item p in hap_freq_dist list
     }
        
 
    ##COVERAGE CALCULATION##
    cat("Calculating Coverage \n")  
    #list to hold frequency distribution of haplotypes for coverage adjustments
     for (p in 1:length(populations)) {
        f1<-length(which(hap_freq_dist[[p]]==1))
        f2<-length(which(hap_freq_dist[[p]]==2))
        n<-sum(hap_freq_dist[[p]])
        coverage<-1-(f1/n)*(((n-1)*f1)/(((n-1)*f1)+2*f2)) ##Chao & Jost 2012
        pop.data[p, "CoverageforActualSampleSize"] <- coverage
        }
    
    ##COVERAGE STANDARDIZED DIVERSITY##
    cat("Calculating Coverage Standardized Diversity \n")
    #find entries in hap_freq_dist with length=1
    singlehap.pops<-which(sapply(hap_freq_dist,length)==1)
    #replace these entries with c(1,1) just as a placeholder so iNEXT will work
    hap_freq_dist[singlehap.pops]<-lapply(1:length(singlehap.pops),function(x) c(1,1))
    #calculate coverage and "species" (haplotype) richness
    coverage <- iNEXT(hap_freq_dist, q=c(0))  #Hill number of 0
    #ERIC - maybe there is a more efficient way to do this loop with an apply function?
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
    }
    
    ##REDUCE TO TRANSVERSIONS ONLY GSL'S CONTAINING AT LEAST ONE POP WITH STANDARDIZED COVERAGE < 0.4 ##  code from LL
    if (any((sapply(SC_standardized_res$SC,max) < 0.4) == TRUE)) {
      #Loop by populaiton - just need to reduce haplotypes within each population
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
        pop.data[p, "CoverageforActualSampleSize"] <- coverage
      }
      
      ##COVERAGE STANDARDIZED DIVERSITY##
      cat("Calculating Coverage Standardized Diversity \n")
      #find entries in hap_freq_dist with length=1
      singlehap.pops<-which(sapply(hap_freq_dist,length)==1)
      #replace these entries with c(1,1) just as a placeholder so iNEXT will work
      hap_freq_dist[singlehap.pops]<-lapply(1:length(singlehap.pops),function(x) c(1,1))
      #calculate coverage and "species" (haplotype) richness
      coverage <- iNEXT(hap_freq_dist, q=c(0))  #Hill number of 0
      #ERIC - maybe there is a more efficient way to do this loop with an apply function?
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
        }  
      }
    
    SC_standardized_res[singlehap.pops,]<-NA  #replace the pops with single haplotypes with NA
    pop.data<-cbind(pop.data, SC_standardized_res$SC)
      
    pop.data$gsl<-gsl
    pop.data$minseqs<-minseqs
 
    coverage_results<-rbind(coverage_results, pop.data)
  }  #end gsl esu_loci
  
  names(coverage_results)[5]<-"SC"
  return(coverage_results)
  
}  #end coverage.db