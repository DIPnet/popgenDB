Title: README - Population Genetic Statistics for DIPnet 
Author: Eric Crandall  
Date: April 20, 2015  

Population Genetic Statistics for the DIPnet dataset
=====================================================

These scripts are intended for use with the DIPnet dataset that will be finished on May 1, 2015. We plan to convert these scripts eventually into an R package that will talk with the DIPnet database as stored in the FIMS framework to calculate multispecies population genetic statistics. The way it currently works is that all the functions are stored in `DIPnet_Stats_Functions.R` and a script for driving those functions is in `DIPnet_Stats_Script.R`. 

There are major functions for looping over the whole database and calculating:
1. Genetic diversity statistics
	1. Unique Haplotype Number
	2. Haplotype Diversity
	3. Shannon-Wiener Diversity
	4. Effective number of haplotypes
	5. Local FST
	6. Nucleotide Diversity (per site and per locus)
	7. Watterson's ThetaS
	8. Tajima's D
	9. Coverage (Chao and Jost 2012)

2. Genetic structure statistics
	1. Nei's Fst
	2. Nei's Gst
	3. Hedrick's G'st
	4. Jost's D
	5. Weir & Cockerham's Theta
	6. Excoffier et al. PhiST
	7. Chi-Square Test of Population Differentiation
	8. Nei's Da
	
3. Hierarchical Structure using AMOVA


May want to break these functions down into component parts (the looping and stats individually for example) to make the package more flexible.

Eric Crandall 