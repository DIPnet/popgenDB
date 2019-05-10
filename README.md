Population Genetic Statistics for the DIPnet dataset
=====================================================
Added April 20, 2015 by Eric Crandall
Updated May 10, 2019 by Eric Crandall

These scripts are intended for use with the DIPnet dataset that is now published in [Global Biogeography and Ecology](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12905). Much of the R workflow for that paper can be viewed [here](http://htmlpreview.github.io/?https://github.com/DIPnet/popgenDB/blob/master/analysis_notebook.nb.html)

The way it currently works is that all the functions are stored in `DIPnet_Stats_Functions.R` and a script for driving those functions is in `DIPnet_Stats_Script.R`. The scripts can be used on any data frame that has the following fields:

1. Genus_species_locus
2. materialSampleID (a unique identifier)
3. locality (just the name)
4. decimalLatitude
5. decimalLongitude
6. sequence (a haploid mtDNA sequence)

You'll want to change the paths in `config_example.R` to values that point to local versions of the DIPnet DB and spatial DB and then save that file as `config.R`. Similarly, you can create an empty directory "/output" save whatever results you get to it. These files won't be synced (pushed or pulled) with GitHub because they are included in the .gitignore file.

There are major functions for looping over the whole database and calculating:

1. Genetic diversity statistics (`genetic.diversity.mtDNA.db()`)
	1. Unique Haplotype Number
	2. Haplotype Diversity
	3. Shannon-Wiener Diversity
	4. Effective number of haplotypes
	5. Local FST
	6. Nucleotide Diversity (per site and per locus)
	7. Watterson's ThetaS
	8. Tajima's D
	9. Coverage (Chao and Jost 2012)
	10. Coverage standardized haplotype diversity

2. Genetic structure statistics (`pairwise.structure.mtDNA.db()`)
	1. Nei's Fst
	2. Nei's Gst
	3. Hedrick's G'st
	4. Jost's D
	5. Weir & Cockerham's Theta
	6. Excoffier et al. PhiST
	7. Chi-Square Test of Population Differentiation
	8. Nei's Da
	
3. Hierarchical Structure using AMOVA (up to 3 levels currently supported; `hierarchical.structure.mtDNA.db()`) 



May want to eventually break these functions down into component parts (the looping and stats individually for example) to make the package more flexible.

