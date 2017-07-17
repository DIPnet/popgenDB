library(gdistance)
library(ncdf4) #requires unix/linux install of netcdf
library(rasterVis)

#futz with the two rasters I downloaded and make a single merged raster, no need to repeat this step
mapTheme <- rasterTheme(region = rev(brewer.pal(10, "RdBu")))

IP_West<-raster("etopo_IP_West.nc") #read in the topographic raster
projection(IP_West)<-"+proj=longlat +datum=WGS84 +ellps=WGS84" #add a projection

IP_East<-raster("etopo_IP_East.nc") #read in the Eastern topographic raster
projection(IP_East)<-"+proj=longlat +datum=WGS84 +ellps=WGS84"

#trim IP_East east to just the eastern IP (it includes more than this), and shift it to positive degrees (+360)
extent2<-extent(IP_East)
extent2@xmax<--109
IP_East2<-crop(IP_East,extent2)
IP_East3<-shift(IP_East2,360)

#create the merged raster and save it for future use 
IPtopo<-raster::merge(IP_West,IP_East3,filename="IPtopo.raster")







# Prelude using already loaded ipdb to get to the same set of population 
# names that were measured for great-circle distance
########################################################################################################################
setwd("~/github/popgenDB/output")
esu_loci <- unique(ipdb$Genus_species_locus)


#read in the topo raster
IPtopo<-raster("IPtopo.grd")
#convert values greater than 0 to 0, values less than 0 to 1
IPtopo[IPtopo>0]<-0 #or the other way?
IPtopo[IPtopo<0]<-1
#add the sites on the map using a 1 for each site
IPtopo<-rasterize(sites2,IPtopo,1,fun="max",update=T)
#create a trarasnsition matrix
IPtrans<-transition(IPtopo,transitionFunction=min,directions=8)
IPtrans.c<-geoCorrection(IPtrans,type="c")


for(gsl in esu_loci){ #gsl<-"Linckia_laevigata_CO1" "Tridacna_crocea_CO1" "Lutjanus_kasmira_CYB" "Acanthaster_planci_CO1"
  
  cat("\n","\n","\n","Now starting", gsl, "\n")
  
  if(any(is.na(diffstats[[gsl]]))){cat("NAs in FST table, No gdm calculated"); next}
  
  if(diffstats[[gsl]]=="Fewer than 3 sampled populations after filtering. No stats calculated"){nostats<-c(nostats,gsl);next}
  
  #pull out the data for this genus-species-locus (gsl)
  sp<-ipdb[which(ipdb$Genus_species_locus==gsl),]
  #clean weird backslashes from names
  sp$locality<-gsub("\"","",sp$locality)
  
  sp$sample<-paste(sp$locality,round(sp$decimalLatitude, digits=0),round(sp$decimalLongitude, digits=0),sep="_")  #sets up a variable that matches the name in Fst table
  sp<-sp[order(sp$sample),]
  # Not all localities are included in Veron's regionalization (e.g. Guam), so get their names and then zap NAs
  nonVeronpops<-unique(sp$sample[is.na(sp$VeronDivis)])
  sp<-sp[!is.na(sp$VeronDivis),]
  
  #subsample Fst 
  gslFST<-diffstats[[gsl]]
  #make a matrix out of gslFST
  gslFSTm<-as.matrix(gslFST)
  
  gslFSTm[which(gslFSTm > 1)] <- 1 #some values that look like 1.0000 are registering as greater than 1
  #gslFSTm[which(gslFSTm < 0)] <- 0 #get rid of artifactual negative values
  gslFSTm<-rescale(gslFSTm)
  #gslFSTm<-gslFSTm/(1-gslFSTm)
  
  #zap weird slashes in the names
  rownames(gslFSTm)<-gsub("\"","",rownames(gslFSTm))
  colnames(gslFSTm)<-rownames(gslFSTm)
  
  #zap the same na populations from the list of non existent pops from VeronDivis
  if(any(rownames(gslFSTm) %in% nonVeronpops)){
    gslFSTm<-gslFSTm[-(which(rownames(gslFSTm) %in% nonVeronpops)),-(which(colnames(gslFSTm) %in% nonVeronpops))]
  }
  
  if(length(rownames(gslFSTm))<5){nostats<-c(nostats,gsl);cat("Fewer than 5 sampled populations");next}
  
  #and filter sp based on the localities that have Fst values
  sp<-sp[sp$sample %in% rownames(gslFSTm),]
  
  #and vice versa
  
  gslFSTm<- gslFSTm[which(rownames(gslFSTm) %in% unique(sp$sample)),which(rownames(gslFSTm) %in% unique(sp$sample))]
  
  
  
  #create a locations data frame that has all the localities plus lats and longs and their Veron region.
  locs<-as.data.frame(unique(sp$sample))
  names(locs)<-"sample"
  #locs$Long<-sp$decimalLongitude[which(locs %in% sp$sample)]
  #can't do a unique on sample, lats and longs because some samples have non-unique lats and longs! So I do a join and take the first match.
  locs<-join(locs,sp[c("sample","decimalLongitude","decimalLatitude"
                       ,"VeronDivis")], by="sample", match="first")
  
  #sort gslFSTm
  gslFSTm<-gslFSTm[order(rownames(gslFSTm)),order(colnames(gslFSTm))]
  # convert to data frame with popsample names as first column
  gslFSTm<-cbind(sample=locs$sample,as.data.frame(gslFSTm))
  
  ######################################################################
  # 3. Calculate Great Circle Distance
  gcdist_km <- pointDistance(locs[,2:3],lonlat=T)/1000
  #symmetricize the matrix
  gcdist_km[upper.tri(gcdist_km)]<-t(gcdist_km)[upper.tri(gcdist_km)]
#Calculate Overwater Distance  
########################################################################################################
  
  #correct negative west-latitude values to positive values
  locs$decimalLongitude[which(locs$decimalLongitude<0)]<-locs$decimalLongitude[which(locs$decimalLongitude<0)] + 360
  
  #create various formats for the sites data
  sites1<-as.data.frame(c(locs[,c(2,3)])) #take the lat longs only
  sites2<-SpatialPoints(sites1) #transform to spatial points
  projection(sites2)<-"+proj=longlat +datum=WGS84 +ellps=WGS84" #add a projection
  sites3<-as.matrix(sites1) #convert to a matrix
  
  #Create a source raster at the same res as your land/water mask, where the GPS point is the uniqueID, everything else is zero 
  gps_src<-rasterize(sites2,IPtopo,1:length(locs$sample),fun="max",background=0) 

  
  #crop the topographic raster to appropriate size, and convert values greater than 0 to 0 and values less than zero to 1
  #IPtopo2<-crop(IPtopo,extent(sites2)) 
  
  
  #NAvalue(IPtopo2)<-0


test<-costDistance(IPtrans.c,sites2,sites2)
  