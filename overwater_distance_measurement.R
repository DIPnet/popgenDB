#1) go to http://maps.ngdc.noaa.gov/ to get ETOPO1 data, download relevant raster
# or http://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.html
library(gdistance)
library(ncdf4) #requires unix/linux install of netcdf
library(rasterVis)

setwd("~/github/popgenDB/output")

sites<- # read in the sites

mapTheme <- rasterTheme(region = rev(brewer.pal(10, "RdBu")))
  
IP_West<-raster("etopo_IP_West.nc") #read in the topographic raster
projection(IP_West)<-"+proj=longlat +datum=WGS84 +ellps=WGS84" #add a projection


#create various formats for the sites data
sites1<-as.data.frame(c(sites[,c(3,2)])) #take the lat longs only
sites2<-SpatialPoints(sites1) #transform to spatial points
projection(sites2)<-"+proj=longlat +datum=WGS84 +ellps=WGS84" #add a projection
sites3<-as.matrix(sites1) #convert to a matrix


#Create a source raster at the same res as your land/water mask, where the GPS point is the uniqueID, everything else is zero 
gps_src<-rasterize(sites2,impatopo,1:28,fun="max",background=0) 
plot(gps_src,xlim=c(100,101),ylim=c(-2,-1))

		
#measure Euclidean distance
eucdist<-pointDistance(sites2,longlat=T)/1000


#crop the topographic raster to appropriate size, and convert values greater than 0 to 0 and values less than zero to 1
impatopo2<-crop(impatopo,c(94,135,-10,17)) # or impatopo2<-crop(impatopo,extent(sites2))

NAvalue(impatopo2)<-0
impatopo2[impatopo2>0]<-0
impatopo2[impatopo2<0]<-1
impatopo2[is.na(impatopo2)]

#add the sites on the map using a 1 for each site
impatopo2<-rasterize(sites2[28,],impatopo2,1,fun="max",update=T)

#plot to check, make sure that all sites are bordering (within 1 pixel) on water!!!
plot(impatopo2)
points(sites2)
text(sites2)

#create a transition matrix with 8-way movements (Queenlike), correct for different distances across the diagonals. This may take awhile (5-10 mins)
tr12<-transition(impatopo2,transitionFunction=min,directions=8)
trgc<-geoCorrection(tr12,type="c")

#Loop through each GPS point, add that point to the converted raster and creating a cost-distance raster with accCost (overwater distances from the point of interest to every point on the raster). Use assign() to give each cost-distance matrix a name. This may take quite awhile (2-3 hours?).
for(i in 1:28){
	assign(paste("gps_src_",i,sep=""),rasterize(sites1[i,],impatopo2,i,background=0))
	assign(paste("gps_dist_",i,sep=""),accCost(trgc,sites3[i,]))
}

#create an empty matrix to populate with values
ow<-matrix(nrow=28,ncol=28)

#loop through each source and destination, using zonal statistics to measure the minimum overwater distance from the destination to the source, put these together in a distance matrix. Use get to call each cost-distance raster by name. This also may take quite awhile (~ 1 hour)
for(src in 1:28){
	for(dest in 1:28){
		z<-zonal(get(paste("gps_dist_",src,sep="")),get(paste("gps_src_",dest,sep="")),min)
		ow[src,dest]<-z[2,2]
	}
}






