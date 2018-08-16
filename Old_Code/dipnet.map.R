### MAPPING FUNCTION FOR METADATAQC ####

dipnet.map<-function(filename){

#set the working directory to the metadata maps directory
#setwd(path.expand(file.path("~","Google Drive","!DIPnet_DataQC","Metadata QC area","Maps")))

#Open the file 
dat<-read.table(file=filename,sep="\t",header=T,stringsAsFactors=F,na.strings=c("","NA","#N/A"),quote="",comment.char="",strip.white=T,fill=T)
#____________________________________________
#Manipulating the data ready for mapping
#____________________________________________


#retaining the name of the mdfasta file to be used in the output file name
dataset_name <- paste(strsplit(strsplit(filename,"/")[[1]][5],"\\.")[[1]][1], Sys.Date(), sep="_")

#add 360 to all negative values in the longitude column; this fits with the map
dat$decimalLongitude<-as.numeric(dat$decimalLongitude)
dat$decimalLongitude <- ifelse(dat$decimalLongitude < 0, dat$decimalLongitude + 360, dat$decimalLongitude)

#making unique IDs for sampled locations, NOTE: locations will be given unique numbers, but they don't mean anything
dat$sites <- as.numeric(interaction(dat$decimalLongitude, dat$decimalLatitude, drop = TRUE))

#____________________________________________
#Plotting maps for each species_locus_PI dataset presented on a separate page, where colours represent sample sizes at each location
#____________________________________________


#putting the relevant data into a list form; samples sizes per location
sample.size <- list()
  for(i in unique(dat$sites)){
    site_i <- dat[dat$sites == i, ]
        sample.size[[i]] <- data.frame(site = site_i$sites[1], long = site_i$decimalLongitude[1], lat = site_i$decimalLatitude[1], size = dim(site_i)[1])
      }

#turning the lists into a dataframe
sample.size.df <- data.frame(do.call("rbind", sample.size), row.names = NULL)
colnames(sample.size.df) <- c("site", "long", "lat", "size")

#defining our colour ramp for the map
ramp.cols<- c("steelblue2", "palegreen3", "yellow", "orange", "red","maroon")
ramp <- colorRampPalette(ramp.cols)(max(sample.size.df$size))

#making a column in our dataframe to hold a sample size per location in the dataset
sample.size.df$col <- NA
  for (i in 1:max(sample.size.df$size)){
    sample.size.df[sample.size.df$size == i, "col"] <- ramp[i]
  }

#set up to print the map to a pdf
pdf(paste(dataset_name,"pdf", sep = "."), paper = "a4r")

layout(matrix(c(1,2), ncol = 1, byrow = TRUE), heights = c(8,2))
  
par(mar = c(0,0,0,0))  
  
map("world2Hires", xlim = c(min(sample.size.df$long-5), max(sample.size.df$long+5)), ylim = c(min(sample.size.df$lat-5), max(sample.size.df$lat+5)), interior = FALSE, col = "grey50", lwd =1.5)  
  
points(sample.size.df$long, sample.size.df$lat, cex = 1, pch = 16, col = sample.size.df$col)

map.axes()

par(mar = c(2,0,0,0))

plot(c(1:max(sample.size.df$size + 1)), seq(0,1, length = max(sample.size.df$size) + 1), type="n",bty="n", xaxt="n", xlab="", yaxt="n", ylab="")
axis(1, at = seq(1, max(sample.size.df$size), length = 4), lab = round(seq(1, max(sample.size.df$size), length = 4)))
  for (j in 1:(length(ramp))) {
    rect(j, 0, j+1, 1, col=ramp[j], border=NA)
  }

text(4,1,"sample size")

#terminates the printing of the pdf into your working directory
dev.off()

#to write a tab delimited file that contains all unique locations within the dataset (named only by random numbers), the corresponding long/lat for that site, and the sample size at that site

#first change the longitudes back to between -180 to 180
sample.size.df$long <- ifelse(sample.size.df$long > 180, sample.size.df$long - 360, sample.size.df$long)

write.table(sample.size.df, file = paste(dataset_name,"txt", sep = "."),quote=F)

#Make sure .pdf and .csv file are in ~/Google Drive/!DIPnet_DataQC/Maps and email relevant PI or contributor

#END________________________________________
}

