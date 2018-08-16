### TURNING OBIS OUTPUT INTO MAPS####
Libby Liggins, 16 June 2015

#____________________________________________
#Accessing packages you will need
#____________________________________________

#the first time you use these packages you must download them from the library
install.packages("maps", dependencies = TRUE)
install.packages("mapdata", dependencies = TRUE)
install.packages("mapproj", dependencies = TRUE)

#after the first time, you will just need to call them from your local repository
require(maps)
require(mapproj)
require(mapdata)

#set the working directory to where the spreadsheet of OBIS data is

#Open the file, substituting "OBIS_Mulloidichthys_vanicolensis.csv" for your file name
OBIS_data<-read.csv(file="OBIS_Mulloidichthys_vanicolensis.csv",header=T)
#exclude any rows that have na for longitude or latitude
OBIS_data <- OBIS_data[!is.na(OBIS_data$longitude),]
OBIS_data <- OBIS_data[!is.na(OBIS_data$latitude),]

#____________________________________________
#Manipulating the data ready for mapping
#____________________________________________

#add 360 to all negative values in the longitude column; this fits with the map
OBIS_data$longitude <- as.numeric(OBIS_data$longitude)
OBIS_data$longitude <- ifelse(OBIS_data$longitude < 0, OBIS_data$longitude + 360, OBIS_data$longitude)

#____________________________________________
#Plotting a map for the OBIS dataset
#____________________________________________

#set up to print the map to a pdf, this map will be placed where your working directory is (i.e. where the dataset that you are using is placed). Substitute "OBIS_Mulloidichthys_vanicolensis" for an appropriate title
pdf(paste("OBIS_Mulloidichthys_vanicolensis","pdf", sep = "."), paper = "a4r")

layout(matrix(c(1,2), ncol = 1, byrow = TRUE), heights = c(8,2))
  
par(mar = c(0,0,0,0))  
  
map("world2Hires", xlim = c(min(OBIS_data$longitude)-5, max(OBIS_data$longitude)+5), ylim = c(min(OBIS_data$latitude)-5, max(OBIS_data$latitude)+5), interior = FALSE, col = "grey50", lwd =1.5)  
  
points(OBIS_data$longitude, OBIS_data$latitude, cex = 0.5, pch = 16, col = "red")

map.axes()

par(mar = c(2,0,0,0))

#terminates the printing of the pdf into your working directory (make sure you have any previsou R generated versions of this map closed)
dev.off()

#END________________________________________


