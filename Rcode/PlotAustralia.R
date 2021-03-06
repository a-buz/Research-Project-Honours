
# directory
setwd("C:/Users/rver4657/ownCloud/working/brucetran/Research-Project-Honours")

# libraries
library(tidyverse)
library(maptools)
library(rgdal)
require(raster)

# read in the data from the hrs stations
stations <- read_csv("data/hrs_station_details.csv", skip = 11)
# read in the Ozdata from this page:http://www.elaliberte.info/code
ozdata <- read_csv("data/ozdata.csv")
# split by state, otherwise plots funny
ozdata_NSW <- ozdata %>% 
  filter(state == "NSW")
ozdata_VIC <- ozdata %>% 
  filter(state == "VIC")
ozdata_SA <- ozdata %>% 
  filter(state == "SA")
ozdata_QLD <- ozdata %>% 
  filter(state == "QLD")
ozdata_WA <- ozdata %>% 
  filter(state == "WA")
ozdata_TAS <- ozdata %>% 
  filter(state == "TAS")
ozdata_NT <- ozdata %>% 
  filter(state == "NT")

# Plot them all together
p <- ggplot(ozdata_NSW, aes(long, lat))
p <- p + geom_polygon(fill="white",colour="black")
p <- p + geom_polygon(data=ozdata_VIC, aes(long, lat),
                      fill="white",colour="black")
p <- p + geom_polygon(data=ozdata_QLD, aes(long, lat),
                      fill="white",colour="black")
p <- p + geom_polygon(data=ozdata_SA, aes(long, lat),
                      fill="white",colour="black")
p <- p + geom_polygon(data=ozdata_NT, aes(long, lat),
                      fill="white",colour="black")
p <- p + geom_polygon(data=ozdata_WA, aes(long, lat),
                      fill="white",colour="black")
p <- p + geom_polygon(data=ozdata_TAS, aes(long, lat),
                      fill="white",colour="black")
p <- p + geom_point(data =stations, aes(Longitude,Latitude), 
                    colour="blue", size=2)
p <- p + coord_equal()
p  <- p + xlab("Longitude") + ylab("Latitude") +
  theme(axis.title.x = element_text(face="bold",  size=16),
        axis.text.x  = element_text(size=12)) +
  theme(axis.title.y = element_text(face="bold",  size=16),
        axis.text.y  = element_text(size=12))
p

tiff("paper/Figure_HRSstations.tif",
     width=8*480, height=8*480,
     res=600, compression="lzw")
print(p)
dev.off()
