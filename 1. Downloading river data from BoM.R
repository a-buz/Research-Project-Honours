#------------------------------------------------#
# Script to download river data files from the   #
# Australian Bureau of Meterology                #
#------------------------------------------------#

# WD set by R Project file

# Libraries
library(RCurl)

# Create rivers directory if it does not exist
if(dir.exists('Data/Rivers') == FALSE) {
  dir.create('Data/Rivers', recursive = TRUE, showWarnings = FALSE)
}

# Read in HRS data and ensure station number and name columns are characters
Data <- read.csv("Data/hrs_station_details.csv", skip = 11)
Data$AWRC.Station.Number <- as.character(Data$AWRC.Station.Number)
Data$Station.Name <- as.character(Data$Station.Name)

# Loop to download each file from BoM
for(i in 1:nrow(Data)) {
  filename <- paste("Data/Rivers/", Data$AWRC.Station.Number[i], ".csv", sep = "")
  
  if(file.exists(filename) == FALSE) { 
    URL <- paste0("http://www.bom.gov.au/water/hrs/content/data/", Data$AWRC.Station.Number[i],
                  "/", Data$AWRC.Station.Number[i], "_daily_ts.csv")

    download.file(URL, 
                  filename,
                  method="libcurl")
  }
}