#------------------------------------------------#
# Script 3
# Determine closest SILO stations
#------------------------------------------------#

# Note this script currently does not assess for quality of SILO station data

library(tidyverse)

# Read river data
rivers <- read_csv("Data/hrs_station_details.csv", skip=11)

# Merge river quality column
riverQuality <- read_csv("Data/CatchmentQuality.csv")
colnames(riverQuality) <- c('AWRC Station Number',
                            'Station Name',
                            'Flow Accuracy')

# Drop station name column
riverQuality <- riverQuality %>% select(-`Station Name`)

# Merge by station number and flow accuracy cut off
riverMerge <- inner_join(rivers, riverQuality, by = 'AWRC Station Number')

# Read in SILO station data
siloLoc <- read_csv('Data/SILOStations.csv')
# Filter where data was being collected until the end of the study period
siloLoc <- siloLoc %>% 
  filter(start_year < 1990) %>% 
  filter(end_year > 2010 | is.na(end_year) == TRUE)

# Find nearest silo station by minimum distance
library(rgeos)
library(sp)

# Create spatial point tables
riversSP <- SpatialPoints(riverMerge[,4:3])
siloSP <- SpatialPoints(siloLoc[,7:6])

# First insert row index
siloIndex <- apply(gDistance(siloSP, riversSP, byid=TRUE), 1, which.min)
# Filter by index using slice
siloStations <- siloLoc %>% 
  slice(siloIndex)

# Rename columns
colnames(siloStations) <- c('SILO Name',
                            'SILO Number',
                            'Start Year',
                            'End Year',
                            'Supplier',
                            'Station Latitude',
                            'Station Longitude',
                            'Station Elevation')

# Bind together
finalData <- bind_cols(riverMerge, siloStations)
# Write
write_csv(finalData, 'Data/FinalStations.csv')
