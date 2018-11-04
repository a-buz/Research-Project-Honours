# 
# Make zoo files
#

library(tidyverse)
library(zoo)

# Read in final stations rivers file
rivers <- read_csv('Data/FinalStations.csv')

# Start loop
for(i in 1:nrow(rivers)) {
  # BoM flow data
  riverFileName <- paste0('Data/Rivers/', rivers$`AWRC Station Number`[i], '.csv')
  river <- read_csv(riverFileName, skip = 26, col_types = list(col_date(), col_double(), col_character()))
  colnames(river) <- c('Date', 'Q', 'QCode')
  
  # Subset river data to within 1990 - 2010 time period
  river <- river %>% 
    select(Date, Q) %>% 
    filter(Date >= '1990-01-01' & Date <= '2010-12-31')
  
  # Get catchment area and convert from ML to mm
  # Read in from header
  areaHeader <- read_csv(riverFileName, col_names = F, skip = 15, n_max = 2)
  
  if(is.double(areaHeader$X3[1]) == FALSE) {
    break # Break if there is an issue with metadata
  }
  
  riverArea <- areaHeader$X3[1]
  print(areaHeader)
  # Convert
  river$Q <- river$Q/riverArea
  
  # SILO data
  silo <- read_csv(paste0('Data/SILO/Patched/', basename(riverFileName)))
  silo$mean_temp <- (silo$max_temp + silo$min_temp)/2
  
  # Select Date, P, Morton pET
  silo <- silo %>%
    select(date, daily_rain, et_morton_potential, mean_temp)
  colnames(silo) <- c('Date', 'P', 'E', 'T')
  
  # Merge tables
  joinedTBL <- left_join(river, silo, by = 'Date') %>% 
    select(Date, P, Q, E, T)
  
  # Create zoo object
  riverZoo <- zoo(joinedTBL[,2:5], order.by = joinedTBL$Date)
  # Create dir to write
  if(dir.exists('Data/Zoo') == FALSE) {
    dir.create('Data/Zoo', recursive = TRUE, showWarnings = FALSE)
  }
  # Write to disk
  write.zoo(riverZoo,
            file = paste0('Data/Zoo/', strsplit(basename(riverFileName), '.csv')[[1]]))
}
