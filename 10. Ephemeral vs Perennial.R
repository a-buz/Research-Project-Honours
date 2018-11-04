#
# Ephemeral vs perennial
# This output will be used in script 9. Only using 2000 classification now to make comparisons easy
# 
library(tidyverse)
library(hydromad)

# Stations
finalStations <- read_csv('Data/FinalStations.csv')

riverType <- finalStations %>% 
  select(`AWRC Station Number`, `Station Name`, `Catchment Area (km2)`) %>% 
  mutate(`Type` = NA)

for (i in 1:nrow(riverType)) {
  river <- read.zoo(paste0('Data/Zoo/', riverType$`AWRC Station Number`[i]),
                    header = TRUE)
  #Getting right dates for the observed flow
  obs <- window(river, start = "1990-04-10", end = "1996-12-31")
  
  proportion <- length(which(obs$Q == 0))/length(obs$Q)*100
  
  riverType$Type[i] <- ifelse(proportion <= 1 , "Perennial", ifelse(proportion > 1, "Ephemeral"))
}

write_csv(riverType, 'Data/RiverType1990-1996.csv')

riverType <- finalStations %>% 
  select(`AWRC Station Number`, `Station Name`, `Catchment Area (km2)`) %>% 
  mutate(`Type` = NA)

for (i in 1:nrow(riverType)) {
  river <- read.zoo(paste0('Data/Zoo/', riverType$`AWRC Station Number`[i]),
                    header = TRUE)
  #Getting right dates for the observed flow
  obs <- window(river, start = "2000-04-10", end = "2006-12-31")
  
  proportion <- length(which(obs$Q == 0))/length(obs$Q)*100
  
  riverType$Type[i] <- ifelse(proportion <= 1 , "Perennial", ifelse(proportion > 1, "Ephemeral"))
}

write_csv(riverType, 'Data/RiverType2000-2006.csv')
