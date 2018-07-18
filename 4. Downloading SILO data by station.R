#------------------------------------------------#
# Script to download SILO data                   #
# This script uses the new SILO API              #
#------------------------------------------------#

library(dplyr)

# Create directories if they do not exist
if(dir.exists('Data/SILO/Patched') == FALSE) {
  dir.create('Data/SILO/Patched', recursive = TRUE)  
}

# Read in final stations rivers file
rivers <- read_csv('Data/FinalStations.csv')

# Set SILO API parameters
apiurl <- 'https://siloapi.longpaddock.qld.gov.au'
apiKey <- 'Akhi03HxgFWclJXwTnzWU4V9oZBf0oO3G2z1F6vL' # Paste SILO API key here
dataFormat <- 'csv'
startDate <- '19900101'
endDate <- '20101231'
variables <- 'max_temp,min_temp,daily_rain,et_morton_potential'

#Only the following variables are available: daily_rain, max_temp, min_temp, evap_pan, vp, mslp, 
# radiation, rh_tmax, rh_tmin, vp_deficit, evap_syn, evap_comb, evap_morton_lake, et_morton_actual,
# et_morton_potential, et_morton_wet, et_short_crop, et_tall_crop.'

for (i in 1:nrow(rivers)) {
  downloadFile <- paste0('Data/SILO/Patched/', rivers$`AWRC Station Number`[i], '.csv')
  if(file.exists(downloadFile) == FALSE) {
    URL <- paste0(apiurl, '/pointdata?',
                  'station=', rivers$`SILO Number`[i],
                  '&apikey=', apiKey,
                  '&start=', startDate,
                  '&finish=', endDate,
                  '&format=', dataFormat,
                  '&variables=', variables)
  
    download.file(URL, destfile = downloadFile, method = 'libcurl')
  }
}
