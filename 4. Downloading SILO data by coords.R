#------------------------------------------------#
# Script to download SILO data                   #
# This script uses the new SILO API              #
#------------------------------------------------#

# Create directories if they do not exist
if(dir.exists('Data/SILO') == FALSE) {
  dir.create('Data/SILO', recursive = TRUE)  
}

# Read in rivers file
Rivers <- read.csv('Data/hrs_station_details.csv', skip = 11, header = TRUE)
Rivers <- Rivers[order(Rivers$AWRC.Station.Number),]

# Set SILO API parameters
apiurl <- 'https://siloapi.longpaddock.qld.gov.au'
apiKey <- '' # Paste SILO API key here
dataFormat <- 'csv'
startDate <- '19900101'
endDate <- '20101231'
variables <- 'max_temp,min_temp,daily_rain,et_morton_potential'
# Only the following variables are available: daily_rain, max_temp, min_temp, evap_pan, vp, mslp, radiation, rh_tmax, rh_tmin, vp_deficit, evap_syn, evap_comb, evap_morton_lake, et_morton_actual, et_morton_potential, et_morton_wet, et_short_crop, et_tall_crop.'

for (i in 1:nrow(Rivers)) {
  downloadFile <- paste0('Data/SILO/', Rivers$AWRC.Station.Number[i], '.csv')
  if(file.exists(downloadFile) == FALSE) {
    URL <- paste0(apiurl, '/pointdata?',
                  'lat=', round(Rivers$Latitude[i]/0.05)*0.05,
                  '&lon=', round(Rivers$Longitude[i]/0.05)*0.05, 
                  '&apikey=', apiKey,
                  '&start=', startDate,
                  '&finish=', endDate,
                  '&format=', dataFormat,
                  '&variables=', variables)
  
    download.file(URL, destfile = downloadFile, method = 'libcurl')
  }
}
