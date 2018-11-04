# Analysis of Results
library(tidyverse)
library(hydromad)

# Final stations
finalStations <- read_csv('Data/FinalStations.csv')
finalStations <- finalStations %>% dplyr::rename(Number = `AWRC Station Number`)

# River type to merge into rows
riverType1990 <- read_csv('Data/RiverType1990-1996.csv')
riverType1990 <- riverType1990 %>% 
  dplyr::rename(Number = `AWRC Station Number`)
riverType2000 <- read_csv('Data/RiverType2000-2006.csv')
riverType2000 <- riverType2000 %>% 
  dplyr::rename(Number = `AWRC Station Number`)

# Koppengeiger numbers to add in
kgRaster <- raster::raster('Data/Aus_KG_1986-2010_5m.grd')
# levels(kgRaster) to see attributes
kgAttr <- raster::levels(kgRaster)[[1]]
# Create points from final stations tables
stationpts <- sp::SpatialPoints(cbind(finalStations$Longitude,
                                      finalStations$Latitude),
                                proj4string = raster::crs(kgRaster))
finalStations$KGext <- kgAttr$climate[raster::extract(kgRaster, stationpts)]


#-------------------------------------------------------------------------------
# Sim

simSummary <- function(model, timePeriod) {
  # Read in summary data
  filePath <- paste('Data/', model, '/Fitted_', timePeriod, '.csv', sep='')
  modelSummary <- read_csv(filePath)
  
  # Remove false convergences and non-convergences
  falseConverge <- 'false convergence (8)'
  evalLimit <- 'function evaluation limit reached without convergence (9)'
  
  # Get converged
  converged <- modelSummary %>% 
    filter(!Message %in% c(falseConverge, evalLimit),
           !Number == 'A0020101') # Remove due to large area
  
  sim.summary <- converged %>% mutate(Cum.flow  = rep(NA, length(converged$Number)),
    FDC.low   = rep(NA, length(converged$Number)),
    FDC.mid   = rep(NA, length(converged$Number)),
    FDC.high  = rep(NA, length(converged$Number)),
    AC        = rep(NA, length(converged$Number)),
    Peaks     = rep(NA, length(converged$Number)),
    Slope     = rep(NA, length(converged$Number)))
  
  sim.summary <- inner_join(sim.summary, 
                            #get(paste0('riverType', substr(timePeriod, 1, 4))) %>% select(Number, Type),
                            riverType2000 %>% select(Number, Type), # Change to 2000 as it shows really if it is peren or eph
                            by = 'Number')
  
  # Merge in KGext also
  sim.summary <- inner_join(sim.summary,
                            finalStations %>% select(Number, KGext),
                            by = 'Number')
  
  # Simulated
  for(i in 1:nrow(sim.summary)) {
    # Record number and objfun
    stationNumber <- sim.summary$Number[i]
    objFun <- sim.summary$objFun[i]
    
    # Read in Rdata file
    load(paste('Data/', model, '/Fitted_', timePeriod, '/', stationNumber, '.Rdata', sep=''))
    
    if(objFun == 'KGE') {
      sim <- fitted$KGE
    } else if(objFun == 'NSE') {
      sim <- fitted$NSE
    } else if(objFun == 'RSquaredLog') {
      sim <- fitted$RSquaredLog
    }
    
    #######cum.flow
    cum.flow <- sum(fitted(sim))
    
    ##FDC.low - Performance metroc: F
    FDC.low <- unname(quantile(fitted(sim), probs = 0.1)) #/quantile(obs$Q, probs = 0.25)))
    
    ####FDC.High - Performance metroc: F
    FDC.high <- unname(quantile(fitted(sim), probs = 0.9))#/quantile(obs$Q, probs = 0.75)))
    
    ######FDC.Mid - Performance metroc: F
    FDC.mid <- unname(quantile(fitted(sim), probs = 0.5))#/quantile(obs$Q, probs = 0.5)))
    
    ######AC - Performance metroc: F
    AC.sim <- acf(fitted(sim),  plot = F)$acf[2]
    
    ##Peak Distribution - Perforamnce Metric: F
    Peak.sim <- (quantile(fitted(sim), probs = 0.9)-(quantile(fitted(sim), probs = 0.5)))/(0.9-0.5)
    
    #slope mid FDC
    slope.sim <- unname(diff(quantile(fitted(sim), probs = c(0.25, 0.75))/0.75-0.25))
    
    ##Putting into d.frame
    #sig.summary[i,2] <- signif(FDC.total, 4)
    sim.summary$Cum.flow[i] <- cum.flow
    sim.summary$FDC.low[i] <- signif(FDC.low, 4)
    sim.summary$FDC.mid[i] <- signif(FDC.mid, 4)
    sim.summary$FDC.high[i] <- signif(FDC.high, 4)
    sim.summary$AC[i] <- signif(AC.sim, 4)
    sim.summary$Peaks[i] <- signif(Peak.sim, 4)
    sim.summary$Slope[i] <- signif(slope.sim, 4)
  }
  
  #Change Inf's to NA
  sim.summary <- do.call(data.frame, lapply(sim.summary, function(x) replace(x, is.infinite(x),NA)))
  
  #Correlation matrix
  modelCor <- cor(sim.summary[,c(17:22)], use = "pairwise.complete.obs", method = "spearman")
  
  outPath <- paste('Data/Analysis/', timePeriod, '/', model, 'Summary.csv', sep = '')
  write_csv(sim.summary, outPath)
  corPath <- paste('Data/Analysis/', timePeriod, '/', model, 'Cor.txt', sep = '')
  write.table(modelCor, corPath)
}

simSummary('GR4J', '1990_1996')
simSummary('GR4J', '2000_2006')

simSummary('SIMHYD', '1990_1996')
simSummary('SIMHYD', '2000_2006')

simSummary('HBV', '1990_1996')
simSummary('HBV', '2000_2006')


# 
# 
# #-------------------------------------------------------------------------------
# # Obs
# # Flow duration Curve Analysis
# obs.summary <- tibble(Number    = finalStations$`AWRC Station Number`,
#                       cum.flow  = rep(NA, nrow(finalStations)),
#                       FDC.low   = rep(NA, nrow(finalStations)),
#                       FDC.mid   = rep(NA, nrow(finalStations)),
#                       FDC.high  = rep(NA, nrow(finalStations)),
#                       AC        = rep(NA, nrow(finalStations)),
#                       Peaks     = rep(NA, nrow(finalStations)),
#                       slope     = rep(NA, nrow(finalStations)))
# 
# 
# for (i in 1:nrow(finalStations)) {
#   # Read in zoo file
#   obs <- read.zoo(paste0('Data/Zoo/', finalStations$`AWRC Station Number`[i]), header = TRUE)
#   obs <- window(obs, start = "1990-04-11", end = "1996-12-31")
#   
#   #######cum.flow
#   cum.flow <- sum(obs$Q)
#   
#   ##FDC.low - Performance metroc: F
#   FDC.low <- unname(quantile(obs$Q, probs = 0.1))
#   
#   ####FDC.High - Performance metroc: F
#   FDC.high <- unname(quantile(obs$Q, probs = 0.99))
#   
#   ######FDC.Mid - Performance metroc: F
#   FDC.mid <- unname(quantile(obs$Q, probs = 0.5))
#   
#   ######AC - Performance metric: F
#   AC.obs <- acf(obs$Q,  plot = F)$acf[2]
#   AC <- AC.obs
#   
#   ##Peak Distribution - Perforamnce Metric: F
#   Peak.obs <- (quantile(obs$Q, probs = 0.9)-(quantile(obs$Q, probs = 0.5)))/(0.9-0.5)
#   Peaks <- Peak.obs
#   
#   #slope mid FDC
#   slope.obs <- unname(diff((quantile(obs$Q, probs = c(0.33, 0.66))/0.66-0.33)))
#   slope <- slope.obs
#   
#   ##Putting into d.frame
#   obs.summary[i,1] <- finalStations$`AWRC Station Number`[i]
#   obs.summary[i,2] <- cum.flow
#   obs.summary[i,3] <- signif(FDC.low, 4)
#   obs.summary[i,4] <- signif(FDC.mid, 4)
#   obs.summary[i,5] <- signif(FDC.high, 4)
#   obs.summary[i,6] <- signif(AC, 4)
#   obs.summary[i,7] <- signif(Peaks, 4)
#   obs.summary[i,8] <- signif(slope, 4)
# }
# 
# ObsCor <- cor(obs.summary[,c(2:8)], use = "pairwise.complete.obs", method = "spearman")
# 
# write_csv(obs.summary, 'Data/Analysis/1990_1996/ObsSummary.csv')
# write.table(ObsCor, 'Data/Analysis/1990_1996/ObsCor.txt')
