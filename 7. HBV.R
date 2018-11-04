#--------------------------------------------------------------------------------------------------#
# HBV model run                                                                                 #
#--------------------------------------------------------------------------------------------------#

library(foreach)
library(doParallel)
library(hydromad)
library(Rcpp)

setwd('~/Dropbox (Sydney Uni)/BTPaper/Research-Project-Honours/')

# Set up monolithic function to simplify unattended runs
parallelRuns <- function(calStart, calEnd) {
  #----------------------------------------------------------------------------#
  # Loop preparation #
  #----------------------------------------------------------------------------#
  # Create KGE objective function
  KGE <- function(Q, X, ...) {
    1 - sqrt(
      (cor(X, Q) - 1)^2 +
        (mean(X)/mean(Q) - 1)^2 +
        (sd(X)/sd(Q) - 1)^2
    )
  }
  
  # River zoo files from previous step containing P, Q, E
  zooFiles <- dir('Data/Zoo', full.names = TRUE)
  # River details look up table
  riverTable <- read.csv("Data/FinalStations.csv")
  riverTable$AWRC.Station.Number <- as.character(riverTable$AWRC.Station.Number)
  riverTable$Station.Name <- as.character(riverTable$Station.Name)
  # Column names for variables collected in loop
  HBVColNames <-  c('Number', 'Name', 'Area', 'objFun', 'Message',
                    'ck2', 'ck1', 'ck0', 'maxbas',
                    'degd', 'degw', 'ttlim', 'perc',
                    'beta', 'lp', 'fcap', 'hl1',
                    'runoff', 'rbias', 'NSE', 'r.sq.sqrt',
                    'r.sq.log', 'KGE')
  
  runName <- paste0('Data/HBV/Fitted_',
                    substr(calStart, 1, 4),
                    '_', substr(calEnd, 1, 4))
  
  # Set up dir for output
  if(dir.exists(runName) == FALSE) {
    dir.create(runName, recursive = TRUE)
  }
  
  #---------------------------------------------------------------------------#
  # Parallel loop                                                             #
  #---------------------------------------------------------------------------#
  # Allocate cores. including logical cores
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)
  
  # Create log file to check progress
  writeLines(c(""), paste0(runName, "_log.txt"))
  
  # Start parloop and assign to collector
  r <- foreach(i=1:length(zooFiles), .combine=rbind, .packages=c('hydromad', 'Rcpp'),
               .verbose=TRUE) %dopar% {
    # Write to logfile
    sink(paste0(runName, "_log.txt"), append=TRUE)
    cat(paste(Sys.time(), "Starting iteration", i, "of", nrow(riverTable), "\n"))
    sink() # End diversion
    
    # Get river data
    riverDetails <- riverTable[riverTable$AWRC.Station.Number == basename(zooFiles[i]),]
    
    # Read in zoo file from csv. read.csv.zoo() will not work for some reason
    river <- read.zoo(zooFiles[i], header = TRUE)
    
    # Set cal period and create hydromad variable
    river.cal <- window(river, start = calStart, end = calEnd)
    
    source('HBV/HBV.R')
    
    # Pars based from http://dx.doi.org/10.1016/j.jhydrol.2014.04.037 
    CMod <- hydromad(DATA = river.cal,
                     sma = 'HBV',
                     routing = NULL,
                     ck2    = c(0.001, 0.15), # Routing: recession coefficient
                     ck1    = c(0.01, 0.4),   # Routing: recession coefficient
                     ck0    = c(0.05, 0.5),   # Routing: recession coefficient
                     maxbas = c(1, 7),        # Routing: weighting function length
                     degd   = c(0, 20),       # Snow: degree day factor
                     degw   = 1,              # Snow: snowmelt threshold
                     ttlim  = -1.41934,       # Snow: snowfall threshold
                     beta   = c(1, 6),        # Soil: shape coefficient
                     lp     = c(0.3, 1),      # Soil: SM threshold for evap reduction
                     fcap   = c(50, 2000),    # Soil: field capacity
                     perc   = c(0, 3),        # GW: max percolation
                     hl1    = c(10, 100),     # GW: quick runoff threshold
                     return_state = TRUE)
    
    # Fit KGE
    fitKGE <- fitByOptim(CMod,
                         objective = KGE,
                         samples = 500,
                         method = "PORT")
    
    resultKGE  <- c(riverDetails$AWRC.Station.Number,
                    riverDetails$Station.Name,
                    riverDetails$Catchment.Area..km2.,
                    'KGE',
                    fitKGE$fit.result$message,
                    fitKGE$parlist$ck2,
                    fitKGE$parlist$ck1,
                    fitKGE$parlist$ck0,
                    fitKGE$parlist$maxbas,
                    fitKGE$parlist$degd,
                    fitKGE$parlist$degw,
                    fitKGE$parlist$ttlim,
                    fitKGE$parlist$perc,
                    fitKGE$parlist$beta,
                    fitKGE$parlist$lp,
                    fitKGE$parlist$fcap,
                    fitKGE$parlist$hl1,
                    summary(fitKGE)$runoff,
                    summary(fitKGE)$rel.bias,
                    summary(fitKGE)$r.squared,
                    summary(fitKGE)$r.sq.sqrt,
                    summary(fitKGE)$r.sq.log,
                    objFunVal(fitKGE, KGE)
    )
    
    # Fit with NSE objective function
    fitNSE <- fitByOptim(CMod,
                         objective = hmadstat('r.squared'),
                         samples=500,
                         method="PORT")
    
    resultNSE  <- c(riverDetails$AWRC.Station.Number,
                    riverDetails$Station.Name,
                    riverDetails$Catchment.Area..km2.,
                    'NSE',
                    fitNSE$fit.result$message,
                    fitNSE$parlist$ck2,
                    fitNSE$parlist$ck1,
                    fitNSE$parlist$ck0,
                    fitNSE$parlist$maxbas,
                    fitNSE$parlist$degd,
                    fitNSE$parlist$degw,
                    fitNSE$parlist$ttlim,
                    fitNSE$parlist$perc,
                    fitNSE$parlist$beta,
                    fitNSE$parlist$lp,
                    fitNSE$parlist$fcap,
                    fitNSE$parlist$hl1,
                    summary(fitNSE)$runoff,
                    summary(fitNSE)$rel.bias,
                    summary(fitNSE)$r.squared,
                    summary(fitNSE)$r.sq.sqrt,
                    summary(fitNSE)$r.sq.log,
                    objFunVal(fitNSE, KGE)
    )
    
    # Fit with NSE objective function
    fitRSQLog <- fitByOptim(CMod,
                            objective = hmadstat('r.sq.log'),
                            samples = 500,
                            method = "PORT")
    
    resultRSQLog  <- c(riverDetails$AWRC.Station.Number,
                       riverDetails$Station.Name,
                       riverDetails$Catchment.Area..km2.,
                       'RSquaredLog',
                       fitRSQLog$fit.result$message,
                       fitRSQLog$parlist$ck2,
                       fitRSQLog$parlist$ck1,
                       fitRSQLog$parlist$ck0,
                       fitRSQLog$parlist$maxbas,
                       fitRSQLog$parlist$degd,
                       fitRSQLog$parlist$degw,
                       fitRSQLog$parlist$ttlim,
                       fitRSQLog$parlist$perc,
                       fitRSQLog$parlist$beta,
                       fitRSQLog$parlist$lp,
                       fitRSQLog$parlist$fcap,
                       fitRSQLog$parlist$hl1,
                       summary(fitRSQLog)$runoff,
                       summary(fitRSQLog)$rel.bias,
                       summary(fitRSQLog)$r.squared,
                       summary(fitRSQLog)$r.sq.sqrt,
                       summary(fitRSQLog)$r.sq.log,
                       objFunVal(fitRSQLog, KGE)
    )
    
    # Create list of fitted mode objects
    fitted <- list('KGE' = fitKGE,
                   'NSE' = fitNSE,
                   'RSquaredLog' = fitRSQLog)
    
    # Write
    save(fitted, file = paste(runName, '/', riverDetails$AWRC.Station.Number, ".Rdata", sep = ""))
    
    # Get summary results
    results <- t(matrix(c(resultKGE,resultNSE,resultRSQLog), ncol = 3))
    return(results)
  }
  
  # Stop clusters
  stopImplicitCluster()
  stopCluster(cl)
  
  HBVSummary <- as.data.frame(r, row.names = FALSE)
  colnames(HBVSummary) <- HBVColNames
  
  write.csv(HBVSummary, paste0(runName, '.csv'), row.names = FALSE)
}

# Runs
parallelRuns(calStart = '1990-01-01',
             calEnd = '1996-12-31')

parallelRuns(calStart = '2000-01-01',
             calEnd = '2006-12-31')
