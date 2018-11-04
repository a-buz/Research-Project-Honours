#----------------------------------------------------------------------------
# SIMHYD model run
#----------------------------------------------------------------------------

library(foreach)
library(doParallel)
library(hydromad)
library(Rcpp)

setwd('~/Dropbox (Sydney Uni)/BTPaper/Research-Project-Honours/')

# Set up monolithic function to simplify unattended runs
parallelRuns <- function(calStart, calEnd) {
  #----------------------------------------------------------------------------
  # Loop preparation
  #----------------------------------------------------------------------------
  # Create KGE objective function
  KGE <- function(Q, X, ...) {
    1 - sqrt(
      (cor(X, Q) - 1)^2 +
        (mean(X)/mean(Q) - 1)^2 +
        (sd(X)/sd(Q) - 1)^2
    )
  }
  
  # Source SIMHYD scripts
  # rcode_dir <- 'SIMHYDCode'
  # source(paste(rcode_dir, 'Simhyd.R', sep = '/'))
  
  # River zoo files from previous step containing P, Q, E
  zooFiles <- dir('Data/Zoo', full.names = TRUE)
  # River details look up table
  riverTable <- read.csv("Data/FinalStations.csv")
  riverTable$AWRC.Station.Number <- as.character(riverTable$AWRC.Station.Number)
  riverTable$Station.Name <- as.character(riverTable$Station.Name)
  # Column names for variables collected in loop
  SIMHYDColNames <-  c('Number', 'Name', 'Area', 'objFun', 'Message',
                       'INSC', 'COEFF', 'SQ', 'SMSC', 'SUB', 'CRAK',
                       'K', 'runoff', 'rbias', 'NSE', 'r.sq.sqrt',
                       'r.sq.log', 'KGE')
  
  runName <- paste0('Data/SIMHYD/Fitted_',
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
    river <- read.table(zooFiles[i], header = TRUE)
    river$Index <- as.Date(river$Index)
    river <- zoo(river[,2:4], order.by = river$Index)
    
    # Set cal period and create hydromad variable
    river.cal <- window(river, start = calStart, end = calEnd)
    
    source('SIMHYD/SIMHYD_C2009.R')
    
    # Pars based from here https://journals.ametsoc.org/doi/pdf/10.1175/2009JHM1061.1
    CMod <- hydromad(DATA = river.cal,
                     sma = "simhyd_C2009",
                     routing = NULL,
                     INSC  = c(0.5, 5),     # Interception store capacity
                     COEFF = c(50, 500),    # Max infiltration loss
                     SQ    = c(0, 6),       # Infiltration loss exponent
                     SMSC  = c(50, 2000),    # Soil moisture store capacity
                     SUB   = c(0, 1),       # Interflow equation constant
                     CRAK  = c(0, 1),       # GW recharge constant
                     K     = c(0.003, 0.3), # Baseflow linear recession par
                     etmult = 1, # Supplied potential ET
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
                    #fitKGE$parlist$DELAY,
                    #fitKGE$parlist$X_m,
                    fitKGE$parlist$INSC,
                    fitKGE$parlist$COEFF,
                    fitKGE$parlist$SQ,
                    fitKGE$parlist$SMSC,
                    fitKGE$parlist$SUB,
                    fitKGE$parlist$CRAK,
                    fitKGE$parlist$K,
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
                    #fitNSE$parlist$DELAY,
                    #fitNSE$parlist$X_m,
                    fitNSE$parlist$INSC,
                    fitNSE$parlist$COEFF,
                    fitNSE$parlist$SQ,
                    fitNSE$parlist$SMSC,
                    fitNSE$parlist$SUB,
                    fitNSE$parlist$CRAK,
                    fitNSE$parlist$K,
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
                       #fitRSQLog$parlist$DELAY,
                       #fitRSQLog$parlist$X_m,
                       fitRSQLog$parlist$INSC,
                       fitRSQLog$parlist$COEFF,
                       fitRSQLog$parlist$SQ,
                       fitRSQLog$parlist$SMSC,
                       fitRSQLog$parlist$SUB,
                       fitRSQLog$parlist$CRAK,
                       fitRSQLog$parlist$K,
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
  
  SIMHYDSummary <- as.data.frame(r, row.names = FALSE)
  colnames(SIMHYDSummary) <- SIMHYDColNames
  
  write.csv(SIMHYDSummary, paste0(runName, '.csv'), row.names = FALSE)
}

# Runs
parallelRuns(calStart = '1990-01-01',
             calEnd = '1996-12-31')

parallelRuns(calStart = '2000-01-01',
             calEnd = '2006-12-31')
