#--------------------------------------------------------------------------------------------------#
# GR4J model run                                                                                   #
#--------------------------------------------------------------------------------------------------#

library(foreach)
library(doParallel)
library(hydromad)

setwd('~/Dropbox (Sydney Uni)/BTPaper/Research-Project-Honours/')

# Set up monolithic function to simplify unattended runs
parallelRuns <- function(calStart, calEnd, x2) {
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
  GR4JColNames <-  c('Number', 
                     'Name',
                     'Area',
                     'objFun',
                     'Message',
                     'x1',
                     'x2',
                     'x3',
                     'x4',
                     'runoff',
                     'rbias',
                     'NSE',
                     'r.sq.sqrt',
                     'r.sq.log',
                     'KGE')
  # Non changing
  x1 <- c(50,2000)
  x3 <- c(5,500)
  x4 <- c(0.5,10)
  
  # Bruce's naming convention
  if(min(x2) == -10 & max(x2) == 10) {
    suffix <- '3' 
  } else if(min(x2) == -30 & max(x2) == 20) {
    suffix <- '2'
  }
  
  runName <- paste0('Data/GR4J/Fitted_', substr(calStart, 1, 4),
                    '_', substr(calEnd, 1, 4), '_', suffix)
  
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
  r <- foreach(i=1:nrow(riverTable), .combine=rbind, .packages=c('hydromad')) %dopar% {
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
    CMod <- hydromad(river.cal,
                     sma="gr4j",
                     routing="gr4jrouting",
                     etmult=0.15,
                     x1 = x1,
                     x2 = x2,
                     x3 = x3,
                     x4 = x4)
    
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
                    fitKGE$parlist$x1,
                    fitKGE$parlist$x2,
                    fitKGE$parlist$x3,
                    fitKGE$parlist$x4,
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
                    fitNSE$parlist$x1,
                    fitNSE$parlist$x2,
                    fitNSE$parlist$x3,
                    fitNSE$parlist$x4,
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
                       fitRSQLog$parlist$x1,
                       fitRSQLog$parlist$x2,
                       fitRSQLog$parlist$x3,
                       fitRSQLog$parlist$x4,
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
  
  GR4JSummary <- as.data.frame(r, row.names = FALSE)
  colnames(GR4JSummary) <- GR4JColNames
  
  write.csv(GR4JSummary, paste0(runName, '.csv'), row.names = FALSE)
}

# Runs
# Bruce's naming convention
# x2 = c(-10, 10) = Fitted_xxxx_xxxx_3
# x2 = c(-30, 20) = Fitted_xxxx_xxxx_2
# I should change this

parallelRuns(calStart = '1990-01-01',
             calEnd = '1996-12-31',
             x2 = c(-10,10))
parallelRuns(calStart = '1990-01-01',
             calEnd = '1996-12-31',
             x2 = c(-30,20))
parallelRuns(calStart = '2000-01-01',
             calEnd = '2006-12-31',
             x2 = c(-10,10))
parallelRuns(calStart = '2000-01-01',
             calEnd = '2006-12-31',
             x2 = c(-30,20))