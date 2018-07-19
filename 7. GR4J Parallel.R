#--------------------------------------------------------------------------------------------------#
# GR4J model run                                                                                   #
#--------------------------------------------------------------------------------------------------#

library(foreach)
library(doParallel)

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
                   'x1',
                   'x2',
                   'x3',
                   'x4',
                   'runoff',
                   'rbias',
                   'NSE',
                   'r.sq.sqrt',
                   'r.sq.log',
                   'objFunVal')

# Set up dir for output
if(dir.exists('Data/GR4J/Fitted_1990_1996_3') == FALSE) {
  dir.create('Data/GR4J/Fitted_1990_1996_3', recursive = TRUE)
}

#--------------------------------------------------------------------------------------------------#
# Parallel loop                                                                                    #
#--------------------------------------------------------------------------------------------------#
# Set cores to 1 minus cores on machine, including logical cores.
# Set cluster to 1 less of the amount of cores 
cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)

# Create log file to check progress
writeLines(c(""), "log.txt")

# Start parloop and assign to collector
r <- foreach(i=1:nrow(riverTable), .combine=rbind, .packages=c('hydromad', 'zoo')) %dopar% {
  # Write to logfile
  sink("log.txt", append=TRUE)
  cat(paste(Sys.time(), "Starting iteration", i, "of", nrow(riverTable), "\n"))
  sink() # End diversion
  
  # Get river data
  riverDetails <- riverTable[riverTable$AWRC.Station.Number == basename(zooFiles[i]),]
  
  # Read in zoo file from csv. read.csv.zoo() will not work for some reason
  river <- read.table(zooFiles[i], header = TRUE)
  river$Index <- as.Date(river$Index)
  river <- zoo(river[,2:4], order.by = river$Index)
  
  # Set cal period and create hydromad variable
  river.cal <- window(river, start = "1990-01-01", end = "1996-12-31")
  CMod <- hydromad(river.cal,
                   sma="gr4j",
                   routing="gr4jrouting",
                   etmult=0.15,
                   x1 = c(50,2000),
                   x2 = c(-10,10),
                   x3 =c(5,500),
                   x4 = c(0.5,10))
  
  # Fit KGE
  fitKGE <- fitByOptim(CMod,
                       objective = function(Q, X, ...) {
                         1 - sqrt(
                           (cor(X, Q) - 1)^2 +
                             (mean(X)/mean(Q) - 1)^2 +
                             (sd(X)/sd(Q) - 1)^2
                         )
                       },
                       samples = 500,
                       method = "PORT")
  
  resultKGE  <- c(riverDetails$AWRC.Station.Number,
                  riverDetails$Station.Name,
                  riverDetails$Catchment.Area..km2.,
                  'KGE',
                  fitKGE$parlist$x1,
                  fitKGE$parlist$x2,
                  fitKGE$parlist$x3,
                  fitKGE$parlist$x4,
                  summary(fitKGE)$runoff,
                  summary(fitKGE)$rel.bias,
                  summary(fitKGE)$r.squared,
                  summary(fitKGE)$r.sq.sqrt,
                  summary(fitKGE)$r.sq.log,
                  objFunVal(fitKGE)
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
                  fitNSE$parlist$x1,
                  fitNSE$parlist$x2,
                  fitNSE$parlist$x3,
                  fitNSE$parlist$x4,
                  summary(fitNSE)$runoff,
                  summary(fitNSE)$rel.bias,
                  summary(fitNSE)$r.squared,
                  summary(fitNSE)$r.sq.sqrt,
                  summary(fitNSE)$r.sq.log,
                  objFunVal(fitNSE)
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
                     fitNSE$parlist$x1,
                     fitNSE$parlist$x2,
                     fitNSE$parlist$x3,
                     fitNSE$parlist$x4,
                     summary(fitNSE)$runoff,
                     summary(fitNSE)$rel.bias,
                     summary(fitNSE)$r.squared,
                     summary(fitNSE)$r.sq.sqrt,
                     summary(fitNSE)$r.sq.log,
                     objFunVal(fitNSE)
  )
  
  # Create list of fitted mode objects
  fitted <- list('KGE' = fitKGE,
                 'NSE' = fitNSE,
                 'RSquaredLog' = fitRSQLog)
  
  # Write
  save(fitted, file = paste("Data/GR4J/Fitted_1990_1996_3/", riverDetails$AWRC.Station.Number, ".Rdata", sep = ""))
  
  # Get summary results
  results <- t(matrix(c(resultKGE,resultNSE,resultRSQLog), ncol = 3))
  return(results)
}

# Stop clusters
stopImplicitCluster()
stopCluster(cl)

GR4JSummary <- as.data.frame(r, row.names = FALSE)
colnames(GR4JSummary) <- GR4JColNames

write.csv(GR4JSummary, 'Data/GR4JSummary.csv', row.names = FALSE)