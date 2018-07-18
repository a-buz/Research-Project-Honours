#
# Script 6
#

# Libraries
library(tidyverse)
library(zoo)
library(hydromad)

# Flow duration curves
zooFiles <- dir('Data/Zoo', full.names = TRUE)
riverTable <- read_csv("Data/FinalStations.csv")

FDC.summary <- data.frame(Number = rep(NA, length(zooFiles)), 
                          Name = rep(NA, length(zooFiles)),
                          P25.FDC1 = rep(NA, length(zooFiles)),
                          P50.FDC1 = rep(NA, length(zooFiles)),
                          P75.FDC1 = rep(NA, length(zooFiles)),
                          P90.FDC1 = rep(NA, length(zooFiles)),
                          P25.FDC2 = rep(NA, length(zooFiles)),
                          P50.FDC2 = rep(NA, length(zooFiles)),
                          P75.FDC2 = rep(NA, length(zooFiles)),
                          P90.FDC2 = rep(NA, length(zooFiles)))

percentiles <- c(0.25, 0.5, 0.75, 0.9)

# Plot FDC
for (i in 1:length(zooFiles)) {
  # River name
  riverDetails <- riverTable %>%
    filter(`AWRC Station Number` == basename(zooFiles[i]))
  
  # Read in zoo file
  river <- read.table(zooFiles[i], header = TRUE)
  river$Index <- as.Date(river$Index)
  river <- zoo(river[,2:4], order.by = river$Index)
  
  #Subset zoo files into 1990-1996 and 2000-2006 time periods
  FDC1 <- window(river, start = "1990-01-01", end = "2000-12-31")
  FDC2 <- window(river, start = "2001-01-01", end = "2010-12-31")
  
  #streamflow duration curve analysis
  n <- nrow(FDC1)
  sort.flow1 <- sort(as.numeric(FDC1$Q), decreasing = TRUE, na.last=FALSE)
  #Sort - just a series of numbers now
  rank.flow1 <- 1:n
  Prob1 <- rank.flow1/(n+1) 
  
  m <- nrow(FDC2)
  sort.flow2 <- sort(as.numeric(FDC2$Q), decreasing = TRUE, na.last=FALSE)
  #Sort - just a series of numbers now
  rank.flow2 <- 1:m
  Prob2 <- rank.flow2/(m+1) 
  
  ##Extracting Quantlies from FDCs
  FDC1.quantiles <- quantile(sort.flow1, percentiles)
  FDC2.quantiles <- quantile(sort.flow2, percentiles)
  
  #Putting results into data.frame
  FDC.summary[i,1] <- riverDetails$`AWRC Station Number`
  FDC.summary[i,2] <- riverDetails$`Station Name`
  FDC.summary[i,3] <- FDC1.quantiles[1]
  FDC.summary[i,4] <- FDC1.quantiles[2]
  FDC.summary[i,5] <- FDC1.quantiles[3]
  FDC.summary[i,6] <- FDC1.quantiles[4]
  FDC.summary[i,7] <- FDC2.quantiles[1]
  FDC.summary[i,8] <- FDC2.quantiles[2]
  FDC.summary[i,9] <- FDC2.quantiles[3]
  FDC.summary[i,10] <- FDC2.quantiles[4]

}

par(mfrow=c(2,2))
plot(FDC.summary[,3], FDC.summary[,7], xlab = "1990-2000", ylab = "2001-2010", col = "red", main = "25th percentile")
abline(0,1)
plot(FDC.summary[,4], FDC.summary[,8], xlab = "1990-2000", ylab = "2001-20010", col = "blue", main = "50th percentile")
abline(0,1)
plot(FDC.summary[,5], FDC.summary[,9], xlab = "1990-2000", ylab = "2001-2010", col = "green", main = "75th percentile")
abline(0,1)
plot(FDC.summary[,6], FDC.summary[,10], xlab = "1990-2000", ylab = "2001-2010", col = "purple", main = "90th percentile")
abline(0,1)

ratio = FDC.summary[,3]/FDC.summary[,7]
length(ratio[which(ratio > 1)])

#-------------------------------------------------------------------------------------------------#
# GR4J
#-------------------------------------------------------------------------------------------------#
# 
# GR4J.summary <-  data.frame(Number = rep(NA, length(zooFiles)), 
#                             Name = rep(NA, length(zooFiles)),
#                             Area = rep(NA, length(zooFiles)),
#                             x1 = rep(NA, length(zooFiles)),
#                             x2 = rep(NA, length(zooFiles)),
#                             x3 = rep(NA, length(zooFiles)),
#                             x4 = rep(NA, length(zooFiles)), 
#                             runoff = rep(NA, length(zooFiles)),
#                             rbias = rep(NA, length(zooFiles)),
#                             NSE = rep(NA, length(zooFiles)),
#                             r.sq.sqrt = rep(NA, length(zooFiles)),
#                             r.sq.log = rep(NA, length(zooFiles)),
#                             objFun = rep(NA, length(zooFiles)))
# 
# 
# #Setting up the table to input information - pretty much everything of importance from the model. 
# #Warning 
# 
# GR4J.summary$Number <- as.character(GR4J.summary$Number)
# GR4J.summary$Name <- as.character(GR4J.summary$Name)
# 
# # Make function
# fitRiver <- function(data) {
#   
# }
# 
# ##Reading in the data
# for (i in 1:5) {#length(zooFiles)) {
#   #load(zooFiles[i])
#   
#   # River name
#   riverDetails <- riverTable %>%
#     filter(`AWRC Station Number` == basename(zooFiles[i]))
#   
#   # Read in zoo file
#   river <- read.table(zooFiles[i], header = TRUE)
#   river$Index <- as.Date(river$Index)
#   river <- zoo(river[,2:4], order.by = river$Index)
#   
#   # Set cal period
#   river.cal <- window(river, start = "1990-01-01", end = "1996-12-31")
#   
#   CMod <- hydromad(river.cal,
#                    sma="gr4j",
#                    routing="gr4jrouting",
#                    etmult=0.15,
#                    x1 = c(50,2000),
#                    x2 = c(-10,10),
#                    x3 =c(5,500),
#                    x4 = c(0.5,10)) 
#   #S_0 = 0.5, R_0 = 0.5)
#   # hydromad.stats("viney" = function(Q, X, ...) {
#   #   hmadstat("r.squared")(Q, X, ...) -
#   #     5*(abs(log(1+hmadstat("rel.bias")(Q,X)))^2.5)})
#   
#   # Fit. This uses viney which needs to be added to your hydromad stats
#   # See bottom of help(hydromad.stats)
#   riverFit <- fitByOptim(CMod,
#                          objective = hmadstat("viney"),
#                          samples=500,
#                          method="PORT") ##Let's start by maximising R2 (NSE))
#   
#   GR4J.summary[i,1] <- riverDetails$`AWRC Station Number`
#   GR4J.summary[i,2] <- riverDetails$`Station Name`
#   GR4J.summary[i,3] <- riverDetails$`Catchment Area (km2)`
#   GR4J.summary[i,4] <- riverFit$parlist$x1
#   GR4J.summary[i,5] <- riverFit$parlist$x2
#   GR4J.summary[i,6] <- riverFit$parlist$x3
#   GR4J.summary[i,7] <- riverFit$parlist$x4
#   GR4J.summary[i,8] <- summary(riverFit)$runoff
#   GR4J.summary[i,9] <- summary(riverFit)$rel.bias
#   GR4J.summary[i,10] <- summary(riverFit)$r.squared
#   GR4J.summary[i,11] <- summary(riverFit)$r.sq.sqrt
#   GR4J.summary[i,12] <- summary(riverFit)$r.sq.log
#   GR4J.summary[i,13] <- objFunVal(riverFit)
#   
#   fitted <- riverFit$fitted.values
#   
#   if(dir.exists('Data/GR4J/Fitted_1990_1996_3') == FALSE) {
#     dir.create('Data/GR4J/Fitted_1990_1996_3', recursive = TRUE)
#   }
#   
#   save(fitted, file = paste("Data/GR4J/Fitted_1990_1996_3/", riverDetails$`AWRC Station Number`, ".rdata", sep = ""))
# }
# 
# # The actual model itself. An rdata file was created for every catchment modelled. A summary file was also created. 
# # At this stage, I was still contemplating using 6-10 years, so 4 runs were required using this code, or rather, 1 run using 4 computers. 
# # With each run, the dates, folder paths, and folder name for results had to be changed. 
# # The "_2" and "_3" refer to different parameter boundaries for x2. I decided to compare it because I was unsure which was the best for it. 
# # GR4J_1990_1996_3 <- x2(-10,10)
# # "           "_2 <- x2(-30,20) 
