##Analysis of Results:

setwd("F:\\University\\Honours\\Research Project\\Modelling\\Results\\Round 5")
GR4J.sum <- dir(pattern = "csv")
GR4J.sum

#Checking COnvergences
Convergences <- read.csv("F:\\University\\Honours\\Research Project\\Modelling\\Results\\Convergence Checks.csv",
                         header = TRUE)
for (i in 1:8) {
  exclude <- which(Convergences[,i+1] == "false convergence (8)"|
                   Convergences[,i+1] == "function evaluation limit reached without convergence (9)")
  print(exclude)
}

exclude_1990 <- which(Convergences$y6.1990.2 == "false convergence (8)"|
                        Convergences$y6.1990.2 == "function evaluation limit reached without convergence (9)")
exclude_2000 <- which(Convergences$y6.2000.2 == "false convergence (8)"|
                        Convergences$y6.2000.2 == "function evaluation limit reached without convergence (9)")

##Updating based on the exclusions

##Flow duration Curve Analysis
library(hydromad)
library(ggplot2)
rm(river)
rm(riverFit)

##Signature Analysis ##Chiew et al. (2013) line
GR4J_1990_1996_2 <- read.csv(GR4J.sum[1], header = TRUE)
deleted_1990 <- GR4J_1990_1996_2[c(exclude_1990),]
deleted_1990_2 <-  GR4J_2000_2006_2[c(exclude_1990),]
GR4J_1990_1996_2 <- GR4J_1990_1996_2[-c(exclude_1990),]

GR4J_2000_2006_2 <- read.csv(GR4J.sum[5], header = TRUE)
deleted_2000 <- GR4J_2000_2006_2[c(exclude_2000),]
deleted_2000_2 <- GR4J_1990_1996_2[c(exclude_2000),]
GR4J_2000_2006_2 <- GR4J_2000_2006_2[-c(exclude_2000),]


#Rounding off these figures to 2 dp
GR4J_1990_1996_2[,c(4:13)] <- round(GR4J_1990_1996_2[,c(4:13)], 2)
GR4J_2000_2006_2[,c(4:13)] <- round(GR4J_2000_2006_2[,c(4:13)], 2)

##Plotting these figures

##Low flows >= 90th  percentile
##High flows <= 25th percentile
##Mid flow 25th < mid flow < 90th

obs_flow <- dir("F:\\University\\Honours\\Research Project\\Modelling\\zoo SILO", 
                pattern = ".Rdata")
obs_flow_1990 <- obs_flow[c(exclude_1990)]
obs_flow_2000 <- obs_flow[c(exclude_2000)]

sim_flow <- dir("E:\\University\\Honours\\Research Project\\Modelling\\Results\\Round 5\\Fitted_1990_1996_2")
sim_flow_1990 <- sim_flow[-c(exclude_1990)]
sim_flow_2000 <- sim_flow[-c(exclude_2000)]

obs.summary2 <- data.frame(Number   = rep(NA, length(obs_flow_2000)), 
                          cum.flow  = rep(NA, length(obs_flow_2000)),
                          FDC.low   = rep(NA, length(obs_flow_2000)),
                          FDC.mid   = rep(NA, length(obs_flow_2000)),
                          FDC.high  = rep(NA, length(obs_flow_2000)),
                          AC        = rep(NA, length(obs_flow_2000)),
                          Peaks     = rep(NA, length(obs_flow_2000)),
                          slope     = rep(NA, length(obs_flow_2000)))

for (i in 1:length(obs_flow_2000)) {
  #load(paste(".\\Fitted_2000_2006_2\\", sim_flow_2000[i], sep = ""))
  load(paste("E:\\University\\Honours\\Data Access\\Modelling\\zoo SILO\\", obs_flow_1990[i], sep = ""))
  
  #Getting right dates for the observed flow
  obs <- window(river, start = "1990-04-10", end = "1996-12-31")

  #Extracting number, name and Area for conversion
  river.number <- as.character(strsplit(obs_flow_1990[i], ".Rdata"))
  Area <- read.csv(paste("E:\\University\\Honours\\Data Access\\Modelling\\Rivers\\",river.number, ".csv", sep = ""),
                   sep = c(" ", ",")) ##Read in river file to extract area
  Area <- as.numeric(substr(Area[16,1], 1, gregexpr("\\,",Area[16,1])[[1]]-1)) ##Extract Area
  obs <- convertFlow(obs, from ="ML", area.km2= Area) #Convert flow from ML to mm
  
  ##Flow duration curve ranking and sorting
  # n <- nrow(obs)
  # sort.flow1 <- sort(as.numeric(obs$Q), decreasing = TRUE, na.last=FALSE)
  # #Sort - just a series of numbers now
  # rank.flow1 <- 1:n
  # Prob1 <- rank.flow1/(n+1)
  
  # m <- length(fitted(riverFit))
  # sort.flow2 <- sort(as.numeric(fitted(riverFit)), decreasing = TRUE, na.last=FALSE)
  # #Sort - just a series of numbers now
  # rank.flow2 <- 1:m
  # Prob2 <- rank.flow2/(m+1)
  
  ##FDC total - Performance metric: NSE
  # FDC.total <- nseStat(sort.flow1)#, sort.flow2)
  
  #######cum.flow
  cum.flow <- sum(obs$Q)
  #cum.flow <- sum(fitted(riverFit))
  
  ##FDC.low - Performance metroc: F
  FDC.low <- unname(quantile(obs$Q, probs = 0.1)) #abs(1 - quantile(fitted(riverFit), probs = 0.25)/))
  #FDC.low <- unname(quantile(fitted(riverFit), probs = 0.1)) #/quantile(obs$Q, probs = 0.25)))
  
  ####FDC.High - Performance metroc: F
  FDC.high <- unname(quantile(obs$Q, probs = 0.99))#abs(1 - quantile(fitted(riverFit), probs = 0.75)/))
  #FDC.high <- unname(quantile(fitted(riverFit), probs = 0.9))#/quantile(obs$Q, probs = 0.75)))
  
  ######FDC.Mid - Performance metroc: F
  FDC.mid <- unname(quantile(obs$Q, probs = 0.5)) #abs(1 - quantile(fitted(riverFit), probs = 0.5)/))
  #FDC.mid <- unname(quantile(fitted(riverFit), probs = 0.5))#/quantile(obs$Q, probs = 0.5)))
  
  ######AC - Performance metroc: F
  AC.obs <- acf(obs$Q,  plot = F)$acf[2]
  #AC.sim <- acf(fitted(riverFit),  plot = F)$acf[2]
  # AC <- abs(1 - AC.sim/AC.obs)
  # AC <- AC.sim/AC.obs
  AC <- AC.obs
  
  ##Peak Distribution - Perforamnce Metric: F
  Peak.obs <- (quantile(obs$Q, probs = 0.9)-(quantile(obs$Q, probs = 0.5)))/(0.9-0.5)
  #Peak.sim <- (quantile(fitted(riverFit), probs = 0.9)-(quantile(fitted(riverFit), probs = 0.5)))/(0.9-0.5)
  #Peaks <- abs(1 - Peak.sim/Peak.obs)
  # Peaks <- Peak.sim/Peak.obs
  Peaks <- Peak.obs
  
  #slope mid FDC
  slope.obs <- unname(diff((quantile(obs$Q, probs = c(0.33, 0.66))/0.66-0.33)))
  #slope.sim <- unname(diff(quantile(fitted(riverFit), probs = c(0.25, 0.75))/0.75-0.25))
  #slope <- abs(1 - slope.sim/slope.obs)
  # slope <- slope.sim/slope.obs
  slope <- slope.obs
  
  ##Putting into d.frame
  obs.summary2[i,1] <- river.number
  #sig.summary[i,2] <- signif(FDC.total, 4)
  obs.summary2[i,2] <- cum.flow
  obs.summary2[i,3] <- signif(FDC.low, 4)
  obs.summary2[i,4] <- signif(FDC.mid, 4)
  obs.summary2[i,5] <- signif(FDC.high, 4)
  obs.summary2[i,6] <- signif(AC, 4)
  obs.summary2[i,7] <- signif(Peaks, 4)
  obs.summary2[i,8] <- signif(slope, 4)
  
  rm(obs)
  #rm(riverFit)
}
xyplot(riverFit)
write.csv(obs.summary2, file = "obs.summary3_2000.csv")

#Change Inf's to NA
sim.summary <- do.call(data.frame,lapply(sim.summary, function(x) replace(x, is.infinite(x),NA)))

#Correlation matrix
cor(sig.summary[,c(2:9)], use = "pairwise.complete.obs", method = "spearman")

##install.packages("hydroTSM") ##For nseStat function
library(hydroTSM)

########-----------Plotting over Australia---------#############
#List of stations to extract coordinates:
Gauges <- read.csv("E:\\University\\Honours\\Data Access\\Modelling\\Final_RiverStations.csv",
                    header = T)
Gauges <- Gauges[order(Gauges[,1]),]
rownames(Gauges) <- NULL

sig.summary$Lat <- Gauges$LAT
sig.summary$Lon <- Gauges$LON

library(ggplot2)

ggplot(sig.summary, aes(x = Lon, y = Lat))  + geom_point(aes(color = GR4J_1990_1996_2$objFun)) + 
  scale_color_continuous(low = "blue", high = "green")

###---------annual variation and random shit using hydromad?acf------##############

ann.sim <- aggregate(fitted(riverFit),by=list(Year=format(time(fitted(riverFit)),"%Y")), sum, na.rm=TRUE)
ann.obs <- aggregate(obs$Q, by=list(Year=format(time(obs),"%Y")), sum, na.rm=TRUE)


plot(as.numeric(time(ann.sim)),ann.sim, ylab="Flow (ML/year)", xlab="Year",type="h", col="blue")
lines(as.numeric(time(ann.obs)),ann.obs, ylab="Flow (ML/year)", xlab="Year",type="h", col="red")
ann.sim-ann.obs

xyplot(riverFit)
stats <- summary(riverFit, breaks = "1 year")
statsSeries <- stats[,c("r.squared", "r.sq.sqrt", "rel.bias", "runoff")]
statsSeries[,c(1,2)] <- pmax(statsSeries[, 1:2], 0)
c(xyplot(statsSeries, type = "s", lwd = 2, ylab = "statistic",
         xlab = NULL), `observed streamflow` = xyplot(observed(riverFit)),
  layout = c(1, 5), x.same = TRUE) + layer_(panel.refline(h = 0,
                                                          v = time(statsSeries)))
xyplot(window(obs, start = "2002-01-01", end = "2002-12-31"))

qqmath(riverFit, scales = list(y = list(log = TRUE)),
       type = c("l", "g"))

qqmath(riverFit, type = c("l", "g"), scales = list(y = list(log = TRUE)),
       xlab = "Standard normal variate", ylab = "Flow (mm/day)",
       f.value = ppoints(100), tails.n = 50, as.table = TRUE)

nseStat(window(fitted(riverFit), start = "2000-01-01", end = "2001-12-31"),window(obs$Q, start = "2000-01-01", end = "2001-12-31"))


