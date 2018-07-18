##Analysis of Results:

setwd("E:\\University\\Honours\\Data Access\\Modelling\\Results\\Round 5")

##Fitted modelled flow (2000_2006)
sim6y1990.2   <- dir(path = ".\\Fitted_1990_1996_2\\", pattern = ".rdata")
sim6y1990.3   <- dir(path = ".\\Fitted_1990_1996_3\\", pattern = ".rdata")
sim10y1990.2  <- dir(path = ".\\Fitted_1990_1996_2\\", pattern = ".rdata")
sim10y1990.3  <- dir(path = ".\\Fitted_1990_1996_3\\", pattern = ".rdata")
sim6y2000.2   <- dir(path = ".\\Fitted_2000_2006_2\\", pattern = ".rdata")
sim6y2000.3   <- dir(path = ".\\Fitted_2000_2006_3\\", pattern = ".rdata")
sim10y2000.2  <- dir(path = ".\\Fitted_2000_2006_2\\", pattern = ".rdata")
sim10y2000.3  <- dir(path = ".\\Fitted_2000_2006_3\\", pattern = ".rdata")

obs_flow <- dir(path = "E:\\University\\Honours\\Data Access\\Modelling\\zoo SILO", pattern = ".Rdata")

##Making a table for the messages
Message <- data.frame(Number = rep(NA, length(obs_flow)), 
                      y6.1990.2  =  rep(NA, length(obs_flow)),
                      y6.1990.3  =  rep(NA, length(obs_flow)),
                      y10.1990.2 =  rep(NA, length(obs_flow)),
                      y10.1990.3 =  rep(NA, length(obs_flow)),
                      y6.2000.2  =  rep(NA, length(obs_flow)),
                      y6.2000.3  =  rep(NA, length(obs_flow)),
                      y10.2000.2 =  rep(NA, length(obs_flow)),
                      y10.2000.3 =  rep(NA, length(obs_flow)))

#Names
for (i in 1:length(obs_flow)) {
  load(paste("E:\\University\\Honours\\Data Access\\Modelling\\zoo SILO\\", obs_flow[i], sep = ""))
  
  ##Extracting number, name and Area for conversion
  river.number <- as.character(strsplit(obs_flow[i], ".Rdata"))
  
  Message[i, 1] <- river.number
}

##Making a list of functions
messages <-   function (x) {
    load(paste(dir()[1], "\\", x, sep = ""))
    
    message <- riverFit$fit.result$message
  }
message <- lapply(sim6y1990.2, FUN = messages)
Message$y6.1990.2 <- message

messages <-   function (x) {
  load(paste(dir()[2], "\\", x, sep = ""))
  
  message <- riverFit$fit.result$message
}
message <- lapply(sim6y1990.3, FUN = messages)
Message$y6.1990.3 <- message

messages <-   function (x) {
  load(paste(dir()[3], "\\", x, sep = ""))
  
  message <- riverFit$fit.result$message
}
message <- lapply(sim10y1990.2, FUN = messages)
Message$y10.1990.2 <- message

messages <-   function (x) {
  load(paste(dir()[4], "\\", x, sep = ""))
  
  message <- riverFit$fit.result$message
}
message <- lapply(sim10y1990.3, FUN = messages)
Message$y10.1990.3 <- message

messages <-   function (x) {
  load(paste(dir()[5], "\\", x, sep = ""))
  
  message <- riverFit$fit.result$message
}
message <- lapply(sim6y1990.2, FUN = messages)
Message$y6.2000.2 <- message

messages <-   function (x) {
  load(paste(dir()[6], "\\", x, sep = ""))
  
  message <- riverFit$fit.result$message
}
message <- lapply(sim6y1990.3, FUN = messages)
Message$y6.2000.3 <- message

messages <-   function (x) {
  load(paste(dir()[7], "\\", x, sep = ""))
  
  message <- riverFit$fit.result$message
}
message<- lapply(sim10y1990.2, FUN = messages)
Message$y10.2000.2 <- message

messages <-   function (x) {
  load(paste(dir()[8], "\\", x, sep = ""))
  
  message <- riverFit$fit.result$message
}
message <- lapply(sim10y1990.3, FUN = messages)
Message$y10.2000.3 <- message

##Taking care of the "unimplemented type 'list' in 'EncodeElement" error
for (i in 1:8) {
  Message[,i+1] <- vapply(Message[,i+1], paste, collapse = ", ", character(1L))
}

write.csv(Message, file = "E:\\University\\Honours\\Data Access\\Modelling\\Results\\Convergence Checks.csv", row.names = FALSE)