
for (i in 1:length(obs_flow_1990)){
  load(paste("G:\\University\\Honours\\Data Access\\Modelling\\zoo SILO\\", obs_flow_1990[i], sep = ""))
  
  #Getting right dates for the observed flow
  obs <- window(river, start = "1990-04-10", end = "1996-12-31")
  #Extracting number, name and Area for conversion
  river.number <- as.character(strsplit(obs_flow_1990[i], ".Rdata"))
  Area <- read.csv(paste("G:\\University\\Honours\\Data Access\\Modelling\\Rivers\\",river.number, ".csv", sep = ""),
                   sep = c(" ", ",")) ##Read in river file to extract area
  Area <- as.numeric(substr(Area[16,1], 1, gregexpr("\\,",Area[16,1])[[1]]-1)) ##Extract Area
  obs <- convertFlow(obs, from ="ML", area.km2= Area) #Convert flow from ML to mm
  
  Proportion <- length(which(obs$Q == 0))/length(obs$Q)*100
  
  Type[i, 1] <- river.number
  Type[i, 2] <- ifelse(Proportion <= 1 , "Perennial", ifelse(Proportion > 1, "Ephemeral"))
  
  rm(river)
}

write.csv(Type, file = "G:\\University\\Honours\\Data Access\\Modelling\\Results\\Type_1990.csv")
