#------------------------------------------------#
# Script 2                                       #
# Extract river accuracy data                    #
#------------------------------------------------#

#Libraries
library(tidyverse)

# Get list of river data files downloaded from BoM
rivers <- list.files("Data/Rivers", full.names = TRUE)
# Set up column names for use in loop
name <- c('Date', 'Flow', 'QCode')
# Empty dataframe to store the data in 
river.accuracy <- data.frame(AWRC.Station.Number = rep(NA, length(rivers)), 
                             Station.Name = rep(NA, length(rivers)),
                             Flow.Accuracy = rep(NA, length(rivers))) 

# Extract river accuracy
for (i in 1:length(rivers)) {
  # Read header file from data to extract AWRC river number and name
  header <- read_csv(rivers[i], col_names = F, skip = 14, n_max = 1)
  # Get full name and AWRC number from header by collapsing rows
  river.name <- paste(header[1,-1], collapse = ' ')
  # Extract number using gsub
  river.accuracy[i,1] <- gsub(".*\\((.*)\\).*", "\\1", river.name)
  # Extract name
  river.accuracy[i,2] <- gsub("\\s*\\([^\\)]+\\)", "", river.name)
  
  # Read river file. Skipping header and default column names
  riverFile <- read_csv(rivers[i], skip = 27, col_names = name)
  # Convert data to POSIXct format and extract study dates
  riverFile$Date <- as.POSIXct(as.character(riverFile$Date), format = "%Y-%m-%d")
  riverFile <- subset(riverFile, Date >= "1990-01-01" & Date <= "2010-12-31")
  
  # Make sure the station operates for all dates
  # Find how much of the data has an accuracy code of A
  if(min(riverFile$Date) == as.Date('1990-01-01') & max(riverFile$Date) == as.Date('2010-12-31')) {
    river.accuracy[i,3] <- length(riverFile$QCode[which(riverFile$QCode == "A")])/length(riverFile$QCode)*100  
  } else {
    river.accuracy[i,3] <- 0 # 0 if not for all dates
  }
}
# End loop

# Order by accuracy and then cut off less than 90
river.accuracy <- river.accuracy[order(-river.accuracy[,3], river.accuracy[,2], river.accuracy[,1]),]
river.accuracy <- river.accuracy[which(river.accuracy$Flow.Accuracy > 90),]

# Write CSV
write_csv(river.accuracy, path = "Data/CatchmentQuality.csv")
