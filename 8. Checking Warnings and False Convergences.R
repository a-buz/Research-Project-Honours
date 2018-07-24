# Analysis of Results
library(tidyverse)

# Read in station list
finalStations <- read_csv('Data/FinalStations.csv')

# Fitted models
fittedModels <- list.dirs('Data/GR4J')
# Remove base dir
fittedModels <- fittedModels[-1]

# Making a table for the messages
messages <- matrix(ncol = length(fittedModels)+1,
                   nrow = nrow(finalStations))
colnames(messages) <- c('Number', basename(fittedModels))
messages <- as.tibble(messages)
messages$Number <- finalStations$`AWRC Station Number`

# Make separate tables for the 3 objFuns
KGE <- messages
NSE <- messages
RSqLog <- messages

# Fill in tables
for (i in 1:nrow(messages)) {
  for(j in 2:length(colnames(messages))) {
    load(paste0('Data/GR4J/', colnames(messages)[j], '/', messages$Number[i], '.Rdata'))
    KGE[i,j] <- fitted$KGE$fit.result$message
    NSE[i,j] <- fitted$NSE$fit.result$message
    RSqLog[i,j] <- fitted$RSquaredLog$fit.result$message
  }
}

# Join tables together
KGE <- KGE %>% 
  gather(Model, Messages, -Number) %>% 
  mutate(ObjFun = 'KGE')
NSE <- NSE %>% 
  gather(Model, Messages, -Number) %>% 
  mutate(ObjFun = 'NSE')
RSqLog <- RSqLog %>% 
  gather(Model, Messages, -Number) %>% 
  mutate(ObjFun = 'RSqLog')

# Sort
combinedFuns <- bind_rows(KGE,NSE,RSqLog) %>% 
  group_by(Model, Number) %>% 
  arrange(.by_group=T)

# Find false convergences and non-convergences
falseConverge <- 'false convergence (8)'
evalLimit <- 'function evaluation limit reached without convergence (9)'

converged <- combinedFuns %>% 
  filter(!Messages %in% c(falseConverge, evalLimit))

# Write
write_csv(converged, 'Data/GR4J/ConvergenceChecked.csv')
