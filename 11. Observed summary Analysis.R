#
# Observed summary analysis
#

library(tidyverse)

analyseMod <- function(model, timePeriod) {
  # Model summary
  modelSummary <- read_csv(paste('Data/Analysis/', timePeriod, '/', model, 'Summary.csv', sep = ''))
  
  # Split into different obj funs
  kgeMod <- modelSummary %>% 
    filter(objFun=='KGE')
  nseMod <- modelSummary %>% 
    filter(objFun=='NSE')
  rsqlMod <- modelSummary %>% 
    filter(objFun=='RSquaredLog')
  
  # PCA vars
  # objfun -> ln(area), ln(cf), ln(fdc.low), ln(fdc.mid), ln(fdc.high), AC, ln(peaks), ln(slope)
  
  # KGE
  kgeVars <- kgeMod %>% 
    dplyr::select(KGE, Area, Cum.flow, FDC.low,
                  FDC.mid, FDC.high, AC, Peaks,
                  Slope, KGext, Type)
  # Log vars
  kgeVarsTf <- kgeVars %>% 
    mutate(Area     = log(Area),
           Cum.flow = log(Cum.flow + 0.1),
           FDC.low  = log(FDC.low + 0.1),
           FDC.mid  = log(FDC.mid + 0.1),
           FDC.high = log(FDC.high + 0.1),
           Peaks    = log(Peaks + 0.1),
           Slope    = log(Slope + 0.1))
  
  kgePCA <- prcomp(kgeVarsTf %>% select(-KGext, -Type), scale = T, retx = T, center = T)
  kgeLM <- lm(KGE ~ ., data = kgeVarsTf %>% select(-KGext))
  
  # NSE
  nseVars <- nseMod %>% 
    dplyr::select(NSE, Area, Cum.flow, FDC.low,
                  FDC.mid, FDC.high, AC, Peaks,
                  Slope, KGext, Type)
  # Log vars
  nseVarsTf <- nseVars %>% 
    mutate(Area     = log(Area),
           Cum.flow = log(Cum.flow + 0.1),
           FDC.low  = log(FDC.low + 0.1),
           FDC.mid  = log(FDC.mid + 0.1),
           FDC.high = log(FDC.high + 0.1),
           Peaks    = log(Peaks + 0.1),
           Slope    = log(Slope + 0.1))
  
  nsePCA <- prcomp(nseVarsTf %>% select(-KGext, -Type), scale = T, retx = T, center = T)
  nseLM <- lm(NSE ~ ., data = nseVarsTf %>% select(-KGext))
  
  # RSqLog
  rsqlVars <- rsqlMod %>% 
    dplyr::select(r.sq.log, Area, Cum.flow, FDC.low,
                  FDC.mid, FDC.high, AC, Peaks,
                  Slope, KGext, Type)
  # Log vars
  rsqlVarsTf <- rsqlVars %>% 
    mutate(Area     = log(Area),
           Cum.flow = log(Cum.flow + 0.1),
           FDC.low  = log(FDC.low + 0.1),
           FDC.mid  = log(FDC.mid + 0.1),
           FDC.high = log(FDC.high + 0.1),
           Peaks    = log(Peaks + 0.1),
           Slope    = log(Slope + 0.1))
  
  rsqlPCA <- prcomp(rsqlVarsTf %>% select(-KGext, -Type), scale = T, retx = T, center = T)
  rsqlLM <- lm(r.sq.log ~ ., data = rsqlVarsTf %>% select(-KGext))
  
  # Regression
  returnList <- list(kgePCA  = kgePCA,  kgeVarsTf  = kgeVarsTf,  kgeLM  = kgeLM,
                     nsePCA  = nsePCA,  nseVarsTf  = nseVarsTf,  nseLM  = nseLM,
                     rsqlPCA = rsqlPCA, rsqlVarsTf = rsqlVarsTf, rsqlLM = rsqlLM)
  return(returnList)
}

gr4j1990 <- analyseMod('GR4J', '1990_1996')
gr4j2000 <- analyseMod('GR4J', '2000_2006')

simhyd1990 <- analyseMod('SIMHYD', '1990_1996')
simhyd2000 <- analyseMod('SIMHYD', '2000_2006')

hbv1990 <- analyseMod('HBV', '1990_1996')
hbv2000 <- analyseMod('HBV', '2000_2006')

# ------------------------------------------------------------------------------
# PCA analysis
# Biplots
# Type
library(car)
library(ggfortify)
autoplot(gr4j1990$kgePCA, data = gr4j1990$kgeVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(gr4j1990$nsePCA, data = gr4j1990$nseVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(gr4j1990$rsqlPCA, data = gr4j1990$rsqlVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

autoplot(gr4j2000$kgePCA, data = gr4j2000$kgeVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(gr4j2000$nsePCA, data = gr4j2000$nseVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(gr4j2000$rsqlPCA, data = gr4j2000$rsqlVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

autoplot(simhyd1990$kgePCA, data = simhyd1990$kgeVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(simhyd1990$nsePCA, data = simhyd1990$nseVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(simhyd1990$rsqlPCA, data = simhyd1990$rsqlVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

autoplot(hbv1990$kgePCA, data = hbv1990$kgeVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(hbv1990$nsePCA, data = hbv1990$nseVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(hbv1990$rsqlPCA, data = hbv1990$rsqlVarsTf, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

# KGext
autoplot(gr4j1990$kgePCA, data = gr4j1990$kgeVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(gr4j1990$nsePCA, data = gr4j1990$nseVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(gr4j1990$rsqlPCA, data = gr4j1990$rsqlVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

autoplot(simhyd1990$kgePCA, data = simhyd1990$kgeVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(simhyd1990$nsePCA, data = simhyd1990$nseVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(simhyd1990$rsqlPCA, data = simhyd1990$rsqlVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

autoplot(hbv1990$kgePCA, data = hbv1990$kgeVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(hbv1990$nsePCA, data = hbv1990$nseVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(hbv1990$rsqlPCA, data = hbv1990$rsqlVarsTf, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

# ------------------------------------------------------------------------------
# Regression
library(car)
# GR4J
vif(gr4j1990$kgeLM)
summary(gr4j1990$kgeLM)
anova(gr4j1990$kgeLM)

vif(gr4j2000$kgeLM)
summary(gr4j2000$kgeLM)
anova(gr4j2000$kgeLM)

cor(gr4j1990$kgeVarsTf[,1:9], method = 'pearson',  use = 'pairwise.complete.obs')
cor(gr4j2000$kgeVarsTf[,1:9], method = 'spearman', use = 'pairwise.complete.obs')

# Histograms of signatures
par(mfrow = c(4,2))
colnames(gr4j1990$kgeVarsTf)
hist(gr4j1990$kgeVarsTf$Area, xlab = "log(Area)", main = NULL)
hist(gr4j1990$kgeVarsTf$Cum.flow, xlab = "log(CF)", main = NULL)
hist(gr4j1990$kgeVarsTf$FDC.low, xlab = "log(FDC.low)", main = NULL)
hist(gr4j1990$kgeVarsTf$FDC.mid, xlab = "log(FDC.mid)", main = NULL)
hist(gr4j1990$kgeVarsTf$FDC.high, xlab = "log(FDC.high)", main = NULL)
hist(gr4j1990$kgeVarsTf$AC, xlab = "AC", main = NULL)
hist(gr4j1990$kgeVarsTf$Peaks, xlab = "log(Peaks)", main = NULL)
hist(gr4j1990$kgeVarsTf$Slope, xlab = "log(Slope)", main = NULL)

# HBV
vif(hbv1990$kgeLM)
summary(hbv1990$kgeLM)
anova(hbv1990$kgeLM)

vif(hbv2000$kgeLM)
summary(hbv2000$kgeLM)
anova(hbv2000$kgeLM)

cor(hbv1990$kgeVarsTf[,1:9], method = 'pearson',  use = 'pairwise.complete.obs')
cor(hbv2000$kgeVarsTf[,1:9], method = 'spearman', use = 'pairwise.complete.obs')

# Histograms of signatures
par(mfrow = c(4,2))
colnames(hbv1990$kgeVarsTf)
hist(hbv1990$kgeVarsTf$Area, xlab = "log(Area)", main = NULL)
hist(hbv1990$kgeVarsTf$Cum.flow, xlab = "log(CF)", main = NULL)
hist(hbv1990$kgeVarsTf$FDC.low, xlab = "log(FDC.low)", main = NULL)
hist(hbv1990$kgeVarsTf$FDC.mid, xlab = "log(FDC.mid)", main = NULL)
hist(hbv1990$kgeVarsTf$FDC.high, xlab = "log(FDC.high)", main = NULL)
hist(hbv1990$kgeVarsTf$AC, xlab = "AC", main = NULL)
hist(hbv1990$kgeVarsTf$Peaks, xlab = "log(Peaks)", main = NULL)
hist(hbv1990$kgeVarsTf$Slope, xlab = "log(Slope)", main = NULL)

# SIMHYD
vif(simhyd1990$kgeLM)
summary(simhyd1990$kgeLM)
anova(simhyd1990$kgeLM)

vif(simhyd2000$kgeLM)
summary(simhyd2000$kgeLM)
anova(simhyd2000$kgeLM)

cor(simhyd1990$kgeVarsTf[,1:9], method = 'pearson',  use = 'pairwise.complete.obs')
cor(simhyd2000$kgeVarsTf[,1:9], method = 'spearman', use = 'pairwise.complete.obs')

# Histograms of signatures
par(mfrow = c(4,2))
colnames(simhyd1990$kgeVarsTf)
hist(simhyd1990$kgeVarsTf$Area, xlab = "log(Area)", main = NULL)
hist(simhyd1990$kgeVarsTf$Cum.flow, xlab = "log(CF)", main = NULL)
hist(simhyd1990$kgeVarsTf$FDC.low, xlab = "log(FDC.low)", main = NULL)
hist(simhyd1990$kgeVarsTf$FDC.mid, xlab = "log(FDC.mid)", main = NULL)
hist(simhyd1990$kgeVarsTf$FDC.high, xlab = "log(FDC.high)", main = NULL)
hist(simhyd1990$kgeVarsTf$AC, xlab = "AC", main = NULL)
hist(simhyd1990$kgeVarsTf$Peaks, xlab = "log(Peaks)", main = NULL)
hist(simhyd1990$kgeVarsTf$Slope, xlab = "log(Slope)", main = NULL)