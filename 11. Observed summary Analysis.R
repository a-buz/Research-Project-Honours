##Observed summary

setwd("F:\\University\\Honours\\Data Access\\Modelling\\Results")

# obs_summary_1990 <- read.csv("obs.summary_1990.csv", header = TRUE)
# obs_summary_2000 <- read.csv("obs.summary_2000.csv", header = TRUE)
# sim_summary_1990 <- read.csv("sim.summary_1990.csv", header = TRUE)
# sim_summary_2000 <- read.csv("sim.summary_2000.csv", header = TRUE)
# 
# par(mfrow = c(1,2))
# boxplot(obs_summary_1990[,18:23])
# boxplot(sim_summary_1990[,18:23])
# 
# range01 <- function(x){1-abs(x)/(max(x)-min(x))}
# range02 <- function(x){(x - mean(x))/sd(x)}
# 
# ##Finding raw differences in signatures
# sig_1990 <- data.frame(obs_summary_1990[,1:18], 
#                        (sim_summary_1990[,19:21] - obs_summary_1990[,19:21]),
#                        AC = abs(1 - sim_summary_1990[,22]/obs_summary_1990[,22]),
#                        sim_summary_1990[,23:24] - obs_summary_1990[,23:24],
#                        KGext = obs_summary_1990[,25])
# sig_2000 <- data.frame(obs_summary_2000[,1:18], 
#                        (sim_summary_2000[,19:21] - obs_summary_2000[,19:21]),
#                        AC = abs(1 - sim_summary_2000[,22]/obs_summary_2000[,22]),
#                        (sim_summary_2000[,23:24] - obs_summary_2000[,23:24]),
#                        KGext = obs_summary_2000[,25])
# 
# hist(sig_1990$FDC.low)


# sig_1990[,19:24] <- lapply(sig_1990[,19:24], range01)
# sig_2000[,19:24] <- lapply(sig_2000[,19:24], range01)

# obs_summary_1990[,17:22] <- lapply(obs_summary_1990[,17:22], range02)
# obs_summary_2000[,17:22] <- lapply(obs_summary_2000[,17:22], range02)
# sim_summary_1990[,17:22] <- lapply(sim_summary_1990[,17:22], range02)
# sim_summary_2000[,17:22] <- lapply(sim_summary_2000[,17:22], range02)

#write.csv(sig_1990, "sig_1990.csv")
#write.csv(sig_2000, "sig_2000.csv")

######-------------Analysis----------------###################
setwd("F:\\University\\Honours\\Data Access\\Modelling\\Results")

obs_1990 <- read.csv("obs.summary2_1990.csv", header = TRUE)
obs_2000 <- read.csv("obs.summary2_2000.csv", header = TRUE)
obs_2000$Number <- as.character(obs_2000$Number)
obs_1990$Number <- as.character(obs_1990$Number)
obs_1990$KGext <- as.factor(obs_1990$KGext)
obs_2000$KGext <- as.factor(obs_2000$KGext)

##PCA of observed
obs1_1990 <- obs_1990[,c(10:11,17:22)]
obs1_2000 <- obs_2000[,c(10:11,17:22)]

obs1_1990[,c(1:2, 8)] <- log(obs1_1990[,c(1:2, 8)])
obs1_1990$FDC.high <- log(obs1_1990$FDC.high)
obs1_1990$FDC.low <- log(obs1_1990$FDC.low+0.1)
obs1_1990$FDC.mid <- log(obs1_1990$FDC.mid+0.1)

obs1_2000[,c(1:2)] <- log(obs1_2000[,c(1:2)])
obs1_2000$Peaks <- log(obs1_2000$Peaks+0.1)
obs1_2000$FDC.high <- log(obs1_2000$FDC.high)
obs1_2000$FDC.low <- log(obs1_2000$FDC.low+0.1)
obs1_2000$FDC.mid <- log(obs1_2000$FDC.mid+0.1)

# obs1_1990 <- do.call(data.frame,lapply(obs1_1990, function(x) replace(x, is.infinite(x),NA)))
# obs1_2000 <- do.call(data.frame,lapply(obs1_2000, function(x) replace(x, is.infinite(x),NA)))

obs1_1990pc <- prcomp(obs1_1990, scale = T, retx = T, center = T)
obs1_2000pc <- prcomp(obs1_2000, scale = T, retx = T, center = T)

##According the the biplots (1990), the FDCs, cum.flow, x3, r.sq.log, AC are correlated to Viney

par(mfrow = c(1,1))
#pca
library(car)
library(ggfortify)
autoplot(obs1_1990pc, data = obs_1990, colour = "KGext", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(obs1_2000pc, data = obs_2000, colour = "KGext",
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(obs1_1990pc, data = obs_1990, colour = "Type", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)
autoplot(obs1_2000pc, data = obs_2000, colour = "Type",
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

summary(obs1_1990pc)
summary(obs1_2000pc)
##So there does seem to be a monotonic relationship between the Viney and the FDC signatures - clearer after logging. 

###Regression time

obs1_1990$KGext <- obs_1990$KGext
obs1_2000$KGext <- obs_2000$KGext
obs1_1990$Type <- obs_1990$Type
obs1_2000$Type <- obs_2000$Type

cor(obs1_1990[,1:8], method = "pearson")
cor(obs1_1990[,1:8], method = "spearman")

# lm_obs1_1990 <- lm(objFun ~ KGext+x3+x4+Area+rbias+r.sq.log
#                    +AC, data = obs1_1990)
lm_obs1_1990 <- lm(Viney ~ ., data = obs1_1990)
#lm_obs1_1990 <- lm(Viney ~ .-CF-FDC.low-Slope, data = obs1_1990)
lm_obs1_1990 <- lm(Viney ~ Area + FDC.mid + Peaks, data = obs1_1990)
lm_obs1_1990 <- step(lm_obs1_1990, method = "both")
vif(lm_obs1_1990)
summary(lm_obs1_1990)
autoplot(lm_obs1_1990)
anova(lm_obs1_1990)

##Here, many of the variables are very highly correlated, such as the signatures 
#(FDC.low, mid, high, Peaks and slope) > 0.8. Cum.flow is also high correlated to runoff (which makes sense)
#Area, runoff, Peaks are significant. 

cor(obs1_2000[,1:9], method = "pearson", use = "pairwise.complete.obs")
cor(obs1_2000[,1:9], method = "spearman", use = "pairwise.complete.obs")

lm_obs1_2000 <- lm(Viney ~ .-FDC.low-FDC.mid,data = obs1_2000)
lm_obs1_2000 <- lm(Viney ~ Type, data = obs1_2000)
lm_obs1_2000 <- step(lm_obs1_2000, method = "backward")
#library(car)
vif(lm_obs1_2000)
summary(lm_obs1_2000)
autoplot(lm_obs1_2000)
anova(lm_obs1_2000)

#In the 1990 linear regression, Area, runoff, FDC,high, AC, and Peaks are significant
#In the 2000 linear regression, Area, FDC.mid, Peaks are significant. 

#Adjust x3 on 2000

##Story --> does memory affect streamflow? We used signatures
##It does depending on the catchment an the climate. Drier catchments, AC has heavier focus compared to wetter. 
##Correlated signatures. 

##FLOW DURATION CURVE
par(mfrow=c(1,1))
n <- nrow(obs_1990)
sort.flow_1990 <- sort(as.numeric(obs_1990$Viney), decreasing = TRUE, na.last=FALSE)

m <- nrow(obs_2000)
sort.flow_2000 <- sort(as.numeric(obs_2000$Viney), decreasing = TRUE, na.last=FALSE)

#Sort - just a series of numbers now
rank.flow_1990 <- 1:n
Prob_1990 <- rank.flow_1990/(n+1) 

rank.flow_2000 <- 1:m
Prob_2000 <- rank.flow_2000/(m+1)

plot(Prob_1990*100, sort.flow_1990, type = "l", xlab="Probability where Viney is exceeded", 
     ylab="Viney", col = "red", ylim = c(0,1))
lines(Prob_2000*100, sort.flow_2000, col="blue")
legend("topright", col = c("red", "blue"), c("Wet", "Dry"), 
       lty = c(1,1))

#boxplots of signatures
par(mfrow = c(4,2))
colnames(obs1_1990)
hist(obs1_1990$Area, xlab = "log(Area)", main = NULL)
hist(obs1_1990$CF, xlab = "log(CF)", main = NULL)
hist(obs1_1990$FDC.low, xlab = "log(FDC.low)", main = NULL)
hist(obs1_1990$FDC.mid, xlab = "log(FDC.mid)", main = NULL)
hist(obs1_1990$FDC.high, xlab = "log(FDC.high)", main = NULL)
hist(obs1_1990$AC, xlab = "AC", main = NULL)
hist(obs1_1990$Peaks, xlab = "log(Peaks)", main = NULL)
hist(obs1_1990$Slope, xlab = "log(Slope)", main = NULL)


