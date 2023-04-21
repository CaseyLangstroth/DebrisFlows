### INSTALL PACKAGES ----
library(DHARMa)
library(glmmTMB)
library(compositions)
library(robustbase)
library(car)
library(effects)
library(MuMIn)

# look at raw distribution of data
hist(master$matureAll_num, breaks=20)
hist(master$matureallCpue, breaks=20)

### % SILT (CONTINUOUS) GLM ----

# subset data just for silt
master_silt <- subset(master, perc_silt != "NA")

options(na.action = "na.fail")
M1 <- glmmTMB(matureAll_num ~ perc_silt + offset(log(effort_hours)),
              zi= ~perc_silt,
          family = nbinom2,
          data = master_silt)
summary(M1)

### TEST MODEL ASSUMPTIONS USING PACKAGE DHARMa (silt only) ----

## set up display
par(mfrow=c(1,1))

# overdispersion
testDispersion(M1)

# calculate residuals, should range from 0.25-0.75 homogeneously, positive results are bad and highlighted in red, indicate heteroskedasticity/outliers
simulationOutput <- simulateResiduals(fittedModel = M1, plot = T)

#can plot residuals against any predictor, will plot categorical or small datasets as boxplots
plotResiduals(simulationOutput, form = master_silt$perc_rock)

# uniformity
testUniformity(M1)

# outliers
testOutliers(M1)

# zero-inflation
testZeroInflation(M1)

### DEPTH (ORDINAL) GLM ----

# the model
M2 <- glmmTMB(matureAll_num ~ ordered(depth_bin) + offset(log(effort_hours)),
              zi= ~ordered(depth_bin),
          family = nbinom2,
          data = master)
summary(M2)

#use package DHARMa to test for overdispersion, takes a simulation approach
testDispersion(M2)

#QQ plot residuals
M2_resids <- simulateResiduals(M2)
par(mfrow=c(2,1))
plot(M2_resids)

#model selection
MuMIn :: dredge(M2)

# check for dispersion manually
E1 <- resid(M2, type = "pearson")
N <- nrow(master)
k <- length(coef(M2))
Disp <- sum(E1^2/(N-k))
Disp

#check for linearity and overdispersion in plots
par(mfrow = c(2,2)) # see all of the diagnostic plots at once
plot(M2)

help(compositions)
??compositions

### SUBSTRATE AND DEPTH GLM, ALL AVAILABLE DATA ----

##the data
master_noNA <- subset(master, master$net_type != "HOP")
#master_2 <- master_noNA[-c(19,25,34,46),], OUTLIERS
sub_ilr <- acomp(master_noNA[,17:20])

table(master_noNA$matureAll_num)
# the model
options(na.action = "na.fail")
M3 <- glmmTMB(matureAll_num ~ ilr(sub_ilr) + as.factor(depth_bin) + offset(log(effort_hours)),
              zi=~ilr(sub_ilr) + as.factor(depth_bin),
              family = poisson(link="log"),
              data=master_noNA)
summary(M3)

M3_resids <- simulateResiduals(M3)
plotResiduals(M3_resids, form= as.factor(master_noNA$depth_bin))
plotResiduals(M3_resids, form= ilr(sub_ilr))
par(mfrow=c(2,1))
plot(M3_resids)

#use package DHARMa for model checking
par(mfrow= c(1,1))
testDispersion(M3)

#confidence intervals
confint(M3)

#diagnoses potential causes of convergence problems
diagnose(M3)

#variance and covariance matrix structure
vcov(M3)

#QQ plot residuals
M3_resids <- simulateResiduals(M3)
par(mfrow=c(2,1))
plot(M3_resids)

#model selection
MuMIn :: dredge(M3)

### SUBSTRATE, PLANT COVER, DEPTH GLM, JULY ONLY ----

# data
master_plant <- subset(master, master$avg_vegcover != "NA")
#master_plant <- master_plant[-c(19,25,34,46),], OUTLIERS
sub_ilr_plant <- acomp(master_plant[,17:20])

# the model
options(na.action = "na.fail")
M4 <- glmmTMB(matureAll_num ~ ilr(sub_ilr_plant) + ordered(depth_bin) + avg_vegcover + offset(log(effort_hours)),
              zi=~ilr(sub_ilr_plant) + ordered(depth_bin) + avg_vegcover,
              family = nbinom2,
              data=master_plant)
summary(M4)

#use package DHARMa for model checking
par(mfrow= c(1,1))
testDispersion(M4)

#confidence intervals
confint(M4)

#diagnoses potential causes of convergence problems
diagnose(M4)

#variance and covariance matrix structure
vcov(M4)

#QQ plot residuals
M4_resids <- simulateResiduals(M4)
par(mfrow=c(2,1))
plot(M4_resids)

#model selection
MuMIn :: dredge(M4)

### SUBSTRATE & DEPTH GLM, JULY ONLY ---- 

# the model
options(na.action = "na.fail")
M5 <- glmmTMB(matureAll_num ~ ilr(sub_ilr_plant) + ordered(depth_bin) + offset(log(effort_hours)),
              zi= ~ilr(sub_ilr_plant) + ordered(depth_bin),
              family = nbinom2,
              data=master_plant)
summary(M5)

#use package DHARMa for model checking
par(mfrow= c(1,1))
testDispersion(M5)

#confidence intervals
confint(M5)

#diagnoses potential causes of convergence problems
diagnose(M5)

#variance and covariance matrix structure
vcov(M5)

#QQ plot residuals
M5_resids <- simulateResiduals(M5)
par(mfrow=c(2,1))
plot(M5_resids)

#model selection
MuMIn :: dredge(M5)

