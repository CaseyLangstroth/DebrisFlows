### INSTALL PACKAGES ----
library(DHARMa)
library(glmmTMB)
library(compositions)
library(robustbase)
library(car)
library(effects)
library(MuMIn)

## USING PACKAGE glmmTMB FOR ALL MODELS ----

# look at raw distribution of data
hist(master$matureAll_num, breaks=20)
hist(master$matureallCpue, breaks=20)


### % SILT (CONTINUOUS) GLM ----

# subset data just for silt
master_silt <- subset(master, perc_silt != "NA")

# the model
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

# can plot residuals against any predictor, will plot categorical or small datasets as boxplots
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

### TEST MODEL ASSUMPTIONS USING PACKAGE DHARMa (depth only) ----

## set up display
par(mfrow=c(1,1))

# over-dispersion
testDispersion(M2)

# calculate residuals, should range from 0.25-0.75 homogeneously, positive results are bad and highlighted in red, indicate heteroskedasticity/outliers
simulationOutput <- simulateResiduals(fittedModel = M2, plot = T)

# can plot residuals against any predictor, will plot categorical or small datasets as boxplots
plotResiduals(simulationOutput, form = master$depth_bin)

# uniformity
testUniformity(M2)

# outliers
testOutliers(M2)

# zero-inflation
testZeroInflation(M2)

### ALL SUBSTRATES AND DEPTH GLMM, ALL MONTHS ----

# data manipulation
master_noNA <- subset(master, master$net_type != "HOP")

### OUTLIERS THAT CAN BE REMOVED ####
#master_2 <- master_noNA[-c(19,25,34,46),]

sub_ilr <- acomp(master_noNA[,17:20])

table(master_noNA$matureAll_num)

# the model
options(na.action = "na.fail")
M3 <- glmmTMB(matureAll_num ~ ilr(sub_ilr) + as.factor(depth_bin) + offset(log(effort_hours)),
              zi = ~ilr(sub_ilr) + as.factor(depth_bin),
              family = poisson(link="log"),
              data = master_noNA)
summary(M3)

### TEST MODEL ASSUMPTIONS USING PACKAGE DHARMa (all substrates with depth) ----

# over-dispersion
testDispersion(M3)

# calculate residuals, should range from 0.25-0.75 homogeneously, positive results are bad and highlighted in red, indicate heteroskedasticity/outliers
simulationOutput <- simulateResiduals(fittedModel = M3, plot = T)

# can plot residuals against any predictor, will plot categorical or small datasets as boxplots
plotResiduals(simulationOutput, form = master$depth_bin)

# uniformity
testUniformity(M3)

# outliers
testOutliers(M3)

# zero-inflation
testZeroInflation(M3)

#diagnoses potential causes of convergence problems
diagnose(M3)

#variance and covariance matrix structure
vcov(M3)

#model selection
MuMIn :: dredge(M3)

### SUBSTRATE, VEGETATION COVER, DEPTH GLMM, JULY ONLY ----

# data MANIPULATION
master_plant <- subset(master, master$avg_vegcover != "NA")

### OUTLIERS THAT CAN BE REMOVED ####
#master_plant <- master_plant[-c(19,25,34,46),]

sub_ilr_plant <- acomp(master_plant[,17:20])

# the model
options(na.action = "na.fail")
M4 <- glmmTMB(matureAll_num ~ ilr(sub_ilr_plant) + ordered(depth_bin) + avg_vegcover + offset(log(effort_hours)),
              zi = ~ilr(sub_ilr_plant) + ordered(depth_bin) + avg_vegcover,
              family = nbinom2,
              data = master_plant)
summary(M4)

### TEST MODEL ASSUMPTIONS USING PACKAGE DHARMa (all substrates with depth) ----

# over-dispersion
testDispersion(M4)

# calculate residuals, should range from 0.25-0.75 homogeneously, positive results are bad and highlighted in red, indicate heteroskedasticity/outliers
simulationOutput <- simulateResiduals(fittedModel = M4, plot = T)

# can plot residuals against any predictor, will plot categorical or small datasets as boxplots
plotResiduals(simulationOutput, form = master$depth_bin)

# uniformity
testUniformity(M4)

# outliers
testOutliers(M4)

# zero-inflation
testZeroInflation(M4)

#diagnoses potential causes of convergence problems
diagnose(M4)

#variance and covariance matrix structure
vcov(M4)

#model selection
MuMIn :: dredge(M4)


### SUBSTRATE & DEPTH GLMM, JULY ONLY ---- 

# data used is from same GLMM with vegetation included for july only

# the model
options(na.action = "na.fail")
M5 <- glmmTMB(matureAll_num ~ ilr(sub_ilr_plant) + ordered(depth_bin) + offset(log(effort_hours)),
              zi= ~ilr(sub_ilr_plant) + ordered(depth_bin),
              family = nbinom2,
              data=master_plant)
summary(M5)

### TEST MODEL ASSUMPTIONS USING PACKAGE DHARMa (all substrates with depth) ----

# over-dispersion
testDispersion(M5)

# calculate residuals, should range from 0.25-0.75 homogeneously, positive results are bad and highlighted in red, indicate heteroskedasticity/outliers
simulationOutput <- simulateResiduals(fittedModel = M5, plot = T)

# can plot residuals against any predictor, will plot categorical or small datasets as boxplots
plotResiduals(simulationOutput, form = master$depth_bin)

# uniformity
testUniformity(M5)

# outliers
testOutliers(M5)

# zero-inflation
testZeroInflation(M5)

#diagnoses potential causes of convergence problems
diagnose(M5)

#variance and covariance matrix structure
vcov(M5)

#model selection
MuMIn :: dredge(M5)

