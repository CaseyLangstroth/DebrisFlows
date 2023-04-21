#### INSTALL PACKAGES ----
library(compositions)
library(robustbase)

#### FOLLOWING THE VIGNETTE FOR THE COMPOSITIONS PACKAGE ----
# http://www.stat.boogaart.de/compositions/UsingCompositions.pdf

# dataset with just compositional data
sub <- master[,c(17:20)]

# view data as compositions
comps <- rcomp(sub) 

# ternary diagrams, can only be used to understand 3 parts, data must be generated using acomp or rcomp (not aplus or rplus)
plot(comps) 
plot(mean(comps), pch=40, add=T, col="red")
ellipses(mean(comps), var(comps), col="blue", r=2)
straight(mean(comps), princomp(comps)$Loadings)

# view 3+ parts in a scatterplot matrix of ternary diagrams that pits 2 components against some margin of the other components
plot(acomp(sub), margin="rcomp")

# view data as amounts
amounts <- aplus(sub) 

# scatterplot matrix in log-scale
plot(amounts) 

mean(comps)
mean(amounts)

# show ratio of amounts
boxplot(amounts)

## relative amounts

# set as acomp to see ratios of amounts in log-geometry (more symmetric)
# rcomp simply displays relative portion itself, susceptible to extreme skew
boxplot(rcomp(amounts), dots=T)

# linear model with CPUE as response and composition as regressors
x <- rcomp(master, c("perc_silt", "perc_sand", "perc_rock", "perc_clay"))
y <- master[, 16]
comp_lm <- lm(y~idt(x), data=list(y=y, x=x))
summary(comp_lm)

#### FOLLOW COMPOSITIONS ANALYSIS BOOK CH 5.1 ON MY DATA ----

### DATA MANIPULATION AND VISUALIZATION ----

# dataset without NAs in substrate data (NAs come from hoop nets that are too shallow to be covered by hydro acoustics)
master_noNA <- subset(master, master$net_type != "HOP")
X <- acomp(master_noNA[-c(19,25,34,46),17:20])
Y <- master_noNA$matureallCpue
Y <- Y[-c(19,25,34,46)]
Y <- log1p(Y)

# visualize dependence of variables in composition
opar <- par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), oma=c(3,3,0,0))
for(i in 1:4){
  for(j in 1:4){
    plot(log(X[,i]/X[,j]), Y, pch=ifelse(i!=j,19,""))
    if(i==j){text (x=0, y=mean(range(Y)),
                  labels=colnames(X)[i],cex=1.5)}
  }
}
  
mtext(text=c("pairwise log-ratio", "dependent variable"),
      side=c(1,2), at=0.5, line=2, outer=TRUE)

### THE MODEL (starting on pg 107 in the book) ----

model <- lm(Y~ilr(X))
model

a <- coef(model)[1]
a

b <- ilrInv(coef(model)[-1], orig=X)
b

# PLOT MODEL ----
## I DO NOT UNDERSTAND THESE ##

plot(X)
straight(mean(X), b, lwd = 2, col = "black", lty=2)

myY <- pretty(Y)
refX <- mean(X) + ((myY - a)/norm(b)^2) * b
plot(refX, add= TRUE, pch = 19)

# parameter confidence regions
varb <- ilrvar2clr(vcov(model) [-1,-1])
varb

par(mar=c(3,3,0,1))
names(b) <- colnames(X)
scaleB <- 1

# all of the compositions with their centroids plotted
plot(scaleB*b) 

# add the dot for the specified substrate types
plot(0*b, add=TRUE, pch=20)

# create 95% CI ellipses
# a black dot outside the 95% centroid indicates the data is significantly different from 0
alpha <- 0.05
rF <- sqrt(qf(1-alpha, nrow(varb)-1, model$df))
ellipses(scaleB*b, scaleB^2*varb, rF)

### CHECK THE MODEL ----

# plot of predicted versus observed
opar <- par(mfrow = c(1,1), mar=c(1,1,1,1))
Predicted=predict(model)
plot(Predicted, Y) #predicted values from the model
abline(0,1) #line of perfect correlation
par(opar)

# Pearson's correlation value
cor(Y, Predicted, method="pearson") 

# centered boxplots of all 3 variables to look at residuals vs. predicted
Residuals <- model$residuals
par(mfrow = c(1,1), mar=c(1,1,1,1))
boxplot(list(Y-mean(Y), Predicted-mean(Y), Residuals))

# model checks are also embedded in the model summary and anova tables
# "the individual ilr coordinates per se have seldom a direct interpretation, and their tests should
# be ignored" - quote from the book
summary(model)
anova(model)

# predicted vs. residuals plot
plot(Predicted, Residuals)

# all the model plots
par(mfrow=c(2,2))
plot(model, add.smooth=FALSE)

# check the normality assumption
par(mfrow=c(1,1))
StudentizedResiduals <- rstudent(model)
qqnorm(StudentizedResiduals) # looks the same as the qqplot above...?
qqline(StudentizedResiduals)
shapiro.test(StudentizedResiduals)

# if p-value below a critical level (say 0.05) suggests non-normality
boxplot(StudentizedResiduals)

# check constant variance
# if this is outside p=0.05 indicates heteroscedasticity (non-constant variance)
cor.test(Predicted, abs(StudentizedResiduals), method="spearman") 


# Robust Regression-study presence of outliers and minimize their influence
modelRob <- lmrob(Y~ilr(X))
summary(modelRob)
