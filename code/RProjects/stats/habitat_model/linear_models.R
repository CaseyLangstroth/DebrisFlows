#### LINEAR MODELS ----

### RAW RESPONSE ----

# linear model with subsetting
model_1 <- lm(matureallCpue~bottom_type_perc, data=hydro_fish, subset=(BottomType == "4"))
summary(model_1)
plot(model_1)

# linear model without subsetting
model_2 <- lm(matureallCpue~avg_vegcover, data=hydro_fish)
summary(model_2)
plot(model_2)

### LOG-TRANSFORMED RESPONSE ----

# LOG(X)
hydro_fish$lnCPUE <- log(hydro_fish$matureallCpue + 1)
fish_datE$lnCPUE <- log(fish_datE$matureallCpue + 1)
scuba_fish$lnCPUE <- log(scuba_fish$matureallCpue + 1)

# LOG(Y)
hydro_fish$lnX <- log(hydro_fish$bottom_type_perc + 1)

# linear model with subsetting
model_1 <- lm(lnCPUE~bottom_type_perc, data=hydro_fish, subset=(BottomType == "1"))
summary(model_1)
plot(model_1)

# linear model without subsetting
model_2 <- lm(matureallCpue~avg_vegcover, data=hydro_fish)
summary(model_2)
plot(model_2)

