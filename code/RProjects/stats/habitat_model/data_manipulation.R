#### INSTALL PACKAGES ----
library(tidyverse)
library(tidyr)
library(reshape2)
library(knitr)
library(dplyr)
library(lubridate)
library(data.table)

#### DATA MANIPULATION 
###--------------VOLUME PERCENTAGES----

# read the csv
datafolder = 'C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\'
setwd(datafolder)
gross_sum_tdf = read.csv('processed_summary_data\\Gross_Summary_Vol_tdf.csv', header = T, na.strings = ' ')

head(gross_sum_tdf)

#edit column names to strings
colnames(gross_sum_tdf) = c('Site','EstInt', '0', '1','2','3','4','5','6','7','8','9','10')

#make new percent dataframe
gross_sum_t_perc = gross_sum_tdf

#change Estimated Initial Column to 100%
gross_sum_t_perc %>%
  mutate(EstInt = 100)

# get percent change using for loop
l1 = list(3:13) # 13 columns, starting at column 3
for (i in l1){
  gross_sum_t_perc[,i] = round((gross_sum_tdf[,i]/gross_sum_tdf[,2]*100),1)
}

head(gross_sum_t_perc)

#Subset fires to their own datasets
brian_head = gross_sum_t_perc %>%
  filter(grepl("Brian Head", Site))
dollar_ridge =   gross_sum_t_perc %>%
  filter(grepl("Dollar Ridge", Site))
pole_creek =   gross_sum_t_perc %>%
  filter(grepl("Lake Fork|POCR", Site))
twitchell = gross_sum_t_perc %>%
  filter(grepl('Twitchell|Fish Creek|Shingle Creek', Site))
seeley=   gross_sum_t_perc %>%
  filter(grepl("Seeley", Site))
shingle=  gross_sum_t_perc %>%
  filter(grepl("Shingle", Site))
trail_mountain =   gross_sum_t_perc %>%
  filter(grepl("Trail Mountain", Site))
clay_springs =   gross_sum_t_perc %>%
  filter(grepl("Clay Springs", Site))


# pivot longer (for plots)
brian_head_t = brian_head %>% pivot_longer(cols = EstInt:'10' , names_to = 'Time', values_to = 'Volume')
dollar_ridge_t = dollar_ridge %>%pivot_longer(cols = EstInt:'10' , names_to = 'Time', values_to = 'Volume')
twitchell_t = twitchell %>% pivot_longer(cols = EstInt:'10', names_to = 'Time', values_to = 'Volume')
pole_creek_t = pole_creek %>% pivot_longer(cols = EstInt:'10' , names_to = 'Time', values_to = 'Volume')
seeley_t = seeley %>% pivot_longer(cols = EstInt:'10' , names_to = 'Time', values_to = 'Volume')
shingle_t = shingle %>% pivot_longer(cols = EstInt:'10' , names_to = 'Time', values_to = 'Volume')
clay_springs_t = clay_springs %>% pivot_longer(cols =  EstInt:'10' , names_to = 'Time', values_to = 'Volume')
trail_mountain_t = trail_mountain %>% pivot_longer(cols =  EstInt:'10' , names_to = 'Time', values_to = 'Volume')
gross_sum_t_perc_t = gross_sum_t_perc %>% pivot_longer(cols =  EstInt:'10' , names_to = 'Time', values_to = 'Volume')
gross_sum_t_perc_t$Fire = NA

gross_sum_t_perc_t[c(1:24),4] = 'Brian Head'
gross_sum_t_perc_t[c(25:120),4] = 'Dollar Ridge'
gross_sum_t_perc_t[c(121:372),4]= 'Twitchell'
gross_sum_t_perc_t[c(373:408),4] = 'Pole Creek'
gross_sum_t_perc_t[c(409:480),4] = 'Seeley'
gross_sum_t_perc_t[c(481:504),4] = 'Shingle'
gross_sum_t_perc_t[c(505:564),4] = 'Trail Mountain'
gross_sum_t_perc_t[c(565:624),4] = 'Clay Springs'

# subset to just 5-30 m
hydro_sites <- hydro_sites %>%
  filter(BottomElevation_m >= -30 & BottomElevation_m <= -5)

# make column that identifies by depth band
hydro_sites <- hydro_sites %>%
  mutate(depth_bin = case_when(
    BottomElevation_m <= -5 & BottomElevation_m >= -10 ~ 1,
    BottomElevation_m < -10 & BottomElevation_m >= -20 ~ 2,
    BottomElevation_m < -20 & BottomElevation_m >= -30 ~3
  ))

# make depth and site factors so you can use summarise
hydro_sites$site <- as.factor(hydro_sites$site)
hydro_sites$depth_bin <- as.factor(hydro_sites$depth_bin)

# average by site and depth
sum_hydro <- hydro_sites %>% group_by(site, depth_bin) %>% summarise(
                                       avg_vegheight = mean(PlantHeight_m, na.rm=T),
                                       se_vegheight = std.error(PlantHeight_m, na.rm=T),
                                       avg_vegcover = mean(PercentCoverage, na.rm=T),
                                       se_vegcover = std.error(PercentCoverage, na.rm=T))

# filter out files where BottomTypeStatus is invalid
bot <- hydro_sites %>%
  filter(BottomStatus == "Valid")

# count by bottom type, depth, and site 
count_bot <- bot %>% group_by(site, depth_bin, BottomType) %>% tally()
count_depth <- bot %>% group_by(site, depth_bin) %>% tally()

# create an intermediate dataframe for calculating percent bottom type by depth at each site
inter_bot <- merge(x=count_bot, y=count_depth, by = c('site', 'depth_bin'))
inter_bot <- inter_bot %>% mutate(bottom_type_perc = (n.x/n.y))

# add the intermediate dataframe to the dataframe with vegetation averages
sum_hydro <- merge(x=inter_bot, y=sum_hydro, by = c('site', 'depth_bin'))

# clean up names
names(sum_hydro)[4:5] <- c("num_type", "num_depth")

# change bottom types from numerical to text
sum_hydro <- sum_hydro %>%
  mutate(BottomType = case_when(
    BottomType == "1" ~ "Clay",
    BottomType == "2" ~ "Sand",
    BottomType == "3" ~ "Silt",
    BottomType == "4" ~ "Rock"
    ))

# delete un-neccesary columns
sum_hydro <- subset(sum_hydro, select = c(-num_type, -num_depth))

# turn bottom type from long form to wide form so for each site/depth_bin combo there is one row with all the percents
sum_hydro <- spread(sum_hydro, key= "BottomType", value= bottom_type_perc)

# change NAs to 0s
sum_hydro[is.na(sum_hydro)] <- 0

# change names of new columns
names(sum_hydro)[8:11] <- c("perc_clay", "perc_rock", "perc_sand", "perc_silt")

# re-order...
sum_hydro <- sum_hydro[ , c(1,7,2,8:11,3:6)]

# force depth_bin to as.numeric for later joins
sum_hydro$depth_bin <- as.numeric(sum_hydro$depth_bin)

## create 3 datasets, one for each month and then conjoin them

# for month = 5, set vegetation to NA
# IMPORTANT= these steps must be done in precise order because they continuously reuse sum_hydro and then reassign it
sum_hydro$month <- 5
s1 <- sum_hydro
s1 <- s1 %>% mutate(avg_vegheight = replace(avg_vegheight, avg_vegheight < 100, NA)) %>% mutate(se_vegheight = replace(se_vegheight, se_vegheight < 100, NA)) %>% mutate(avg_vegcover = replace(avg_vegcover, avg_vegcover < 101, NA)) %>% mutate(se_vegcover = replace(se_vegcover, se_vegcover < 100, NA))

sum_hydro$month <- 6
s2 <- sum_hydro

sum_hydro$month <- 7
s3 <- sum_hydro

sum_hydro <- rbind(s1, s2, s3)

#### DATA MANIPULATION SCUBA ----

# recode depth so it matches sum_hydro
scuba_R <- scuba_R %>% mutate(depth_bin = case_when(depth == "20" ~ "2", 
                                                depth == "10"~"1",
                                                depth == "5" ~ "1",
                                                depth == "0" ~ "0"))

#delete unnecessary things
scuba_R <- subset(scuba_R, select = c(-date, -divers))

#name clean up
names(scuba_R)[11:14] <- c("veg_low", "veg_mid", "veg_high", "clad")

#make a rock column
scuba_R <- scuba_R %>% mutate(perc_rock = ((bedrock + boulder + cobble + gravel + smallgravel)/100))

#### DATA MANIPULATION LIMNOLOGY DATA ----

names(limno)[c(1:7)] <- c("site", "surface_temp", "clarity", "conductivity", "weather", "month", "date")
limno <- subset(limno, select = -weather)
limno <- limno %>% drop_na(site)

# aggregate by month and site to get to 27 obs = 9 sites x 3 months (less messy), larval tow samples currently messing it up
limno <- aggregate(cbind(surface_temp, clarity, conductivity)~ month + site, data = limno, FUN= mean)

# add depth_bin column always equal to 1 for the master dataset
limno$depth_bin <- as.numeric(1)

#### DATA MANIPULATION GEAR DATA ----

#subset gear_shp to pertinent columns
gear_shp <- subset(gear_shp, select = c(-surface_temp, -clarity, -conductivity, -lat, -long))

#filter gear_shp data to just gill and hoop nets
gear_shp <- filter(gear_shp, net == "GIL10" | net == "GIL20" | net == "GIL30" | net == "HOPS")

# turn date time cells from characters into POSIX
gear_shp$date_time_set <- mdy_hm(gear_shp$date_time_set)
gear_shp$date_time_pull <- mdy_hm(gear_shp$date_time_pull)

#### DATA MANIPULATION FISH DATA ----

fecund_fish <- fecund_fish %>% 
  mutate(depth_bin = case_when(
    net == "GIL10" ~ "1",
    net == "HOPS" ~ "0",
    net == "GIL20" ~ "2",
    net == "GIL30" ~ "3"))

# make a matureall_cpue and matureall_num column
names(fecund_fish)[4:5] <- c("matureF_num", "matureM_num")
fecund_fish <- fecund_fish %>% mutate(matureallCpue = matureFCpue + matureMCpue, matureAll_num = matureM_num + matureF_num)

# get rid of HOPS for now
# fish_dat <- fish_dat %>% filter(net != "HOPS")

# make a fecund_fish with a net type column and a column with net type as a binary
# fecund_fish <- fecund_fish%>% mutate(net_type = substring(fecund_fish$net, first = 1, last = 3))
# fecund_fish <- fecund_fish %>% mutate(net_bin = case_when(
                      #net_type == "GIL" ~ "1",
                      #net_type == "HOP" ~ "2"))

#### DATA MANIPULATION ADD GEAR TO FISH DATA ----

# add fish data to gear data (only fecund fish) by site month and net
fish_gear <- merge(x=gear_shp, y=fecund_fish, by = c('site', 'month', 'net'))

# add effort hours column
fish_gear$effort_hours <- as.numeric(fish_gear$date_time_pull - fish_gear$date_time_set)

# re-order...
fish_gear <- fish_gear[, c(1,2,12,3,5,6,4,7,15,8,9,14,10,11,13)]

# force depth bin to numeric for future joins
fish_gear$depth_bin <- as.numeric(fish_gear$depth_bin)

#### DATA MANIPULATION MY TEMP LOGGER DATA ----

# delete weird extra columns
temp_loggers_master <- temp_loggers_master[,1:7]

# turn date time from character to POSIXct
temp_loggers_master$date_time <- mdy_hm(temp_loggers_master$date_time)

# force depth_bin to be numeric for the join
temp_loggers_master$depth_bin <- as.numeric(temp_loggers_master$depth_bin)

# use a non-equi join in data.table to filter temp_loggers master by month, depth_bin, site and month for the time
# intervals specified in fecundonly_habitat
setDT(fish_gear)
setDT(temp_loggers_master)

test <- fish_gear[temp_loggers_master, on=.(month, site, depth_bin), nomatch=0, allow.cartesian = TRUE]
test2 <- test[date_time_set <= date_time & date_time <= date_time_pull]
test3 <- aggregate(data= test2, temp ~ depth_bin + month + site, FUN=mean)
test4 <- merge(x=test3, y=test2, by = c('site', 'month', 'depth_bin'))
test5 <- subset(test4, select = c(-temp.y, -date_time)) # it is important to delete date_time to reduce the dataset to just unique values, you will
                                                        # retain date_time_set and pull
test6 <-test5 %>% distinct(.keep_all= TRUE)
names(test6)[4] <- c("temp_avg")

## HOLD
# for the following analysis, change depth_bin of 0 to 1
# temp_loggers_master['depth_bin'][temp_loggers_master['depth_bin'] == 0] <- 1
## HOLD

#### DATA MANIPULATION TEMP LOGGER GARY'S DATA ----

# turn date time from character to POSIXct
master_fromGPT$date_time <- mdy_hm(master_fromGPT$date_time)

# force depth_bin to be numeric for the join
master_fromGPT$depth_bin <- as.numeric(master_fromGPT$depth_bin)

# use a non-equi join in data.table to filter master_fromGPT by month and depth_bin for the time
# intervals specified in fecundonly_habitat
setDT(fish_gear)
setDT(master_fromGPT)

int <- fish_gear[master_fromGPT, on=.(month, depth_bin), nomatch=0, allow.cartesian = TRUE]
int2 <- int[date_time_set <= date_time & date_time <= date_time_pull]
int3 <- aggregate(data= int2, temp ~ depth_bin + month + site, FUN=mean)
int4 <- merge(x=int3, y=int2, by = c('month', 'depth_bin', 'site'))
int5 <- subset(int4, select = c(-temp.y, -date_time, -depth)) # it is important to delete date_time and depth to reduce the dataset 
                                                              # to just unique values, you will retain date_time_set and pull
                                                              # unfortunately you lose the depths of the temp loggers (which is what that depth is)
                                                              # so it is important to note for GPT temp data: bin 0 = pendant @ 3 m
                                                                                                          #  bin 1 = pendant @ 5 m
                                                                                                          #  bin 2 = pendant @ 13 and 15 m
                                                                                                          #  bin 3 = pendant @ 25, 27 and 30 m
int6 <-int5 %>% distinct(.keep_all= TRUE)
names(int6)[4] <- c("temp_avg")

#### ATTACH FISH DATA TO ALL OF THE DIFFERENT TYPES OF HABITAT DATA ----
hydro_fish <- merge(x=sum_hydro, y=fish_gear, by = c("site", "month", "depth_bin"))

# not really using this so I have largely not "maintained" the dataset or checked this
# if you want to use this need to re-examine data manipulation section and re-run this
scuba_fish <- merge(x=scuba_R, y=fish_gear, by = c("site", "depth_bin")) 

# only have fish from GIL10 here, these variables are highly correlated with better
# variables so I don't think I will use them
# if I do use them may be worth trying to average catch or replicate catch like in the 
# hydro acoustics data across the depth bins
limno_fish <- merge(x=limno, y=fish_gear, by = c("site", "month", "depth_bin")) 
limno_fish <- limno_fish[,c(1:3, 7:12, 4:6, 13:18)]

temp_fish <- int6[,c(3,1,2,5,4,6:16)]

#### MAKE A MASTER DATASET ----
# join hydro acoustics [substrate and vegetation], limnology, temperature and fish data
df_list <- list(sum_hydro, limno, temp_fish)
master <- df_list %>% reduce(full_join, by=c("site", "depth_bin", "month"))
master <- master %>% drop_na(net)

# make a net type column
master <- master %>% mutate(net_type = case_when(net == "GIL10" ~ "GIL", 
                                                    net == "GIL20"~"GIL",
                                                    net == "GIL30" ~ "GIL",
                                                    net == "HOPS" ~ "HOP"))
master <- master [, c(1:3,15,28, 17:27, 4:14, 16)]

#### LINEAR MODELS ----

# linear model with subsetting
model_1 <- lm(matureallCpue~bottom_type_perc, data=hydro_fish, subset=(BottomType == "4"))
summary(model_1)
plot(model_1)

#linear model without
model_2 <- lm(matureallCpue~avg_vegcover, data=hydro_fish)
summary(model_2)
plot(model_2)

### log transforming
# X
hydro_fish$lnCPUE <- log(hydro_fish$matureallCpue + 1)
fish_datE$lnCPUE <- log(fish_datE$matureallCpue + 1)
scuba_fish$lnCPUE <- log(scuba_fish$matureallCpue + 1)

#Y
hydro_fish$lnX <- log(hydro_fish$bottom_type_perc + 1)

# linear model with subsetting
model_1 <- lm(lnCPUE~bottom_type_perc, data=hydro_fish, subset=(BottomType == "1"))
summary(model_1)
plot(model_1)

#linear model without
model_2 <- lm(matureallCpue~avg_vegcover, data=hydro_fish)
summary(model_2)
plot(model_2)

### log transforming
# X
hydro_fish$lnCPUE <- log(hydro_fish$matureallCpue)

#Y
hydro_fish$lnX <- log(hydro_fish$avg_vegcover)
