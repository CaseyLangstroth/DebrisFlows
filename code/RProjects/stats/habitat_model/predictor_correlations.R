### EXPLORING CORRELATIONS BETWEEN PREDICTORS----

### INSTALL PACKAGES ----
library(ggplot2)
library(GGally)

### SCATTER PLOT MATRICES (continuous) WITH BOTTOM TYPE COLORED (using ggpairs) ----
cor_matrix_scat <- ggpairs(master[,c(8, 17:21, 23, 25:28)]) + theme_bw()
cor_matrix_scat

### BOXPLOTS FOR CATEGORICAL VS CONTINUOUS (bottom types in separate plots) ----
#note that substrate and vegetation can not be done for month because they were only collected in 1 month, and substrate doesn't change by month
#surface temp, conductivity and clarity can not be done with depth bin because they were only measured at one depth

## CALCULATE CORRELATION VALUES ----

# bottom types (depth bins)
silt_dep_cor <- cor(master$depth_bin, master$perc_silt, method="spearman", use="complete.obs")
sand_dep_cor <- cor(master$depth_bin, master$perc_sand, method="spearman", use="complete.obs")
clay_dep_cor <- cor(master$depth_bin, master$perc_clay, method="spearman", use="complete.obs")
rock_dep_cor <- cor(master$depth_bin, master$perc_rock, method="spearman", use="complete.obs")

# average vegetation height (depth bins)
vegH_dep_cor <- cor(master$depth_bin, master$avg_vegheight, method = "spearman", use="complete.obs")

# average vegetation cover (depth bins)
vegC_dep_cor <- cor(master$depth_bin, master$avg_vegcover , method = "spearman", use="complete.obs")

# limnology variables (month)
surfT_mon_cor <- cor(master$month, master$surface_temp, method = "spearman", use="complete.obs")
clar_mon_cor <- cor(master$month, master$clarity, method = "spearman", use="complete.obs")
cond_mon_cor <- cor(master$month, master$conductivity, method = "spearman", use="complete.obs")

# temperature (depth bin and month)
tempA_mon_cor <- cor(master$month, master$temp_avg, method = "spearman", use="complete.obs")
tempA_dep_cor <- cor(master$depth_bin, master$temp_avg, method = "spearman", use="complete.obs")

comp_cor <- list(silt_dep_cor, sand_dep_cor, clay_dep_cor, rock_dep_cor, vegH_dep_cor, vegC_dep_cor, 
                 surfT_mon_cor, clar_mon_cor, cond_mon_cor, tempA_mon_cor, tempA_dep_cor)
comp_cor

## BOTTOM TYPES VS. DEPTH BINS ----

# silt
silt_dep_box <- ggplot(master, mapping = aes(x = as.factor(depth_bin), y = perc_silt)) + 
  geom_boxplot() +
  labs( x= "Depth Bin", y= "% Silt") +
  annotate("text", x=1.5, y=0.65, size=4, label="Sp. Cor = 0.75" ) +
  theme_bw()
silt_dep_box

# sand
sand_dep_box <- ggplot(master, mapping = aes(x = as.factor(depth_bin), y = perc_sand)) + 
  geom_boxplot() +
  labs( x= "Depth Bin", y= "% Sand") +
  annotate("text", x=2, y=0.9, size=4, label="Sp. Cor = 0.26" ) +
  theme_bw()
sand_dep_box

# clay
clay_dep_box <- ggplot(master, mapping = aes(x = as.factor(depth_bin), y = perc_clay)) + 
  geom_boxplot() +
  labs( x= "Depth Bin", y= "% Clay") +
  annotate("text", x=2.5, y=0.75, size=4, label="Sp. Cor = -0.71" ) +
  theme_bw()
clay_dep_box

# rock
rock_dep_box <- ggplot(master, mapping = aes(x = as.factor(depth_bin), y = perc_rock)) + 
  geom_boxplot() +
  labs( x= "Depth Bin", y= "% Rock") +
  annotate("text", x=1, y=0.15, size=4, label="Sp. Cor = 0.76" ) +
  theme_bw()
rock_dep_box

## AVERAGE VEGETATION HEIGHT VS. DEPTH BINS ----
veg_height_box <- ggplot(subset(master, depth_bin != 0), mapping = aes(x = as.factor(depth_bin), y = avg_vegheight)) + 
  geom_boxplot() +
  labs( x= "Depth Bin", y= "Average Vegetation Height (m)") +
  annotate("text", x=2.5, y=0.4, size=4, label="Sp. Cor = -0.43" ) +
    theme_bw()
veg_height_box

## AVERAGE VEGETATION % COVER VS. DEPTH BINS ----
veg_cover_box <- ggplot(subset(master, depth_bin != 0), mapping = aes(x = as.factor(depth_bin), y = avg_vegcover)) + 
  geom_boxplot() +
  labs( x= "Depth Bin", y= "Average % Cover") +
  annotate("text", x=2.5, y=60, size=4, label="Sp. Cor = -0.49" ) +
  theme_bw()
veg_cover_box

## SURFACE TEMP VS. MONTH ----
srf_temp_box <- ggplot(master, mapping = aes(x = as.factor(month), y = surface_temp)) + 
  geom_boxplot() +
  labs( x= "Month", y= "Surface Temp") +
  annotate("text", x=1.5, y=25, size=4, label="Sp. Cor = 0.92" ) +
  theme_bw()
srf_temp_box

## CLARITY VS. MONTH ----
clar_box <- ggplot(master, mapping = aes(x = as.factor(month), y = clarity)) + 
  geom_boxplot() +
  labs( x= "Month", y= "Clarity") +
  annotate("text", x=2.5, y=8, size=4, label="Sp. Cor = 0.33" ) +
  theme_bw()
clar_box

## CONDUCTIVITY VS. MONTH ----
cond_box <- ggplot(master, mapping = aes(x = as.factor(month), y = conductivity)) + 
  geom_boxplot() +
  labs( x= "Month", y= "Conductivity") +
  annotate("text", x=1.5, y=8.4, size=4, label="Sp. Cor = 0.91" ) +
  theme_bw()
cond_box

## TEMPERATURE VS. MONTH ----
# note here that there are 6 less deep sets because the temperatures are cut based on net sets and we only set the net at 30 m 3x instead
# of 9x in the late July set
temp_avg_mon_box <- ggplot(master, mapping = aes(x = as.factor(month), y = temp_avg)) + 
  geom_boxplot() +
  labs( x= "Month", y= "Average Logger Temperature") +
  annotate("text", x=1.5, y=22, size=4, label="Sp. Cor = 0.84" ) +
  theme_bw()
temp_avg_mon_box

## TEMPERATURE VS. DEPTH BIN ----
# same here
temp_avg_dep_box <- ggplot(master, mapping = aes(x = as.factor(depth_bin), y = temp_avg)) + 
  geom_boxplot() +
  labs( x= "Depth Bin", y= "Average Logger Temperature") +
  annotate("text", x=3.5, y=25, size=4, label="Sp. Cor = -0.73" ) +
  theme_bw()
temp_avg_dep_box

## ALL TOGETHER ----

# make a blank for aesthitics
blank <- ggplot() +
  geom_blank()
blank

comp_box_cor <- grid.arrange(clay_dep_box, sand_dep_box, silt_dep_box, rock_dep_box, 
                             veg_height_box, veg_cover_box,temp_avg_dep_box, blank, srf_temp_box, clar_box, cond_box,
                             temp_avg_mon_box,
                             nrow=3, ncol=4)
comp_box_cor
