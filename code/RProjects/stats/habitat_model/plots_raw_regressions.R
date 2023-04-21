#### INSTALL PACKAGES ----
library(tidyverse)
library(tidyr)
library(plotrix)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(knitr)
library(dplyr)
theme_set(theme_cowplot())

#### PLOTS ----

### HYDROACOUSTICS ----

# substrate

clay <- ggplot(master, aes(x=perc_clay, y=matureallCpue)) + 
  geom_point()+
  geom_smooth(method=lm) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  ylim(0,8) +
  scale_color_brewer(palette="Dark2") +
  labs(x= "% Clay (<.0005 mm)", y= "CPUE of Fecund Fish (fish/hour)", color = "Month")
clay 

sand <- ggplot(master, aes(perc_sand, y=matureallCpue)) + 
  geom_point()+
  geom_smooth(method=lm) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  ylim(0,8) +
  scale_color_brewer(palette="Dark2") +
  labs( x= "% Sand (0.05-2 mm)", y= "CPUE of Fecund Fish (fish/hour)")
sand

silt<- ggplot(master, aes(x=perc_silt, y=matureallCpue)) + 
  geom_point()+
  geom_smooth(method=lm) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  ylim(0,8) +
  scale_color_brewer(palette="Dark2") +
  labs( x= "% Silt (.0005 - 0.05 mm)", y= "CPUE of Fecund Fish (fish/hour)")
silt

rock <- ggplot(master, aes(x=perc_rock, y=matureallCpue)) + 
  geom_point()+
  geom_smooth(method=lm) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  ylim(0,8) +
  scale_color_brewer(palette="Dark2") +
  labs( x= "% Rock (>2 mm)", y= "CPUE of Fecund Fish (fish/hour)")
rock

grid.arrange(clay, sand, silt, rock, nrow=2, ncol=2)

# vegetation

veg_height <- ggplot(master, aes(x=avg_vegheight, y=matureallCpue)) + 
  geom_point()+
  geom_smooth(method=loess, color="turquoise4") +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Average Vegetation Height (m)", y= "CPUE of Fecund Fish (fish/hour)")
veg_height

veg_cover <- ggplot(master, aes(x=avg_vegcover, y=matureallCpue)) + 
  geom_point()+
  #scale_y_continuous(trans="log2") +
  geom_smooth(method=lm, color="turquoise4") +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Average Vegetation Cover (%)", y= "CPUE of Fecund Fish (fish/hour)")
veg_cover

### DEPTH BIN ----
depth_bin <- ggplot(master, aes(x=as.factor(depth_bin), y=matureallCpue)) + 
  geom_boxplot()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Depth Bin", y= "CPUE of Fecund Fish (fish/hour)")
depth_bin

### MONTH ----
month <- ggplot(master, aes(x=as.factor(month), y=matureallCpue)) + 
  geom_boxplot() +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12, angle = 90),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Month", y= "CPUE of Fecund Fish (fish/hour)")
month

### DATE ----
date <- ggplot(master, aes(x=date_time_set, y=matureallCpue)) + 
  geom_point()+
  geom_jitter() +
  geom_smooth(method=loess, color = "turquoise4", se=F)+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12, angle = 90),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Date", y= "CPUE of Fecund Fish (fish/hour)") +
  facet_wrap(~depth_bin)
date

### NET TYPE ----
net_type <- ggplot(master, aes(x=net_type, y=matureallCpue)) + 
  geom_boxplot()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12, angle = 90),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Net Type", y= "CPUE of Fecund Fish (fish/hour)")
net_type

### SCUBA ----

# substrate
bedrock <- ggplot(scuba_fish, aes(x=bedrock, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% bedrock", y= "CPUE of Fecund Fish (fish/hour)")
bedrock

boulder <- ggplot(scuba_fish, aes(x=boulder, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% boulder", y= "CPUE of Fecund Fish (fish/hour)")
boulder

cobble <- ggplot(scuba_fish, aes(x=cobble, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% cobble", y= "CPUE of Fecund Fish (fish/hour)")
cobble

gravel <- ggplot(scuba_fish, aes(x=gravel, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% gravel", y= "CPUE of Fecund Fish (fish/hour)")
gravel

smallgravel <- ggplot(scuba_fish, aes(x=smallgravel, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% smallgravel", y= "CPUE of Fecund Fish (fish/hour)")
smallgravel

sand <- ggplot(scuba_fish, aes(x=sand, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% sand", y= "CPUE of Fecund Fish (fish/hour)")
sand

clay_fines <- ggplot(scuba_fish, aes(x=clay_fines, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% Silt", y= "CPUE of Fecund Fish (fish/hour)")
clay_fines

hard_clay <- ggplot(scuba_fish, aes(x=hard_clay, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% Clay", y= "CPUE of Fecund Fish (fish/hour)")
hard_clay

perc_rock <- ggplot(scuba_fish, aes(x=perc_rock, y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "% Rock", y= "log(CPUE of Fecund Fish +1) (fish/hour)")
perc_rock

# vegetation/cladophora (algae)
veg_low <- ggplot(data=scuba_fish, aes(x= veg_low , y= lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "<.75 mm vegetation", y= "CPUE of Fecund Fish (fish/hour)")
veg_low

veg_mid <- ggplot(data=scuba_fish, aes(x= veg_mid , y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= ".75-100 mm vegetation", y= "CPUE of Fecund Fish (fish/hour)")
veg_mid

veg_high <- ggplot(data=scuba_fish, aes(x= veg_high , y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "100 mm+ vegetation", y= "CPUE of Fecund Fish (fish/hour)")
veg_high

clad <- ggplot(data=scuba_fish, aes(x= clad , y=lnCPUE)) + 
  geom_point()+
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Cladophora", y= "CPUE of Fecund Fish (fish/hour)")
clad

### LIMNOLOGY ----

# clarity
clarity <- ggplot(data=master, aes(x= clarity , y= matureallCpue)) + 
  geom_point()+
  #geom_smooth(method=loess) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Water Transparency (m)", y= "CPUE of Fecund Fish (fish/hour)")
clarity

# conductivity
conduc <- ggplot(master, aes(x= conductivity , y= matureallCpue)) + 
  geom_point()+
  #geom_smooth(method=lm) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Conductivity (mg/L)", y= "CPUE of Fecund Fish (fish/hour)")
conduc

# temperature-from Gary's data
temp <- ggplot(master, aes(x= temp_avg , y= matureallCpue)) + 
  geom_point()+
  #geom_smooth(method=lm) +
  theme(plot.title = element_text(size=16, hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) +
  labs( x= "Average Temperature (C)", y= "CPUE of Fecund Fish (fish/hour)")
temp
