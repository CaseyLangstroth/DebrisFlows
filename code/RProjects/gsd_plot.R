gsd = read.csv('..//..//data//data_sheets//processed_summary_data//DF_gsd.csv', header = T, na.strings = '')

library(ggplot2)

d16phi = gsd[,c(1, 3:5)]
d50phi = gsd[,c(1, 6)]
colnames(d50phi) = c('Site', 'size')
d50phi$perc = 50
d84phi = gsd[, c(1, 7)]
colnames(d84phi) = c('Site', 'size')
d84phi$perc = 84

library(tidyverse)
library(magrittr)
d16phi_t = d16phi %>% pivot_longer(cols = 2:4, names_to = 'yr', values_to = 'size')
d16phi_t$perc = 16
d16phi_t = d16phi_t[,-c(2)]

all = rbind(d16phi_t, d50phi, d84phi)
all = mutate(all, mm = 10**(-1*size)/3.322)

gsdplot = ggplot(data =na.omit(all), aes(size,perc, group= Site))+ geom_line() + theme(legend.position = 'none') + scale_x_reverse()
gsdplot + scale_y_continuous(breaks = c(16,50,84))
