library(ggplot2)
library(tidyverse)
df_morph = read.csv('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\processed_summary_data\\df_morphology.csv', header = T, na.strings = '')

confinement = df_morph[, c(2,9, 19)]
colnames(confinement) = c('Site', 'inst_r_perc', 'confinement')
confinementplot = ggplot(data = na.omit(confinement), aes(confinement, inst_r_perc)) + geom_point()
confinementplot + ylab('Percent Instantly Removed') + xlab(' ')
ggsave('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\figures\\plots\\confinement_scatter.png', plot = confinementplot, device = 'png')

acc_space = df_morph[,c(2, 9, 22)]
acc_space = cbind(acc_space, confinement[3])
colnames(acc_space) = c('Site', 'inst_r_perc', 'acc_space', 'confinement')
accplot = ggplot(data = na.omit(acc_space), aes(acc_space, inst_r_perc, color = confinement)) + geom_point()
accplot+ ylab('Percent Instantly Removed')
ggsave('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\figures\\plots\\accspace_scatter.png', plot = accplot, device = 'png')

runout = df_morph[,c(2,9, 21)]
runout = cbind(runout, confinement[3])
colnames(runout) = c('Site', 'inst_r_perc', 'runout_perc','confinement')
runout$runout_perc[runout$runout_perc > 100] = 100
runoutplot = ggplot(data = na.omit(runout), aes(runout_perc, inst_r_perc, color = confinement)) + geom_point()
runoutplot
ggsave('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\figures\\plots\\runout_scatter.png', plot = runoutplot, device = 'png')


RF_angle = df_morph[,c(2,9,6)]
RF_angle = cbind(RF_angle, confinement[3])
angleplot = ggplot(data = na.omit(RF_angle), aes(RF_angle, inst_r_perc, color = confinement)) + geom_point()
angleplot
ggsave('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\figures\\plots\\RFangle_scatter.png', plot = angleplot, device = 'png')


sinuosity_pre = df_morph[,c(2,9,13)]
sinuosity_pre = cbind(sinuosity_pre, confinement[3])
sinuosity_preplot = ggplot(data = na.omit(sinuosity_pre), aes(sin_pre, inst_r_perc, color = confinement)) + geom_point()
sinuosity_preplot
#sinuosity_preplot + geom_point(data = na.omit(sinuosity_post), aes(sin_post, inst_r_perc))

sinuosity_post = df_morph[,c(2,9,14)]
sinuosity_post = cbind(sinuosity_post, confinement[3])
sinuosity_postplot = ggplot(data = na.omit(sinuosity_post), aes(sin_post, inst_r_perc, color = confinement)) + geom_point()
sinuosity_postplot

sinuosity_change = cbind(sinuosity_pre, sinuosity_post[3])
sinuosity_change$time = NA
sinchangeplot = ggplot(data = na.omit(sinuosity_change), aes(sin_pre, sin_post)) + geom_point()
sinchangeplot

dabf = df_morph[,c(2,9,24)]
dabf = cbind(dabf, confinement[3])
dabfplot = ggplot(data = na.omit(dabf), aes(dabf_postdep, inst_r_perc, color = confinement)) + geom_point()
dabfplot
