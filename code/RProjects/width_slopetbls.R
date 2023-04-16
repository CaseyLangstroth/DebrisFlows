library(tidyverse)
library(magrittr)
library(data.table)
twitchell_widths = list.files(path = 'C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Twitchell\\channel\\pts', pattern = '*_width.csv') %>%
  map_df(~read_csv(.))

twitchell_slope = list.files(path = 'C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Twitchell\\channel\\pts', pattern = '*_slope.csv') %>%
  map_df(~read_csv(.))


twitchell_longprof = list.files(path ='C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Twitchell\\channel\\long_profiles', pattern = '.csv') %>%
  map_df(~read_csv(.))


twitch_widthplot = ggplot(data = twitchell_widths, aes(ID, Width, color =DF_ID ))+geom_point() + ylim(0, 11) + geom_point(data = widthtest2, color = 'red') + geom_point(data = widthtest3, color = 'blue')
twitch_widthplot

twitch_slopeplot = ggplot(data = twitchell_slope, aes(ID, Slope, color = DF_ID)) + geom_line() + ylim(0,30)
twitch_slopeplot

dr_prewidths = list.files(path = 'C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Dollar_Ridge\\channel\\pts', pattern = '_pts_width.csv') %>%
  map_df(~read_csv(.))

dr_postwidths = list.files(path =  'C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Dollar_Ridge\\channel\\pts', pattern = '_post_width.csv') %>%
  map_df(~read_csv(.))

library(readxl)

Twitch_morph_2011 = read_excel('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Twitchell\\channel\\channel_morphology.xlsx', sheet = '2011')
Twitch_morph_2014 = read_excel('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Twitchell\\channel\\channel_morphology.xlsx', sheet = '2014')
Twitch_morph_2021 = read_excel('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Twitchell\\channel\\channel_morphology.xlsx', sheet = '2021')

T_wid_2011 = Twitch_morph_2011[-c(3, 6,9,12,15,18,21,24),c(1:2, 6)]
T_wid_2014 = Twitch_morph_2014[-c(3, 6,9,12,15,18,21,24),c(1:2, 6)]
T_wid_2021 = Twitch_morph_2021[-c(3, 6,9,12,15,18,21,24),c(1:2, 6)]

T_wid_2011$Year = '2011'
T_wid_2014$Year = '2014'
T_wid_2021$Year = '2021'

T_wid_all = rbind(T_wid_2011, T_wid_2014, T_wid_2021)

T_widthsplot = ggplot(data = na.omit(T_wid_all), aes(ID, avg_width, color = Year, group= DF_ID)) +geom_point() + geom_line()
T_widthsplot


DR_morph_pre = read_excel('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Dollar_Ridge\\channel\\channel.xlsx', sheet = 'pre')
DR_morph_post = read_excel('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets\\raw_site_data\\Dollar_Ridge\\channel\\channel.xlsx', sheet = 'post')

DR_wid_pre = DR_morph_pre[-c(3, 6, 9 , 12), c(1,2,6)]
DR_wid_pre$Time = 'Pre-Fire'
DR_wid_post = DR_morph_post[-c(3, 6, 9 , 12), c(1,2,6)]
DR_wid_post$Time = 'Post-Fire'

DR_wid_pre_full = DR_morph_pre[c(3, 6, 9 , 12), c(1,2,6)]
DR_wid_pre_full$Time = 'Pre-Fire'
DR_wid_post_full = DR_morph_post[c(3, 6, 9 , 12), c(1,2,6)]
DR_wid_post_full$Time = 'Post-Fire'

DR_wid_all = rbind(DR_wid_pre, DR_wid_post)
DR_wid_all_full = rbind(DR_wid_pre_full, DR_wid_post_full)

DR_sin_pre_full =  DR_morph_pre[c(3, 6, 9 , 12), c(1,2,5)]
DR_sin_pre_full$Time = 'Pre-Fire'
DR_sin_post_full = DR_morph_post[c(3, 6, 9 , 12), c(1,2,5)]
DR_sin_post_full$Time = 'Post-Fire'
DR_sin_all_full = rbind(DR_sin_pre_full, DR_sin_post_full)

level_order = c('Pre-Fire', 'Post-Fire')

DR_wid_allplot = ggplot(data = DR_wid_pre, aes(ID, avg_width, color = Time, group = DF_ID)) +geom_point() +geom_point(data = DR_wid_post, aes(color = Time)) + ylim(5, 11) 
DR_wid_allplot

DR_wid_full_plot = ggplot(data = DR_wid_all_full, aes(factor(Time, level = level_order), avg_width, color = DF_ID)) + geom_point(aes(group = DF_ID)) + geom_line(aes(group = DF_ID)) + xlab('Time') + ylab('Average Channel Width (m3)')
DR_wid_full_plot
ggsave(file = 'C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\figures\\plots\\DR_width_fullreach.png', device = 'png', plot = DR_wid_full_plot)
DR_sin_full_plot = ggplot(data = DR_sin_all_full, aes(factor(Time, level = level_order), sinuosity, color = DF_ID)) + geom_point(aes(group = DF_ID)) + geom_line(aes(group = DF_ID)) + xlab('Time') + ylab('Reach Sinuosity')
DR_sin_full_plot
ggsave(file = 'C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\figures\\plots\\DR_sinuosity_fullreach.png', device = 'png', plot = DR_sin_full_plot)
