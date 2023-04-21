#### Building a database
setwd('C:\\Users\\clang\\Documents\\USU\\DF_Areas\\debris_flow_stats\\data\\data_sheets')

#### Connecting and building the db

library(DBI)
df_db = dbConnect(RSQLite::SQLite(), 'df_db')

### Creating tables in db

#DF ID table
dbExecute(df_db, "create table df_ID (
          Site varchar(50) NOT NULL primary key,
          Fire varchar(50) NOT NULL,
          Lat double,
          Long double);")
#DF geospatial  volume estimates
dbExecute(df_db, "create table df_volume(
          DF_ID varchar(5) primary key,
          Site varchar(50) NOT NULL,
          EstInt double,
          t_0 double,
          t_1 double,
          t_2 double,
          t_3 double,
          t_4 double,
          t_5 double,
          t_6 double,
          t_7 double,
          t_8 double,
          t_9 double,
          t_10 double,
          foreign key (Site) references df_ID(Site));")
#df grain size distributions (modelled and observed)
dbExecute(df_db, "create table df_gsd(
          DF_ID varchar(5) primary key,
          Site varchar(50) NOT NULL,
          initialsub_D16phi double,
          initialsub_D50phi double,
          initialsub_D84phi double,
          initialsub_D16mm double,
          initialsub_D50mm double,
          initialsub_D84mm double,
          obssur_D16mm double,
          obssur_D50mm double,
          obssur_D84mm double,
          obssub_D16mm double,
          obssub_D50mm double,
          obssub_D84mm double,
          foreign key (Site) references df_ID(Site));
          ")
#geospatial debris flow morphology
dbExecute(df_db, "create table df_morphology(
          DF_ID varchar(5) primary key,
          site_name varchar(50) NOT NULL,
          runout_L double,
          RF_angle double,
          int_vol double,
          foreign key (site_name) references df_ID(site_name));")
#geospatial reach morphology
dbExecute(df_db, "create table reach_morphology(
          reach_ID integer check(reach_ID in (1,2,3)),
          DF_name varchar(50) primary key,
          DF_ID varchar(5), 
          site_ID varchar(10) generated always as (DF_ID + '_' + reach_ID) stored,
          stream_class char(1) check(stream_class in ('P', 'I')),
          reach_lengthm integer,
          valley_widthm double,
          foreign key (DF_ID) references df_morphology(DF_ID));")
#geospatial/modeled channel morphology
dbExecute(df_db, "create table channel_morphology(
          DF_ID varchar(5) primary key,
          channel_ID varchar(10),
          IMP_BFWm double,
          IMP_BFDm double,
          flow_depthm double,
          D50phi_Snyder double,
          D50mm_Snyder double,
          river_lengthm double,
          sinuosity double,
          usdakm2 double,
          gradient double,
          slope_deg double,
          q2m3 double,
          q5m3 double,
          q10m3 double,
          q50m3 double,
          q100m3 double,
          q500m3 double,
          foreign key (channel_ID) references reach_morphology(site_ID));")
#geospatial subcatchment morphology
dbExecute(df_db, "create table subcatchment_morphology(
          catch_id integer primary key autoincrement,
          DF_ID varchar(5),
          drainage_areakm2 double,
          relief double,
          Cp2008 double,
          Cp2011 double,
          Cp2016 double,
          Sp23 double,
          Mean_Elevation double,
          MgO double,
          Rd double,
          WI double,
          CS double,
          Om double,
          Ro double,
          AnnP double,
          T double,
          K double,
          HC double,
          S23_Areakm2 double,
          Bmh_Areakm2 double,
          Channel_Devlp char(1) check(Channel_Devlp in ('Y','N')),
          foreign key (DF_ID) references reach_morphology(DF_ID));")
### Import csv to dbtables

#first, read csvs
sites = read.csv('Site_Locations.csv', header = T)
df_volumes = read.csv('processed_summary_data//Gross_Summary_Vol_tdf.csv', header = T, na.strings = ' ')
df_gsd = read.csv('processed_summary_data//df_gsd.csv', header = T, na.strings = ' ')
df_morphology = read.csv('processed_summary_data//df_morphology.csv', header = T, na.strings = ' ')
reach_morphology = read.csv('processed_summary_data//reach_morphology.csv', header = T, na.strings = ' ')
channel_morphology = read.csv('processed_summary_data//channel_morphology.csv', header = T, na.strings = ' ')
subcatch_morphology = read.csv('processed_summary_data//subcatchment_morphology.csv', header = T, na.strings = ' ')

# import csv into tables (for completed tables)
dbWriteTable(df_db, 'df_ID', sites, append  = TRUE)
dbWriteTable(df_db, 'df_volume', df_volumes, append = TRUE)
dbWriteTable(df_db, 'df_gsd', df_gsd, append = TRUE)
