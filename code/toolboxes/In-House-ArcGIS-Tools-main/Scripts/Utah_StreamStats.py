"""

##ungaged utah

#origional publication
#https://pubs.usgs.gov/sir/2007/5158/pdf/SIR2007_5158_v4.pdf

all equations used in this code can be found in the above publication
 Drainage area (DRNAREA) was
used in all equations, mean basin elevation (ELEV) was
included in three regions, and the variables mean annual
precipitation (PRECIP), average basin slope (BSLDEM10M),
and percent of basin defined as herbaceous upland
(HERBUPLND) were each used in the equations for one
region

"""
#%% import required modules
import numpy as np
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *

#%%------------------Booleans & License Checkout----------------------------
arcpy.env.overwriteOutput=True
arcpy.CheckOutExtension("Spatial")
usegui=True

#%% user inputs
if not usegui:
    indir1=r'D:\Box Sync\Wildfires\PLI_watersheds\Input_Data'
    indir2=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data'
    basenamein='a000'
    
    riv_net=r"D:\Box Sync\Wildfires\PLI_watersheds\Output_Data\a000\a000_network.shp"# river network shapefile
    DAfield="usarea_km2"
    wtrshed_extent=os.path.join(indir2,basenamein,basenamein+"_watershed.shp")
    DEM=os.path.join(indir2,basenamein,basenamein+"demf.tif")
    
    regionshp=r'C:\Users\Scott\Desktop\PKregions\UT_PK_regions.shp'
    precip_ras=r"D:\Box Sync\Wildfires\PLI_watersheds\Input_Data\Annual_Rainfall\PRISM_ppt_30yr_normal_800mM3_annualrainfall.tif"
    Landcover=os.path.join(indir1,'Landcover','nlcd_2019_land_cover.tif')
    RIin=[2,5,10,25,50,100,200,500] # Recurrence interval of interest options are 2,10,25,50,100,200,500
    basename='a000'
    out_tmp_dir=r"D:\test\Qcalcs"

if usegui:

    riv_net=arcpy.GetParameterAsText(0) # river network shapefile
    DAfield=arcpy.GetParameterAsText(1)# drainage area field
    wtrshed_extent=arcpy.GetParameterAsText(2)
    DEM=arcpy.GetParameterAsText(3)
    
    regionshp=arcpy.GetParameterAsText(4)# geohydro regions shapefile
    precip_ras=arcpy.GetParameterAsText(5)#
    Landcover=arcpy.GetParameterAsText(6)
    RIin=[2,5,10,25,50,100,200,500]  # Recurrence interval of interest options are 2,10,25,50,100,200,500
    basename=arcpy.GetParameterAsText(7)
    out_tmp_dir=arcpy.GetParameterAsText(8)

#%% Internal functions 
# toggle between arcgis interface and command prompt printing
def print2(instr,usegui): 
    if usegui==True:
        arcpy.AddMessage(instr)
    else:
        print(instr)

# Check if a directory exists if not make it    
def chk_mk_dir(input_dir):
    try:
        os.mkdir(input_dir)
    except:
        pass

# cause an error and report a msg
def reporterror(errmsg,usegui):
    if usegui==True:
        arcpy.AddError(errmsg)
    else:
        raise ValueError(errmsg)

# check if a field already exists
def chk_field_exists(infc,field):
    lstFields=arcpy.ListFields(infc)
    if field in lstFields:
        return True
    else:
        return False


#%% Directory managament
chk_mk_dir(out_tmp_dir)

# make sure everything snaps to DEM
arcpy.env.extent=DEM
arcpy.env.snapRaster=DEM

#%% SOLVE FOR CORRECT REGION IN UTAH
tmpwtrshd=basename+'_wtrshdtmp.shp'
wtrshd=os.path.join(out_tmp_dir,tmpwtrshd)

# make a scrap copy of the watershed to append data to
arcpy.FeatureClassToFeatureClass_conversion(wtrshed_extent, out_tmp_dir, 
                                            tmpwtrshd)
# spatial join watershed and regions to find most overlap aka region id
wtr_region=os.path.join(out_tmp_dir,basename+"_wtr_region.shp")
arcpy.analysis.SpatialJoin(wtrshd,regionshp,wtr_region, "JOIN_ONE_TO_ONE",
                           "KEEP_ALL", None, "LARGEST_OVERLAP", None, '')

# Extract the region id from the join
try:
    arr = arcpy.da.FeatureClassToNumPyArray(wtr_region, "Name" ,skip_nulls=True)
except:
    errmsg='Name field likely missing in the regions shapefile'
    reporterror(errmsg,usegui)

Regionid=arr[0][0]# extract sting of interest
Regionstr=Regionid.split('_')# split string to seperate region form id number
Region=int(Regionstr[1])# final value of the region id
del arr, Regionid, Regionstr

# add this back to the watershed shapefile to attribute for error checking
RID="RegionID"
arcpy.AddField_management(wtrshd,RID)
arcpy.CalculateField_management(wtrshd,RID,str(Region))


#%% Calculate zonal stats to parameterize required constants in the equations

# add a field cald UID to use for table joining 
uid_field="UID"
arcpy.AddField_management(wtrshd,uid_field)
arcpy.CalculateField_management(wtrshd,uid_field,"!FID!")

#----------------Compute average elevation------------------------------------
# Zonal stats to get mean elev in the watershed
elevbtble=os.path.join(out_tmp_dir,basename+'elev.dbf')
ZonalStatisticsAsTable(wtrshd, uid_field, DEM,elevbtble, "DATA", "MEAN")
                          
outfieldE='Avg_Elev'
arcpy.AddField_management(elevbtble, outfieldE, "DOUBLE")
arcpy.CalculateField_management(elevbtble, outfieldE, "!MEAN!","PYTHON","")
arcpy.JoinField_management(wtrshd,uid_field,elevbtble,uid_field, outfieldE)

#---------------------Average slope in %---------------------------------------
slp=os.path.join(out_tmp_dir,basename+"_slp.tif")
outSlope = Slope(DEM, "PERCENT_RISE")
outSlope.save(slp)
del outSlope

# Zonal stats to get mean elev in the watershed
slptble=os.path.join(out_tmp_dir,basename+'_slp.dbf')
ZonalStatisticsAsTable(wtrshd, uid_field,slp,slptble, "DATA", "MEAN")
                          
outfieldslp='Avg_Slp'
arcpy.AddField_management(slptble, outfieldslp, "DOUBLE")
arcpy.CalculateField_management(slptble, outfieldslp, "!MEAN!","PYTHON","")
arcpy.JoinField_management(wtrshd,uid_field,slptble,uid_field, outfieldslp)


#----------------Compute average annual precipitation--------------------------
# Zonal stats to get mean elev in the watershed
preciptble=os.path.join(out_tmp_dir,basename+'precip.dbf')
ZonalStatisticsAsTable(wtrshd, uid_field,precip_ras,preciptble, "DATA", "MEAN")
                          
outfieldP='Avg_Precip'
arcpy.AddField_management(preciptble, outfieldP, "DOUBLE")
arcpy.CalculateField_management(preciptble, outfieldP, "!MEAN!","PYTHON","")
arcpy.JoinField_management(wtrshd,uid_field,preciptble,uid_field, outfieldP)

#----------------Compute % of basin defined as hebacious upland----------------
LC_wtshd=os.path.join(out_tmp_dir,basename+"LC.tif")
arcpy.gp.ExtractByMask_sa(Landcover, DEM, LC_wtshd)

herb_rc=os.path.join(out_tmp_dir,basename+"herbaceous_rc.tif")
LCfield="NLCD_Land"
LCherb='\'Herbaceous\''
lcexpress=LCfield+' = '+LCherb
arcpy.gp.Con_sa(LC_wtshd, "1", herb_rc, "0", lcexpress)

outherb="percherb"
herb_tbl=os.path.join(out_tmp_dir,basename+'_prec_conifer.dbf')
ZonalStatisticsAsTable(wtrshd,uid_field, herb_rc, herb_tbl, "NODATA", "MEAN")
arcpy.AddField_management(herb_tbl, outherb, "DOUBLE")
arcpy.CalculateField_management(herb_tbl, outherb, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(wtrshd,uid_field, herb_tbl,uid_field, outherb)


#%%------------------Prepare inputs for equations------------------------------
# Extract the region id from the join
flds=[outfieldE,outfieldslp,outfieldP,outherb]
try:
    arrW = arcpy.da.FeatureClassToNumPyArray(wtrshd, flds ,skip_nulls=True)
except:
    errmsg='Fields are missing after spatial averaging'
    reporterror(errmsg,usegui)

try:
    arrR = arcpy.da.FeatureClassToNumPyArray(riv_net, DAfield ,skip_nulls=True)
except:
    errmsg='Fields are missing after spatial averaging'
    reporterror(errmsg,usegui)

#Extract values to arrays and convert the units

# from the above analysis units are km2 for DA, m for E, mm for precip 
# constants below are to convert them
DA=arrR[DAfield]/2.59# drainage area sq miles
E=arrW[outfieldE]*3.281# average basin elevation in feet
P=arrW[outfieldP]/25.4#*0.0393701# mean annual precip in INCHES 
BSLP=arrW[outfieldslp]# AVERAGE BASIN SLOPE unitis in % 
AHERB=arrW[outherb]# % area of HERBACIOUS UPLAND $nlcd

for RI in RIin:
    #  flag the index of the recurrence intervale of interst
    if RI==2:
        idx=0
    elif RI==5:
        idx=1
    elif RI==10:
        idx=2
    elif RI==25:
        idx=3
    elif RI==50:
        idx=4
    elif RI==100:
        idx=5
    elif RI==200:
        idx=6
    elif RI==500:
        idx=7
    else:
        print2("invalid RI ="+str(RI),usegui)
    
    
    if Region==1:
        print2('Region 1',usegui)
        a=[1.52,5.49,10.3,19.7,29.4,40.4,58.3,85.4]
        b=[0.677,0.614,0.581,0.547,0.524,0.512,0.483,0.457]
        c=[1.39,1.30,1.25,1.21,1.19,1.17,1.15,1.13]
        Qout=a[idx]*(np.power(DA,b[idx])*(c[idx]**(E/1000)))
        
    elif Region==2:
        print2('Region 2',usegui)
        a=[0.585,1.56,2.51,4.0,5.36,6.92,8.79,12.0]
        b=[0.847,0.747,0.703,0.661,0.635,0.613,0.592,0.555]
        c=[1.07,1.07,1.06,1.06,1.06,1.06,1.055,1.05]
        Qout=a[idx]*np.power(DA,b[idx])*(c[idx]**P)
        
    elif Region==3:
        print2('Region 3',usegui)
        a=[14.5,47.6,83.7,148,215,300,411,599]
        b=[0.328,0.287,0.289,0.289,0.302,0.303,0.301,0.299]
        Qout=a[idx]*np.power(DA,b[idx])
        
    elif Region==4:
        print2('Region 4',usegui)
        a=[0.083,0.539,0.753,1.64,2.68,4.18,6.29,10.5]
        b=[0.822,0.816,0.811,0.804,0.798,0.792,0.786,0.778]
        c=[0.656,0.537,0.5,0.414,0.373,0.334,0.299,0.256]
        d=[0.039,0.035,0.032,0.03,0.028,0.023,0.021,0.018]
        Qout=a[idx]*np.power(DA,b[idx])*(2.72**((c[idx]*(E/1000))-(d[idx]*BSLP)))
        
    elif Region==5:
        print2('Region 5',usegui)
        a=[4.32,11.7,18.4,28.8,38.4,50.2,64.7,88.3]
        b=[0.623,0.575,0.555,0.538,0.536,0.515,0.504,0.489]
        c=[0.503,0.425,0.388,0.352,0.331,0.316,0.3,0.285]
        Qout=a[idx]*np.power(DA,b[idx])*np.power(AHERB+1,c[idx])
        
    elif Region==6:
        print2('Region 6',usegui)
        a=[4150,13100,24700,49500,77400,115000,166000,258000]
        b=[0.553,0.479,0.444,0.411,0.391,0.391,0.361,0.344]
        c=[-2.45,-2.44,-2.47,-2.51,-2.54,-2.58,-2.61,-2.65]
        Qout=a[idx]*np.power(DA,b[idx])*np.power(E/1000,c[idx])
        
    elif Region==7:
        print2('Region 7',usegui)
        a=[18.4,67.4,134,278,446,683,1010,1620]
        b=[0.63,0.539,0.487,0.429,0.39,0.355,0.321,0.280]
        Qout=a[idx]*np.power(DA,b[idx])
    else:
        print2("invalid region ",usegui)
    
    # convert Qout from cfs to cms
    Qout=Qout*0.0283168
    
    # Write the output back to the shapefile 
    Qoutfield='Qcms'+str(int(RI))+'yr'
    arcpy.AddField_management(riv_net,Qoutfield,"DOUBLE")
    
    with arcpy.da.UpdateCursor(riv_net, Qoutfield) as cursor: #Add data
        i = 0
        for row in cursor:
            row[0]=Qout[i]
            cursor.updateRow(row)
            i+=1

#Calculate bankfull depth and discharge
# convert DA back to km2
DA=DA*2.59# drainage area sq miles

# Bank Full Width/Depth
#Rocky Mountain West
RMS_W=1.24*np.power(DA,0.435)
RMS_D=0.23*np.power(DA,0.225)

#Intermontain Plateau
IMP_W=1.11*np.power(DA,0.415)
IMP_D=0.07*np.power(DA,0.329)

# # convert ft to meters
# RMS_W=RMS_W*0.3048
# RMS_D=RMS_D*0.3048
# IMP_W=IMP_W*0.3048
# IMP_D=IMP_D*0.3048



flds=["RMS_BFw_m","RMS_BFd_m","IMP_BFw_m","IMP_BFd_m"]


for fld in flds:
    arcpy.AddField_management(riv_net,fld,"DOUBLE")

# build 2D array to pull data from
outdata = np.vstack((RMS_W,RMS_D,IMP_W,IMP_D)).T

# write the data back to the shapefile
jj=0
for fieldname in flds:
    with arcpy.da.UpdateCursor(riv_net,fieldname) as cursor: #Add data
        ii = 0
        for row in cursor:
            row[0]=outdata[ii,jj]
            cursor.updateRow(row)
            ii+=1
    jj+=1

