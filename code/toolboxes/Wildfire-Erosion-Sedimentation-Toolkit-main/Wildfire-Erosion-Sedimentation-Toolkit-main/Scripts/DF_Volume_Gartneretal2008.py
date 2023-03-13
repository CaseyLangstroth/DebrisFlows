# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 23:22:51 2022

@author: Scott
"""

#--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import math
import numpy as np


arcpy.env.overwriteOutput=True
usegui=False

#--------------User Inputs----------------------------------
if usegui==False:
    indir1=r'D:\Box Sync\Wildfires\PLI_watersheds\Input_Data'
    indir2=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data'
    basename='a000'
    
    
    # input filled DEM
    inDEM=os.path.join(indir2,basename,basename+"demf.tif")
    subcatchin=os.path.join(indir2,'west2',"i152yrmodburn_subcatchments2.shp")
    uid_field="sub_ID"
    useProb=True
    P_field="df_Prob"
    P_thresh=0.5
    in_wtrshd=os.path.join(indir2,basename,basename+"_watershed.shp")
    fire_perimeter=in_wtrshd
    statsgo=os.path.join(indir1,'Soil_Data','STATSGO_westernUS.shp')

    RI=os.path.join(indir1,'Precipitation','precip2yr15ma.tif')
    
    out_temp_dir=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data\west2\gart08'
    makenewoutput=True
    subcatchout=os.path.join(indir2,'west2',basename+"_subcatchmentsG08.shp")
    slope_thresh=30

#--------------End User Inputs--------------------------------
if usegui==True:
    inDEM=arcpy.GetParameterAsText(0)
    subcatchin=arcpy.GetParameterAsText(1)
    uid_field=arcpy.GetParameterAsText(2)
    
    useProb=arcpy.GetParameterAsText(3)
    try:
        P_field=arcpy.GetParameterAsText(4)
    except:
        P_field=[]
    
    try:
        P_thresh=arcpy.GetParameter(5)    
    except:
        P_thresh=[]
    
    
    in_wtrshd=arcpy.GetParameterAsText(6)
    fire_perimeter=arcpy.GetParameterAsText(7)
    statsgo=arcpy.GetParameterAsText(8)
    sthicfield="THICK"

    RI=arcpy.GetParameterAsText(9)
    
    basename=arcpy.GetParameterAsText(10)    
    out_temp_dir=arcpy.GetParameterAsText(11)

    makenewoutput=arcpy.GetParameterAsText(12)
    try:
        subcatchout=arcpy.GetParameterAsText(13)
    except:
        subcatchout=[]
    slope_thresh=30

#%%------------------------- Basic Functions-----------------------------------
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

# check if a value is float
def isfloat(RI):
    try:
        float(RI)
        return True
    except ValueError:
        return False

# check if a field already exists
def chk_field_exists(infc,field):
    lstFields=arcpy.ListFields(infc)
    if field in lstFields:
        return True
    else:
        return False

#this corrects for arcpy.getparameterastext
def string2boolean(x):
        if x == 'true':
             x = True
        else:
             x = False
        return x
#%% ..................setup environmental parameters..............................
#Extract input DEM Raster Resolution
# cell size x direction
dx_temp=arcpy.GetRasterProperties_management(inDEM, "CELLSIZEX")
dx_indem=float(dx_temp.getOutput(0))
# cell size y direction
dy_temp=arcpy.GetRasterProperties_management(inDEM, "CELLSIZEY")
dy_indem=float(dy_temp.getOutput(0))

arcpy.env.extent=inDEM
arcpy.env.snapRaster=inDEM
arcpy.env.cellSize=dx_indem


SR_wtrshd = arcpy.Describe(inDEM).spatialReference
SR=str(SR_wtrshd.Name)
SR=SR.replace("_"," ")
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(SR)

chk_mk_dir(out_temp_dir)
print2('All intermediate files will be output to: \n'+out_temp_dir,usegui)


# Change 'true' or 'false' to True or False
if usegui:
    makenewoutput=string2boolean(makenewoutput)
    useProb=string2boolean(useProb)

if useProb:
    #check all probability parameters are provided if using that analysis
    if P_thresh>0 and not P_field:
        errmsg="You have specified a probability threshold and not a  field "\
                "containing probabilities from the Staley et al. code./n"\
                "Please set probability theshold to zero or specify a probability field"
        reporterror(errmsg,usegui)
        
    if P_thresh==0 and  P_field:
        warmsg="WARNING: You have specified a probability field and a threshold"\
                "of zero. The probability field is not being used."
        print2(warmsg,usegui)


    
# error checking on input RI
RITF = isfloat(RI)
if not RITF:
    if not arcpy.Exists(RI): 
        errmsg=('Your input rainfall intensity raster does not exist')

#%%----------------------------Begin analysis----------------------------------
# make a copy of input file and write to it or write to the input file
if makenewoutput:
    subcatchoutdir,subcatchoutfid=os.path.split(subcatchout)
    arcpy.conversion.FeatureClassToFeatureClass(subcatchin, subcatchoutdir, 
                                            subcatchoutfid)
    subcatch=subcatchout
else:
    subcatch=subcatchin


#........................Flag the catchments that burned.......................
arcpy.AddField_management(subcatch, "Burned", "Double")
arcpy.MakeFeatureLayer_management(subcatch, 'lyrSC')
arcpy.management.SelectLayerByLocation('lyrSC', "INTERSECT", fire_perimeter, 
                                        None, "NEW_SELECTION", "NOT_INVERT")
arcpy.CalculateField_management('lyrSC', "Burned", "1", "PYTHON")
#arcpy.DeleteFeatures_management('lyrSC')

# export those basins to a temp file to run calculations on and join back
arcpy.conversion.FeatureClassToFeatureClass(subcatch, out_temp_dir, 
                                            basename+"_subcatchburned.shp",
                                            "Burned = 1")
burn_basin=os.path.join(out_temp_dir,basename+"_subcatchburned.shp")

#.............compute the average i 10 for subcatchment........................
if not RITF:
    print2("Using Input Rainfall Intensity Raster",usegui)
    RI_tbl=os.path.join(out_temp_dir,basename+'_prec_i10')
    ZonalStatisticsAsTable(burn_basin,uid_field, RI, RI_tbl, "DATA", "MEAN")
    outfieldi10="avg_i10"
    arcpy.AddField_management(RI_tbl, outfieldi10, "DOUBLE")
    arcpy.CalculateField_management(RI_tbl, outfieldi10, "!MEAN!", "PYTHON", "")
    arcpy.JoinField_management(subcatch,uid_field, RI_tbl,uid_field, outfieldi10)

if RITF:
    print2("Using Constant Rainfall Intensity")
    outfieldi10="avg_i10"
    arcpy.AddField_management(subcatch, outfieldi10, "DOUBLE")
    arcpy.CalculateField_management(subcatch, outfieldi10, RI, "PYTHON", "")
#..............basin area with slopes greater than or equal to 30%.............

#Compute a cellular slope.
outslope=arcpy.sa.Slope(inDEM, "DEGREE", "1", "PLANAR", "METER")
slp=os.path.join(out_temp_dir,"slp.tif")
outslope.save(slp)
# make raster a binary where 1 exceeds some threshodld
slp_rc=os.path.join(out_temp_dir,"slp_rc.tif")
arcpy.gp.Con_sa(slp, "1", slp_rc, "0", "VALUE>"+str(slope_thresh))

slp_tbl=os.path.join(out_temp_dir,basename+'_slp30')
ZonalStatisticsAsTable(burn_basin,uid_field, slp_rc, slp_tbl, "DATA", "SUM")
slpfield='area_'+str(slope_thresh)
arcpy.AddField_management(slp_tbl, slpfield, "DOUBLE")
areaexp="!SUM!*0.000001*"+str(dx_indem)+"*"+str(dy_indem)
arcpy.CalculateField_management(slp_tbl, slpfield, areaexp, "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, slp_tbl,uid_field, slpfield)

#.............................Soil Volume......................................
arcpy.env.extent=inDEM
arcpy.env.snapRaster=inDEM
statgoaoi=os.path.join(out_temp_dir,basename+"_statsgo.shp")
arcpy.Clip_analysis(statsgo, in_wtrshd, statgoaoi,"")

statgoaoiras_ST=os.path.join(out_temp_dir,basename+"_soilthick.tif")
arcpy.PolygonToRaster_conversion(statgoaoi,sthicfield,statgoaoiras_ST, "CELL_CENTER", "NONE", str(dx_indem))
thick_tbl=os.path.join(out_temp_dir,basename+"_soilthick")
arcpy.gp.ZonalStatisticsAsTable_sa(subcatch, uid_field, statgoaoiras_ST,thick_tbl, "DATA", "SUM")
soilvol='soilvol_m3'
arcpy.AddField_management(thick_tbl, soilvol, "DOUBLE")
volume="!SUM!*"+str(dx_indem)+"*"+str(dx_indem)+'*0.01'# 0.01 CONVERTS CM TO M
arcpy.CalculateField_management(thick_tbl, soilvol, volume, "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field,thick_tbl, uid_field, soilvol)


#............................... Calculate Volumes............................ 

if useProb:
    fields=[slpfield,outfieldi10,soilvol,P_field]
else:
    fields=[slpfield,outfieldi10,soilvol]
arr_db =  arcpy.da.FeatureClassToNumPyArray(subcatch, fields, skip_nulls=True)

# Volume Model Inputs
S = arr_db[fields[0]]#S=basin area >30deg 
R = arr_db[fields[1]]#P=peak i10 rainfall mm/hr

SoilVol=arr_db[fields[2]]# volume of soil available


# Gartner et al., (2008) Volumetric Estimates
lnV = 0.72*np.log(S)-(0.02*R)+8.54 # Modeled ln(Volume)
Vol = np.exp(lnV) #Debris Flow Volume, m3

# set those with a probability below the threshold to nan/0 m3
if useProb:
    if P_thresh>0:
        P=arr_db[fields[4]]
        Vol[P<P_thresh] = np.nan; #Did not generate debris flow

# final Vol Gartner output now make and adj version capping max at soil vol
Vol_Gartner = Vol.copy()
# total volume of sediment eroded
#Total_Vol_Gartner = np.nansum(Vol_Gartner)

# flag if erroded more sediment than available making it supply limited
supply_lim=np.zeros(Vol.size)
supply_lim[Vol>SoilVol]=1

#Correction Factor
#If predicted volume, Vol, is > total soil volume, SoilVol, set Vol = SoilVol
Vol[Vol>SoilVol]=SoilVol[Vol>SoilVol]

#print2(Frac_Input)
SoilLossFrac = Vol/SoilVol

# change all nan to zero
Vol_Gartner[np.isnan(Vol_Gartner)]=0
Vol[np.isnan(Vol)]=0
SoilLossFrac[np.isnan(SoilLossFrac)]=0


# add in a binary 1 is bad 0 is good
# then add a supplylimit (sl_vol) field that is Soiladj


#-----------------------Write Data back to shapefile---------------------------
# Add fields of interest to the shape file
flds=['Vol','supp_lim','sl_vol','FrSoilLos']

for fieldname in flds:
    arcpy.AddField_management(subcatch,fieldname, "Double")

# build 2D array to pull data from
outdata = np.vstack((Vol_Gartner,supply_lim,Vol,SoilLossFrac)).T

jj=0
for fieldname in flds:
    with arcpy.da.UpdateCursor(subcatch,fieldname) as cursor: #Add data
        ii = 0
        for row in cursor:
            row[0]=outdata[ii,jj]
            cursor.updateRow(row)
            ii+=1
    jj+=1

# Ensure non burned basins have a value of zero
arcpy.MakeFeatureLayer_management(subcatch, 'lyrSC')
s_exp="Burned=0"
arcpy.management.SelectLayerByAttribute("lyrSC", "NEW_SELECTION", s_exp)
for field in flds:
    arcpy.CalculateField_management("lyrSC", field, "0","PYTHON","")
    