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
##if usegui==False:
##    indir1=r'D:\Box Sync\Wildfires\PLI_watersheds\Input_Data'
##    indir2=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data'
##    basenamein='a000'
##    
##    #USUAL INPUTS
##    # input filled DEM
##    inDEM=os.path.join(indir2,basenamein,basenamein+"demf.tif")
##    #input subcatchments
##    subcatchin=os.path.join(indir2,'west2',"i152yrmodburn_subcatchments2.shp")
##    useProb=True
##    P_field="df_Prob"
##    P_thresh=0.5
##
##    uid_field="sub_ID"
##    in_wtrshd=os.path.join(indir2,basenamein,basenamein+"_watershed.shp")
##    
##    
##    # REQUIRED DATA
##    #STATSGO
##    statsgo=os.path.join(indir1,'Soil_Data','STATSGO_westernUS.shp')
##    sthicfield="THICK"
##    #FIRE EXTENT
##    fire_perimeter=in_wtrshd
##    #BURN SEVERITY RASTER
##    ftype=os.path.join(indir1,'PLI_FireSeverity','PLI_class_50.tif')
##    #RAINFALL INTENSITY RASTER
##    RI=os.path.join(indir1,'Precipitation','precip2yr15ma.tif')
##    #OUPUT DATA
##    basename='i152yrmodburn'
##    out_temp_dir=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data\west2\gart14'
##    #output subcatments option set to [] and it will update the input
##    makenewoutput=True
##    subcatchout=os.path.join(indir2,'west2',basename+"_subcatchmentsG14v3.shp")
##    
    

    
#--------------End User Inputs--------------------------------

#USUAL INPUTS
# input filled DEM
inDEM=arcpy.GetParameterAsText(0)
#input subcatchments
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


# REQUIRED DATA
#FIRE EXTENT
fire_perimeter=arcpy.GetParameterAsText(7)

#STATSGO
statsgo=arcpy.GetParameterAsText(8)
sthicfield="THICK"
#BURN SEVERITY RASTER
ftype=arcpy.GetParameterAsText(9)
#RAINFALL INTENSITY RASTER/Constant value
RI=arcpy.GetParameterAsText(10)
#OUPUT DATA

basename=arcpy.GetParameterAsText(11)
out_temp_dir=arcpy.GetParameterAsText(12)
#output subcatments option set to [] and it will update the input
makenewoutput=arcpy.GetParameterAsText(13)
try:
    subcatchout=arcpy.GetParameterAsText(14)
except:
    subcatchout=[]
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

#.....................sqrt of i15 rainfall intensity...........................
#first check if an i15 field has already been computed in staley
outfieldi15="avg_i15"
fc_tf=chk_field_exists(subcatch,outfieldi15)

if not fc_tf: 
    if not RITF:
        # if i15 mean data doesn't exist compute it
        RI_tbl=os.path.join(out_temp_dir,basename+'_prec_i15')
        RI_ZStat = ZonalStatisticsAsTable(burn_basin,uid_field, RI, RI_tbl, "DATA", "MEAN")
        arcpy.AddField_management(RI_tbl, outfieldi15, "DOUBLE")
        arcpy.CalculateField_management(RI_tbl, outfieldi15, "!MEAN!", "PYTHON", "")
        arcpy.JoinField_management(subcatch,uid_field, RI_tbl,uid_field, outfieldi15)
    if RITF:
        print2("Using Constant Rainfall Intensity")
        arcpy.AddField_management(subcatch, outfieldi15, "DOUBLE")
        arcpy.CalculateField_management(subcatch, outfieldi15, RI, "PYTHON", "")

#............. area burned at moderate to high intensity.......................
# get file names for output naming
# Get fire type raster data
ftype_wtshd=os.path.join(out_temp_dir,basename+"ftype.tif")
arcpy.gp.ExtractByMask_sa(ftype, inDEM, ftype_wtshd)

#Make a binary raster where 1 is a crown burn
ftype_rc=os.path.join(out_temp_dir,basename+"ftype_rc.tif")
arcpy.gp.Con_sa(ftype_wtshd, "1", ftype_rc, "0", '"VALUE" > 2')


#Debris Flow Volume Term 2-Natural log of area burned at mod to high severity
ftypelogtbl=os.path.join(out_temp_dir,basename+"_area_modhigh")
arcpy.gp.ZonalStatisticsAsTable_sa(burn_basin, uid_field, ftype_rc, ftypelogtbl, "DATA", "MEAN")
outfieldmh='mh_ln_area'
arcpy.AddField_management(ftypelogtbl, outfieldmh, "DOUBLE")
expression = "lnBMH(!MEAN!,!AREA!)"
codeblock = """def lnBMH(MEAN,AREA):
    if MEAN == 0:
        return math.log(9/10000)
    else:
        return math.log(MEAN*AREA/1000000 )"""
arcpy.CalculateField_management(ftypelogtbl, outfieldmh,expression,"PYTHON_9.3",codeblock)
arcpy.JoinField_management(subcatch,uid_field,ftypelogtbl, uid_field,outfieldmh)


#........................ sqrt of basin releif.................................

# Square root of basin relief
relief_tbl=os.path.join(out_temp_dir,basename+"_relief")
arcpy.gp.ZonalStatisticsAsTable_sa(subcatch, uid_field, inDEM, relief_tbl, "DATA", "RANGE")
relief='relief'
arcpy.AddField_management(relief_tbl, relief, "DOUBLE")
arcpy.CalculateField_management(relief_tbl, relief, "math.sqrt( !RANGE! )", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field,relief_tbl, uid_field, relief)

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


#............................ Compute Volumes..................................
#.....................Calculate debris flow probability........................
if useProb:
    fields=[relief,outfieldmh,outfieldi15,soilvol,P_field]
else:
    fields=[relief,outfieldmh,outfieldi15,soilvol]
arr_db =  arcpy.da.FeatureClassToNumPyArray(subcatch, fields, skip_nulls=True)

# Volume Model Inputs
V_X1 = arr_db[fields[0]] # sqrt of basin releif
V_X2 = arr_db[fields[1]]#natural log of burned area at mod to high severity 
V_X3 = np.sqrt(arr_db[fields[2]]) # sqrt of i15 rainfall
SoilVol=arr_db[fields[3]]# volume of soil available


# Gartner et al., (2014) Volumetric Estimates
lnV = 4.22+(0.39*V_X3)+(0.36*V_X2)+(0.13*V_X1) # Modeled ln(Volume)
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
