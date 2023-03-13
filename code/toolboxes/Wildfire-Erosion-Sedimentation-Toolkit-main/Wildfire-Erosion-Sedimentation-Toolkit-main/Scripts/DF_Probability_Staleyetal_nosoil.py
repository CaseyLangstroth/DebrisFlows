#--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import numpy as np


arcpy.env.overwriteOutput=True
usegui=True

#--------------User Inputs----------------------------------
if usegui==False:
    indir1=r'D:\Box Sync\Wildfires\PLI_watersheds\Input_Data'
    indir2=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data'
    basenamein='a000'
    #USUAL INPUTS
    # input filled DEM
    inDEM=os.path.join(indir2,basenamein,basenamein+"demf.tif")
    #input subcatchments
    subcatchin=os.path.join(indir2,basenamein,basenamein+"_subcatchments.shp")
    uid_field="sub_ID"
    in_wtrshd=os.path.join(indir2,basenamein,basenamein+"_watershed.shp")
    
    # REQUIRED DATA
    #STATSGO
    # statsgo=os.path.join(indir1,'Soil_Data','STATSGO_westernUS.shp')
    #FIRE EXTENT
    fire_perimeter=in_wtrshd
    #BURN SEVERITY RASTER
    ftype=os.path.join(indir1,'PLI_FireSeverity','PLI_class_50.tif')
    #dnbr RASTER
    dNBR=os.path.join(indir1,'PLI_dnbr','PLI_dnbr_50.tif')
    #RAINFALL INTENSITY RASTER
    RI=os.path.join(indir1,'Precipitation','precip2yr15ma.tif')
    #OUPUT DATA
    basename='nosoil '
    out_temp_dir=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data\west2\staley'
    
    #output subcatments option set to [] and it will update the input
    makenewoutput=True
    subcatchout=os.path.join(indir2,'west2',basename+"_subcatchments2.shp")
    slope_thresh=23
    

#--------------End User Inputs--------------------------------
if usegui==True:
    
    #USUAL INPUTS
    # input filled DEM
    inDEM=arcpy.GetParameterAsText(0)
    #input subcatchments
    subcatchin=arcpy.GetParameterAsText(1)
    uid_field=arcpy.GetParameterAsText(2)
    in_wtrshd=arcpy.GetParameterAsText(3)
    
    # REQUIRED DATA
    #STATSGO
    # statsgo=arcpy.GetParameterAsText(4)
    #FIRE EXTENT
    fire_perimeter=arcpy.GetParameterAsText(4)
    #BURN SEVERITY RASTER
    ftype=arcpy.GetParameterAsText(5)
    #dnbr RASTER
    dNBR=arcpy.GetParameterAsText(6)
    #RAINFALL INTENSITY RASTER
    RI=arcpy.GetParameterAsText(7)
    #OUPUT DATA
    basename=arcpy.GetParameterAsText(8)
    out_temp_dir=arcpy.GetParameterAsText(9)
    makenewoutput=arcpy.GetParameterAsText(10)
    #output subcatments option set to [] and it will update the input
    try:
        subcatchout=arcpy.GetParameterAsText(11)
    except:
        subcatout=[]
    slope_thresh=23

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


# error checking on input RI
RITF = isfloat(RI)
if not RITF:
    if not arcpy.Exists(RI): 
        errmsg=('Your input rainfall intensity raster does not exist')
if usegui:
    makenewoutput=string2boolean(makenewoutput)
#%%------------------------------------Begin Analysis--------------------------
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

#.........................Soil Data Anlysis....................................
# #STATSGO Kffactor
# # Clip soil data to watershed
# statgoaoi=os.path.join(out_temp_dir,basename+"_statsgo.shp")
# arcpy.Clip_analysis(statsgo, in_wtrshd, statgoaoi,"")


# # Rasterize K factor from statsgo shape file
# statgoaoiras_k=os.path.join(out_temp_dir,basename+"_kffact.tif")
# arcpy.PolygonToRaster_conversion(statgoaoi,"KFFACT",statgoaoiras_k,
#                                  "CELL_CENTER", "NONE", str(dx_indem))

# # Zonal stats to get mean k in each debris flow basin
# statsgotble=os.path.join(out_temp_dir,basename+'kfactor')
# ZonalStatisticsAsTable(burn_basin, uid_field, statgoaoiras_k,
#                           statsgotble, "DATA", "MEAN")

# outfieldk='k_factor'
# arcpy.AddField_management(statsgotble, outfieldk, "DOUBLE")
# arcpy.CalculateField_management(statsgotble, outfieldk, "!MEAN!","PYTHON","")
# arcpy.JoinField_management(subcatch,uid_field,statsgotble,uid_field, outfieldk)

#.........Percent Area Burned at Mod-to-High Severity @ >23 degree slope.......
#Compute a cellular slope.
outslope=arcpy.sa.Slope(inDEM, "DEGREE", "1", "PLANAR", "METER")
slp=os.path.join(out_temp_dir,"slp.tif")
outslope.save(slp)

# make raster a binary where 1 exceeds some threshodld
slp_rc=os.path.join(out_temp_dir,"slp_rc.tif")
arcpy.gp.Con_sa(slp, "1", slp_rc, "0", "VALUE>"+str(slope_thresh))


#Percent Area Burned at Mod-to-High Severity

# get file names for output naming
# Get fire type raster data
ftype_wtshd=os.path.join(out_temp_dir,basename+"ftype.tif")
arcpy.gp.ExtractByMask_sa(ftype, inDEM, ftype_wtshd)

#Make a binary raster where 1 is a crown burn
ftype_rc=os.path.join(out_temp_dir,basename+"ftype_rc.tif")
arcpy.gp.Con_sa(ftype_wtshd, "1", ftype_rc, "0", '"VALUE" > 2')

# make a single mask raster where slope and had a crown burn
ftype_sevslp=os.path.join(out_temp_dir,basename+"ftype_sevslp.tif")
arcpy.gp.Times_sa(ftype_rc, slp_rc, ftype_sevslp)

# Calculate area of high slope and crown burn
ftype_tbl=os.path.join(out_temp_dir,basename+"_percCBslp")
arcpy.gp.ZonalStatisticsAsTable_sa(burn_basin, uid_field, ftype_sevslp, 
                                   ftype_tbl, "DATA", "MEAN")
outfieldmh='area_MH_'+str(slope_thresh) # area mod high >23
arcpy.AddField_management(ftype_tbl, outfieldmh, "DOUBLE")
arcpy.CalculateField_management(ftype_tbl, outfieldmh, "!MEAN!","PYTHON","")

# Add the table back to the subcatchments
arcpy.JoinField_management(subcatch,uid_field,ftype_tbl, uid_field, outfieldmh)

#.........................Average dNBR/1000...................................
dNBR_tbl=os.path.join(out_temp_dir,'dNBRtbl')
# get average dnbr under each burned basin
arcpy.gp.ZonalStatisticsAsTable_sa(burn_basin, uid_field, dNBR, dNBR_tbl, "DATA", "MEAN")
outfielddnbr='avgdNBR'
arcpy.AddField_management(dNBR_tbl, outfielddnbr, "DOUBLE")
arcpy.CalculateField_management(dNBR_tbl, outfielddnbr, "!MEAN!/1000","PYTHON","")
arcpy.JoinField_management(subcatch,uid_field,dNBR_tbl, uid_field, outfielddnbr)


#.........................Rain Fall Intensities................................

if not RITF:
    print2("Using Input Rainfall Intensity Raster",usegui)
    RI_tbl=os.path.join(out_temp_dir,basename+'_prec_i15')
    RI_ZStat = ZonalStatisticsAsTable(subcatch,uid_field, RI, RI_tbl, "DATA", "MEAN")
    outfieldi15="avg_i15"
    arcpy.AddField_management(RI_tbl, outfieldi15, "DOUBLE")
    arcpy.CalculateField_management(RI_tbl, outfieldi15, "!MEAN!", "PYTHON", "")
    arcpy.JoinField_management(subcatch,uid_field, RI_tbl,uid_field, outfieldi15)

if RITF:
    print2("Using Constant Rainfall Intensity",usegui)
    outfieldi15="avg_i15"
    arcpy.AddField_management(subcatch, outfieldi15, "DOUBLE")
    arcpy.CalculateField_management(subcatch, outfieldi15, RI, "PYTHON", "")

#.....................Calculate debris flow probability........................
fields=[outfieldmh,outfielddnbr,outfieldi15]

arr_db =  arcpy.da.FeatureClassToNumPyArray(subcatch, fields, skip_nulls=True)

# Probability Model Inputs
M1_X1 = arr_db[fields[0]] # Percent area burned at mod to high severity
M1_X2 = arr_db[fields[1]]# average dNBR
Rain_mag=arr_db[fields[3]]/4 # 15-min rainfall magnitude (mm)

#Link function
chi= -3.36+(0.54*M1_X1*Rain_mag)+(0.77*M1_X2*Rain_mag)
#Debris Flow Probability
P = np.exp(chi)/(1+np.exp(chi)) 

# Write the probability back to the output
arcpy.AddField_management(subcatch, "df_prob", "Double")

with arcpy.da.UpdateCursor(subcatch, "df_prob") as cursor: #Add data
    i = 0
    for row in cursor:
        row[0]=P[i]
        cursor.updateRow(row)
        i+=1

# Ensure non burned basins have a value of zero
arcpy.MakeFeatureLayer_management(subcatch, 'lyrSC')
s_exp="Burned=0"
arcpy.management.SelectLayerByAttribute("lyrSC", "NEW_SELECTION", s_exp)
arcpy.CalculateField_management("lyrSC", 'df_prob', "0","PYTHON","")
