#--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import math
import numpy as np


arcpy.env.overwriteOutput=True
usegui=True

#%%---------------------------User Inputs----------------------------------
if usegui==False:
    indir1=r'D:\Box Sync\Wildfires\PLI_watersheds\Input_Data'
    indir2=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data'
    basenamein='a000'
    
    #USUAL INPUTS
    # input filled DEM
    inDEM=os.path.join(indir2,basenamein,basenamein+"demf.tif")
    fac=os.path.join(indir2,basenamein,basenamein+"fac.tif")
    #input subcatchments
    subcatchin=os.path.join(indir2,'west2',"i152yrmodburn_subcatchments2.shp")
    uid_field="sub_ID"
    useProb=True
    P_field="df_Prob"
    P_thresh=0.5
    
    
    in_wtrshd=os.path.join(indir2,basenamein,basenamein+"_watershed.shp")
    #FIRE EXTENT
    fire_perimeter=in_wtrshd
    #BURN SEVERITY RASTER
    ftype=os.path.join(indir1,'PLI_FireSeverity','PLI_class_50.tif')
    
    Temp=os.path.join(indir1,'Temperature','PRISM_tmean_30yr_normal_800mM3_annual.tif')
    LCdata=os.path.join(indir1,'Landcover','nlcd_2019_land_cover.tif')
    statsgo=os.path.join(indir1,'Soil_Data','STATSGO_westernUS.shp')
    
    comstrR=os.path.join(indir1,'streamcat_underlying','CompressStrength.tif')
    mgoR=os.path.join(indir1,'streamcat_underlying','MgO_raster.tif')
    runoffR=os.path.join(indir1,'streamcat_underlying','runoff.tif')
    
    basename='i152yrmodburn'
    makenewoutput=True

    out_temp_dir=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data\west2\wall5'
    subcatchout=os.path.join(indir2,'west2',basename+"_subcatchmentswall7.shp")

if usegui==True:
    #USUAL INPUTS
    # input filled DEM
    inDEM=arcpy.GetParameterAsText(0)
    fac=arcpy.GetParameterAsText(1)
    #input subcatchments
    subcatchin=arcpy.GetParameterAsText(2)
    uid_field=arcpy.GetParameterAsText(3)
    useProb=arcpy.GetParameterAsText(4)
    
    try:
        P_field=arcpy.GetParameterAsText(5)
    except:
        P_field=[]
    try:
        P_thresh=arcpy.GetParameter(6)
    except:
        P_thresh=[]
    
    
    in_wtrshd=arcpy.GetParameterAsText(7)
    #FIRE EXTENT
    fire_perimeter=arcpy.GetParameterAsText(8)
    #BURN SEVERITY RASTER
    ftype=arcpy.GetParameterAsText(9)
    
    Temp=arcpy.GetParameterAsText(10)
    LCdata=arcpy.GetParameterAsText(11)
    statsgo=arcpy.GetParameterAsText(12)
    
    comstrR=arcpy.GetParameterAsText(13)
    mgoR=arcpy.GetParameterAsText(14)
    runoffR=arcpy.GetParameterAsText(15)
    
    basename=arcpy.GetParameterAsText(16)

    out_temp_dir=arcpy.GetParameterAsText(17)
    makenewoutput=arcpy.GetParameterAsText(18)

    try:
        subcatchout=arcpy.GetParameterAsText(19)
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
    
def computeTWI(slp,fac,TWIrastername):
    dx_temp=arcpy.GetRasterProperties_management(fac, "CELLSIZEX")
    dx=float(dx_temp.getOutput(0))
    
    # get array of slope values in degrees
    S_deg=arcpy.RasterToNumPyArray(slp,nodata_to_value=np.nan)
    # Convert slope from degrees to radieans
    S_rad=np.radians(S_deg)
    del (S_deg)
    
    # get flow accumulation as an array
    faccum=arcpy.RasterToNumPyArray(fac,nodata_to_value=np.nan)
    #uplslope contributing area
    uca=faccum*dx*dx
    del faccum
    S_rad[S_rad==0]=0.0001
    # comput topographic wetness index
    TWI=np.log(uca/np.tan(S_rad))
    
    # create output raster
    TWI[np.isnan(TWI)] = -9999
    # Get coordinate system for output rasters
    dscRas = arcpy.Raster(slp)
    lowerLeft = arcpy.Point(dscRas.extent.XMin,dscRas.extent.YMin)
    dsc=arcpy.Describe(dscRas)
    cellSize=dscRas.meanCellWidth
    coord_sys=dsc.spatialReference
    outRaster = arcpy.NumPyArrayToRaster(TWI,lowerLeft,cellSize, value_to_nodata=-9999)
    outRaster.save(TWIrastername)
    del outRaster
    arcpy.DefineProjection_management(TWIrastername,coord_sys) 

#this corrects for arcpy.getparameterastext
def string2boolean(x):
        if x == 'true':
             x = True
        else:
             x = False
        return x


print2('Add in usual code for getting cell area in meters',usegui)
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
#%%----------------------------Begin analysis----------------------------------
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



#.......................Compute required statistics............................

#..............................Average Temperature.............................
outtemp="avg_Temp"
temp_tbl=os.path.join(out_temp_dir,basename+'_prec_i15.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, Temp, temp_tbl, "NODATA", "MEAN")
arcpy.AddField_management(temp_tbl, outtemp, "DOUBLE")
arcpy.CalculateField_management(temp_tbl, outtemp, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, temp_tbl,uid_field, outtemp)


#..............................Percent Conifer.................................
LC_wtshd=os.path.join(out_temp_dir,basename+"LC.tif")
arcpy.gp.ExtractByMask_sa(LCdata, inDEM, LC_wtshd)

#Make a binary raster where 1 is a crown burn
conifer_rc=os.path.join(out_temp_dir,basename+"conifer_rc.tif")
LCfield="NLCD_Land"
LCevergreen='\'Evergreen Forest\''
lcexpress=LCfield+' = '+LCevergreen
arcpy.gp.Con_sa(LC_wtshd, "1", conifer_rc, "0", lcexpress)

outconifier="percconif"
conifer_tbl=os.path.join(out_temp_dir,basename+'_prec_conifer.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, conifer_rc, conifer_tbl, "NODATA", "MEAN")
arcpy.AddField_management(conifer_tbl, outconifier, "DOUBLE")
arcpy.CalculateField_management(conifer_tbl, outconifier, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, conifer_tbl,uid_field, outconifier)

#..........................Area burned mod-high--------------------------------
ftype_wtshd=os.path.join(out_temp_dir,basename+"ftype.tif")
arcpy.gp.ExtractByMask_sa(ftype, inDEM, ftype_wtshd)

#Make a binary raster where 1 is a crown burn
ftype_rc=os.path.join(out_temp_dir,basename+"ftype_rc.tif")
arcpy.gp.Con_sa(ftype_wtshd, "1", ftype_rc, "0", '"VALUE" > 2')


#Debris Flow Volume Term 2-Natural log of area burned at mod to high severity
ftypelogtbl=os.path.join(out_temp_dir,basename+"_area_modhigh.dbf")
arcpy.gp.ZonalStatisticsAsTable_sa(burn_basin, uid_field, ftype_rc, 
                                   ftypelogtbl, "DATA", "SUM")
outfieldmh='mh_area'
arcpy.AddField_management(ftypelogtbl, outfieldmh, "DOUBLE")
MHAexp="!SUM!*0.000001*"+str(dx_indem)+"*"+str(dx_indem)#in km2
arcpy.CalculateField_management(ftypelogtbl, outfieldmh, MHAexp, "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, ftypelogtbl,uid_field, outfieldmh)

#..........................percent/ area slope >23.............................
#Compute a cellular slope.
outslope=arcpy.sa.Slope(inDEM, "DEGREE", "1", "PLANAR", "METER")
#outcon= arcpy.sa.Con(in_conditional_raster = outslope, in_true_raster_or_constant= 0.01, where_clause= 'VALUE = 0')
slp=os.path.join(out_temp_dir,"slp.tif")
outslope.save(slp)

# make raster a binary where 1 exceeds some threshodld
slope_thresh=23
slp_rc=os.path.join(out_temp_dir,"slp_rc.tif")
arcpy.gp.Con_sa(slp, "1", slp_rc, "0", "VALUE>"+str(slope_thresh))

# percent of area with slope >23
slp_tbl=os.path.join(out_temp_dir,basename+'_prcslp23.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, slp_rc, slp_tbl, "DATA", "MEAN")
slpfield='percslp_'+str(slope_thresh)
arcpy.AddField_management(slp_tbl, slpfield, "DOUBLE")
arcpy.CalculateField_management(slp_tbl, slpfield, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, slp_tbl,uid_field, slpfield)

# total area of slope >23 km2
slpA_tbl=os.path.join(out_temp_dir,basename+'_areaslp23.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, slp_rc, slpA_tbl, "DATA", "SUM")
slpAfield='areaslp_'+str(slope_thresh)
arcpy.AddField_management(slpA_tbl, slpAfield, "DOUBLE")
slpAexp="!SUM!*0.000001*"+str(dx_indem)+"*"+str(dx_indem)#in km2
arcpy.CalculateField_management(slpA_tbl, slpAfield, slpAexp, "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, slpA_tbl,uid_field, slpAfield)


#..............................STATSGO Processing..............................
#Extract STATSGO
# Clip soil data to watershed
statgoaoi=os.path.join(out_temp_dir,basename+"_statsgo.shp")
arcpy.Clip_analysis(statsgo, in_wtrshd, statgoaoi,"")

# thickness
statgoaoiras_ST=os.path.join(out_temp_dir,basename+"_sthick.tif")
arcpy.PolygonToRaster_conversion(statgoaoi,"THICK",statgoaoiras_ST,  "CELL_CENTER", "NONE", str(dx_indem))


thick_tbl=os.path.join(out_temp_dir,basename+"_soilthick.dbf")
arcpy.gp.ZonalStatisticsAsTable_sa(burn_basin, uid_field, statgoaoiras_ST,
                                   thick_tbl, "DATA", "MEAN")

# average thickness
soilthick='soilthick'
arcpy.AddField_management(thick_tbl, soilthick, "DOUBLE")
thick="!MEAN!"# 0.01 CONVERTS CM TO M
arcpy.CalculateField_management(thick_tbl, soilthick, thick, "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field,thick_tbl, uid_field, soilthick)

vol_tbl=os.path.join(out_temp_dir,basename+"_soilthick.dbf")
arcpy.gp.ZonalStatisticsAsTable_sa(burn_basin, uid_field, statgoaoiras_ST,
                                   vol_tbl, "DATA", "SUM")

# total volume
soilvol='soilvol_m3'
arcpy.AddField_management(vol_tbl, soilvol, "DOUBLE")
volume="!SUM!*"+str(dx_indem)+"*"+str(dx_indem)+'*0.01'# 0.01 CONVERTS CM TO M
arcpy.CalculateField_management(vol_tbl, soilvol, volume, "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field,vol_tbl, uid_field, soilvol)

# permeability
statgo_perm=os.path.join(out_temp_dir,basename+"_perm.tif")
arcpy.PolygonToRaster_conversion(statgoaoi,"PERM",statgo_perm, 
                                 "CELL_CENTER", "NONE", str(dx_indem))

permtbl=os.path.join(out_temp_dir,basename+'perm.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin, uid_field, statgo_perm,
                          permtbl, "DATA", "MEAN")
outfieldperm='avg_perm'
arcpy.AddField_management(permtbl, outfieldperm, "DOUBLE")
arcpy.CalculateField_management(permtbl, outfieldperm, "!MEAN!","PYTHON","")
arcpy.JoinField_management(subcatch,uid_field,permtbl,uid_field, outfieldperm)

# percent organic matter
statgo_OM=os.path.join(out_temp_dir,basename+"_orgmatter.tif")
arcpy.PolygonToRaster_conversion(statgoaoi,"OrgMatter",statgo_OM, 
                                 "CELL_CENTER", "NONE", str(dx_indem))

OMtbl=os.path.join(out_temp_dir,basename+'OrgMatter.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin, uid_field, statgo_OM,
                          OMtbl, "DATA", "MEAN")
outfieldom='perc_org'
arcpy.AddField_management(OMtbl, outfieldom, "DOUBLE")
arcpy.CalculateField_management(OMtbl, outfieldom, "!MEAN!","PYTHON","")
arcpy.JoinField_management(subcatch,uid_field,OMtbl,uid_field, outfieldom)

#rock comp strengh
LC_wtshd=os.path.join(out_temp_dir,basename+"LC.tif")
arcpy.gp.ExtractByMask_sa(LCdata, inDEM, LC_wtshd)

#Percent Conifier
conifer_rc=os.path.join(out_temp_dir,basename+"conifer_rc.tif")
print2('make these a user input',usegui)
LCfield="NLCD_Land"
LCevergreen='\'Evergreen Forest\''
lcexpress=LCfield+' = '+LCevergreen
arcpy.gp.Con_sa(LC_wtshd, "1", conifer_rc, "0", lcexpress)

outconifier="percconif"
conifer_tbl=os.path.join(out_temp_dir,basename+'_prec_conifer.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, conifer_rc, conifer_tbl, "NODATA", "MEAN")
arcpy.AddField_management(conifer_tbl, outconifier, "DOUBLE")
arcpy.CalculateField_management(conifer_tbl, outconifier, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, conifer_tbl,uid_field, outconifier)

#.....................StreamCat underlying layers..............................
#mgo
MGO_wtshd=os.path.join(out_temp_dir,basename+"MGO.tif")
arcpy.gp.ExtractByMask_sa(mgoR, inDEM, MGO_wtshd)

outMGO="avg_MgO"
MGO_tbl=os.path.join(out_temp_dir,basename+'_avg_MgO.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, MGO_wtshd, MGO_tbl, "NODATA", "MEAN")
arcpy.AddField_management(MGO_tbl, outMGO, "DOUBLE")
arcpy.CalculateField_management(MGO_tbl, outMGO, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, MGO_tbl,uid_field, outMGO)

#compressive strengh
CST_wtshd=os.path.join(out_temp_dir,basename+"comp_strength.tif")
arcpy.gp.ExtractByMask_sa(comstrR, inDEM, CST_wtshd)

outCS="avg_ComStr"
CS_tbl=os.path.join(out_temp_dir,basename+'_avg_comstr.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, CST_wtshd, CS_tbl, "NODATA", "MEAN")
arcpy.AddField_management(CS_tbl, outCS, "DOUBLE")
arcpy.CalculateField_management(CS_tbl, outCS, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, CS_tbl,uid_field, outCS)

# #RUNOFF
runoffaoi=os.path.join(out_temp_dir,basename+"_runoff.tif")
arcpy.gp.ExtractByMask_sa(runoffR, inDEM, runoffaoi)

outRO="avg_runoff"
runoff_tbl=os.path.join(out_temp_dir,basename+'_avg_runoff.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, runoffaoi, runoff_tbl, "NODATA", "MEAN")
arcpy.AddField_management(runoff_tbl, outRO, "DOUBLE")
arcpy.CalculateField_management(runoff_tbl, outRO, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, runoff_tbl,uid_field, outRO)

#.......................Average Elevation......................................
outE="avg_Elev"
E_tbl=os.path.join(out_temp_dir,basename+'_avg_Elev.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, inDEM, E_tbl, "NODATA", "MEAN")
arcpy.AddField_management(E_tbl, outE, "DOUBLE")
arcpy.CalculateField_management(E_tbl, outE, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, E_tbl,uid_field, outE)

#--------------------------Wetness Index---------------------------------------
#TWI
TWIR=os.path.join(out_temp_dir,basename+'_TWI.tif')
computeTWI(slp,fac,TWIR)

outTWI="avg_Wetind"
TWI_tbl=os.path.join(out_temp_dir,basename+'_wetness.dbf')
arcpy.sa.ZonalStatisticsAsTable(burn_basin,uid_field, TWIR, TWI_tbl, "DATA", "MEAN")
arcpy.AddField_management(TWI_tbl, outTWI, "DOUBLE")
arcpy.CalculateField_management(TWI_tbl, outTWI, "!MEAN!", "PYTHON", "")
arcpy.JoinField_management(subcatch,uid_field, TWI_tbl,uid_field, outTWI)


#------------------------Run the calculations----------------------------------
if useProb:
    fields=[outtemp,outfieldmh,outconifier,slpfield,slpAfield,soilthick,soilvol,
            outfieldperm,outfieldom,outMGO,outCS,outRO,outE,outTWI,
            P_field]
else:
    fields=[outtemp,outfieldmh,outconifier,slpfield,slpAfield,soilthick,soilvol,
            outfieldperm,outfieldom,outMGO,outCS,outRO,outE,outTWI]
  
arr_db =arcpy.da.FeatureClassToNumPyArray(subcatch, fields, skip_nulls=True)



# assign values to parameters
Om = arr_db[outfieldom]#organic matter
arcpy.AddMessage('Om:')
arcpy.AddMessage(Om)
Sa=arr_db[slpAfield]# slope area >23deg
arcpy.AddMessage('Sa:')
arcpy.AddMessage(Sa)
Bmh=arr_db[outfieldmh]# area burned moderate high
arcpy.AddMessage('mharea:')
arcpy.AddMessage(Bmh)
Ro=arr_db[outRO]#average runoff
arcpy.AddMessage('Ro:')
arcpy.AddMessage(Ro)
T=arr_db[outtemp]# average temp
arcpy.AddMessage('T:')
arcpy.AddMessage(T)
Cp=arr_db[outconifier]# % conifer
arcpy.AddMessage('Cp:')
arcpy.AddMessage(Cp)
Sp=arr_db[slpfield] #% slope >23
arcpy.AddMessage('Sp:')
arcpy.AddMessage(Sp)
Rd=arr_db[soilthick]#average soil thickness
arcpy.AddMessage('Rd:')
arcpy.AddMessage(Rd)
CS=arr_db[outCS]#avg compressive strench
arcpy.AddMessage('Cs:')
arcpy.AddMessage(CS)
K=arr_db[outfieldperm]#permeability
arcpy.AddMessage('K:')
arcpy.AddMessage(K)
WI=arr_db[outTWI]# wetness index
arcpy.AddMessage('WI:')
arcpy.AddMessage(WI)
E=arr_db[outE]#avg elevation
arcpy.AddMessage('E:')
arcpy.AddMessage(E)
MgO=arr_db[outMGO]#MgO content
arcpy.AddMessage('MgO:')
arcpy.AddMessage(MgO)
SoilVol=arr_db[soilvol]# soil available to erode
arcpy.AddMessage('soilvol:')
arcpy.AddMessage(SoilVol)

#................................ D16..........................................
D16=-7.21+(2.83*np.log(T))-(0.03*Cp)-(0.03*Sp)+(0.04*Rd)
D16=2**D16 
#................................ D50..........................................
D50=5.52-(0.036*CS)+(1.17*np.log(T))-(0.35*K)
D50=2**D50

#................................ D84..........................................
D84=4.83+(0.66*np.log(T))-(0.004*Ro)-(0.27*K)+(0.16*np.sqrt(Sp))
D84=2**D84

#................................D84B..........................................
D84B=3.75-(0.002*WI)+(0.022*np.sqrt(Sp))+(0.093*MgO)+(0.75*np.log(E))
D84B=2**D84B

#...............................volume.........................................
# soil orgain matter; Slope # area burned at mod to high; runoff
lnV=2.7+(1.9*Om)+(0.175*np.sqrt(Sa))+(0.8*np.sqrt(Bmh))+(0.003*Ro)
Vol = np.exp(lnV) #Debris Flow Volume, m3
arcpy.AddMessage('lnV:')
arcpy.AddMessage(lnV)
#...............................join calculations back.........................
# set those with a probability below the threshold to nan/0 m3
if useProb:
    if P_thresh>0:
        P=arr_db[P_field]
        Vol[P<P_thresh] = np.nan; #Did not generate debris flow
        D16[P<P_thresh] = np.nan; #Did not generate debris flow
        D50[P<P_thresh] = np.nan; #Did not generate debris flow
        D84[P<P_thresh] = np.nan; #Did not generate debris flow
        D84B[P<P_thresh] = np.nan; #Did not generate debris flow

# final Vol Gartner output now make and adj version capping max at soil vol
Vol_Wall = Vol
# total volume of sediment eroded
#Total_Vol_Gartner = np.nansum(Vol_Gartner)


# flag if erroded more sediment than available making it supply limited
supply_lim=np.zeros(Vol.size)
supply_lim[Vol>SoilVol]=1

#Correction Factor
#If predicted volume, Vol, is > total soil volume, SoilVol, set Vol = SoilVol
#Vol[Vol>SoilVol]=SoilVol[Vol>SoilVol]

#print2(Frac_Input)
SoilLossFrac = Vol_Wall/SoilVol

# change all nan to zero
Vol_Wall[np.isnan(Vol_Wall)]=0
D16[np.isnan(D16)]=0
D50[np.isnan(D50)]=0
D84[np.isnan(D84)]=0
D84B[np.isnan(D84B)]=0

Vol[np.isnan(Vol)]=0
SoilLossFrac[np.isnan(SoilLossFrac)]=0



# Write Data back to shapefile
# Add fields of interest to the shape file
flds=['Vol','supp_lim','sl_vol','FrSoilLos','D16','D50','D84','D84B']

for fieldname in flds:
    arcpy.AddField_management(subcatch,fieldname, "Double")

# build 2D array to pull data from
outdata = np.vstack((Vol_Wall,supply_lim,Vol,SoilLossFrac,D16,D50,D84,D84B)).T

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
    