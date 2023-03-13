# RUSLE Post Wildfire
# Will-perform a RUSLE analysis in all basins that did not generate a debris flow
# A=R*K*LS*C*P

#--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import math
import numpy as np
from scipy import signal
# import matplotlib.pyplot as plt
# import glob

#--------------Booleans & License Checkout----------------------------
arcpy.env.overwriteOutput=True
arcpy.CheckOutExtension("Spatial")
usegui=True
usethresholds=True
P_thresh=0.5
#--------------User Inputs--------------------------------
if usegui==False:
    indir1=r'D:\PLI_Analysis\Input_Data'
    indir2=r'D:\PLI_Analysis\TestOutput\a060N\a060'
    indir3=r'D:\PLI_Analysis\TestOutput\a060N\a060\Sev05i2yr'
    basenamein='a060'
    
    #USUAL INPUTS
    # input filled DEM
    inDEM=os.path.join(indir2,basenamein+"_demf.tif")
    fac=os.path.join(indir2,basenamein+"_fac.tif")
    fdr=os.path.join(indir2,basenamein+"_fdr.tif")
    in_wtrshd=os.path.join(indir2,basenamein+"_watershed.shp")
    fire_perimeter=in_wtrshd
    
    interfluves=os.path.join(indir3,basenamein+"_f_interfluves.shp")
    fluveid="influv_ID"
    
    subcatchments=os.path.join(indir3,basenamein+"_subcatchments_S05i2yr.shp")
    subcatchid="sub_ID"
    probfield="Prob"
    statsgo=os.path.join(indir1,'Soil_Data','STATSGO_westernUS.shp')
    kfactfield="KFFACT"
    sthicfield="THICK"
    BDfield="Db3rdbar"
    psandfield="Sand"
    LCdata=os.path.join(indir1,'Landcover','nlcd_2019_land_cover.tif')
    ftype=os.path.join(indir1,'fireseverity','ftype','dnbr_05_class.tif')

    #Rainfalltype="Annual" #or "Event Rainfall" ; "Annual Erosivity"
    Rainfalltype="Event Rainfall"
    
    RI=os.path.join(indir1,'i15_mmhr','precip2yr15ma.tif')
    # R parameters
    duration=60 # min
    
    basename='i152yrLburn' 
    out_temp_dir=r'D:\PLI_Analysis\TestOutput\a060N\a060\temp'


if usegui==True:
    inDEM=arcpy.GetParameterAsText(0)
    fac=arcpy.GetParameterAsText(1)
    fdr=arcpy.GetParameterAsText(2)
    in_wtrshd=arcpy.GetParameterAsText(3)
    fire_perimeter=arcpy.GetParameterAsText(4)
    
    interfluvesin=arcpy.GetParameterAsText(5)
    fluveid=arcpy.GetParameterAsText(6)
    subcatchmentsin=arcpy.GetParameterAsText(7)
    subcatchid=arcpy.GetParameterAsText(8)
    probfield=arcpy.GetParameterAsText(9)   
    
    
    statsgo=arcpy.GetParameterAsText(10)
    kfactfield=arcpy.GetParameterAsText(11)
    sthicfield="THICK"
    BDfield="Db3rdbar"
    psandfield="Sand"
    

    LCdata=arcpy.GetParameterAsText(12)
    ftype=arcpy.GetParameterAsText(13)
    
    #Rainfalltype="Annual" #or "Event Rainfall" ; "Annual Erosivity"
    Rainfalltype=arcpy.GetParameterAsText(14)
    
    RI=arcpy.GetParameterAsText(15)
    # R parameters
    duration=arcpy.GetParameter(16) # min
    
    try:
        interfluvesout=arcpy.GetParameterAsText(17)
    except:
        interfluvesout=[]
    try:
        subcatchmentsout=arcpy.GetParameterAsText(18)
    except:
        subcatchmentsout=[]
    basename=arcpy.GetParameterAsText(19)
    out_temp_dir=arcpy.GetParameterAsText(20)
    
#----------------Functions-------------------------------------
# toggle between arcgis interface and command prompt printing
def print2(instr): 
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

def reporterror(errmsg):
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
#---------------------------Directory Managment--------------------------------
# Make a temporary directory to store all intermediate outputs
# Check if output directory exists-- if not make the directory
chk_mk_dir(out_temp_dir)
# print2('Final outputs will be written to: \n'+out_dir)
fid=basename

if interfluvesout:
    interfluves=interfluvesout
    [Itempdir,Itempfid]=os.path.split(interfluves)
    arcpy.FeatureClassToFeatureClass_conversion(interfluvesin,Itempdir,Itempfid)
else:
    interfluves=interfluvesin

if subcatchmentsout:
    subcatchments=subcatchmentsout
    [catchtempdir,catchtempfid]=os.path.split(subcatchmentsout)
    arcpy.FeatureClassToFeatureClass_conversion(subcatchmentsin,catchtempdir,catchtempfid)
    print2('not working')
else:
    subcatchments=subcatchmentsin
    print2('working')


#----------------------GET INPUT SIZES AND SET ENV PARMETERS------------------

#Extract input DEM Raster Resolution
# cell size x direction
dx_temp=arcpy.GetRasterProperties_management(inDEM, "CELLSIZEX")
dx_indem=float(dx_temp.getOutput(0))

# cell size y direction
dy_temp=arcpy.GetRasterProperties_management(inDEM, "CELLSIZEY")
dy_indem=float(dy_temp.getOutput(0))

#Set up spatial reference external data is likely differnt than the watershed
SR_wtrshd = arcpy.Describe(in_wtrshd).spatialReference
SR=str(SR_wtrshd.Name)
SR=SR.replace("_"," ")
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(SR)
arcpy.env.extent=inDEM
arcpy.env.snapRaster=inDEM
arcpy.env.cellSize=dx_indem

# Get coordinate system for output rasters
dscRas = arcpy.Raster(inDEM)
lowerLeft = arcpy.Point(dscRas.extent.XMin,dscRas.extent.YMin)
dsc=arcpy.Describe(dscRas)
cellSize=dscRas.meanCellWidth
coord_sys=dsc.spatialReference

# error checking on input RI
RITF = isfloat(RI)
if not RITF:
    if not arcpy.Exists(RI): 
        errmsg=('Your input rainfall intensity raster does not exist')
        reporterror(errmsg)

#------------------------------Begin Analysis----------------------------------

# get fire type raster
ftype_aoi=os.path.join(out_temp_dir,basename+"ftype.tif")
arcpy.gp.ExtractByMask_sa(ftype, inDEM, ftype_aoi)

# make 2D numpy array from the raster
fras=arcpy.RasterToNumPyArray(ftype_aoi,nodata_to_value=-9999)
fras=fras.astype(np.float32)
fras[fras==-9999]=np.nan

# fire type ids
nof=1.0 #no fire
lowf=2.0 #low severity fire
medf=3.0 #moderate severity
highf=4.0 # high severity
#extf=5.0 # extreme severity
    
    
    
#--------Soil Erodibility (K)-------------------------
#STATSGO Kffactor
# Clip soil data to watershed
statgoaoi=os.path.join(out_temp_dir,fid+"_statsgo.shp")
arcpy.Clip_analysis(statsgo, in_wtrshd, statgoaoi,"")

# Rasterize K factor from statsgo shape file
statgoaoiras_k=os.path.join(out_temp_dir,fid+"_kRUSLEraw.tif")
arcpy.PolygonToRaster_conversion(statgoaoi,kfactfield,statgoaoiras_k, 
                                 "CELL_CENTER", "NONE", str(dx_indem))

# Convert to numpy arraay
Ki=arcpy.RasterToNumPyArray(statgoaoiras_k,nodata_to_value=np.nan)
#Ki=Ki/7.59 # convert to SI uniits t ha h/(ha MJ mm) Renard et al., 1997


# create matrix of scalars
frask=fras.copy() # make matrix equal to fire type raster
frask[frask==nof]=1.0
frask[frask==lowf]=1.5
frask[frask==medf]=1.75
frask[frask>=highf]=2

# Scale K based on fire severity
K=Ki*frask
    
del frask, Ki


#-------------Slope-Length Factor (LS)-----------------
#Get the slope steepness of topogrpahy
outslope=arcpy.sa.Slope(inDEM, "DEGREE", "1", "PLANAR", "METER")
slp=os.path.join(out_temp_dir,fid+"slp.tif")
outslope.save(slp)
theta_d=arcpy.RasterToNumPyArray(slp,nodata_to_value=np.nan)
theta=np.deg2rad(theta_d)# convert degrees to radians

# Get the aspect of the topography
aspect=arcpy.sa.Aspect(inDEM)
asp=os.path.join(out_temp_dir,fid+"asp.tif")
aspect.save(asp)
alpha_d=arcpy.RasterToNumPyArray(asp,nodata_to_value=np.nan)
alpha=np.deg2rad(alpha_d)
# Get flow accumulation/upstream drainage area
dxfactmp=arcpy.GetRasterProperties_management(fac, "CELLSIZEX")
facdx=int(dxfactmp.getOutput(0))
A=arcpy.RasterToNumPyArray(fac,nodata_to_value=np.nan)*facdx*facdx

# Threshold data following Gannon et al 2019
#We capped theta at 55% slope to avoid extrapolating beyond the
# range of Nearings data (Theobald et al. 2010; Litschert et al. 2014)
percslope=np.tan(theta)*100
slopethreshold=np.arctan(55/100.0)
if usethresholds==True:
    theta[theta>slopethreshold]=slopethreshold

# area threshold
#When calculating LS, we capped A at 0.9 ha to approximate the 
#maximum hillslope length of 300 m as suggested in Renard et al. (1997).
if usethresholds==True:
    A[A>9000]=9000 # 1ha = 10000m2

# print2('test thresholds binary on and off over multiple basins')
# slope calculation
S=-1.5+(17/(1+np.exp(2.3-(6.1*np.sin(theta)))))

B=(np.sin(theta)/0.0896)/(3*np.sin(theta**0.8)+0.56)
X=np.abs(np.sin(alpha))+np.abs(np.cos(alpha))
m=B/(1+B)
D=facdx
LS=S*(((np.power((A+(D**2)),m+1)-np.power(A,m+1)))/
    (np.power(D,m+2)*np.power(X,m)*np.power(22.13,m)))

#We also limited LS values to the maximum value of 72.15 listed in Renard et al. (1997).
if usethresholds==True:
    LS[LS>72.15]=72.15
    
#----------------landcover (C)-------------------------
#naming for land cover data in area of interest
LCaoi=os.path.join(out_temp_dir,fid+"_LC.tif")
# extract and resample raster to demf resoltion
arcpy.gp.ExtractByMask_sa(LCdata, inDEM, LCaoi)

#print2('Double check the classes and defining values')
#print2('REVISIT AND REFINE THE CLASSIFICATION')
LC=arcpy.RasterToNumPyArray(LCaoi,nodata_to_value=-9999)
LC=LC.astype(float)
LC[LC==-9999]=np.nan


'''
Baseline C factor
values ranged from 0.001 to 0.003 for forests, 0.025 to 0.029 for
shrublands, 0.012 to 0.080 for grasslands, and up to 1.0 for the
rare barren areas disturbed by agriculture or mining. Barren
alpine areas above 2900 m elevation were assigned a C value of
0.002 owing to high rock cover.

'''
# Basline C values from Gannon et al. 2019
C=LC.copy()# make the matrix to house data
LCRC=C.copy()#matrix to house reclassified data

# consider making this a user input
# Reclassify landcover into groups
LCRC[(LC<35)]=1# OTHER NEED TO SUBDIVIDE BARREN OUT
LCRC[LC>=83]=1# wetlands 
LCRC[(LC>=2)*(LC<50)]=2# FOREST
LCRC[LC==52]=3# SHURBLAND
LCRC[(LC>70)*(LC<=81)]=4# GRASSLAND and pastures
LCRC[(LC>81)*(LC<=82)]=5# cultivated crops 


# Basline C values from Gannon et al. 2019
C[LCRC==1]=0.002# OTHER NEED TO SUBDIVIDE BARREN OUT
C[LCRC==2]=0.002# FOREST
C[LCRC==3]=0.027# SHURBLAND
C[LCRC==4]=0.046 # GRASSLAND 
C[LCRC==5]=1 # crops

# Further reclasify based on fire
# forested and burned
C[(fras==lowf)*(LC==2)]=0.01 
C[(fras==medf)*(LC==2)]=0.05
C[(fras>=highf)*(LC==2)]=0.2

# non forested and burned
C[(fras==lowf)*(LC!=2)]=C[(fras==lowf)*(LC!=2)]*1.2 
C[(fras==medf)*(LC!=2)]=C[(fras==medf)*(LC!=2)]*1.5
C[(fras>=highf)*(LC!=2)]=C[(fras>=highf)*(LC!=2)]*2

#------------------landuse (P)-------------------------
#Assuming no anti-erosion practice
P=fras*0.0+1.0 # make a matrix of ones 

# here is where we could add that in land practice
# also if we do maybe look for a better eqn 
#P = 0.2 + (0.03*S) Look in renard 1997 if we incorporate this 
# S= slope (%); maximum value that P can reach is 1.0.

#P[P>1]=1


#------------------soil thickness---------------------------------------------
# #e.g. can't erode more than you have
# # soil thickness raster
# statgoaoiras_ST=os.path.join(out_temp_dir,fid+"_soilthick.tif")
# arcpy.PolygonToRaster_conversion(statgoaoi,"THICK",statgoaoiras_ST,
#     "CELL_CENTER", "NONE", str(int(round(dx_indem))))
    
# # make numpy array of soil thickness
# Sthick=arcpy.RasterToNumPyArray(statgoaoiras_ST,nodata_to_value=-9999)
# Sthick=Sthick*dx_indem*dx_indem
# print2('check these units after fixing wall et al')

# #Get a bulk density to convert A to m^3
# statgoaoiras_den=os.path.join(out_temp_dir,fid+"_density.tif")
# arcpy.PolygonToRaster_conversion(statgoaoi,"Db3rdbar",statgoaoiras_den,
#     "CELL_CENTER", "NONE", str(int(round(dx_indem))))

# #Get a bulk density to convert A to m^3
# statgoaoiras_sand=os.path.join(out_temp_dir,fid+"_percSand.tif")
# arcpy.PolygonToRaster_conversion(statgoaoi,"Sand",statgoaoiras_sand,
#     "CELL_CENTER", "NONE", str(int(round(dx_indem))))

# rho=arcpy.RasterToNumPyArray(statgoaoiras_den,nodata_to_value=-9999)
# print2("Change this")
# rho[rho==0]=1e6
    
#--------------------Rainfall Runoff Erosivity (R)-----------------------------

# if using a input raster clip it to the dem extent
if not RITF: 
    #naming for rainfall intensity raster in area of interest
    [RIdir,RIid]=os.path.split(RI)
    RIaoi=os.path.join(out_temp_dir,fid+"_"+RIid+".tif")
    # extract and resample raster to demf resoltion
    arcpy.gp.ExtractByMask_sa(RI, inDEM, RIaoi)
    i=arcpy.RasterToNumPyArray(RIaoi,nodata_to_value=np.nan)

if RITF:
    i=Raster(inDEM)*0+RI
    
if Rainfalltype=="Event Rainfall":
    # convert to numpy array
    # i=i*10 # convert cm to mm
    # print2('Code is converting rainfall units from cm to mm!')
    # erositivity calculation
    #e = 0.119 + 0.0873log10(i)
    #e = 0.29[1 - 0.72 exp(-.05*i)] #Brown and Foster (1987)
    
    
    E=np.zeros(i.shape) # output array
    R=np.zeros(i.shape) # output array
    
    n,m=i.shape
    sigma=(duration/2)*0.34 # one standard deviation
    gprof=signal.gaussian(duration,sigma)
    dt=1/duration # duration in min.
    
    for x in range(0,n):
        for y in range(0,m):
            i15g=gprof*i[x,y]
            e = 0.29*(1-0.72*np.exp(-0.082*(i15g))) #McGregor et al. (1995)
            E[x,y]=np.nansum(e*i15g*dt)# Energy per storm MJ/ha
            R[x,y]=E[x,y]*np.max(i15g) # MJ mm /(ha hr)
            del i15g
    del E, e, gprof, i, sigma, n, m
    
# Renard and Freimund, 1994
#R2=0.0483*np.power(AP,1.61)
elif Rainfalltype=="Annual Rainfall": 
    R=587.8-(1.219*i+(0.004105*np.power(i,2)))
    #R=-823.8+(5.213*AP)
    #R=C*0.0+R2
    
elif Rainfalltype=="Annual Erosivity":
    print('clean up based on new user inputs')
    R=arcpy.RasterToNumPyArray(R_values,nodata_to_value=-9999)
else:
    errmsg='Input type (Annual, Event, or Annual Erosivity) was not specified'
    reporterror(errmsg)

#R=E*numstorms

# future maybe add in an option to input a 

#----------------------------soil data-----------------------------------------
#e.g. can't erode more than you have
# soil thickness raster
statgoaoiras_ST=os.path.join(out_temp_dir,fid+"_soilthick.tif")
arcpy.PolygonToRaster_conversion(statgoaoi,sthicfield,statgoaoiras_ST,
    "CELL_CENTER", "NONE", str(int(round(dx_indem))))
    
# make numpy array of soil thickness
Sthick=arcpy.RasterToNumPyArray(statgoaoiras_ST,nodata_to_value=-9999)
Sthick=Sthick*dx_indem*dx_indem
print2('check these units after fixing wall et al')

#Get a bulk density to convert A to m^3
statgoaoiras_den=os.path.join(out_temp_dir,fid+"_density.tif")
arcpy.PolygonToRaster_conversion(statgoaoi,BDfield,statgoaoiras_den,
    "CELL_CENTER", "NONE", str(int(round(dx_indem))))

#Get a bulk density to convert A to m^3
statgoaoiras_sand=os.path.join(out_temp_dir,fid+"_percSand.tif")
arcpy.PolygonToRaster_conversion(statgoaoi,psandfield,statgoaoiras_sand,
    "CELL_CENTER", "NONE", str(int(round(dx_indem))))

rho=arcpy.RasterToNumPyArray(statgoaoiras_den,nodata_to_value=-9999)
print2("Change this")
rho[rho==0]=1e6


#------------------Erosion (A)-------------------------

A=R*K*LS*C*P #Mg/ha-->Change to volume

A[A<0]=0
# Compute the Volume Eroded
vol_eroded=A*dx_indem*dx_indem /(1e4*rho)

#Threshold to soil thickness ----------------
# if eroding more soil than available threshold it
vol_eroded[vol_eroded>Sthick]=Sthick[vol_eroded>Sthick]
    
#------------------output rasters-------------------------

# Output Rasters
A[np.isnan(A)] = -9999
outRaster = arcpy.NumPyArrayToRaster(A,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=fid+'_A_RUSLE.tif'
outputDEM=os.path.join(out_temp_dir,outputid)
outRaster.save(outputDEM)
del outRaster
arcpy.DefineProjection_management(outputDEM,coord_sys) 

vol_eroded[np.isnan(vol_eroded)] = -9999
outRaster = arcpy.NumPyArrayToRaster(vol_eroded,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=fid+'_RUSLE_volume_m3.tif'
outputDEM=os.path.join(out_temp_dir,outputid)
RvolRas=outputDEM
outRaster.save(outputDEM)
del outRaster
arcpy.DefineProjection_management(outputDEM,coord_sys) 

R[np.isnan(R)] = -9999
outRaster = arcpy.NumPyArrayToRaster(R,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=fid+'_R_RUSLE.tif'
outputDEM=os.path.join(out_temp_dir,outputid)
outRaster.save(outputDEM)
del outRaster
arcpy.DefineProjection_management(outputDEM,coord_sys)

del A, R, vol_eroded
    
# End for loop here
P[np.isnan(P)] = -9999
outRaster = arcpy.NumPyArrayToRaster(P,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=fid+'_P_RUSLE.tif'
outputDEM=os.path.join(out_temp_dir,outputid)
outRaster.save(outputDEM)
del outRaster
arcpy.DefineProjection_management(outputDEM,coord_sys)

C[np.isnan(C)] = -9999
outRaster = arcpy.NumPyArrayToRaster(C,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=fid+'_C_RUSLE.tif'
outputDEM=os.path.join(out_temp_dir,outputid)
outRaster.save(outputDEM)
del outRaster
arcpy.DefineProjection_management(outputDEM,coord_sys) 

LS[np.isnan(LS)] = -9999
outRaster = arcpy.NumPyArrayToRaster(LS,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=fid+'_LS_RUSLE.tif'
outputDEM=os.path.join(out_temp_dir,outputid)
outRaster.save(outputDEM)
del outRaster
arcpy.DefineProjection_management(outputDEM,coord_sys) 

K[np.isnan(K)] = -9999
outRaster = arcpy.NumPyArrayToRaster(K,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=fid+'_K_RUSLE.tif'
outputDEM=os.path.join(out_temp_dir,outputid)
outRaster.save(outputDEM)
del outRaster
arcpy.DefineProjection_management(outputDEM,coord_sys) 
    
#%%-------------------Compute hsd----------------------------------------------

# merge subcatchments and interfluves
norivshp=os.path.join(out_temp_dir,basename+"noriv.shp")
arcpy.Merge_management(interfluves+";"+subcatchments, norivshp)


# get the flow direction around that excludes the river
hsfdr=os.path.join(out_temp_dir,basename+"hsfdr.tif")
fdrnoriv=arcpy.sa.ExtractByMask(fdr, norivshp)
fdrnoriv.save(hsfdr)

#Get cell size
dx_temp=arcpy.GetRasterProperties_management(hsfdr, "CELLSIZEX")# cell size x direction
dx_inFL=float(dx_temp.getOutput(0))
cell_diag=(dx_inFL**2+dx_inFL**2)**0.5

# # calculate flow accross each pixel
FLpix=arcpy.RasterToNumPyArray(hsfdr,nodata_to_value=-9999)
FLpix=FLpix.astype(float)
# get flow length accross each pixel
FLpix[(FLpix==2)+(FLpix==8)+(FLpix==32)+(FLpix==128)]=cell_diag
FLpix[(FLpix==1)+(FLpix==4)+(FLpix==16)+(FLpix==64)]=dx_inFL

# # Get the flow length for cells to the river
flowlen=os.path.join(out_temp_dir,basename+"flen.tif")
FL=arcpy.sa.FlowLength(hsfdr,"DOWNSTREAM", "")
FL.save(flowlen)
totalFL=arcpy.RasterToNumPyArray(flowlen,nodata_to_value=-9999)
totalFL=totalFL.astype(float)

totalFL[totalFL==-9999]=np.nan
FLpix[FLpix==-9999]=np.nan
LR=totalFL/FLpix
# hillslope delivery ratio Gannon et al. 2019
hSDR=10**(-0.56-(0.0094*LR))

# output georeferencing information
dscRas = arcpy.Raster(fdr)
lowerLeft = arcpy.Point(dscRas.extent.XMin,dscRas.extent.YMin)
dsc=arcpy.Describe(dscRas)
cellSize=dscRas.meanCellWidth
coord_sys=dsc.spatialReference

# ouput hillslope delivery ratio  to raster
hSDR[np.isnan(hSDR)] = -9999
outRaster = arcpy.NumPyArrayToRaster(hSDR,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=basename+'_hsdr.tif'
hSDRdem=os.path.join(out_temp_dir,outputid)
outRaster.save(hSDRdem)
del outRaster
arcpy.DefineProjection_management(hSDRdem,coord_sys)


#Adjust for a delivery ratio
delout=basename+'_sed_del.tif'
delivery=arcpy.sa.Times(hSDRdem,RvolRas)
sed_delivered_dem=os.path.join(out_temp_dir,delout)
delivery.save(sed_delivered_dem)


#%%------Run zonal stats to get the total amounts erroded and deliverd---------

#------------get percent sand in each polygon-------------------

# # this is created in the first RUSLE tool-->raster of % sand
# statgoaoiras_sand=os.path.join(out_temp_dir,fid+"_percSand.tif")
#statgoaoiras_sand
# # Zonal stats to get total erosion in each debris flow & hillslope basin
rusledebtbl=os.path.join(out_temp_dir,fid+'rusledeb.dbf')
rusledirtbl=os.path.join(out_temp_dir,fid+'ruslehill.dbf')

#Zonal stats to get sediment delivery to river in each debris flow & hillslope basin
# compute average sand concentration in pecent
arcpy.gp.ZonalStatisticsAsTable_sa(subcatchments,subcatchid, statgoaoiras_sand,
                                   rusledebtbl, "DATA", "Mean")

arcpy.gp.ZonalStatisticsAsTable_sa(interfluves,fluveid, statgoaoiras_sand,
                                   rusledirtbl, "DATA", "Mean")

fld='percSand'
arcpy.AddField_management(rusledirtbl, fld, "Double")
arcpy.AddField_management(rusledebtbl, fld, "Double")

# compute sum aka just copy it over
arcpy.CalculateField_management(rusledebtbl, fld, "!MEAN!","PYTHON","")
arcpy.CalculateField_management(rusledirtbl, fld, "!MEAN!","PYTHON","")

# Join the field back to the shapefiles
arcpy.JoinField_management(subcatchments,subcatchid,rusledebtbl,subcatchid,fld)
arcpy.JoinField_management(interfluves,fluveid,rusledirtbl,fluveid, fld)


#--------------------Caluculate the total eroded volumes-----------------------
vcatchtbl=os.path.join(out_temp_dir,fid+'ruslevolcatch.dbf')
vintertbl=os.path.join(out_temp_dir,fid+'ruslevolinter.dbf')

arcpy.gp.ZonalStatisticsAsTable_sa(subcatchments,subcatchid, RvolRas,
                                   vcatchtbl, "DATA", "SUM")
arcpy.gp.ZonalStatisticsAsTable_sa(interfluves,fluveid, RvolRas,
                                   vintertbl, "DATA", "SUM")

sumfld='RtotVolm3'
arcpy.AddField_management(vintertbl, sumfld, "Double")
arcpy.AddField_management(vcatchtbl, sumfld, "Double")

# compute sum aka just copy it over
arcpy.CalculateField_management(vcatchtbl, sumfld, "!SUM!","PYTHON","")
arcpy.CalculateField_management(vintertbl, sumfld, "!SUM!","PYTHON","")

# Join the field back to the shapefiles
arcpy.JoinField_management(subcatchments,subcatchid,vcatchtbl,subcatchid,sumfld)
arcpy.JoinField_management(interfluves,fluveid,vintertbl,fluveid,sumfld)

#-----------------Caluculate the total delivered volumes-----------------------
# compute fraction delivered to stream
dcatchtbl=os.path.join(out_temp_dir,fid+'rusledelcatch.dbf')
dintertbl=os.path.join(out_temp_dir,fid+'rusledelinter.dbf')
arcpy.gp.ZonalStatisticsAsTable_sa(subcatchments,subcatchid, sed_delivered_dem,dcatchtbl, "DATA", "SUM")
arcpy.gp.ZonalStatisticsAsTable_sa(interfluves,fluveid,sed_delivered_dem,dintertbl, "DATA", "SUM")

hsdfld='RdelVolm3'
arcpy.AddField_management(dintertbl, hsdfld, "Double")
arcpy.AddField_management(dcatchtbl, hsdfld, "Double")

# # compute sum aka just copy it over
arcpy.CalculateField_management(dcatchtbl, hsdfld, "!SUM!","PYTHON","")
arcpy.CalculateField_management(dintertbl, hsdfld, "!SUM!","PYTHON","")

# # Join the field back to the shapefiles
arcpy.JoinField_management(subcatchments,subcatchid,dcatchtbl, subcatchid, hsdfld)
arcpy.JoinField_management(interfluves,fluveid,dintertbl,fluveid, hsdfld)

#%%---------Flag unburned basins and force their values to zero----------------
#interfluves
outburn1=os.path.join(out_temp_dir,basename+"burnedinterfluves.shp")
inputshps=[interfluves,fire_perimeter]
arcpy.Intersect_analysis(inputshps, outburn1, "ALL", "", "INPUT")
arcpy.CalculateField_management(outburn1, "!Burned!=", "1","PYTHON","")
arcpy.JoinField_management(interfluves, fluveid, outburn1, fluveid, "Burned")

# subacatchamnets
outburn2=os.path.join(out_temp_dir,basename+"burnedsubcatch.shp")
inputshps=[subcatchments,fire_perimeter]
arcpy.Intersect_analysis(inputshps, outburn2, "ALL", "", "INPUT")
arcpy.CalculateField_management(outburn2, "!Burned!=", "1","PYTHON","")
arcpy.JoinField_management(subcatchments, subcatchid, outburn2, subcatchid,"Burned")

#force unburned values to zero 
# need to make a layer file to add the data to intermediately 
# lyrid=fid+"_RUSLEhill_lyr"
# arcpy.MakeFeatureLayer_management(interfluves, lyrid)
# arcpy.SelectLayerByAttribute_management(lyrid, "NEW_SELECTION", '"Burned"=0')
# arcpy.CalculateField_management(lyrid, sumfld, 0,"PYTHON","")
# arcpy.CalculateField_management(lyrid, hsdfld, 0,"PYTHON","")
# arcpy.SelectLayerByAttribute_management(lyrid, "CLEAR_SELECTION")

#filter out direct inputs that didnt burn
# need to make a layer file to add the data to intermediately 
lyrid2=fid+"_RUSLEdeb_lyr"
arcpy.MakeFeatureLayer_management(subcatchments, lyrid2)
arcpy.SelectLayerByAttribute_management(lyrid2, "NEW_SELECTION", '"Burned"=0')
arcpy.SelectLayerByAttribute_management(lyrid2, "ADD_TO_SELECTION", 
                                        probfield+'>'+str(P_thresh))
arcpy.CalculateField_management(lyrid2, sumfld, 0,"PYTHON","")
arcpy.CalculateField_management(lyrid2, hsdfld, 0,"PYTHON","")
arcpy.SelectLayerByAttribute_management(lyrid2, "CLEAR_SELECTION")