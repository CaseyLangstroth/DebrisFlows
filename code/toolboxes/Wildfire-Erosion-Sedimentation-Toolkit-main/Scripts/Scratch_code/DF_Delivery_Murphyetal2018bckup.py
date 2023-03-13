#--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import math
import numpy as np
import pandas as pd

#------------------Booleans & License Checkout----------------------------
arcpy.env.overwriteOutput=True
arcpy.CheckOutExtension("Spatial")
usegui=False

#--------------User Inputs----------------------------------
if usegui==False:
    
    indir2=r'D:\Box Sync\Wildfires\PLI_watersheds\Output_Data'
    # debris flow basins shapefile
    basename='i152yrmodburn'

    subcatchin=os.path.join(indir2,'west2',basename+"_subcatchmentsG14.shp")
    #"D:\\Box Sync\\Wildfires\\Watershed_preprocessing\\Output_Data\\a000\\a000_debrisflow_basins.shp"
    # Attributed debris flow pour points
    sub_pp=r"D:\Box Sync\Wildfires\PLI_watersheds\Output_Data\a000\a000_subcatchments_pp_snap_att.shp"
    #"D:\\Box Sync\\Wildfires\\Watershed_preprocessing\\Output_Data\\a000\\a000_debrisflow_pp_snap_netatt.shp"
    RivNet=r"D:\Box Sync\Wildfires\PLI_watersheds\Output_Data\a000\a000_network.shp"
    fdr=r"D:\Box Sync\Wildfires\PLI_watersheds\Output_Data\a000\a000fdr.tif"
    Rwidthfield="Riv_width"
    VBwidthfield="VB_width"

    basename='murphtest4'
    out_temp_dir=r'D:\test\murph4'
    alpha = 60 # Average angle of Fan Splay, degrees
    #theta = 9 # Average Surface Slope of Fan, degrees
    K = 3.1 # Runout Scalar Fitting Parameter 
if usegui==True:

    # debris flow basins shapefile
    subcatchin=arcpy.GetParameterAsText(0)

    # Attributed debris flow pour points
    sub_pp=arcpy.GetParameterAsText(1)

    RivNet=arcpy.GetParameterAsText(2)
    
    fdr=arcpy.GetParameterAsText(3)
   
    Rwidthfield=arcpy.GetParameterAsText(4)
    VBwidthfield=arcpy.GetParameterAsText(5)

    basename=arcpy.GetParameterAsText(6)
    out_temp_dir=arcpy.GetParameterAsText(7)
    
    alpha = arcpy.GetParameterAsText(8) # Average angle of Fan Splay, degrees
    #theta = 9 # Average Surface Slope of Fan, degrees
    K = arcpy.GetParameterAsText(9) # Runout Scalar Fitting Parameter 
    
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

def extractflowpath(XSind,YSind,FD,maxdist):
    
    # setup intial starting point
    Xind=XSind.copy()
    Yind=YSind.copy()
    tdist=0# initial total distance 
    # output indices, initalize with the startig point
    iout=Xind.copy()
    jout=Yind.copy()
    FDi=0.0
    mFD=FD[Xind][Yind]# get the current flow dist at the starting point
    sFD=mFD.copy()# log the starting flow distance point
    
    outlet=np.nanmin(FD)
    while tdist<maxdist and mFD>outlet:
        # build the window to look over # D8 walk
        mrange=np.arange(Xind-1,Xind+2)
        nrange=np.arange(Yind-1,Yind+2)
        
        # intialize temp variables to hold the data
        vtemp=[]# flow distance value at x,y index
        itemp=[]# x index
        jtemp=[]# y index
        # loop through combos/ do D8 walk and extract the values
        for i in mrange:
            for j in nrange:
                vtemp.append(FD[i,j])
                itemp.append(i)
                jtemp.append(j)
        
        # find the minimum flow distance value that will be the next downstream cell
        mFD=np.nanmin(vtemp)
        
        # get location/index of minimum value
        mloc=np.where(mFD==vtemp)
        istep=itemp[mloc[0][0]]
        jstep=jtemp[mloc[0][0]]
        
        # log the indices of interest for a output
        iout=np.append(iout,istep)
        jout=np.append(jout,jstep)
        
        # log the flow distance 
        FDi=np.append(FDi,sFD-mFD)
        tdist=np.sum(np.diff(FDi))
        
        #update the indices
        Xind=istep
        Yind=jstep
        
    # get distance between each cell
    dist_step=np.diff(FDi)
    return iout,jout,FDi,dist_step,tdist


#%% Murphy et al scenarios   
#....................... .....Scenario 1 ......................................
def do_scenario1(L_toe,runout,H_max,w,b,slope,L_fanbody,dx,alpha,theta,FD,
                 DownValley_Input,VBw_ras):
    cond = 1
    
    # step 1 solve for unobstructed fan geometry
    
    # Toe Volume
    t = 1-(L_toe/runout)
    Vol_toe = ((H_max*runout**2.0)/3.0)*((np.pi/2.0)-(2.0*t*np.sqrt(1.0-t**2.0))-(np.arcsin(t))+(t**3*np.arccosh(1.0/t)))

    #Fan Body Overrun Volume
    #discritize fan body
    L_step = np.flipud(np.arange(w,L_fanbody,dx))
    
    x = L_fanbody-L_step
    h_step = (x/np.cos(np.deg2rad(alpha/2.0)))*np.tan(np.deg2rad(theta))
    radius_step = L_step/np.cos(np.deg2rad(alpha/2))
    k = 1-(2.0*np.sin(np.deg2rad(alpha/4.0))**2.0)
    c = 2.0*L_step*np.tan(np.deg2rad(alpha/2.0))
    
    # step 2 solve for over run body volumes 
    Vol_hyperbola = (H_max-h_step)*radius_step*(np.sqrt(1-(k**2))-(k**2)*np.arccosh(1./k))*dx;
    Vol_base = h_step*c*dx
    Vol_step = Vol_hyperbola+Vol_base
    Vol_body = np.max(np.cumsum(Vol_step))

    Vol_Overrun = Vol_toe + Vol_body

    #Direct Fan Input Riverbed (pre aggradation calc)
    L_step = np.flipud(np.arange(w-b,w,dx))
    x = L_fanbody-L_step
    h_step = (x/np.cos(np.deg2rad(alpha/2)))*np.tan(np.deg2rad(theta))
    radius_step = L_step/np.cos(np.deg2rad(alpha/2))
    k = 1-(2*np.sin(np.deg2rad(alpha/4))**2)
    c = 2*L_step*np.tan(np.deg2rad(alpha/2))

    Vol_hyperbola = (H_max-h_step)*radius_step*(np.sqrt(1-(k**2))-(k**2)*np.arccosh(1/k))*dx
    Vol_base = h_step*c*dx
    Vol_step = Vol_hyperbola+Vol_base
    Vol_body = np.max(np.cumsum(Vol_step))

    DirectInput_Vol = Vol_body

    # Fan Length & Height on Opposite Valley Wall
    L_river = 2*w*np.tan(np.deg2rad(alpha/2))
    
    #step3 in adj workflow
    # z0 is sed thickness at valley wall
    h_vw = L_toe*np.tan(slope)
    
    # z1 is thickness is thickness of the wedge
    z1 = L_river*np.tan(slope)
    
    #total thickness of sediment 
    z=h_vw+z1
    
    
    #Max Possible Aggradation Wedge Below Fan
    PlanFanArea = 0.5*w*L_river
    MaxWedgeVol = PlanFanArea*z1*0.5
    AggWedgeVol = 0
    h_temp = 0.01

    #If Still Volume Remaining - Fill Completely 
    if MaxWedgeVol <= Vol_Overrun:
        AggWedgeVol = MaxWedgeVol.copy()
        wedge_h = z.copy()
        River_Agg_h = np.median(h_step)+z
        fillcont = 1 #continue filling 
    #If Volume Insufficient - Fill & Stop 
    else:
        while AggWedgeVol < Vol_Overrun:
            h_temp = h_temp + 0.01
            AggWedgeVol = PlanFanArea*h_temp*0.5
        AggWedgeVol = PlanFanArea*(h_temp-0.01)*0.5 #take one step back down
        wedge_h = h_temp-0.01
        River_Agg_h = np.median(h_step)+h_temp-0.01
        fillcont = 0 #stop filling
        
    River_Wedge_Vol = 0.5*wedge_h*b*L_river

    #If Volume Remaining, Fill Down Valley Adjacent Wedge
    if fillcont == 1:
        adj_base_vol = 0.5 * (L_river/2)*w*z
        adj_top_vol = 0.25 * w*(L_river/2.0)*(H_max-z)
        MaxAdjVol = adj_base_vol+adj_top_vol
       
        AdjVol = 0
        h_temp = 0.01
        
        #If Still Volume Remaining - Fill Completely
        if MaxAdjVol <= (Vol_Overrun-AggWedgeVol):
            AdjVol = MaxAdjVol.copy()
            adj_wedge_h = z.copy()
            fillcont = 1 #continue filling
        
        #If Volume Insufficient - Fill & Stop 
        else:
            while AdjVol < (Vol_Overrun-AggWedgeVol):
                h_temp = h_temp + 0.01
                adj_base_vol = 0.5 * (L_river/2)*w*h_temp
                adj_top_vol = 0.25 * w*(L_river/2)*(H_max-h_temp)
                AdjVol = adj_base_vol+adj_top_vol

            adj_base_vol = 0.5 * (L_river/2)*w*(h_temp-0.01)
            adj_top_vol = 0.25 * w*(L_river/2)*(H_max-(h_temp-0.01)) #take one step back down
            AdjVol = adj_base_vol+adj_top_vol
            adj_wedge_h = h_temp-0.01
            fillcont = 0 #stop filling 

        #Input to River
        AdjVol_base_Input = 0.5*(b**2)*np.tan(np.deg2rad(alpha/2))*adj_wedge_h # Adjacent Wedge Base in River
        AdjVol_top_Input = 0.25*(b**3)*np.tan(np.deg2rad(theta))*np.tan(np.deg2rad(alpha/2.0)) #Adjacent Wedge Top in River
        AdjVol_Input = AdjVol_base_Input+AdjVol_top_Input
       
    else:
        AdjVol_Input = 0

    #If Volume Remaining, Fill Down Valley Bottom
    if fillcont == 1:
        # no downvalley routing
            # Vol_remain = Vol_Overrun-AggWedgeVol-AdjVol
            # dist_downstream = Vol_remain/(0.5*w*z)
            # DownValley_Input = 0.5*z*dist_downstream*b
            # Vol_input = DirectInput_Vol+River_Wedge_Vol+AdjVol_Input+DownValley_Input

            
        # downvalley routing V1  
        Vol_remain = Vol_Overrun-AggWedgeVol-AdjVol
        # distance dowstream to route the sediment
        dist_downstream = Vol_remain/(0.5*b*z)
        
        
        # find the indices of where the pour point falls in the 2D array
        # set true where flow distance and trib id are both equal to pp data
        Sarray=[[FD==FDid][0]*[trib==tribid][0]][0]
        
        # extract the x and y indices of the location to start the routing
        Sloc=np.where(Sarray==True)
        XSind=Sloc[0][0]
        YSind=Sloc[1][0]
        #extract indices and distances of path down stream
        xi,yi,dst_down,dst,tdst=extractflowpath(XSind,YSind,FD,dist_downstream)
        
        # get the downvalley wedge angle beta using full downstream distance
        beta=np.arctan(z/dist_downstream)
        
       #  # if hits the end of the domain append the full distance traveled to 
       #  # cummulative celluslar distance 
       #  if tdst<dist_downstream:
       #      dst_down=np.append(dst_down,dist_downstream)
        
       #  # build the profile of the wedge e.g. how thick is the sediment
       #  sed_thick=np.tan(beta)*dst_down
        
       #  dl=len(sed_profile)
       # # Ased=np.zeros(dl)
       #  Ased=np.zeros(dl)
       #  print('change to river width')
       #  VBi=VBw_ras[xi,yi]
        
       #  for k in np.arange(0,dl-1):  
       #    sed_thick[k+1]=np.trapz([sed_profile[k],sed_profile[k+1]],[dst_down[k],dst_down[k+1]])
            
        # old version to compute area and convert to sediment volume
        #     Ased[k+1]=np.trapz([sed_profile[k],sed_profile[k+1]],[dst_down[k],dst_down[k+1]])

        #     #print(k)
        #     #print(sed_profile[k])
        # vol_sed_down=Ased*VBi
        # #Input to River*\
        # DownValley_Input[xi,yi] = DownValley_Input[xi,yi]+vol_sed_down
        #DVi=DownValley_Input.copy()
    #Total Volume Input To River
    Vol_input = DirectInput_Vol+River_Wedge_Vol+AdjVol_Input#+DownValley_Input
    
    return Vol_input, cond, DownValley_Input#, DownValley_Input

#%%..........................Scenario 2.1......................................
#[Vol_input[i], cond[i]]=do_scenario2_1(L_toe,runout,H_max,L_overrun[i],L_fanbody,theta,alpha,dx):
def do_scenario2_1(L_toe,runout,H_max,L_overrun,L_fanbody,theta,alpha,dx):
    cond= 2.1
    # Toe Volume
    t = 1-(L_toe/runout)
    Vol_toe = ((H_max*runout**2)/3)*((np.pi/2)-(2*t*np.sqrt(1-t**2))-(np.arcsin(t))+(t**3*np.arccosh(1.0/t)))

    # Body Volume
    x_L = L_overrun-L_toe
    L_step = np.flipud(np.arange(L_fanbody-x_L,L_fanbody,dx))
    x = L_fanbody-L_step
    h_step = (x/np.cos(np.deg2rad(alpha/2)))*np.tan(np.deg2rad(theta));
    radius_step = L_step/np.cos(np.deg2rad(alpha/2))
    k = 1-(2*np.sin(np.deg2rad(alpha/4))**2)
    c = 2*L_step*np.tan(np.deg2rad(alpha/2))

    Vol_hyperbola = (H_max-h_step)*radius_step*(np.sqrt(1-(k**2))-(k**2)*np.arccosh(1/k))*dx
    Vol_base = h_step*c*dx
    Vol_step = Vol_hyperbola+Vol_base
    Vol_body = np.max(np.cumsum(Vol_step))
    Vol_input=Vol_body+Vol_toe
    return Vol_input, cond

#%%..........................Scenario 2.2......................................
#[Vol_input[i], cond[i]]=do_scenario2_2(H_max,L_overrun[i],runout)  
def do_scenario2_2(H_max,L_overrun,runout):
    cond = 2.2
    t=1-(L_overrun/runout)
    Vol_toe = ((H_max*runout**2)/3)*((np.pi/2)-(2*t*np.sqrt(1-t**2))-(np.arcsin(t))+(t**3*np.arccosh(1/t)))
    #Total Volume Input To River
    return Vol_toe, cond #Vol_toe=Vol_input

#%%...........................Scenario 3.......................................
#[Vol_input[i], cond[i]]=do_scenario3()
def do_scenario3():
    #No calculations needed
    cond = 3.0
    #Total Volume Input To River
    Vol_input= 0.0
    return Vol_input, cond



#%%---------------------------Begin Analysis-----------------------------------
chk_mk_dir(out_temp_dir)

out_tmp_id=os.path.join(out_temp_dir,basename)
#---------------- prepare down river routing data------------------------------

# build raster mask of river channel and assign a value to each tributary
# cell size x direction
dx_temp=arcpy.GetRasterProperties_management(fdr, "CELLSIZEX")
dx=float(dx_temp.getOutput(0))
arcpy.env.extent=fdr
arcpy.env.snapRaster=fdr
arcpy.env.cellSize=fdr




# dissovle to get each tributary as an individual feature
tribshp=out_tmp_id+'test_trib_diss.shp'
arcpy.management.Dissolve(RivNet, tribshp, None, None, 
                          "SINGLE_PART", "DISSOLVE_LINES", '')
flagfield="limbID"
arcpy.AddField_management(tribshp, flagfield, "DOUBLE")
arcpy.management.CalculateField(tribshp, flagfield, "!FID!", "PYTHON3")

# convert river shp to raster
chan_ras=out_tmp_id+"_chanras.tif"
arcpy.conversion.PolylineToRaster(tribshp, flagfield,chan_ras, 
                                  "MAXIMUM_LENGTH", "NONE", dx, "BUILD")

# convert river shp to raster
VBWidth_ras=out_tmp_id+"_VB_Widthras.tif"
arcpy.conversion.PolylineToRaster(RivNet, VBwidthfield,VBWidth_ras, 
                                  "MAXIMUM_LENGTH", "NONE", dx, "BUILD")


chan_points=out_tmp_id+"_chanpoints.shp"
arcpy.conversion.RasterToPoint(chan_ras, chan_points)

FL = arcpy.sa.FlowLength(fdr, "DOWNSTREAM", None); 
#FL.save(out_tmp_id+'allflowlenth.tif')
FLrd = arcpy.ia.RoundDown(FL); 
#FLrd.save(out_tmp_id+'allflowlenthRD.tif')

outras = arcpy.sa.ExtractByMask(FLrd, chan_ras, "INSIDE"); 
FL_chan=out_tmp_id+'_flow_length.tif'
outras.save(FL_chan)

del FL, FLrd, outras

# get the channel id and flow length field appended to the subcatch pour points

newppfid=basename+"_ppin.shp"
arcpy.FeatureClassToFeatureClass_conversion(sub_pp, out_temp_dir, 
    newppfid)
tmp_pp=os.path.join(out_temp_dir,newppfid)
snapenv=[chan_points,"VERTEX",str(dx*2)+' Unknown']# unknown sets it to map units
arcpy.edit.Snap(tmp_pp, [snapenv])

# add tributarty id 
arcpy.sa.AddSurfaceInformation(tmp_pp,chan_ras, "Z", "BILINEAR", None, 1, 0,'')
tribfield="trib_id"
arcpy.AddField_management(tmp_pp, tribfield, "Double")
arcpy.CalculateField_management(tmp_pp, tribfield, "!Z!","PYTHON","")
arcpy.management.DeleteField(tmp_pp, "Z")

#add flow distance field
arcpy.sa.AddSurfaceInformation(tmp_pp,FL_chan, "Z", "BILINEAR", None, 1, 0,'')
FDfield="flowdist"
arcpy.AddField_management(tmp_pp, FDfield, "Double")
arcpy.CalculateField_management(tmp_pp, FDfield, "!Z!","PYTHON","")
arcpy.management.DeleteField(tmp_pp, "Z")



#---------------------------Extract variables of interest----------------------
#.......................Get the data from subcatchments........................
scflds = ["Vol","sub_ID"]
try:
    arr_sc = arcpy.da.FeatureClassToNumPyArray(subcatchin,scflds, skip_nulls=True)
except:
    errmsg='Check all attributs are in the attributed river network shapefile'
    reporterror(errmsg,usegui)

# Extract the values to arrays
sc_id=arr_sc["sub_ID"]
VolV=arr_sc["Vol"]

#....................... Get the data from pour point nodes....................
ppflds=['sub_ID','GridID',tribfield,FDfield]
try:
    arr_pp = arcpy.da.FeatureClassToNumPyArray( tmp_pp, ppflds, skip_nulls=True)
except:
    errmsg='Check all attributs are in the attributed pour point shapefile'
    reporterror(errmsg,usegui)

# Extract the values to arrays
ppgidtemp=arr_pp['GridID']
ppSidtemp=arr_pp['sub_ID']
pptribtemp=arr_pp[tribfield]
ppFDtemp=arr_pp[FDfield]


#check the order of the pp vs subcatch data and adjust pp if needed
inds=np.where(ppSidtemp==sc_id,ppSidtemp,sc_id)

# reorder to pour points to match subcatchments
pp_gid=ppgidtemp[inds]
pp_Sid=ppSidtemp[inds]
pp_tid=pptribtemp[inds]
pp_FDid=ppFDtemp[inds]

# datachk = np.vstack((sc_id,VolV,pp_Sid,pp_gid,pp_FDid)).T

#............... laod in data from the river network array.....................
#'usarea_km2',
print('make user inputs with defualt correct values')
scflds=['GridID','ToLink','Slope','RES',Rwidthfield,VBwidthfield]
try:
    arr_RN = arcpy.da.FeatureClassToNumPyArray(RivNet, scflds ,skip_nulls=True)
except:
    errmsg='Check all attributs are in the attributed pour point shapefile'
    reporterror(errmsg,usegui)
    
slopeV=arr_RN['Slope']
gidV=arr_RN['GridID']
bV=arr_RN[Rwidthfield]
wV=arr_RN[VBwidthfield]
RESV=arr_RN['RES']
To_linkV=arr_RN['ToLink']


#--------------------------Set Parameters--------------------------------------
# Where valley widths, w, are < river widths, b, set w = b
wV[wV<bV]=bV[wV<bV]

#Unobstructed Fan Geometry
F = alpha/360.0 # Fraction of Cone [degrees/degrees]

#Rickenmann (1999) Runout Model:
R_V = K*(VolV**(1/3.)) #Debris Flow Runout Length [m] 

theta_temp_radians=np.arctan(VolV/(F*np.pi*(R_V**3)))
theta_V=3*np.rad2deg(theta_temp_radians) # vector of thetas


#%%------------------------Run Calculations------------------------------------
#Volumetric Sediment Delivery Model
# Modified from Murphy et al., 2018
# initialize arrays for data writing

# debris flow parameters
L_overrun=np.zeros(R_V.size)
cond=np.zeros(R_V.size)
Vol_input=np.zeros(R_V.size)
dist_downstream=np.zeros(R_V.size)
#sub_id_out=np.zeros(R_V.size)*np.nan
tdist=np.zeros(R_V.size)

# get 2D array of tributary ids
trib=arcpy.RasterToNumPyArray(chan_ras,nodata_to_value=-9999)*1.0
trib[trib==-9999]=np.nan

# get 2d array of valley bottom widths
VBw_ras=arcpy.RasterToNumPyArray(VBWidth_ras, nodata_to_value=np.nan)

# get flow distance 2D array
FD=arcpy.RasterToNumPyArray(FL_chan,nodata_to_value=np.nan)

# 2D array to house ouput volumes
DownValley_Input=FD*0

# Route the debris flow to the channel and solve for how much is input to the river
#For-loop of debris flow input nodes

DFrange=np.arange(0,R_V.size)
for i in DFrange:
    print(sc_id[i],'/',pp_FDid[i])
    theta=theta_V[i]
    Vol=VolV[i]
    R=R_V[i]
    
    
    
    # find GID start to extract relevant parameters from river network
    slope=slopeV[gidV==pp_gid[i]][0]
    gid=gidV[gidV==pp_gid[i]][0]
    b=bV[gidV==pp_gid[i]][0]
    w=wV[gidV==pp_gid[i]][0]
    RES=RESV[gidV==pp_gid[i]][0]
    tribid=pp_tid[i]
    FDid=pp_FDid[i]

    # print(RES)
    #To_link=To_linkV[i]
    
    # start the scenerio based caculations
    if RES==1:
        #Total Volume Input To Reservoir
        Vol_input[i] = Vol
    else:
        runout = R
        L_fanbody = runout*np.cos(np.deg2rad(alpha/2))
        H_max = runout*np.tan(np.deg2rad(theta))
        L_toe = runout-L_fanbody
        dx = 0.001 #[m]
        dp = 0.001 #[m]
    
        #If Fan Exists & Reaches River, Determine Length of Overrun
        if np.isnan(runout)==1 or runout<(w-b):
            L_overrun[i]= 0
        else:
            L_overrun[i] = runout-(w-b)
        
        #If Fan Reaches River
        if L_overrun[i] > 0: 
        
            #Scenario 2.2: Toe Reaches River
            if L_overrun[i] <= L_toe:
                [Vol_input[i], cond[i]]=do_scenario2_2(H_max,L_overrun[i],runout)   
            
            #Scenario 2.1: Fan  Reaches River But Fan Body Length < Valley Width
            elif L_overrun[i] > L_toe and L_fanbody <= w:
                #Total Volume Input To River
                [Vol_input[i], cond[i]]=do_scenario2_1(L_toe,runout,H_max,L_overrun[i],
                    L_fanbody,theta,alpha,dx)


            #Scenario 1: If Runout >> Valley Width
            elif L_fanbody > w:
                Vol_input[i], cond[i],DownValley_Input=do_scenario1(L_toe,runout,H_max,w,b,slope,L_fanbody,dx,alpha,theta,FD,DownValley_Input,VBw_ras)
                
          # Scenario 3: If Fan Does Not Reach River
        else:
            #No calculations needed
            # just writes 0 to vol and 3 to cond
            [Vol_input[i], cond[i]]=do_scenario3()
            


#------------------- Manage ouputs back to shapefile------------------------------
# # Downvalley Inputs Raster Exporting
# # Get coordinate system for output rasters
# dscRas = arcpy.Raster(FL_chan)
# lowerLeft = arcpy.Point(dscRas.extent.XMin,dscRas.extent.YMin)
# dsc=arcpy.Describe(dscRas)
# cellSize=dscRas.meanCellWidth
# coord_sys=dsc.spatialReference

# # write the raster
# DownValley_Input[np.isnan(DownValley_Input)] = -9999
# outRaster = arcpy.NumPyArrayToRaster(DownValley_Input,lowerLeft,cellSize, value_to_nodata=-9999)
# outputid=basename+'_DownvalleySed.tif'
# outputDEM=os.path.join(out_temp_dir,outputid)
# outRaster.save(outputDEM)
# del outRaster
# arcpy.DefineProjection_management(outputDEM,coord_sys)


# # add the downvalley inputs to the river link
# outtemp="DVSedm3"
# temp_tbl=os.path.join(out_temp_dir,basename+'_totalDVsed.dbf')
# ZonalStatisticsAsTable(RivNet,"FID", outputDEM, temp_tbl, "NODATA", "SUM")
# arcpy.AddField_management(temp_tbl, outtemp, "DOUBLE")
# arcpy.CalculateField_management(temp_tbl, outtemp, "!SUM!", "PYTHON", "")
# arcpy.JoinField_management(RivNet,"FID", temp_tbl,"FID_", outtemp)




# # Add fields of interest to the shape file
# flds=['VolInput','Scenerio']

# for fieldname in flds:
#     arcpy.AddField_management(subcatchin,fieldname, "Double")

# # build 2D array to pull data from
# outdata = np.vstack((Vol_input,cond)).T

# jj=0
# for fieldname in flds:
#     with arcpy.da.UpdateCursor(subcatchin,fieldname) as cursor: #Add data
#         ii = 0
#         for row in cursor:
#             row[0]=outdata[ii,jj]
#             cursor.updateRow(row)
#             ii+=1
#     jj+=1
