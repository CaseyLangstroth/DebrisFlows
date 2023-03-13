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
    outfolders=['Sev05i2Y','Sev50i5Y','Sev95i25Y']
    indir1=r'D:\PLI_Analysis\Outputs\a076N\a076'
    indir2t=r'D:\PLI_Analysis\Outputs\a076N\a076'
    indir2=os.path.join(indir2t,outfolders[2])
    # debris flow basins shapefile
    basename='a076'

    subcatchin=os.path.join(indir2,basename+'_subcatchments.shp')#os.path.join(indir2,basename+"_subcatchments.shp")
    # Attributed debris flow pour points
    sub_pp=os.path.join(indir2,basename+'_subcatchments_pp_snap_att.shp')
    #"D:\\Box Sync\\Wildfires\\Watershed_preprocessing\\Output_Data\\a000\\a000_debrisflow_pp_snap_netatt.shp"
    RivNetin=os.path.join(indir1,basename+"_network.shp")
    fdr=os.path.join(indir1,basename+"_fdr.tif")
    chan_ras=os.path.join(indir1,'temp','Subbasin',basename+"_chanras.tif")
    Rwidthfield="Riv_width"
    Rdepthfield="RMS_BFw_m"
    VBwidthfield="VB_width"

    out_temp_dir=os.path.join(indir2,'temp')
    
    makenewoutput=True
    subcatchout=[]#os.path.join(indir2,basename+'_subcatchmentsMurph.shp')    
    RivNetout=os.path.join(indir2,basename+'_netork.shp')    
    
    alpha = 60 # Average angle of Fan Splay, degrees
    #theta = 9 # Average Surface Slope of Fan, degrees
    K = 3.1 # Runout Scalar Fitting Parameter 
    levang=45 # levee angle
    bypass=False
    
if usegui==True:

    # debris flow basins shapefile
    subcatchin=arcpy.GetParameterAsText(0)

    # Attributed debris flow pour points
    sub_pp=arcpy.GetParameterAsText(1)

    RivNetin=arcpy.GetParameterAsText(2)
    
    fdr=arcpy.GetParameterAsText(3)
   
    Rwidthfield=arcpy.GetParameterAsText(4)
    VBwidthfield=arcpy.GetParameterAsText(5)
    Rdepthfield=arcpy.GetParameterAsText(6)

    basename=arcpy.GetParameterAsText(7)
    out_temp_dir=arcpy.GetParameterAsText(8)
    
    try:
        subcatchout=arcpy.GetParameterAsText(9)
    except:
        subcatchout=[]

    try:
        RivNetout=arcpy.GetParameterAsText(10)
    except:
        RivNetout=[]
    
    alpha = arcpy.GetParameter(11) # Average angle of Fan Splay, degrees
    #theta = 9 # Average Surface Slope of Fan, degrees
    K = arcpy.GetParameter(12) # Runout Scalar Fitting Parameter 
    levang=arcpy.GetParameter(13)  # levee angle
    bypass=False
    chan_ras=arcpy.GetParameterAsText(14)
    
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


def extractflowpath(XSind,YSind,FD,maxdist):
    # setup intial starting point
    Xind=XSind.copy()
    Yind=YSind.copy()
    totaldist=0# initial total distance 
    # output indices, initalize with the startig point
    iout=Xind.copy()
    jout=Yind.copy()
    FDi=0.0
    mFD=FD[Xind][Yind]# get the current flow dist at the starting point
    sFD=mFD.copy()# log the starting flow distance point
    
    outlet=np.nanmin(FD)
    while totaldist<maxdist and mFD>outlet:
        #print(totaldist)
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
        
        #update the indices
        Xind=istep
        Yind=jstep
    
        totaldist=np.sum(np.diff(FDi))
  
    return iout,jout,FDi,totaldist

# levee builder for when channnels completely aggrade
def buildlevees(b,Hd, sedthick,levang):
    
    # b = 3; %channel width
    # Hd = 2; % bankfull depth
    # sedthick = 10; %total depoisiton estimated at reach-----sed thickness
    # Dd = EstDeposition - Hd; %deposit depth relative to bankfull
    # A_avail = b*Dd; %Area available for levees + channel deposition
    # levang = 45; %levee slope angle (in degrees)
    
    # b = 3; %channel width
    # Hd = 2; % bankfull depth
    # EstDeposition = 10; %total depoisiton estimated at reach-----sed thickness
    levangrad=np.deg2rad(levang)
    Dd = sedthick - Hd#deposit depth relative to bankfull
    A_avail = b*Dd #Area available for levees + channel deposition    
    
    B = b*np.tan(levangrad) #%constant 1
    A = B*(Dd+Hd) #%constant 2
    Lh = 0.5*(np.sqrt((4*A)+B**2) - B)# %Levee height
    
    Lb = Lh/np.tan(levangrad)# %Levee width
    
    Ddf = Lh-Hd# Final channel deposit depth on top of bankfull
    
    Afinal = Ddf*b + (Lb*Lh)
    
    FinalAgg = Hd+Ddf
    LeveeFrac = (Lh*Lb)/(A_avail+(Hd*b))
    ChanFrac = (b*FinalAgg)/(A_avail+(Hd*b))
    return FinalAgg, LeveeFrac, Lb

#%% Murphy et al scenarios   
#....................... .....Scenario 1 ......................................
def do_scenario1(L_toe,runout,H_max,w,b,h,slope,L_fanbody,dx,alpha,theta,FD,RW_ras,
                 DownValley_Input,bypass):
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
    
    #print("I think this is wrong its a volume not a thickness but this is what I had in my notes")
    h_bf=h#h*b*L_river#river xsectional area *L_river--green bit step 3
    
    # z1; h_fbw is thickness of the wedge this is slide 4
    h_fbw = L_river*np.tan(slope)#h_fbw

    #total thickness of sediment 
    z=h_vw+h_fbw+h_bf
    
    #Max Possible Aggradation Wedge Below Fan
    PlanFanArea = 0.5*w*L_river
    MaxWedgeVol = PlanFanArea*h_fbw*0.5
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
        # downvalley routing 
        # calculate the volume we have to route downstream
        Vol_remain = Vol_Overrun-AggWedgeVol-AdjVol
        
        # assuming uniform geometry find unobstructed distance down stream it could go
        dist_downstream = Vol_remain/(0.5*b*z)
        
        # get the downvalley wedge angle beta using full downstream distance
        # that is could possibly run out assuming uniform channel width
        beta=np.arctan(z/dist_downstream)
        
        
        # find the indices of where the pour point falls in the 2D array
        # set true where flow distance and trib id are both equal to pp data
        Sarray=[[FD==FDid][0]*[trib==tribid][0]][0]
        
        # extract the x and y indices of the location to start the routing
        Sloc=np.where(Sarray==True)
        XSind=Sloc[0][0]
        YSind=Sloc[1][0]
        
        #extract indices and distances of path down stream,overdomain
        xi,yi,dst_down,totdst=extractflowpath(XSind,YSind,FD,dist_downstream)
        
        # determine if the debris flow over runs the domain
        #check if the total distance it could go is longer than the domain
        if totdst<dist_downstream:
            overdomain=1
            dst_down=np.append(dst_down,dist_downstream)
        else:
            overdomain=0
        
        # reverse the array so largest value is upstream and smallest is down
        # needed for wedge calculations so 0 thickness is downstream most which 
        # is calculated when d=0
        dst_up=np.max(dst_down)-dst_down
        dist_step=np.abs(np.diff(dst_up))
        
        # solve for midpoint distances between cells
        half_step=dist_step/2
        midpointdist=dst_up[0:-1]-half_step
        
        # build array at midpoints between cell centers and first and last values
        distprof=np.append(midpointdist,0)
        distprof=np.insert(distprof,0,dst_up[0])
        profdist=np.abs(np.diff(distprof))
        
        #construct the sediment profile
        sed_profile=np.tan(beta)*distprof
        sed_area=np.zeros(len(dst_up))
        dl=len(sed_area)
        
        for k in np.arange(0,dl):
            sed_area[k]=np.trapz([sed_profile[k+1],sed_profile[k]],x=[distprof[k+1],distprof[k]])
        
        # adjust where the sediment is if we are by passing 
        if overdomain:
            if bypass:
                # if letting sediment out of the domain dump the last cell
                sed_area=sed_area[:-1]
            else:
                # move the overrun sediment into the last cell in the domain
                sedoverflow=sed_area[-1]
                sed_area=sed_area[:-1]
                sed_area[-1]=sed_area[-1]+sedoverflow
            
            
        # conservation of mass checking
        tvol_init=np.max(np.cumsum(sed_area*b))
        
        # now account for varying channel widths
        # get channel widths in area of interest
        chanwidth=RW_ras[xi,yi]
        # if overdomain:
        #     chanwidth=np.append(chanwidth,chanwidth[-1])
        
        
        # account for variable channel width
        sed_vol=sed_area*chanwidth
        
        # recalculate new total volume
        adj_tot_vol=np.max(np.cumsum(sed_vol))
        
        # given channel width likely increase we are now likely predicting more
        # volume than what is available; need to truncate the debris flow
        #if  new vol > vol available find where sediment would end and truncate
        if adj_tot_vol>Vol_remain:
            cum_sed_vol=np.cumsum(sed_vol)
            idx_a=np.asarray(np.where(cum_sed_vol>Vol_remain))
            idx_a=idx_a[0]
            idx=np.min(idx_a)
            # replace last cell with any remaining sediment
            vt=Vol_remain-cum_sed_vol[idx-1]
            if vt>0:
                sed_vol[idx]=vt
                if len(idx_a)>1:
                    sed_vol=np.delete(sed_vol,idx_a[1:])
                    xi=np.delete(xi,idx_a[1:])
                    yi=np.delete(yi,idx_a[1:])
            else:
                sed_vol=np.delete(sed_vol,idx_a)
                xi=np.delete(xi,idx_a)
                yi=np.delete(yi,idx_a)
                
            # pick up here need to del all values past idx
            # idxlen=len(idx_a)
            # sed_vol[idx!=idx_a]
            
        # recalculate new total volume # conservation of mass check
        adj_tot_vol2=np.max(np.cumsum(sed_vol))
        # sediment routed down the river channel
        DownValley_Input[xi,yi]=DownValley_Input[xi,yi]+sed_vol
        
    # generate outputs
    # dumped locally into the river
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

if RivNetout:
    RivNet=RivNetout
    [Rtempdir,Rtempfid]=os.path.split(RivNet)
    arcpy.FeatureClassToFeatureClass_conversion(RivNetin,Rtempdir,Rtempfid)
else:
    RivNet=RivNetin


if subcatchout:
    subcatch=subcatchout
    [catchtempdir,catchtempfid]=os.path.split(subcatchout)
    arcpy.FeatureClassToFeatureClass_conversion(subcatchin,catchtempdir,catchtempfid)
else:
    subcatch=subcatchin


# dissovle to get each tributary as an individual feature
# tribshp=out_tmp_id+'test_trib_diss.shp'
# arcpy.management.Dissolve(RivNet, tribshp, None, None, 
#                           "SINGLE_PART", "DISSOLVE_LINES", '')
# flagfield="limbID"
# arcpy.AddField_management(tribshp, flagfield, "DOUBLE")
# arcpy.management.CalculateField(tribshp, flagfield, "!FID!", "PYTHON3")

# # convert river shp to raster
# chan_ras=out_tmp_id+"_chanras.tif"
# arcpy.conversion.PolylineToRaster(tribshp, flagfield,chan_ras, 
#                                   "MAXIMUM_LENGTH", "NONE", dx, "BUILD")

# convert river shp to raster
VBWidth_ras=out_tmp_id+"_VB_Widthras.tif"
arcpy.conversion.PolylineToRaster(RivNet, VBwidthfield,VBWidth_ras, 
                                  "MAXIMUM_LENGTH", "NONE", dx, "BUILD")
# convert river shp to raster
RivWidth_ras=out_tmp_id+"_riv_Widthras.tif"
arcpy.conversion.PolylineToRaster(RivNet, Rwidthfield,RivWidth_ras, 
                                  "MAXIMUM_LENGTH", "NONE", dx, "BUILD")
# convert river shp to raster
RivDepth_ras=out_tmp_id+"_Riv_depthras.tif"
arcpy.conversion.PolylineToRaster(RivNet, Rdepthfield,RivDepth_ras, 
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
    arr_sc = arcpy.da.FeatureClassToNumPyArray(subcatch,scflds, skip_nulls=True)
except:
    errmsg='Check all attributs are in the attributed subcatchment shapefile'
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
scflds=['GridID','ToLink','Slope','in_RES',Rwidthfield,VBwidthfield,Rdepthfield]
try:
    arr_RN = arcpy.da.FeatureClassToNumPyArray(RivNet, scflds ,skip_nulls=True)
except:
    errmsg='Check all attributs are in the river network shapefile'
    reporterror(errmsg,usegui)
    
slopeV=arr_RN['Slope']
gidV=arr_RN['GridID']
bV=arr_RN[Rwidthfield]
wV=arr_RN[VBwidthfield]
RESV=arr_RN['in_RES']
To_linkV=arr_RN['ToLink']
hV=arr_RN[Rdepthfield]


#--------------------------Set Parameters--------------------------------------
# Where valley widths, w, are < river widths, b, set w = b
wV[wV<bV]=bV[wV<bV]

#Unobstructed Fan Geometry
F = alpha/360.0 # Fraction of Cone [degrees/degrees]

#Rickenmann (1999) Runout Model:
R_V = K*(VolV**(1/3.)) #Debris Flow Runout Length [m] 
R_V[R_V==0]=np.nan
theta_temp_radians=np.arctan(VolV/(F*np.pi*(R_V**3)))
theta_temp_radians[np.isnan(theta_temp_radians)]=0
theta_V=3*np.rad2deg(theta_temp_radians) # vector of thetas


#%% ------------------------Run Calculations------------------------------------
#Volumetric Sediment Delivery Model
# Modified from Murphy et al., 2018
# initialize arrays for data writing

# debris flow parameters
L_overrun=np.zeros(R_V.size)
cond=np.zeros(R_V.size)
Vol_input=np.zeros(R_V.size)
dist_downstream=np.zeros(R_V.size)
#sub_id_out=np.zeros(R_V.size)*np.nan
#tdist=np.zeros(R_V.size)

# get 2D array of tributary ids
trib=arcpy.RasterToNumPyArray(chan_ras,nodata_to_value=-9999)*1.0
trib[trib==-9999]=np.nan

# get 2d array of valley bottom widths
VBw_ras=arcpy.RasterToNumPyArray(VBWidth_ras, nodata_to_value=np.nan)
RW_ras=arcpy.RasterToNumPyArray(RivWidth_ras, nodata_to_value=np.nan)
RD_ras=arcpy.RasterToNumPyArray(RivDepth_ras, nodata_to_value=np.nan)

# get flow distance 2D array
FD=arcpy.RasterToNumPyArray(FL_chan,nodata_to_value=np.nan)

# 2D array to house ouput volumes
DownValley_Input=FD*0

# Route the debris flow to the channel and solve for how much is input to the river
#For-loop of debris flow input nodes

DFrange=np.arange(0,R_V.size)
for i in DFrange:
    #print(i)
    #print(sc_id[i],'/',pp_FDid[i])
    theta=theta_V[i]
    Vol=VolV[i]
    R=R_V[i]
    
    
    
    # find GID start to extract relevant parameters from river network
    slope=slopeV[gidV==pp_gid[i]][0]
    gid=gidV[gidV==pp_gid[i]][0]
    b=bV[gidV==pp_gid[i]][0]
    w=wV[gidV==pp_gid[i]][0]
    h=hV[gidV==pp_gid[i]][0]
    RES=RESV[gidV==pp_gid[i]][0]
    tribid=pp_tid[i]
    FDid=pp_FDid[i]

    # print(RES)
    #To_link=To_linkV[i]
    
    # start the scenerio based caculations
    if RES==1:
        #Total Volume Input To Reservoir
        Vol_input[i] = Vol
        
    elif Vol==0:
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
                Vol_input[i], cond[i],DownValley_Input=do_scenario1(L_toe,runout, H_max,w,b,h,slope,L_fanbody,dx,alpha,theta,FD,RW_ras,DownValley_Input,bypass)
                
          # Scenario 3: If Fan Does Not Reach River
        else:
            #No calculations needed
            # just writes 0 to vol and 3 to cond
            [Vol_input[i], cond[i]]=do_scenario3()
            

#------------------- Manage ouputs back to shapefile------------------------------

# Downvalley Inputs Raster Exporting

# Get coordinate system for output rasters
dscRas = arcpy.Raster(FL_chan)
lowerLeft = arcpy.Point(dscRas.extent.XMin,dscRas.extent.YMin)
dsc=arcpy.Describe(dscRas)
cellSize=dscRas.meanCellWidth
coord_sys=dsc.spatialReference

# # write the raster
DownValley_Input[np.isnan(DownValley_Input)] = -9999
outRaster = arcpy.NumPyArrayToRaster(DownValley_Input,lowerLeft,cellSize, value_to_nodata=-9999)
outputid=basename+'_DownvalleySed.tif'
outputDEM=os.path.join(out_temp_dir,outputid)
outRaster.save(outputDEM)
del outRaster
arcpy.DefineProjection_management(outputDEM,coord_sys)

print('computing valley sums')
# # add the downvalley inputs to the river link
outtemp="DVSedm3"
temp_tbl=os.path.join(out_temp_dir,basename+'_totalDVsed.dbf')
ZonalStatisticsAsTable(RivNet,"FID", outputDEM, temp_tbl, "NODATA", "SUM")
arcpy.AddField_management(temp_tbl, outtemp, "DOUBLE")
arcpy.CalculateField_management(temp_tbl, outtemp, "!SUM!", "PYTHON", "")
arcpy.JoinField_management(RivNet,"FID", temp_tbl,"FID_", outtemp)


#....................... Get the data from river depostion....................
rivsedflds=[outtemp,Rwidthfield,Rdepthfield,"Length_m"]

try:
    arr_riv = arcpy.da.FeatureClassToNumPyArray( RivNet, rivsedflds, skip_nulls=True)
except:
    errmsg='Check all attributs are in the attributedrivernetwork'
    reporterror(errmsg,usegui)

# # Extract the values to arrays
r_w=arr_riv[Rwidthfield]
r_d=arr_riv[Rdepthfield]
r_l=arr_riv["Length_m"]
r_vol=arr_riv[outtemp]

r_sed_thick=r_vol/(r_w*r_l)

# levee aggradation and final volume calculation
LinkVol=np.zeros(len(r_w))
LeveeFrac=np.zeros(len(r_w))
Lb=np.zeros(len(r_w)) # levee height
ST=np.zeros(len(r_w)) # sed thickness in river
for i in np.arange(0,len(r_sed_thick)):
    if r_sed_thick[i]>r_d[i]:
        [ST[i], LeveeFrac[i], Lb[i]]=buildlevees(r_w[i], r_d[i], r_sed_thick[i], levang)
    else:
        LeveeFrac[i]=0
        Lb[i]=0
        ST[i]=r_sed_thick[i]

#remove the volume of sediment dumped into levees
r_vol_adj=ST*r_w*r_l
#add these in: levee height; final bed agg after levee 
#------------------- Write  ouputs back to shapefile---------------------------

# # Add fields of interest to the shape 
# adding data to the sub-catchments
flds=['VolInput','Scenerio']

for fieldname in flds:
    arcpy.AddField_management(subcatch,fieldname, "Double")

# build 2D array to pull data from
outdata = np.vstack((Vol_input,cond)).T

jj=0
for fieldname in flds:
    with arcpy.da.UpdateCursor(subcatch,fieldname) as cursor: #Add data
        ii = 0
        for row in cursor:
            row[0]=outdata[ii,jj]
            cursor.updateRow(row)
            ii+=1
    jj+=1


# add data to the river network
#add these in: levee height; final bed agg after levee 


flds=['BedVol','BedAgg','LeveeW','LeveeFrac']

for fieldname in flds:
    arcpy.AddField_management(RivNet,fieldname, "Double")

# build 2D array to pull data from
outdata = np.vstack((r_vol_adj,ST,Lb,LeveeFrac)).T

jj=0
for fieldname in flds:
    with arcpy.da.UpdateCursor(RivNet,fieldname) as cursor: #Add data
        ii = 0
        for row in cursor:
            row[0]=outdata[ii,jj]
            cursor.updateRow(row)
            ii+=1
    jj+=1
