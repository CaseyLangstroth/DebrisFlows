
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
    Qfield="Q_2yr"
    CWfield="Riv_width"#ChannelWidthField
    Sfield="Slope"#Channel Slope
    
    n = 0.02 # Manning's n
    tau_crit = 0.04# Critical shields stress
    rho_w = 1000 #density of water, kg/m3
    rho_s = 2650 #density of sediment, kg/m3
    g=9.81# gravity
    phi = 90 # bank angle, degrees
    
    #outputs
    # outptus
    D50field="riv_D50"
    Hfield="BF_depth0"
    

if usegui:

    riv_net=arcpy.GetParameterAsText(0) # river network shapefile
    Qfield=arcpy.GetParameterAsText(1)
    CWfield=arcpy.GetParameterAsText(2)
    Sfield=arcpy.GetParameterAsText(3)
    
    n = arcpy.GetParameter(4) # Manning's n
    
    rho_w = arcpy.GetParameter(5) #density of water, kg/m3
    rho_s = arcpy.GetParameter(6) #density of sediment, kg/m3
    tau_crit =arcpy.GetParameter(7) # Critical shields stress
    phi = arcpy.GetParameter(8) # bank angle, degrees
    g=9.81# gravity

    #outputs
    # outptus
    D50field=arcpy.GetParameterAsText(9)
    Hfield=arcpy.GetParameterAsText(10)

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


#%% Extract data to arrays

fields=[Qfield,CWfield,Sfield]
arr =  arcpy.da.FeatureClassToNumPyArray(riv_net, fields, skip_nulls=True)

B = arr[CWfield] # bottom channel width, meters
S = arr[Sfield] # bed slope, -/-
Q = arr[Qfield] # discharge, cms

## Finding the minimum for a sol'n
SSG = (rho_s-rho_w)/rho_w# sumberged specific weight of sed

# Manning Eqn solver
x = 1/np.tan(np.deg2rad(phi))
y=np.arange(0.01,15,0.01)

h=np.zeros(np.size(B))
for i in np.arange(len(B)):
    # Solve for a range of Q
    Q_vec=((np.sqrt(S[i])*(((B[i] + x*y)*y)/(B[i]+2*y*np.sqrt(1+x**2)))**(2/3))/n) * ((B[i] + x*y)*y)
    
    # Find minimum difference between input Q and range of Q 
    Qs=np.abs(Q_vec-Q[i])
    mQ=min(Qs)
    # extract that depth values
    h[i]=y[Qs==mQ]

# Snyder model for D50
D50 = (n**(3/5) * Q**(3/5) * S**(7/10))/(SSG * tau_crit * B**(3/5)) # meters


#Ahammad 2021 gravel depth solver
# ks=2*D50
# H=((Q*(ks**(1/6)))/(8.1*B*np.sqrt(g*S)))**(3/5);

#%% write outputs back to the shapefile
# Write the output back to the shapefile 

flds=[D50field,Hfield]

for fld in flds:
    arcpy.AddField_management(riv_net,fld,"DOUBLE")

with arcpy.da.UpdateCursor(riv_net, D50field) as cursor: #Add data
    i = 0
    for row in cursor:
        row[0]=D50[i]
        cursor.updateRow(row)
        i+=1

with arcpy.da.UpdateCursor(riv_net, Hfield) as cursor: #Add data
    i = 0
    for row in cursor:
        row[0]=h[i]
        cursor.updateRow(row)
        i+=1