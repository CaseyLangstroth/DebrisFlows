
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

#--------------User Inputs--------------------------------
if usegui==False:

    indir3=r'D:\PLI_Analysis\TestOutput\a060N\a060\Sev05i2yr'
    basenamein='a060'
    
    #USUAL INPUTS
    
    interfluves=os.path.join(indir3,basenamein+"_f_interfluves.shp")
    interfluvespp=[]
    fluveid="influv_ID"
    
    subcatchments=os.path.join(indir3,basenamein+"_subcatchments_S05i2yr.shp")
    subcatchmentspp=[]
    subcatchid="sub_ID"
   


if usegui==True:
    inDEM=arcpy.GetParameterAsText(0)
    fac=arcpy.GetParameterAsText(1)
    fdr=arcpy.GetParameterAsText(2)
    in_wtrshd=arcpy.GetParameterAsText(3)

arcpy.management.JoinField(scpp, "sub_ID", sc, "sub_ID", None, "USE_FM","")

arcpy.management.JoinField(fluv, "influv_ID", sc, "influv_ID", None, "USE_FM","")

