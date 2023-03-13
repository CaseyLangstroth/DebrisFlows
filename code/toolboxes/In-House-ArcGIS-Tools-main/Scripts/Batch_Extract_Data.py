# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 00:03:40 2022

@author: Scott
"""


import os
#from osgeo import gdal
import arcpy
from arcpy import env

# demfid=r'D:\Scott\PLI_Analysis\Input_Data\DEM\DEM.tif'

# # input aoi and directory
inputdatapath=r"C:\Users\Scott\Desktop\AOI_PP"
inpoly=os.path.join(inputdatapath,'PLI_AOIV2.shp')

indem=r'D:\Research\Projects\Wildfires\Westwide_data\WesternUS_DEM\Western_US_DEM.tif'

inputdir=r'D:\Research\Projects\Wildfires\Westwide_data\Western_US_Data'
outdir=r'D:\Research\Projects\Wildfires\Westwide_data\PLI_Inputs'
dx=10



#workspace=r'C:\Users\Scott\Desktop\test\trash'
# Recursively extracts all files of a particular extension
def getallfiles(extension,pdir):
    filepaths=[]
    for subdir, dirs, files in os.walk(pdir):
        for file in files:
            filepath = subdir + os.sep + file
            if filepath.endswith(extension):
                filepaths.append(filepath)
    return filepaths

# checks if directy exits if not it makes it
def chk_mk_dir(dir_path):
    try:
        os.mkdir(dir_path)
        print("All outputs will be written to: \n", dir_path)
    except:
        print("All outputs will be written to: \n", dir_path)

# #------------ Directory Management---------------------------------------------
#Check if a directory exists if not make it   
chk_mk_dir(outdir) 
demdir=os.path.join(outdir, 'DEM')
chk_mk_dir(demdir)

# #------------------Setup environmental parameters------------------------------

# # do not project on the first extraction this will cause striations in the dem
# # extract raster data
demfid=os.path.join(demdir,'DEMtemp.tif')
out_raster = arcpy.sa.ExtractByMask(indem, inpoly); 
out_raster.save(demfid)
del out_raster


# # set up the output projections based on input polygon
SR_aoi = arcpy.Describe(inpoly).spatialReference
SR=str(SR_aoi.Name)
SR=SR.replace("_"," ")
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(SR)

arcpy.env.cellSize = dx
arcpy.env.resamplingMethod = "BILINEAR"

# first extract dem and project it
demfidout=os.path.join(demdir,'DEM.tif')
# demfidout=demfid
arcpy.ProjectRaster_management(demfid, demfidout, inpoly, "BILINEAR", "10")

# make sure everything snaps to DEM
arcpy.env.extent=demfidout
arcpy.env.snapRaster=demfidout



# get all shape and tif files housed inside a folder/subfolders
shpfilepaths = getallfiles(".shp",inputdir)
tiffilepaths= getallfiles(".tif",inputdir)

for tiffile in tiffilepaths:
    print(tiffile)
    fid1=tiffile[len(inputdir):]
    fid=os.path.split(fid1)
    outsubdir=os.path.join(outdir,fid[0][1:])
    chk_mk_dir(outsubdir)
    outfid=os.path.join(outsubdir,fid[1])
    # extract raster data
    out_raster = arcpy.sa.ExtractByMask(tiffile, demfid); 
    out_raster.save(outfid)
    del fid,outfid,out_raster,outsubdir,fid1

for shpfile in shpfilepaths:
    print(shpfile)
    fid1=shpfile[len(inputdir):]
    fid=os.path.split(fid1)
    outsubdir=os.path.join(outdir,fid[0][1:])
    chk_mk_dir(outsubdir)
    outfid=os.path.join(outsubdir,fid[1])
    # extract raster data
    arcpy.Clip_analysis(shpfile,inpoly, outfid)
    del fid,outfid,outsubdir,fid1
    
    
