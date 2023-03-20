# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:33:41 2023

@author: clang
"""

import arcpy
from arcpy import env
import os


arcpy.env.overwriteOutput = True #allow overwrite

def print2(instr): 
    arcpy.AddMessage(instr)
  
        
def chk_mk_dir(in_dir):
    try:
        os.mkdir(in_dir)
    except:
        pass
# =============================================================================
# Begin Tool    
# =============================================================================
def ReachDisc (DFGroupLayer, streamshp, subcatch, buffdis, out_dir):
    
    path = os.path.join(out_dir, 'temp\\')
    chk_mk_dir(path)
    fpath = os.path.join(out_dir,'final_outputs\\')
    chk_mk_dir(fpath)
        
# =============================================================================
# Get the Debris Flows from the group layer ##############
# =============================================================================
    aprx = arcpy.mp.ArcGISProject('CURRENT')
    mps =aprx.listMaps()
    DF_lyrs =[]
    DF_Originals = []
    for i in mps:
        DF_lyrs = i.listLayers()
    for i in DF_lyrs:
        if DFGroupLayer in i.longName:
            if i.isGroupLayer== False:
                if i.name.endswith('Assumed_Original'):
                    DF_Originals.append(i)                
# =============================================================================
# Find DF IDs     
# =============================================================================
    for x in DF_Originals:
       print2('Processing:\n' + x.name)
       end = len(x.name) 
       ID1 = x.name[0] #first letter
       spindex = x.name.find('_') #find underscore
       ID2ind = spindex + 1 
       ID2 = x.name[ID2ind] #second letter
       ns1 = x.name[ID2ind+1:end] #look after second letter
       IDnum = ns1.find('_') # find second underscore 
       ID3 = ns1[IDnum+1] #site number
       DF_ID = ID1 + ID2 + ID3 
       if ns1[IDnum+2] == 'A': #if site is added
           ID4 = 'A'
           DF_ID = ID1 + ID2 + ID3 + ID4
           if ns1[IDnum+4] == '2': #if site is secondary/third df and added
               ID5 = '_' + ns1[IDnum+4:IDnum+8]
               DF_ID = ID1 + ID2 + ID3 + ID4 + ID5
       if ns1[IDnum+3] == '2': #if site is secondary/third df (no A)
           ID6 = '_' + ns1[IDnum+3:IDnum+7]
           DF_ID = ID1 + ID2 + ID3 + ID6
# =============================================================================
# Create buffer
# =============================================================================
       ########## full buffer and dissolve ###########################
       buff = buffdis + ' Meters'
       outputbuff = path + DF_ID + 'fullbuffer'
       arcpy.analysis.Buffer(x, outputbuff, buff , '', '', 'ALL')
       
       ########## ring buffer and dissolve ###########################
       outputringbuff = path + DF_ID + 'bufferring' 
       arcpy.analysis.Buffer(x, outputringbuff, buff , 'OUTSIDE_ONLY','', 'ALL')
# =============================================================================
# Delineate stream near DF        
# =============================================================================
       ################ clip streamshp to buffer ############################# 
       outputstream = path + DF_ID + 'outstream' 
       arcpy.analysis.Intersect([streamshp, outputbuff], outputstream)
       # buffer selects lines within DF fan but not on top of fan
       
       ############### delete tributary streams if applicable ################
       select = arcpy.management.SelectLayerByLocation(outputstream, 'INTERSECT', subcatch, '', 'NEW_SELECTION')
       arcpy.management.DeleteFeatures(select)
       dfall = path + DF_ID+ 'dflineall' 
       arcpy.management.Dissolve(outputstream, dfall)
     
       ############## split stream into polygon overlap ######################
       dfpoly = path + DF_ID + 'dfpolystream' 
       arcpy.analysis.Intersect([outputstream,x], dfpoly)
       
       ############## split stream at buffer overlap #########################
       dfbuff = path + DF_ID + 'buffstream' 
       arcpy.analysis.Intersect([outputstream, outputringbuff], dfbuff)
       dfbuffparts = path + DF_ID + 'buffparts'
       dfbuffpts = path + DF_ID+ 'dfbuffpts'
       # parts = []
       # with arcpy.da.SearchCursor(dfbuff, '*') as cursor:
       #     for row in cursor:
       #         select = arcpy.management.SelectLayerByAttribute(dfbuff)
       #         arcpy.conversion.ExportFeatures(select, dfbuffpts)
       #         parts.append(dfbuffpts)
       # arcpy.management.Merge(parts,dfbuffparts)
               
       arcpy.management.FeatureVerticesToPoints(dfbuff, dfbuffpts, 'BOTH_ENDS')
       arcpy.management.SplitLineAtPoint(dfbuff, dfbuffpts, dfbuffparts)
# ============================================================================
# Find Lines that are artifically long due to buffer        
# =============================================================================
       ############ add end points to full line ##############################
       dfendpts = path + DF_ID + 'dfendpts'
       arcpy.management.FeatureVerticesToPoints(dfall, dfendpts, 'BOTH_ENDS')
       
       ############ Find and delete line segment that touches edge of buffer #
       arcpy.analysis.Near(dfbuffparts, dfendpts)
       with arcpy.da.UpdateCursor(dfbuffparts, 'NEAR_DIST') as cursor:
           for row in cursor:
               if row[0] < 1:
                   cursor.deleteRow()
# =============================================================================
# Finally, get df reach line
# =============================================================================
       dfreachp = path + DF_ID + 'dfreachprts'
       arcpy.management.Merge([dfbuffparts, dfpoly], dfreachp)
       dfreach = path + DF_ID + 'dfreach'
       arcpy.management.Dissolve(dfreachp, dfreach)
       
# =============================================================================
# Get the details of the reach
# =============================================================================
       arcpy.management.AddGeometryAttributes(dfreach, 'LENGTH_GEODESIC', 'METERS')
       desc = arcpy.Describe(dfreach)
       fields = desc.fields
       fieldsname = []
       for i in fields:
           fieldsname.append(i.name)
       if 'Id' in fieldsname:
           arcpy.management.CalculateField(dfreach, 'Id', 2)
       else:
            arcpy.management.AddField(dfreach, 'Id', 'SHORT')
            arcpy.management.CalculateField(dfreach, 'Id', 2)
       dfreachl = []
       with arcpy.da.SearchCursor(dfreach, 'LENGTH_GEO') as cursor:
            for row in cursor:
                dfreachl.append(row[0])                 
# =============================================================================
# Begin Reach Discretization
# =============================================================================
       # dissolve original stream because it is split at tribs via USUAL code #
       dstream = path + 'dissolvedstr'
       arcpy.management.Dissolve(streamshp, dstream)
       # get up and downstream start points
       rchpts = path + DF_ID + 'rchpts'
       arcpy.management.FeatureVerticesToPoints(dfreach,rchpts, 'BOTH_ENDS')
       # split dissolved stream to reach points
       splitstr = path + DF_ID+ 'splitstr'
       arcpy.management.SplitLineAtPoint(dstream, rchpts, splitstr, '100 Meters')
       ## delete duplicate df reach
       select = arcpy.management.SelectLayerByLocation(splitstr, 'INTERSECT', x)
       arcpy.management.DeleteFeatures(select)
       # ### generate points to split stream at dfreach length
       distpts = path + DF_ID + 'distpts'
       for i in dfreachl:
           arcpy.management.GeneratePointsAlongLines(splitstr, distpts, 'DISTANCE', i, '', 'END_POINTS')
       # get start/end points of reach
       # arcpy.analysis.Near(distpts, rchpts)
       # with arcpy.da.UpdateCursor(distpts, 'NEAR_DIST') as cursor:
       #     for row in cursor:
       #         if row[0] < 1:
       #             cursor.deleteRow()
       spreach = path + DF_ID + 'spreaches'
       arcpy.management.SplitLineAtPoint(splitstr, distpts, spreach, '100 Meters')
       
      
    return
# This is used to execute code if the file was run but not imported
if __name__ == '__main__':
    
    # Tool parameter accessed with GetParameter or GetParameterAsText
    DFGroupLayer = arcpy.GetParameterAsText(0)
    streamshp = arcpy.GetParameterAsText(1)
    subcatch = arcpy.GetParameterAsText(2)
    buffdis = arcpy.GetParameterAsText(3)
    out_dir = arcpy.GetParameterAsText(4)
    
    ReachDisc(DFGroupLayer, streamshp, subcatch, buffdis, out_dir)
    