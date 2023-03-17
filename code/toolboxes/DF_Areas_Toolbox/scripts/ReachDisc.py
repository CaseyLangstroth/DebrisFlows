# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:33:41 2023

@author: clang
"""

import arcpy
from arcpy import env
import os


arcpy.env.overwriteOutput = True #allow overwrite

def ReachDisc (DFGroupLayer, streamshp):
        
### Get the Debris Flows from the group layer ##############
    aprx = arcpy.mp.ArcGISProject('CURRENT')
    mps =aprx.listMaps()
    DF_lyrs =[]
    for i in mps:
        DF_lyrs = i.listLayers()
    for i in DF_lyrs:
        if DFGroupLayer in i.longName:
            if i.isGroupLayer== False:
                DF_lyrs.append(i)
    for i in DFlyrs:
        if i.name.endswith('Assumed_Original'):
            DF_Originals.append(i)
            arcpy.addMessage(i.name)
            
    IDs = len(DF_Originals) ## get sequential ID's based on how many DF sites there are
    
################################ Dissolve Stream to seperate at tributaries #######################

    arcpy.management.Dissolve(streamshp, 'splitstream', '', '', 'FALSE')
     
########################### Make Buffer around DF polygons and get Stream lines #########################   
    
    for i in range(IDs):
        arcpy.AddMessage(i)
        for x in Df_Originals:
            outputbuff = str(i) + 'buffer'
            outputbuffdiss = str(i) + 'buffdiss'
            arcpy.analysis.Buffer(x, outputbuff, '8 Meters')
            arcpy.analysis.Dissolve(outputbuff, outputbuffdiss)
            ######################## Clip stream at buffer edge #########################
            outputline = str(i) + 'buffline'
            arcpy.analysis.Intersect(['splitstream', outputbuffdiss], outputline)
            ######################### get the lengths of lines at dfs ###########################
            arcpy.management.Intersect(['splitstream', x], 'dflines') 
            ######################## Add fields ##########################
            arcpy.management.AddField('dflines', 'DF_ID', 'SHORT')
            ###################### select the lines that are in the same df to give them the same df_ID#################
            select = arcpy.management.SelectLayerByLocation('dflines', 'INTERSECT', x)
            arcpy.management.CalculateField(select, 'DF_ID', i)
            with arcpy.da.SearchCursor('dflines', ['Shape_Length', 'DF_ID']) as searcher:
                linelengths = []
                for row in searcher:
                    if row[1] == i:
                        linelengths.append(row[0])
                if len(linelengths) > 1: #more than one line in DF removes trib channel from reach count
                    lines_sorted = sorted(linelengths) #sort by lowest to highest
                    line_short = lines_sorted[0] #first entry is shortest length
                    with arcpy.da.UpdateCursor('dflines', 'Shape_Length') as cursor:
                        for row in cursor:
                            if row[0] == line_short:
                                cursor.deleteRow() #remove these tribs from shapefile
            
    return
# This is used to execute code if the file was run but not imported
if __name__ == '__main__':
    
    # Tool parameter accessed with GetParameter or GetParameterAsText
    DFGroupLayer = arcpy.GetParameterAsText(0)
    streamshp = arcpy.GetParameterAsText(1)
    
    ReachDisc(DFGroupLayer, streamshp)
    