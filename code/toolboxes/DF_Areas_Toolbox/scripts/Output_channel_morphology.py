# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 12:33:48 2023

@author: clang
"""

import arcpy
import os
import csv
arcpy.env.overwriteOutput = True

def OutputData(param0, param1, param2, param3, param4, param5):
    ############################################### Get the Layers #########################################################
    reach_layers_all = []
    reach_layers = []
    aprx = arcpy.mp.ArcGISProject('CURRENT')
    
    m1 = aprx.listMaps()
    for i in m1:
            l1 = i.listLayers()
    reachtbls = []
    
    ################################################ Get IDs ##################################################################
    for i in l1:
        if i.name.endswith('reach'): #get df_id
            reachtbls.append(i)
    sitenames = []
    channel_id = []
    for i in reachtbls:
        end = len(i.name)
        eA = end - 5
        sitenames.append(i.name[0:eA])
    reachno = ['1','2','3']
    for i in sitenames:
        for y in reachno:
            channel_id.append(i + '_' + y) #get channel_id
    
    ################################################### Create Table #################################################################
    
    tbl = arcpy.management.CreateTable(r'C:\Users\clang\Documents\USU\DF_Areas\GIS\Project_Files\GDB\DataTables.gdb', 'channel_morphology')
    ## name fields according to SQL database fields
    text_fields = ['DF_ID', 'channel_ID']
    for i in text_fields:
        arcpy.management.AddField(tbl, i, 'TEXT')
    flt_fields = ['IMP_BFWm', 'IMP_BFDm', 'flow_depthm', 'D50phi', 'D50mm', 'river_lengthm', 'straight_lengthm', 'sinuosity', 'usdakm2', 'gradient', 'slope', 'q2m3', 'q5m3', 'q10m3','q50m3','q100m3','q500m3']
    for i in flt_fields:
        arcpy.management.AddField(tbl, i, 'FLOAT')