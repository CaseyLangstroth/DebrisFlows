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