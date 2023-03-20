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
def ReachDisc (DFGroupLayer, streamshp, pp, subcatch, buffdis, out_dir):
    
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
                    

    ################ Copy input stream to avoid editing it ################
    copystream = path + 'copystream'
    arcpy.conversion.ExportFeatures(streamshp, copystream)
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
       ################ clip copystream to buffer ############################# 
       outputstream = path + DF_ID + 'outstream' 
       arcpy.analysis.Intersect([outputbuff, copystream], outputstream)
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
       arcpy.management.MultipartToSinglepart(dfbuff, dfbuffparts)
       arcpy.management.FeatureVerticesToPoints(dfbuffparts, dfbuffpts, 'BOTH_ENDS')
# ============================================================================
# Find Lines that are artifically long due to buffer        
# =============================================================================
       ############ add end points to full line ##############################
       dfendpts = path + DF_ID + 'dfendpts'
       arcpy.management.FeatureVerticesToPoints(dfall, dfendpts, 'BOTH_ENDS')
       ############ Find and delete line segment that touches edge of buffer #
       ptbuff = path + DF_ID + 'ptbuffer'
       arcpy.analysis.Buffer(dfendpts, ptbuff, '1 Meter')
       select = arcpy.management.SelectLayerByLocation(dfbuffparts, 'INTERSECT', ptbuff, '', 'NEW_SELECTION')
       arcpy.management.DeleteFeatures(select)
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
       ##########remove all tribs and dissolve entire stream
       select = arcpy.management.SelectLayerByLocation(copystream, 'INTERSECT', subcatch, '', 'NEW_SELECTION')
       arcpy.management.DeleteFeatures(select)
       dstream = path + 'dissolvedstr'
       arcpy.management.Dissolve(copystream, dstream)
       # get up and downstream start points
       rchpts = path + DF_ID + 'rchpts'
       arcpy.management.FeatureVerticesToPoints(dfreach, rchpts, 'BOTH_ENDS')
       # split dissolved stream to reach points
       splitstr = path + DF_ID+ 'splitstr'
       arcpy.management.SplitLineAtPoint(dstream, rchpts, splitstr, '100 Meters')
       # delete duplicate df reach
       select = arcpy.management.SelectLayerByLocation(splitstr, 'INTERSECT', ptbuff, '', 'NEW_SELECTION', 'INVERT')
       arcpy.management.DeleteFeatures(select)
       # get start/end points of reach
       arcpy.analysis.Near(splitstr, pp)
       dists = []
       with arcpy.da.SearchCursor(splitstr, 'NEAR_DIST') as cursor:
           for row in cursor:
               dists.append(row[0])
       distsort = sorted(dists)
       downstreamr = path + DF_ID + 'dwnstrm'
       upstreamr = path + DF_ID + 'upstrm'
       with arcpy.da.SearchCursor(splitstr, 'NEAR_DIST') as cursor:
           for row in cursor:
               expressiondw = '"NEAR_DIST" ={}'.format(distsort[0])
               select = arcpy.management.SelectLayerByAttribute(splitstr, 'NEW_SELECTION', expressiondw)
               arcpy.conversion.ExportFeatures(select, downstreamr)
       with arcpy.da.SearchCursor(splitstr, 'NEAR_DIST') as cursor:
           for row in cursor:
               expressionup = '"NEAR_DIST" ={}'.format(distsort[1])
               select = arcpy.management.SelectLayerByAttribute(splitstr, 'NEW_SELECTION', expressionup)
               arcpy.conversion.ExportFeatures(select, upstreamr)
       # generate points to split downstream at dfreach length
       dwnpts = path + DF_ID + 'dwnpts'
       for i in dfreachl:
           arcpy.management.GeneratePointsAlongLines(downstreamr, dwnpts, 'DISTANCE', i)
       downstream = path + DF_ID + 'dwnreach'
       arcpy.management.SplitLineAtPoint(downstreamr, dwnpts, downstream)                               
       desc = arcpy.Describe(downstream)
       fields = desc.fields
       fieldsname = []
       for i in fields:
           fieldsname.append(i.name)
       if 'Id' in fieldsname:
           arcpy.management.CalculateField(downstream, 'Id', 3)
       else:
            arcpy.management.AddField(downstream, 'Id', 'SHORT')
            arcpy.management.CalculateField(downstream, 'Id', 3)
       select = arcpy.management.SelectLayerByLocation(downstream, 'WITHIN_A_DISTANCE', x, '1 Meter', '', 'INVERT')
       arcpy.management.DeleteFeatures(select)
       arcpy.management.AddGeometryAttributes(downstream, 'LENGTH_GEODESIC', 'METERS')
       # split upstream portion
       arcpy.edit.FlipLine(upstreamr)
       uppts = path + DF_ID + 'uppts'
       for i in dfreachl:
           arcpy.management.GeneratePointsAlongLines(upstreamr, uppts, 'DISTANCE', i)
       upstream = path + DF_ID + 'upreach'
       arcpy.management.SplitLineAtPoint(upstreamr, uppts, upstream, '100 Meters')                               
       desc = arcpy.Describe(upstream)
       fields = desc.fields
       fieldsname = []
       for i in fields:
           fieldsname.append(i.name)
       if 'Id' in fieldsname:
           arcpy.management.CalculateField(upstream, 'Id', 1)
       else:
            arcpy.management.AddField(upstream, 'Id', 'SHORT')
            arcpy.management.CalculateField(upstream, 'Id', 1)
       select = arcpy.management.SelectLayerByLocation(upstream, 'WITHIN_A_DISTANCE', x, '1 meter', 'NEW_SELECTION', 'INVERT')
       arcpy.management.DeleteFeatures(select)
       arcpy.management.AddGeometryAttributes(upstream, 'LENGTH_GEODESIC', 'METERS')
       ######### get end points of up and downstream reaches
       dpts = path + DF_ID + 'dpts'
       arcpy.management.FeatureVerticesToPoints(downstream, dpts, 'BOTH_ENDS')
       upts = path + DF_ID + 'upts'
       arcpy.management.FeatureVerticesToPoints(upstream, upts, 'BOTH_ENDS')
# ==============================================================================
# Finally merge final lines and points   
# =============================================================================
       reaches = fpath + DF_ID + 'reaches' 
       arcpy.management.Merge([dfreach, upstream, downstream], reaches)
       reachpts = fpath + DF_ID + 'reachpts'
       arcpy.management.Merge([upts, dpts], reachpts)
 
       
    return
# This is used to execute code if the file was run but not imported
if __name__ == '__main__':
    
    # Tool parameter accessed with GetParameter or GetParameterAsText
    DFGroupLayer = arcpy.GetParameterAsText(0)
    streamshp = arcpy.GetParameterAsText(1)
    pp = arcpy.GetParameterAsText(2)
    subcatch = arcpy.GetParameterAsText(3)
    buffdis = arcpy.GetParameterAsText(4)
    out_dir = arcpy.GetParameterAsText(5)
    
    ReachDisc(DFGroupLayer, streamshp, pp, subcatch, buffdis, out_dir)
    