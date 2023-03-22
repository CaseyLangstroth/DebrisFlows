# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:33:41 2023

@author: clang
"""

import arcpy
from arcpy import env
import os
import numpy as np
import math
from arcpy.sa import *


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
def ReachDisc (DFGroupLayer, streamshp, pp, subcatch, buffdis, demf, fac, out_dir):
    
    path = os.path.join(out_dir, 'temp\\')
    chk_mk_dir(path)
    rpath = os.path.join(out_dir,'final_reaches\\')
    chk_mk_dir(rpath)
    ppath = os.path.join(out_dir, 'final_points\\')
    chk_mk_dir(ppath)
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
       arcpy.management.SplitLineAtPoint(downstreamr, dwnpts, downstream, '100 Meters')                               
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
       dpts = path + DF_ID + 'dpts'
       arcpy.management.FeatureVerticesToPoints(downstream, dpts, 'BOTH_ENDS')
       # arcpy.management.AddGeometryAttributes(downstream, 'LENGTH_GEODESIC', 'METERS')
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
       # arcpy.management.AddGeometryAttributes(upstream, 'LENGTH_GEODESIC', 'METERS')
       ######### get end points of up and downstream reaches
       # dpts = path + DF_ID + 'dpts'
       # arcpy.management.FeatureVerticesToPoints(downstream, dpts, 'END')
       
# =============================================================================
# Get and label all four points
# =============================================================================
       upt1 = path + DF_ID + 'upt1'
       arcpy.management.FeatureVerticesToPoints(upstream, upt1, 'END')
       #arcpy.management.AddField(upt1, 'Id', 'DOUBLE')
       arcpy.management.CalculateField(upt1, 'Id', 1)
       upt2 = path + DF_ID + 'upt2'
       arcpy.management.FeatureVerticesToPoints(upstream, upt2, 'START')
       #arcpy.management.AddField(upt2, 'Id', 'DOUBLE')
       arcpy.management.CalculateField(upt2, 'Id', 2)
       
       dpt1 = path + DF_ID + 'dpt1'
       arcpy.management.FeatureVerticesToPoints(downstream, dpt1, 'START')
       #arcpy.management.AddField(dpt1, 'Id', 'DOUBLE')
       arcpy.management.CalculateField(dpt1, 'Id', 3)
       dpt2 = path + DF_ID + 'dpt2'
       arcpy.management.FeatureVerticesToPoints(downstream, dpt2, 'END')
       #arcpy.management.AddField(dpt2, 'Id', 'DOUBLE')
       arcpy.management.CalculateField(dpt2, 'Id', 4)
# ==============================================================================
# Finally merge final lines and points   
# =============================================================================
       reaches = rpath + DF_ID + 'reaches' 
       arcpy.management.Merge([dfreach, upstream, downstream], reaches)
       reachpts = ppath+ DF_ID + 'reachpts'
       arcpy.management.Merge([upt1, upt2, dpt1, dpt2], reachpts)
 # =============================================================================
 #  Reproject reaches and points to UTM for USUAL code
 # =============================================================================
       rutm = rpath + DF_ID + 'reachesUTM'
       uputm = path + DF_ID + 'upreachUTM'
       dwnutm = path + DF_ID + 'downreachUTM'
       putm = ppath + DF_ID + 'reachptsUTM.shp'
       arcpy.management.Project(reaches, rutm, '26912')
       arcpy.management.Project(upstream, uputm, '26912')
       arcpy.management.Project(downstream, dwnutm, '26912')
       arcpy.management.Project(reachpts, putm, '26912')
# =============================================================================
# Get the DATA (modified from USUAL code)
# =============================================================================
       #.........Get Raster Cellsize............
       # cell size x direction
       dx_temp=arcpy.GetRasterProperties_management(fac, "CELLSIZEX")
       dx=float(dx_temp.getOutput(0))
       # cell size y direction
       dy_temp=arcpy.GetRasterProperties_management(fac, "CELLSIZEY")
       dy=float(dy_temp.getOutput(0))
       # cell diagonal length
       cell_diag=float(math.ceil(math.sqrt(dx**2 + dy**2))) 
       # get area of a cell
       cell_area=dx*dy
       #toloerance to to find points to split the line
       tol=cell_diag

       # #Add GridID to attribute table
       # arcpy.AddField_management(rutm, "GridID", "Double")
       # arcpy.CalculateField_management(rutm, "GridID", "!FID!+1","PYTHON","")
       # Add lengths to the table and delete old data
       arcpy.AddGeometryAttributes_management(rutm, "LENGTH_GEODESIC", "METERS")
       arcpy.AddField_management(rutm, "Length_m", "Double")
       arcpy.CalculateField_management(rutm, "Length_m", "!LENGTH_GEO!","PYTHON","")
       arcpy.DeleteField_management(rutm, ["LENGTH_GEO"])
        
       # Get the starting points of the line segments
       # ##flip this line back to original direction
       # arcpy.edit.FlipLine(uputm)
       # uptsutm = path + DF_ID + 'uptsutm'
       # arcpy.management.FeatureVerticesToPoints(uputm, uptsutm, 'START')
       # # # get all the end pointss of the line segment
       # dptsutm = path + DF_ID + 'dptsutm'
       # arcpy.management.FeatureVerticesToPoints(dwnutm, dptsutm, 'END')
       # #project to UTM
       # enutm = path + DF_ID + 'dptsUTM'
       # stutm = path + DF_ID + 'upptsUTM'
       # arcpy.management.Project(dpts, stutm, '26912')
       # arcpy.management.Project(upts, enutm, '26912')
       # get the max GridID value of start points
       # with arcpy.da.SearchCursor(putm, "Id") as cursor:
       #     for row in cursor:
       #         if row[0] == 1
       #             MaxValue = int(row[0])
       #         else:
       #             #MinValue = min(int(row[0]),MinValue)
       #             MaxValue = max(int(row[0]),MaxValue)
        
       # Create a point at the end of the stream
       # delete all endpoints that overlap with start points
       # arcpy.Erase_analysis(enutm, stutm, path + DF_ID + "endpt.shp", "2 Meters")
        
        
        
       # set grid id to +1 largest grid id of 
       # arcpy.AddField_management(dptsutm, "GridID", "Double")
       # arcpy.CalculateField_management(dptsutm, "GridID", MaxValue+1,"PYTHON","")
        
       # # merge start and end points and clean up the data
       # arcpy.Merge_management([dptsutm,uptsutm],path + DF_ID +"_nodes.shp")
       # arcpy.DeleteField_management(path + DF_ID+ "_nodes.shp", ["ORIG_FID"])#,"arcid","grid_code","from_node","to_node"])
       # arcpy.management.DefineProjection(path + DF_ID + '_nodes.shp', '26912') 
       
       # snap pour points 
       arcpy.gp.SnapPourPoint_sa(putm, fac, path +"_nodes_pp.tif", str(cell_diag), "Id")
       
       # convert raster to points
       arcpy.RasterToPoint_conversion(path +"_nodes_pp.tif", putm, raster_field="Value")
        
       #add flow accumulation and elevation data to the points
       arcpy.gp.ExtractMultiValuesToPoints_sa(putm, [[fac,"fac1"],[demf,"elev_m1"]], "NONE")
        
       # Join the two point shapefiles *** Maybe change this to spatial join
       # arcpy.JoinField_management(path + DF_ID+"_nodes.shp","GridID", path + DF_ID+"_nodes_pp.shp", "grid_code", "FID;pointid;grid_code;fac1;elev_m1")
        
       # export to final shape file of nodes
       # arcpy.FeatureClassToFeatureClass_conversion(path + DF_ID+"_nodes.shp", putm)
       # clean up the attribute table 
       # arcpy.DeleteField_management(putm, ["Length_m","pointid","grid_code"])

       
      
       arcpy.gp.ExtractMultiValuesToPoints_sa(putm, [[fac,"fac2"],[demf,"elev_m2"]], "NONE")
       arcpy.AddField_management(putm,"fac", "Double")
       arcpy.AddField_management(putm,"elev_m", "Double")
       calcDA="def calcDA( fac1 , fac2 , fac ):\n    if fac1==0:\n       return(fac2)\n    else:\n        return(fac1)"
       arcpy.CalculateField_management(putm,"fac","calcDA( !fac1! , !fac2! , !fac! )", "PYTHON_9.3",calcDA)
    
       calcE="def calcelev( elev_m1 , elev_m2 , elev_m ):\n    if elev_m1==0:\n       return(elev_m2)\n    else:\n        return(elev_m1)"
       arcpy.CalculateField_management(putm,"elev_m","calcDA( !elev_m1! , !elev_m2! , !elev_m! )", "PYTHON_9.3",calcE)
    
       arcpy.DeleteField_management(putm, ["elev_m1","elev_m2","fac1","fac2"])
    
    
       arcpy.AddField_management(putm, "usarea_m2", "Double")
       arcpy.CalculateField_management(putm, "usarea_m2","!fac!*"+str(cell_area),"PYTHON","")
    
       # join all the data
       arcpy.SpatialJoin_analysis(rutm, putm, path+"_split_join.shp", "JOIN_ONE_TO_MANY", "KEEP_ALL",'', "INTERSECT", "", "")
    
    
       # Extract attributes to arrays
       arr = arcpy.da.FeatureClassToNumPyArray(path+"_split_join.shp", ['Id','Length_m','elev_m','usarea_m2'], skip_nulls=True)
       gid=arr['Id']# line id
       l=arr['Length_m']
       elev=arr['elev_m']
       #nodeid=arr['GridID_1']#['Join_FID']
       usarea=arr['usarea_m2']
    
    
       # Caluculate slopes and to node--- this explicitly assume drainage area increases downstream rather than Jon's approach
       step=np.unique(gid)# get number of network segments
       step=step.astype(int)
    
       # allocate memory to write outputs to
       l_i=np.zeros(len(step))
       dz_i=np.zeros(len(step))
       us_elev=np.zeros(len(step))
       ds_elev=np.zeros(len(step))
       to_link=np.zeros(len(step))
       usarea_km2=np.zeros(len(step))
       maxda=np.nanmax(usarea)
       for i in step:
           tf=gid==i
           #print(i)
           #nid_temp=nodeid[tf]
           elev_temp=elev[tf]
           l_temp=l[tf]
           #nodeid_temp=nodeid[tf]
           usarea_temp=usarea[tf]
           #print2(i)
           # adding unique to help with resevoir data where values are the same
           us_elev[i-1]=np.unique(elev_temp[elev_temp==max(elev_temp)])
           ds_elev[i-1]=np.unique(elev_temp[elev_temp==min(elev_temp)])    
           # if np.max(usarea_temp)==maxda:
           #     to_link[i-1]=0
           # else:
           #     to_link[i-1]=np.min(nodeid_temp[usarea_temp==max(usarea_temp)])
           usarea_km2[i-1]=np.unique(min(usarea_temp)*1e-6)
           dz_i[i-1]=us_elev[i-1]-ds_elev[i-1]
           l_i[i-1]=l_temp[0]
             
       slope=dz_i/l_i
    
      # remove slopes that are likely unrealistically small
       slope[slope<0.001]=0.001
    
      # build 2D array to pull data from
       outdata = np.vstack((to_link,usarea_km2,us_elev,ds_elev,slope)).T
    
      #............................Data Ouputing.....................................
       print2('Writing all attributs to final ouput')
       arcpy.FeatureClassToFeatureClass_conversion(rutm,rpath,DF_ID+"_reaches.shp")
    
       rivnetwork=DF_ID+"_reaches.shp"
       # Add fields of interest to the shape file
       flds=['ToLink','usarea_km2','uselev_m','dselev_m','Slope']
       for fieldname in flds:
            arcpy.AddField_management(rivnetwork, fieldname, "Double")
             
       jj=0
       for fieldname in flds:
           with arcpy.da.UpdateCursor(rivnetwork,fieldname) as cursor: #Add data
               ii = 0
               for row in cursor:
                   row[0]=outdata[ii,jj]
                   cursor.updateRow(row)
                   ii+=1
           jj+=1
    
         # print2("Finished! Final outputs are located: \n"+out_dir,usegui)
    
             
         # if disp_outputs and usegui:
         #     aprx = arcpy.mp.ArcGISProject('current')
         #     cmap = aprx.listMaps()[0]  #data to be added to first map listed
         #     cmap.addDataFromPath(rivnetwork)  
    return
# This is used to execute code if the file was run but not imported
if __name__ == '__main__':
    
    # Tool parameter accessed with GetParameter or GetParameterAsText
    DFGroupLayer = arcpy.GetParameterAsText(0)
    streamshp = arcpy.GetParameterAsText(1)
    pp = arcpy.GetParameterAsText(2)
    subcatch = arcpy.GetParameterAsText(3)
    buffdis = arcpy.GetParameterAsText(4)
    demf = arcpy.GetParameterAsText(5)
    fac = arcpy.GetParameterAsText(6)
    out_dir = arcpy.GetParameterAsText(7)   
    
    ReachDisc(DFGroupLayer, streamshp, pp, subcatch, buffdis, demf, fac, out_dir)
    