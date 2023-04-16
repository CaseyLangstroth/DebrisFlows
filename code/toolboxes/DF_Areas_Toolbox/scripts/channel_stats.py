# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 11:22:29 2023

@author: clang
"""

import arcpy
from arcpy import env
import os
import numpy as np
import math
from arcpy.sa import *
from scipy import spatial


arcpy.env.overwriteOutput = True #allow overwrite

def print2(instr): 
    arcpy.AddMessage(instr)
  
        
def chk_mk_dir(in_dir):
    try:
        os.mkdir(in_dir)
    except:
        pass
 # function ot find closest points
def do_kdtree(combined_x_y_arrays,points):
    mytree = spatial.cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)
    return indexes, dist    

# =============================================================================
# Begin Tool    
# =============================================================================
def ChannelWidth (ChannelGroupLayer, IDfield, cl, demf, tranden, dx, dxr, out_dir, tname):
    
    path = os.path.join(out_dir, 'temp\\')
    chk_mk_dir(path)
    fpath = os.path.join(out_dir,'outputs\\')
    chk_mk_dir(fpath)
    tpath = os.path.join(out_dir, 'tables\\')
    chk_mk_dir(tpath)
   
# =============================================================================
# Get the Channels from the group layer (for testing it is only a single layer)
# =============================================================================S
    aprx = arcpy.mp.ArcGISProject('CURRENT')
    mps =aprx.listMaps()
    layers_ALL =[]
    Channel_lyrs =[]
    aoiex = []
    for i in mps:
        layers_ALL = i.listLayers()
    for i in layers_ALL:
        if ChannelGroupLayer in i.longName:
            if i.isGroupLayer== False:
                Channel_lyrs.append(i)
        # if excludeaoi in i.longName:
        #     if i.isGroupLayer == False:
        #        aoiex.append(i)
# =============================================================================
# Make table
# =============================================================================
    table = arcpy.management.CreateTable(tpath, tname)
    arcpy.management.AddField(table, 'Site', 'TEXT')
    for i in Channel_lyrs:
        out_fid = i.name
        #.......Set up directory with base file id for outputs.......
        #out_tmp_id=path +"\\"+out_fid
        out_id=out_dir+"\\"+out_fid
        # name1 = i.name + '_1'
        # name2 = i.name + '_2'
        # arcpy.management.CalculateField(table, 'Site', name1)
        # arcpy.management.CalculateField(table, 'Site', name2)
        arcpy.management.CalculateField(table, 'Site', out_fid)
        dblfields = ['Reach_Length_m','Avg_Width_m', 'Avg_Slope']
    for i in dblfields:
        arcpy.management.AddField(table, i, 'Double')   
# =============================================================================
# project polygon to UTM for USUAL and change name for centerline tool
# =============================================================================
    for channel in Channel_lyrs:
       print2('Working: \n' + channel.name)
       # print2('creating centerline')
       # utmpoly = path + 'fc.shp'
       # arcpy.management.Project(i, utmpoly, WKID)
       # utmcl = path + 'fcl.shp'
       # arcpy.topographic.PolygonToCenterline(utmpoly, utmcl)
       # #project centerline back to WGS
       # SR = arcpy.Describe(cl).spatialReference
       # if not SR == '26912':
       polyUTM = path  + 'UTM.shp'
       arcpy.management.Project(channel, polyUTM, '26912')
       UTMcl = path + channel.name + '_UTMcl.shp'
       arcpy.management.Project(cl, UTMcl, '26912')
       # if len(excludeaoi) > 0:
       #     for y in aoiex:
       #         arcpy.management.Project(y, path + 'aoiex', WKID)
# =============================================================================
# Run Scott's Fluvial Polygon Transect and Width tool        
# =============================================================================
# =============================================================================
# transect tool
# =============================================================================
       print2('Beginning fluvial polygon transect tool')
       # preprocess data to get edge points split by polygon side
       pfc=[polyUTM,UTMcl]
       outpfc=path+'splitpoly.shp'
       # seperate left and right polygon
       arcpy.management.FeatureToPolygon(pfc, outpfc, None, "NO_ATTRIBUTES", None)
    
      #apply unique line identifier to clean up data later
       arcpy.AddField_management(outpfc, "LID", "LONG")
     # write the origional ids to origid field
       arcpy.CalculateField_management(outpfc, "LID", "!FID!","PYTHON","")
    
       nodeletefield=["FID","LID","Shape"] # fields not to delete
      # clean up all the unnecessary attributes
       for f in arcpy.ListFields(outpfc):
           if f.name not in nodeletefield:
               try:
                    arcpy.DeleteField_management(outpfc,f.name)
               except:
                    print(f.name)
     #convert to polygon edges to lines
       polyperm=path+'LRpolyedge.shp'
       arcpy.management.FeatureToLine(outpfc, polyperm, None, "ATTRIBUTES")
    
       #delete the centerline from the polygon to line
       plineedge=path+'edge.shp'
       arcpy.analysis.Erase(polyperm, UTMcl, plineedge, None)
    
       # dissolve by UID to ensure not breaks along an edge
       edgediss=path+'edgediss.shp'
       arcpy.management.Dissolve(plineedge,edgediss, "LID", None,
                                 "MULTI_PART", "DISSOLVE_LINES")
    
    
       # genererate points along each edge
       den_dist=str(tranden)#mapunits can forace a unit via +" Meters"
       edgepoints=path+'points.shp'       
       arcpy.management.GeneratePointsAlongLines(edgediss, edgepoints, "DISTANCE", 
                                                 den_dist, None, None)
       #add a point unique id 
       #apply unique identifier to clean up data later
       arcpy.AddField_management(edgepoints, "UID", "LONG")
      # write the origional ids to origid field
       arcpy.CalculateField_management(edgepoints, "UID", "!FID!","PYTHON","")
    
    
       #%%------------------ Optional filtering of points-----------------------------
       # remove any points that fall under input polygon(s)
       # if len(excludeaoi) > 0:  
       #     for i in aoiex:
       #         pclip=path+'points_aoiclip.shp'
       #         arcpy.analysis.Erase(edgepoints, exludeaoi, pclip, None)
       #         edgepoints=pclip
           
      
           
       #%%----------------------------Match up points---------------------------------
       print2("Finding point pairs")       
       # fields to get from shape file
       fields=["SHAPE@X","SHAPE@Y","UID","LID"]
       # pull data from shape files and write them to arrays
       arr = arcpy.da.FeatureClassToNumPyArray(edgepoints, fields, skip_nulls=True)
       xp=arr["SHAPE@X"]
       yp=arr["SHAPE@Y"]
       uid=arr["UID"]
       lid=arr["LID"]
    
       # get unique line identifiers to loop over
       stepper=np.unique(lid)
    
       #xy = np.dstack([xp.ravel(),yp.ravel()])[0]
       # build initial stack of fid xy data [2D array]
       xyorig = np.dstack([xp,yp,lid,uid])[0]
       xy=xyorig.copy()# make a copy to edit
       c=0
       out=[0,0,0]      
       for i in stepper:
           # xy coorinates
           #idx_q=np.logical_and(xy[:,2]==i,xy[:,4]==-9999)
           #idx_i=np.logical_and(xy[:,2]!=i,xy[:,4]==-9999)
           xyq=xy[xy[:,2]==i,0:2]#points along the line
           xyi=xy[xy[:,2]!=i,0:2]#other points
          #xyq=xy[idx_q,0:2]#points along the line
           #xyi=xy[idx_i,0:2]#other point
           xyq_uid=xy[xy[:,2]==i,3]
           xyi_uid=xy[xy[:,2]!=i,3]
           #xyq_uid=xy[idx_q,3]
           #    xyi_uid=xy[idx_i,3]
       
           [idx,dist] = do_kdtree(xyi,xyq)
           cxy=xyi[idx]
           
           # this is gross find elegant sol'n in future
           for j in np.arange(len(cxy)):
               cxytemp=np.hstack((cxy[j,:],c))
               xytemp=np.hstack((xyq[j,:],c))
               out=np.vstack((out,cxytemp,xytemp))
               c+=1
       out=np.delete(out,0,0)
    
    
    
       oid=np.arange(np.size(out,0))
       X=out[:,0]
       Y=out[:,1]
       cid=out[:,2]
       SR = arcpy.Describe(polyUTM).spatialReference
    
    
       array = np.array([(oid[0], (X[0], Y[0]),cid[0])],
                          np.dtype([('idfield',np.int32),('XY', '<f8', 2),('cid',np.int32)]))
    
       for i in range(1,len(X)):
           array = np.concatenate((array,np.array([(oid[i], (X[i], Y[i]),cid[i])],
                                  np.dtype([('idfield',np.int32),('XY', '<f8', 2),('cid',np.int32)]))))
    
    
       #%%-----------------------------Generate final outputs-------------------------
       print2('Building lines and final outputs')
       outfc=out_id+"transect_nodes.shp"
       arcpy.da.NumPyArrayToFeatureClass(array, outfc, ['XY'], SR)

       arcpy.env.outputCoordinateSystem = SR

       tranfc=out_id+"transects.shp"
       arcpy.management.PointsToLine(outfc,tranfc , "cid", None, "NO_CLOSE")

       arcpy.CalculateGeometryAttributes_management(tranfc, [["Length_m", "LENGTH"]],"METERS")
       
    # =============================================================================
    # width tool       
    # =============================================================================
       print2('Beginning width tool')
       # check for a Length_m field created in previous code.... if not create it
       if not arcpy.ListFields(tranfc, "Length_m"):
           arcpy.CalculateGeometryAttributes_management(tranfc,
                                                    [["Length_m", "LENGTH"]],"METERS")
    
       # #densify transects for point generation
       # trandensity=path+"densetransect.shp"
       # #make duplicate line to not modify origional data
       # arcpy.conversion.FeatureClassToFeatureClass(tranfc, path, trandensity)
       # trandensity=path+"\\"+trandensity#add path to file name
    
    
      # densify line to have sufficent number of point
       arcpy.edit.Densify(tranfc, "DISTANCE", str(dx))#+" Meters")
    
       #convert line to points
       tranpoints=path+"transectpoints.shp"
       arcpy.management.FeatureVerticesToPoints(tranfc, tranpoints, "ALL")
    
       # get the data projection
       SR = arcpy.Describe(tranpoints).spatialReference
    
       #generate tin where "elevation" is the running width
       indata=tranpoints+" Length_m Mass_Points <None>;"+polyUTM+" <None> Hard_Clip <None>"
       outTIN=path+"widthTIN"
       arcpy.ddd.CreateTin(outTIN, SR, indata, "CONSTRAINED_DELAUNAY")
    
       #convert the tin to a raster       
       widthRaster=path+"_widthraster.tif"
       arcpy.ddd.TinRaster(outTIN, widthRaster, "FLOAT", "LINEAR", "CELLSIZE", 1, dxr)
    # =============================================================================
    # Back to my code       
    # =============================================================================
       print2('Getting attributes')
       #reproject and save raster
       width = fpath + channel.name + '_widthraster.tif'
       arcpy.management.ProjectRaster(widthRaster, width)
       
       # get zonal statistics as table for the entire polygon, not link
       wtbl = fpath + channel.name + '_avgwidth.dbf'
       arcpy.sa.ZonalStatisticsAsTable(channel, IDfield, width, wtbl)
                                       
      # get width along points
      # get slope for latest channel polygon
       
       end = len(channel.name)
       yr = int(channel.name[end-3:end])
       if yr > 2018 or yr == 'post':
           outslope = arcpy.sa.Slope(demf, 'DEGREE')
           pts = fpath + channel.name + '_pts.shp'
           arcpy.management.GeneratePointsAlongLines(UTMcl, pts, 'DISTANCE', '10 Meters','', 'NO_END_POINTS', 'ADD_CHAINAGE')
           arcpy.sa.ExtractMultiValuesToPoints(pts,[width, outslope])
           stbl = fpath + channel.name + '_slope.dbf'
           arcpy.sa.ZonalStatisticsAsTable(channel, IDfield, outslope, stbl)
       else:
          pts = fpath + channel.name + '_pts.shp'
          arcpy.management.GeneratePointsAlongLines(UTMcl, pts, 'DISTANCE', '10 Meters','','NO_END_POINTS', 'ADD_CHAINAGE')
          arcpy.sa.ExtractMultiValuesToPoints(pts,width)
       
       #get straight length distance      
       # arcpy.management.AddGeometryAttributes(UTMcl, 'LENGTH_GEODESIC', 'METERS')
       # sinuosity = []
       # outtbl = path + channel.name + 'strlength'
       # with arcpy.da.SearchCursor(UTMcl, [IDfield, 'LENGTH_GEO']) as cursor:
       #     strlen = []
       #     for row in cursor:
       #         stpt = arcpy.management.FeatureVerticesToPoints(UTMcl, path + channel.name + 'st.shp', 'START')
       #         endpt = arcpy.management.FeatureVerticesToPoints(UTMcl, path + channel.name + 'end.shp', 'END')
       #         arcpy.analysis.PointDistance(stpt, endpt, outtbl)
       #         with arcpy.da.SearchCursor(outtble, 'DISTANCE') as cursor:                  
       #             dist = []
       #             for x in searcher:
       #                 dist.append(x[0])
       #                 distset = set(dist)
       #                 strlen.append(distset)
       #         for i in strlen:
       #             sin = row[1]/i
       #             sinuosity.append(row[0], sin)
    # =============================================================================
    # Gather data
    # =============================================================================
       # print2(sinuosity)
      #reach length
       reachlen = []
       arcpy.management.AddGeometryAttributes(UTMcl, 'LENGTH_GEODESIC', 'METERS')
       with arcpy.da.SearchCursor(UTMcl, 'LENGTH_GEO') as cursor:
          for row in cursor:
              reachlen.append(row)
       print2(reachlen)
      # average width
       avgwidth = []
       with arcpy.da.SearchCursor(wtbl, 'MEAN') as cursor:
           for row in cursor:
               avgwidth.append(row)    
         #average slope
       avgslope = []
       with arcpy.da.SearchCursor(stbl, 'MEAN') as cursor:
            for row in cursor:                
                avgslope.append(row)
    # =============================================================================
    # Populate Table
    # =============================================================================
       rows = []
       with arcpy.da.SearchCursor(table, 'Site') as cursor:
           for row in cursor:
               if row[0] == name1:
                   for a in reachlen:
                       if a[0] == channel.name:
                           rows.append(a[1])
                   for b in avgwidth:
                       if b[0] == channel.name:
                           rows.append(b[1])
                   for c in avgslope:
                       if c[0] == channel.name:
                           rows.append(c[1])
                   # for d in sinuosity:
                   #     if d[0] == channel.name:
                   #         rows.append(d[1])
       print2(rows)
        
       with arcpy.da.InsertCursor(table, dblfields) as cursor:
           cursor.insertRows(rows)

    return
# This is used to execute code if the file was run but not imported
if __name__ == '__main__':
    
    # Tool parameter accessed with GetParameter or GetParameterAsText
    ChannelGroupLayer = arcpy.GetParameterAsText(0)
    IDfield = arcpy.GetParameterAsText(1)
    cl = arcpy.GetParameterAsText(2)
    demf = arcpy.GetParameterAsText(3)
    tranden = arcpy.GetParameter(4)
    dx= arcpy.GetParameter(5)
    dxr = arcpy.GetParameter(6)
    out_dir = arcpy.GetParameterAsText(7) 
    tname = arcpy.GetParameterAsText(8)
    
    ChannelWidth(ChannelGroupLayer, IDfield, cl, demf, tranden, dx, dxr, out_dir, tname)                      
                             
                     
           
     