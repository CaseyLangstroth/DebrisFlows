# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:26:13 2022

@author: BelmontLab
"""
import arcpy
import pandas as pd
import numpy as np
import os
import math
from arcpy.sa import *
arcpy.env.overwriteOutput=True

inData=r"D:\Alec\WEST\inData.xlsx"

#Out Folder
out_folder = r'D:\Alec\WEST\Output'
#toolbox to load in
tbox=r"D:\Alec\WEST\Utah-State-University-Applied-Watershed-Tools-main\Utah-State-University-Applied-Watershed-Tools.tbx"
tbox2=r"D:\Kevin\GIS Tools\Fluvial_corridor\ArcGIS-Pro-2_1_1-2_5_1\FluvialCorridor10_1.tbx"

#load in toolbox
arcpy.ImportToolbox(tbox)
arcpy.ImportToolbox(tbox2)

def chk_mk_dir(input_dir):
    try:
        os.mkdir(input_dir)
    except:
        pass
 
# Load in excel input data
df = pd.read_excel(inData)
c = df.columns
for col in c:
   globals()[col]=np.array(df[col])

indem=r'D:\GIS\Utah_DEM_utm_fixed.tif'
pp=r'D:\Alec\WEST\pp.shp'
aoi=r'D:\Alec\WEST\wtshd_buff.shp'
# try:
#     arcpy.CreateFileGDB_management(out_folder, "centerline.gdb")
# except:
#     pass

for i in np.arange(len(ppfid)):
    print(outid[i])
    outpath=os.path.join(out_folder,outid[i])
    chk_mk_dir(outpath)
    
    #%% Run the USUAL Tools
    # arcpy.USUAL.WatershedRiverDelineation(indem, 5, pp, "FID", str(ppfid[i]), 50, aoi, "FID", str(wtshdfid[i]), outpath, outid[i], False, '', False)
    # print ("Watershed Delineation Done")
    # indemf=fr"{outpath}\{outid[i]}demf.tif"
    # fdr=fr"{outpath}\{outid[i]}fdr.tif"
    # fac=fr"{outpath}\{outid[i]}fac.tif"
    # stream_5k_ras=fr"{outpath}\{outid[i]}_stream_5km.tif"
    # stream_5k_shp=fr"{outpath}\{outid[i]}_stream_5km.shp"
    # wtshd=fr"{outpath}\{outid[i]}_watershed.shp"
    # arcpy.USUAL.SubCatchmentInterfluveDelineation(indemf, fdr, fac, stream_5k_shp, wtshd, 0.1, outpath, outid[i],'',True,True,False,False)
    # print ("SubCatchment Delineation Done")
    # arcpy.USUAL.NetworkDiscretizationAttribution(stream_5k_shp, fac, indemf, 500, 0.001, outpath, outid[i], False)
    # print ("Network Discretization Done")
    # sn=fr"{outpath}\{outid[i]}_network.shp"
    # vb= fr'D:\Alec\Input\VBs\{outid[i]}.shp' #Valley Bottom
    # arcpy.FluvialCorridor.ValleyBottom(False, sn, indemf, 100, 100, 10, 100, -4, 4, None, 100, 20, 1000, 100, vb, False)   
    # print ("Valley Bottom Delineation Done")
    # centerline=fr'D:\Alec\WEST\Output\centerline.gdb\{outid[i]}'  
    # arcpy.topographic.PolygonToCenterline(vb, centerline, None)  
    # arcpy.USUAL.FluvialPolygonTransects(vb, centerline, 10, outpath, outid[i], None, 30, False)
    # transects=fr"{outpath}\{outid[i]}transects.shp"
    # arcpy.USUAL.FluvialLinkAverageWidth(sn, "VB_width", vb, transects, 10, outpath, outid[i], 10, False)
    # print ("Width Calculations Done")
    
    #%% Metrics ####
   
    points=fr'D:\Alec\WEST\Output\{outid[i]}\temp\RiverDisc\{outid[i]}_nodes.shp'
    srb_points=r'D:\Alec\WEST\Input\SRB_model_points_2022.shp'
   
    rdnbr=r'D:\Alec\Input\mtbs8_18_reclassify.tif'
    perimeters='D:\Alec\Input\Final_fire_perimeters.shp'
    modhigh=r'D:\Alec\Input\mod_high.tif'
    cover=r'D:\Alec\Input\US_130EVT_reclass.tif'
    conifer=r'D:\Alec\Input\conifer.tif'
    lw_file=fr'D:\Alec\Input\Wood\{wood[i]}.shp'  ##need to set up
    debrisflow=fr'D:\Al
    
    ec\Input\DF_volume_all2.shp'
    tpi=r'D:\Alec\Input\tpi_Ut.tif'
    
    #input variables ####
    theta=concavity[i]
    a=0.72 #channel width power function variables
    b=0.514
    
    out_temp_dir= os.path.join(outpath,'temp','metrics') # make directory name
    chk_mk_dir(out_temp_dir)
    print(out_temp_dir)
   
    #.......Set up directory with base file id for outputs.......
    out_fid=outid[i]
    out_tmp_id=out_temp_dir+"\\"+out_fid
    out_id=outpath+"\\"+out_fid
    
    disc_poly=outpath+"\\"+outid[i]+"disc_poly.shp"
    wtshds=outpath+"\\"+outid[i]+"wtshds.shp"
    sn=fr"{outpath}\{outid[i]}_network.shp"
    
    # %% Geomorphic Metrics
      ##Normalized Steepness
      #Source- The resilience of logjams to floods- Wohl 2021
      #Calculation- a=Floodplain width b=
      #Importance
    # arcpy.DeleteField_management(sn, ["NormSteep"])
    # arcpy.management.AddField(sn, "NormSteep", "FLOAT")
    # with arcpy.da.UpdateCursor(in_table=sn, field_names=['Slope','usarea_km2','NormSteep']) as updater:
    #     for row in updater:
    #         row[2]=row[0]/math.pow(abs(row[1])*1000000,-theta)
    #         updater.updateRow(row)
                    
    #   ##Stream Power
    #   #Source- The resilience of logjams to floods- Wohl 2021
    #   #Calculation- a=Floodplain width b=
    #   #Importance
    # arcpy.DeleteField_management(sn, ["StreamPow"])
    # arcpy.management.AddField(sn, "StreamPow", "FLOAT")
    # with arcpy.da.UpdateCursor(in_table=sn, field_names=['Slope','usarea_km2','StreamPow']) as updater:
    #     for row in updater:
    #         row[2]=row[0]*row[1]
    #         updater.updateRow(row)      
    
    
    #   ##Finite Differencing
    #   #Source- The resilience of logjams to floods- Wohl 2021
    #   #Calculation- a=Floodplain width b=
    #   #Importance
    # arcpy.DeleteField_management(sn, ["S_Change","VB_Change"])
    # arcpy.management.AddField(sn, "S_Change", "FLOAT")
    # arcpy.management.AddField(sn, "VB_Change", "FLOAT")
    # updater=arcpy.da.UpdateCursor(sn, field_names=['Slope', 'VB_width','S_Change', 'VB_Change','GRIDID'])
    # for row in updater:
    #     drainsfrom=row[4]
    #     us_slope=[]
    #     us_width=[]
    #     searcher= arcpy.da.SearchCursor(sn, field_names=['ToLink','Slope', 'VB_width'],where_clause=f'ToLink={drainsfrom}')
    #     for bow in searcher:
    #         us_slope.append(bow[1])
    #         us_width.append(bow[2])
    #     if len(us_slope)==0:
    #         row[2]=0
    #     else:
    #         avg_slope=min(us_slope)
    #         ds_slope=row[0]
    #         row[2]=(ds_slope-avg_slope)/ds_slope
    #     if len(us_width)==0:
    #         row[3]=0
    #     else:
    #         avg_width=max(us_width)
    #         ds_width=row[1]
    #         if ds_width==0:
    #             row[3]=0
    #         else:
    #             row[3]=(ds_width-avg_width)/ds_width
    #     updater.updateRow(row)
    
    #   ##Floodplain width vs Channel Width
    # #Source- The resilience of logjams to floods- Wohl 2021
    # #Calculation- a=Floodplain width b=
    # #Importance  
    # arcpy.DeleteField_management(sn, ["FpChanRat","ChanWidth"])
    # arcpy.management.AddField(sn, "FpChanRat", "FLOAT")
    # arcpy.management.AddField(sn, "ChanWidth", "FLOAT")
    # with arcpy.da.UpdateCursor(in_table=sn, field_names=['usarea_km2','VB_width','FpChanRat','ChanWidth']) as updater:
    #     for row in updater:
    #         row[2]=row[1]/(a*math.pow(abs(row[0]),b))
    #         row[3]=a*math.pow(abs(row[0]),b)
    #         updater.updateRow(row)
    
    #   ##Confinement
    #   #Source- The resilience of logjams to floods- Wohl 2021
    #   #Calculation- a=Floodplain width b=
    #   #Importance
    # arcpy.DeleteField_management(sn, ["Confine"])
    # arcpy.management.AddField(sn, "Confine", "FLOAT")
    # arcpy.management.AddField(sn, "chan_buf", "FLOAT")
    # arcpy.CalculateField_management(sn, "chan_buf", "!ChanWidth!+15", "PYTHON_9.3")
    # sn_buffer=out_temp_dir+"\sn_buffer_"+outid[i]+".shp"
    # arcpy.analysis.Buffer(sn, sn_buffer, "chan_buf", "FULL", "FLAT", "NONE", None, "PLANAR")
    # boundary_file=out_tmp_id+"boundary_file.shp"
    # arcpy.PolygonToLine_management (vb, boundary_file) # turn the valley bottom polygon into a polyline boundary
    # sn_intersect=out_temp_dir+"\sn_intersect_"+outid[i]+".shp"
    # arcpy.analysis.Intersect(f"{boundary_file};{sn_buffer}", sn_intersect, "ONLY_FID", None, "INPUT")
    # arcpy.management.AddGeometryAttributes(sn_intersect, "LENGTH_GEODESIC")
    # confine_stats=out_tmp_id+"confine_stats.dbf"
    # arcpy.analysis.Statistics(sn_intersect, confine_stats, [["LENGTH_GEO","SUM"]], "FID_sn_buf")
    # arcpy.management.JoinField(sn, "FID", confine_stats, "FID_sn_buf", ['SUM_LENGTH'])
    # arcpy.CalculateField_management(sn, "Confine", "!SUM_LENGTH!/!Length_m!","PYTHON","")
    # arcpy.DeleteField_management(sn, "SUM_LENGTH;chan_buf")

    #   ##Sinuosity
    #   #Source- The resilience of logjams to floods- Wohl 2021
    #   #Calculation- a=Floodplain width b=
    #   #Importance
    # arcpy.DeleteField_management(sn, ["Sinuosity"])
    # arcpy.management.AddField(sn, "Sinuosity", "FLOAT")
    # sin_line=out_tmp_id+"sin_line.shp"
    # sin_points=out_tmp_id+"sin_points.shp"
    # updater=arcpy.da.UpdateCursor(sn, field_names=["Length_m",'shape@','Sinuosity'])
    # for row in updater:
    #     sin_line=row[1]
    #     arcpy.management.FeatureVerticesToPoints(sin_line,sin_points,  "BOTH_ENDS")
    #     searcher= arcpy.da.SearchCursor(in_table=sin_points, field_names=['shape@'])
    #     for bow in searcher:
    #         item=bow[0]
    #         item2=next(searcher)
    #         item3=item2[0]
    #         new_distance=item.distanceTo(item3)
    #     row[2]=row[0]/new_distance
    #     updater.updateRow(row)
    
     # %% Local Landcover Variables 
    
    # ## Develop Discretized Valley Bottom
    # arcpy.CalculateField_management(sn, "vb_buf", "!VB_width!+25", "PYTHON_9.3")
    # sn_vb_buffer=out_tmp_id+"sn_vb_buffer.shp"
    # arcpy.analysis.Buffer(sn, sn_vb_buffer, "vb_buf", "FULL", "FLAT", "NONE", None, "PLANAR")
    # disc_poly=outpath+"\\"+outid[i]+"disc_poly.shp"
    # arcpy.analysis.Clip(sn_vb_buffer, vb, disc_poly, None)
    # arcpy.DeleteField_management(sn, "vb_buf")
     
    # ##Local Dominant Cover Type
    # #Source- The resilience of logjams to floods- Wohl 2021
    # #Calculation- a=Floodplain width b=
    # #Importance
    # arcpy.DeleteField_management(sn, ["Cover"])
    # arcpy.management.AddField(sn, "Cover", "FLOAT")
    # covertable=out_tmp_id+"covertable.dbf"
    # arcpy.sa.ZonalStatisticsAsTable(disc_poly, "FID", cover, covertable, "DATA", "MAJORITY", "CURRENT_SLICE", 90, "AUTO_DETECT")
    # arcpy.management.JoinField(disc_poly, "FID", covertable, "FID_", "MAJORITY")
    # arcpy.management.JoinField(sn, "GridID", disc_poly, "GridID","MAJORITY")
    # arcpy.DeleteField_management(disc_poly, ["MAJORITY"])
    # arcpy.CalculateField_management(sn, "Cover", "!MAJORITY!","PYTHON","")
    # arcpy.DeleteField_management(sn, ["MAJORITY"])
    
    ##Topographic Position Analysis
    #Source- The resilience of logjams to floods- Wohl 2021
    #Calculation- a=Floodplain width b=
    #Importance
    arcpy.DeleteField_management(sn, ["TPI"])
    arcpy.management.AddField(sn, "TPI", "FLOAT")
    tpitable=out_tmp_id+"tpitable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(disc_poly, "FID", tpi, tpitable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.DeleteField_management(disc_poly, ["MEAN"])
    arcpy.management.JoinField(disc_poly, "FID", tpitable, "FID_", "MEAN")
    arcpy.management.JoinField(sn, "GridID", disc_poly, "GridID","MEAN")
    arcpy.CalculateField_management(sn, "TPI", "!MEAN!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MEAN"])
    
    ##Local Burn Severity
    #Source- The resilience of logjams to floods- Wohl 2021
    #Calculation- a=Floodplain width b=
    #Importance
    arcpy.DeleteField_management(sn, ["Burn"])
    arcpy.management.AddField(sn, "Burn", "FLOAT")
    burnsevtable=out_tmp_id+"burnsevtable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(disc_poly, "FID", rdnbr, burnsevtable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.DeleteField_management(disc_poly, ["MEAN"])
    arcpy.management.JoinField(disc_poly, "FID", burnsevtable, "FID_", "MEAN")
    #arcpy.SpatialJoin_analysis(rivnetmidp, disc_poly, joindata)
    arcpy.management.JoinField(sn, "GridID", disc_poly, "GridID","MEAN")
    arcpy.CalculateField_management(sn, "Burn", "!MEAN!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MEAN"])

    
     # %% Upstream Landcover Variables
    # wtshds=outpath+"\\"+outid[i]+"wtshds.shp"
    #     ##Upstream 
    # srs=arcpy.SpatialReference(26912)
    # arcpy.CreateFeatureclass_management(out_path=outpath, out_name=f"{outid[i]}wtshds.shp", geometry_type='POLYGON',spatial_reference=srs)
    # arcpy.management.AddField(wtshds, "GridID", "FLOAT")
    # # Collapse stream network to points
    # rivnetmidp=out_tmp_id+"networkmidpoints.shp"
    # arcpy.FeatureVerticesToPoints_management(sn,rivnetmidp,"MID")
    # searcher= arcpy.da.SearchCursor(in_table=rivnetmidp, field_names=['GridID','shape@'])
    # for row in searcher:
    #       id=row[0]
    #       pourpoint=row[1]
    #       flowdir=out_id+"Fdr.tif"
    #       in_raster=arcpy.gp.Watershed_sa(flowdir, pourpoint)
    #       wtshd=out_tmp_id+"wtshd.shp"
    #       polygon=arcpy.conversion.RasterToPolygon(in_raster, wtshd)
    #       searcher2= arcpy.da.SearchCursor(in_table=wtshd, field_names=['shape@'])
    #       for row in searcher2:
    #           polygon=row[0]
    #       icur=arcpy.da.InsertCursor(in_table=wtshds,field_names=['GridID','shape@'])
    #       icur.insertRow([id,polygon])
    # arcpy.CalculateField_management(wtshds, "Id", "!GridID!", "PYTHON_9.3")
    # print('Done with Upstream!')
    
    #   ##Cover Type    
    #   #Source- The resilience of logjams to floods- Wohl 2021
    #   #Calculation- a=Floodplain width b=
    #   #Importance
    # arcpy.DeleteField_management(sn, ["Cover_shd"])
    # arcpy.management.AddField(sn, "Cover_shd", "FLOAT")
    # covershdtable=out_tmp_id+"covershdtable.dbf"
    # arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", cover, covershdtable, "DATA", "MAJORITY", "CURRENT_SLICE", 90, "AUTO_DETECT")
    # arcpy.management.JoinField(sn, "GridID", covershdtable, "Id", "MAJORITY")
    # arcpy.CalculateField_management(sn, "Cover_Shd", "!MAJORITY!","PYTHON","")
    # arcpy.DeleteField_management(sn, ["MAJORITY"])
    
    #   ##Conifer
    #   #Source- The resilience of logjams to floods- Wohl 2021
    #   #Calculation- a=Floodplain width b=
    #   #Importance
    # arcpy.DeleteField_management(sn, ["Conifer"])
    # arcpy.management.AddField(sn, "Conifer", "FLOAT")
    # conifertable=out_tmp_id+"conifertable.dbf"
    # arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", conifer, conifertable, "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
    # arcpy.management.AddField(conifertable, "conifer_up", "FLOAT")
    # arcpy.CalculateField_management(conifertable, "conifer_up", "!SUM!*900/!AREA!*100", "PYTHON_9.3")
    # arcpy.management.JoinField(sn, "GridID",conifertable , "Id", "conifer_up")
    # arcpy.CalculateField_management(sn, "Conifer", "!conifer_up!","PYTHON","")
    # arcpy.DeleteField_management(sn, ["conifer_up"])
    
      #TPI Upstream
      #Source- The resilience of logjams to floods- Wohl 2021
      #Calculation- a=Floodplain width b=
      #Importance
    # arcpy.DeleteField_management(sn, ["tpi_up"])
    # arcpy.management.AddField(sn, "tpi_up", "FLOAT")
    # wtshds_tpi=out_tmp_id+"wtshds_tpi.shp"
    # arcpy.analysis.Intersect(f"{wtshds};{sn}", wtshds_tpi, "ALL", None, "INPUT")
    # #arcpy.SpatialJoin_analysis(wtshds, sn, wtshds_tpi, "JOIN_ONE_TO_MANY")
    # tpi_dissolve=out_tmp_id+"tpi_dissolve"
    # arcpy.management.Dissolve(wtshds_tpi,tpi_dissolve , "GridID", "TPI MEAN", "MULTI_PART", "DISSOLVE_LINES")
    # arcpy.management.JoinField(sn, "GridID",tpi_dissolve , "GridID", "MEAN_TPI")
    # arcpy.CalculateField_management(sn, "tpi_up", "!MEAN_TPI!","PYTHON","")
    # arcpy.DeleteField_management(sn, ["MEAN_TPI"])
    #arcpy.analysis.SpatialJoin("sf_watershed", "sf_network", r"C:\Users\BelmontLab\Documents\ArcGIS\Projects\MyProject7\MyProject7.gdb\sf_watershed_SpatialJoin", "JOIN_ONE_TO_MANY", "KEEP_ALL", 'Id "Id" true true false 10 Long 0 10,First,#,sf_watershed,Id,-1,-1;gridcode "gridcode" true true false 10 Long 0 10,First,#,sf_watershed,gridcode,-1,-1;Id_1 "Id" true true false 6 Long 0 6,First,#,sf_network,Id,-1,-1;ORIG_FID "ORIG_FID" true true false 10 Long 0 10,First,#,sf_network,ORIG_FID,-1,-1;ORIG_SEQ "ORIG_SEQ" true true false 10 Long 0 10,First,#,sf_network,ORIG_SEQ,-1,-1;GridID "GridID" true true false 19 Double 0 0,First,#,sf_network,GridID,-1,-1;Length_m "Length_m" true true false 19 Double 0 0,First,#,sf_network,Length_m,-1,-1;ToLink "ToLink" true true false 19 Double 0 0,First,#,sf_network,ToLink,-1,-1;usarea_km2 "usarea_km2" true true false 19 Double 0 0,First,#,sf_network,usarea_km2,-1,-1;uselev_m "uselev_m" true true false 19 Double 0 0,First,#,sf_network,uselev_m,-1,-1;dselev_m "dselev_m" true true false 19 Double 0 0,First,#,sf_network,dselev_m,-1,-1;Slope "Slope" true true false 19 Double 0 0,First,#,sf_network,Slope,-1,-1;VB_width "VB_width" true true false 19 Double 0 0,First,#,sf_network,VB_width,-1,-1;StreamPow "StreamPow" true true false 13 Float 0 0,First,#,sf_network,StreamPow,-1,-1;S_Change "S_Change" true true false 13 Float 0 0,First,#,sf_network,S_Change,-1,-1;VB_Change "VB_Change" true true false 13 Float 0 0,First,#,sf_network,VB_Change,-1,-1;FpChanRat "FpChanRat" true true false 13 Float 0 0,First,#,sf_network,FpChanRat,-1,-1;ChanWidth "ChanWidth" true true false 13 Float 0 0,First,#,sf_network,ChanWidth,-1,-1;Confine "Confine" true true false 13 Float 0 0,First,#,sf_network,Confine,-1,-1;Sinuosity "Sinuosity" true true false 13 Float 0 0,First,#,sf_network,Sinuosity,-1,-1;Cover "Cover" true true false 13 Float 0 0,First,#,sf_network,Cover,-1,-1;TPI "TPI" true true false 13 Float 0 0,First,#,sf_network,TPI,-1,-1;Burn "Burn" true true false 13 Float 0 0,First,#,sf_network,Burn,-1,-1;Cover_shd "Cover_shd" true true false 13 Float 0 0,First,#,sf_network,Cover_shd,-1,-1;Conifer "Conifer" true true false 13 Float 0 0,First,#,sf_network,Conifer,-1,-1;Mod_High "Mod_High" true true false 13 Float 0 0,First,#,sf_network,Mod_High,-1,-1;mh_sqm "mh_sqm" true true false 13 Float 0 0,First,#,sf_network,mh_sqm,-1,-1;per_burn "per_burn" true true false 13 Float 0 0,First,#,sf_network,per_burn,-1,-1;srb_up "srb_up" true true false 13 Float 0 0,First,#,sf_network,srb_up,-1,-1;srb_vol_up "srb_vol_up" true true false 13 Float 0 0,First,#,sf_network,srb_vol_up,-1,-1;df_vol "df_vol" true true false 13 Float 0 0,First,#,sf_network,df_vol,-1,-1;df_vol_up "df_vol_up" true true false 13 Float 0 0,First,#,sf_network,df_vol_up,-1,-1;df_up "df_up" true true false 254 Text 0 0,First,#,sf_network,df_up,0,254;Sed_Trans "Sed_Trans" true true false 13 Float 0 0,First,#,sf_network,Sed_Trans,-1,-1;Fire_name "Fire_name" true true false 254 Text 0 0,First,#,sf_network,Fire_name,0,254;Wtshd "Wtshd" true true false 254 Text 0 0,First,#,sf_network,Wtshd,0,254;tpi_up "tpi_up" true true false 13 Float 0 0,First,#,sf_network,tpi_up,-1,-1;NormSteep "NormSteep" true true false 13 Float 0 0,First,#,sf_network,NormSteep,-1,-1;srb_vol "srb_vol" true true false 13 Float 0 0,First,#,sf_network,srb_vol,-1,-1;srb "srb" true true false 13 Float 0 0,First,#,sf_network,srb,-1,-1;VBM_SF "VBM_SF" true true false 254 Text 0 0,First,#,sf_network,VBM_SF,0,254;Cause "Cause" true true false 254 Text 0 0,First,#,sf_network,Cause,0,254', "INTERSECT", None, '')
    #arcpy.analysis.SpatialJoin("sfwtshds", r"Networks\sf_network", r"D:\Alec\WEST\WEST.gdb\sfwtshds_SpatialJoin", "JOIN_ONE_TO_MANY", "KEEP_ALL", r'Id "Id" true true false 6 Long 0 6,First,#,sfwtshds,Id,-1,-1;GridID "GridID" true true false 13 Float 0 0,First,#,sfwtshds,GridID,-1,-1;AREA_GEO "AREA_GEO" true true false 19 Double 0 0,First,#,sfwtshds,AREA_GEO,-1,-1;perburn "perburn" true true false 13 Float 0 0,First,#,sfwtshds,perburn,-1,-1;Id_1 "Id" true true false 6 Long 0 6,First,#,Networks\sf_network,Id,-1,-1;ORIG_FID "ORIG_FID" true true false 10 Long 0 10,First,#,Networks\sf_network,ORIG_FID,-1,-1;ORIG_SEQ "ORIG_SEQ" true true false 10 Long 0 10,First,#,Networks\sf_network,ORIG_SEQ,-1,-1;GridID_1 "GridID" true true false 19 Double 0 0,First,#,Networks\sf_network,GridID,-1,-1;Length_m "Length_m" true true false 19 Double 0 0,First,#,Networks\sf_network,Length_m,-1,-1;ToLink "ToLink" true true false 19 Double 0 0,First,#,Networks\sf_network,ToLink,-1,-1;usarea_km2 "usarea_km2" true true false 19 Double 0 0,First,#,Networks\sf_network,usarea_km2,-1,-1;uselev_m "uselev_m" true true false 19 Double 0 0,First,#,Networks\sf_network,uselev_m,-1,-1;dselev_m "dselev_m" true true false 19 Double 0 0,First,#,Networks\sf_network,dselev_m,-1,-1;Slope "Slope" true true false 19 Double 0 0,First,#,Networks\sf_network,Slope,-1,-1;VB_width "VB_width" true true false 19 Double 0 0,First,#,Networks\sf_network,VB_width,-1,-1;StreamPow "StreamPow" true true false 13 Float 0 0,First,#,Networks\sf_network,StreamPow,-1,-1;S_Change "S_Change" true true false 13 Float 0 0,First,#,Networks\sf_network,S_Change,-1,-1;VB_Change "VB_Change" true true false 13 Float 0 0,First,#,Networks\sf_network,VB_Change,-1,-1;FpChanRat "FpChanRat" true true false 13 Float 0 0,First,#,Networks\sf_network,FpChanRat,-1,-1;ChanWidth "ChanWidth" true true false 13 Float 0 0,First,#,Networks\sf_network,ChanWidth,-1,-1;Confine "Confine" true true false 13 Float 0 0,First,#,Networks\sf_network,Confine,-1,-1;Sinuosity "Sinuosity" true true false 13 Float 0 0,First,#,Networks\sf_network,Sinuosity,-1,-1;Cover "Cover" true true false 13 Float 0 0,First,#,Networks\sf_network,Cover,-1,-1;TPI "TPI" true true false 13 Float 0 0,First,#,Networks\sf_network,TPI,-1,-1;Cover_shd "Cover_shd" true true false 13 Float 0 0,First,#,Networks\sf_network,Cover_shd,-1,-1;Conifer "Conifer" true true false 13 Float 0 0,First,#,Networks\sf_network,Conifer,-1,-1;Mod_High "Mod_High" true true false 13 Float 0 0,First,#,Networks\sf_network,Mod_High,-1,-1;mh_sqm "mh_sqm" true true false 13 Float 0 0,First,#,Networks\sf_network,mh_sqm,-1,-1;per_burn "per_burn" true true false 13 Float 0 0,First,#,Networks\sf_network,per_burn,-1,-1;srb_up "srb_up" true true false 13 Float 0 0,First,#,Networks\sf_network,srb_up,-1,-1;srb_vol_up "srb_vol_up" true true false 13 Float 0 0,First,#,Networks\sf_network,srb_vol_up,-1,-1;df_vol "df_vol" true true false 13 Float 0 0,First,#,Networks\sf_network,df_vol,-1,-1;df_vol_up "df_vol_up" true true false 13 Float 0 0,First,#,Networks\sf_network,df_vol_up,-1,-1;df_up "df_up" true true false 254 Text 0 0,First,#,Networks\sf_network,df_up,0,254;Sed_Trans "Sed_Trans" true true false 13 Float 0 0,First,#,Networks\sf_network,Sed_Trans,-1,-1;Fire_name "Fire_name" true true false 254 Text 0 0,First,#,Networks\sf_network,Fire_name,0,254;Wtshd "Wtshd" true true false 254 Text 0 0,First,#,Networks\sf_network,Wtshd,0,254;srb_vol "srb_vol" true true false 13 Float 0 0,First,#,Networks\sf_network,srb_vol,-1,-1;srb "srb" true true false 13 Float 0 0,First,#,Networks\sf_network,srb,-1,-1;VBM_SF "VBM_SF" true true false 254 Text 0 0,First,#,Networks\sf_network,VBM_SF,0,254;Cause "Cause" true true false 254 Text 0 0,First,#,Networks\sf_network,Cause,0,254;NormSteep "NormSteep" true true false 13 Float 0 0,First,#,Networks\sf_network,NormSteep,-1,-1;Burn "Burn" true true false 13 Float 0 0,First,#,Networks\sf_network,Burn,-1,-1;tpi_up "tpi_up" true true false 13 Float 0 0,First,#,Networks\sf_network,tpi_up,-1,-1', "INTERSECT", None, '')
    #     ##Mod High
    #     #Source- The resilience of logjams to floods- Wohl 2021
    #     #Calculation- a=Floodplain width b=
    #     #Importance
    # arcpy.DeleteField_management(sn, ["Mod_High"])
    # arcpy.management.AddField(sn, "Mod_High", "FLOAT")
    # modhightable=out_tmp_id+"modhightable.dbf"
    # arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", modhigh, modhightable, "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
    # arcpy.management.AddField(modhightable, "modhigh", "FLOAT")
    # arcpy.CalculateField_management(modhightable, "modhigh", "!SUM!/!AREA!*100", "PYTHON_9.3")
    # arcpy.management.JoinField(sn, "GridID", modhightable, "Id", "modhigh")
    # arcpy.CalculateField_management(sn, "Mod_High", "!modhigh!","PYTHON","")
    # arcpy.DeleteField_management(sn, ["modhigh"])
    
    #     ##Mod High Square Meters
    #     #Source- The resilience of logjams to floods- Wohl 2021
    #     #Calculation- a=Floodplain width b=
    #     #Importance
    # arcpy.DeleteField_management(sn, ["mh_sqm"])
    # arcpy.management.AddField(sn, "mh_sqm", "FLOAT")
    # arcpy.CalculateField_management(modhightable, "mhsqm", "!SUM!", "PYTHON_9.3")
    # arcpy.management.JoinField(sn, "GridID", modhightable, "Id", "mhsqm")
    # arcpy.CalculateField_management(sn, "mh_sqm", "!mhsqm!","PYTHON","")
    # arcpy.DeleteField_management(sn, ["mhsqm"])
   
    # #     ##PRISM Data
    # # arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", prism, prism_shd, "DATA", "MAXIMUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
    # # arcpy.management.AddField(prism_shd, "prism_shd", "FLOAT")
    # # arcpy.CalculateField_management(prism_shd, "prism_shd", "!MAX!", "PYTHON_9.3")
    # # arcpy.management.JoinField(sn, "GridID", prism_shd, "Id", "prism_shd")
    
    #     ##Total Burn Percentage
    # arcpy.DeleteField_management(sn, ["per_burn"])
    # arcpy.management.AddField(sn, "per_burn", "FLOAT")
    # arcpy.management.AddGeometryAttributes(wtshds, "AREA_GEODESIC")
    # arcpy.management.AddField(wtshds, "perburn", "FLOAT")
    # perintersect=out_tmp_id+"perintersect.shp"
    # water_feature_name=out_fid+"water_feature"
    # water_feature=out_tmp_id+"water_feature"
    # updater=arcpy.da.UpdateCursor(wtshds, field_names=['AREA_GEO', 'perburn','SHAPE@'])
    # for row in updater:
    #     waterarea=row[0]
    #     water=row[2]
    #     arcpy.conversion.FeatureClassToFeatureClass(water,out_temp_dir, water_feature_name)
    #     arcpy.analysis.Intersect(f"{water_feature};{perimeters}", perintersect)
    #     arcpy.management.AddGeometryAttributes(perintersect, "AREA_GEODESIC")
    #     searcher= arcpy.da.SearchCursor(perintersect, field_names=['AREA_GEO'])
    #     for bow in searcher:
    #         burnarea=bow[0]
    #     row[1]=burnarea/waterarea*100
    #     #print(row[1])
    #     updater.updateRow(row)
    # arcpy.management.JoinField(sn, "GridID", wtshds, "GridID", ['perburn'])
    # arcpy.CalculateField_management(sn, "per_burn", "!perburn!","PYTHON","")
    # arcpy.DeleteField_management(sn, "perburn")
    
    #%% Bottlenecks, Debris Flows and Wood
    
      ##Bottlenecks Present 
    # arcpy.DeleteField_management(sn, ["srb_vol","srb","VBM_SF","Cause","srb_vbm_sf","srb_cause"])
    # arcpy.management.AddField(sn, "srb_vol", "FLOAT")
    # arcpy.management.AddField(sn, "srb", "FLOAT")
    # arcpy.management.AddField(sn, "VBM_SF", "TEXT")
    # arcpy.management.AddField(sn, "Cause", "TEXT")
    # disc_poly2=out_tmp_id+"disc_poly2"
    # arcpy.SpatialJoin_analysis(disc_poly, srb_points, disc_poly2, "JOIN_ONE_TO_MANY")
    # disc_poly2_dissolve=out_tmp_id+"disc_poly2_dissolve"
    # arcpy.management.Dissolve(disc_poly2,disc_poly2_dissolve , "GridID", "Join_Count SUM;Volume SUM;VBM_SF FIRST;Cause FIRST", "MULTI_PART", "DISSOLVE_LINES")
    # arcpy.management.JoinField(sn, "GridID", disc_poly2_dissolve, "GridID",['SUM_Join_C','FIRST_VBM_','FIRST_Caus','SUM_Volume'])
    # arcpy.CalculateField_management(sn, "srb_vol", "!SUM_Volume!","PYTHON","")
    # arcpy.CalculateField_management(sn, "srb", "!SUM_Join_C!","PYTHON","")
    # arcpy.CalculateField_management(sn, "VBM_SF", "!FIRST_VBM_!","PYTHON","")
    # arcpy.CalculateField_management(sn, "Cause", "!FIRST_Caus!","PYTHON","")
    # arcpy.DeleteField_management(sn,'SUM_Join_C;SUM_Volume;FIRST_VBM_;FIRST_Caus')
    
    #   ##Bottleneck Upstream
    # arcpy.DeleteField_management(sn, ["srb_up","srb_vol_up"])
    # arcpy.management.AddField(sn, "srb_up", "FLOAT")
    # arcpy.management.AddField(sn, "srb_vol_up", "FLOAT")
    # srb_up=out_tmp_id+"srb_up"
    # arcpy.SpatialJoin_analysis(wtshds, srb_points, srb_up, "JOIN_ONE_TO_MANY")
    # srb_up_dissolve=out_tmp_id+"srb_up_dissolve"
    # arcpy.management.Dissolve(srb_up,srb_up_dissolve , "GriDID", "Join_Count SUM;Volume SUM", "MULTI_PART", "DISSOLVE_LINES")
    # arcpy.management.JoinField(sn, "GridID", srb_up_dissolve, "GridID", ['SUM_Join_C','SUM_Volume'])
    # arcpy.CalculateField_management(sn, "srb_up", "!SUM_Join_C!","PYTHON","")
    # arcpy.CalculateField_management(sn, "srb_vol_up", "!SUM_Volume!","PYTHON","")
    # arcpy.DeleteField_management(sn,'SUM_Join_C;SUM_Volume')
    
    #   ##Debris Flow Delivery
    # arcpy.DeleteField_management(sn, ["df_vol"])
    # arcpy.management.AddField(sn, "df_vol", "FLOAT")
    # sn_buffer2=out_tmp_id+"sn_buffer2.shp"
    # arcpy.analysis.Buffer(sn, sn_buffer2, "Chanwidth", "FULL", "FLAT", "NONE", None, "PLANAR")
    # sn_intersect2=out_tmp_id+"sn_intersect2.shp"
    # arcpy.analysis.PairwiseIntersect(f"{debrisflow};{sn_buffer2}", sn_intersect2, "ALL", None, "INPUT")
    # arcpy.management.AddField(sn_intersect2, "Vol_input", "FLOAT")
    # arcpy.management.AddGeometryAttributes(sn_intersect2, "AREA_GEODESIC")
    # arcpy.CalculateField_management(sn_intersect2, "Vol_input", "!AREA_GEO!/!m2!*!m3!", "PYTHON_9.3")
    # arcpy.management.JoinField(sn, "FID", sn_intersect2, "ORIG_FID", ['Vol_input'])
    # arcpy.CalculateField_management(sn, "df_vol", "!Vol_input!","PYTHON","")
    # arcpy.DeleteField_management(sn,'Vol_input')
    
    #   ##Debris Flow Delivery Upstream
    # arcpy.DeleteField_management(sn, ["df_up","df_vol_up"])
    # arcpy.management.AddField(sn, "df_vol_up", "FLOAT")
    # df_up=out_tmp_id+"df_up.shp"
    # arcpy.SpatialJoin_analysis(wtshds, sn_intersect2 , df_up, "JOIN_ONE_TO_MANY")
    # df_up_dissolve=out_tmp_id+"df_up_dissolve.shp"
    # arcpy.management.Dissolve(df_up,df_up_dissolve , "GridID", "Join_Count SUM;Vol_input SUM", "MULTI_PART", "DISSOLVE_LINES")
    # arcpy.management.JoinField(sn, "GridID", df_up_dissolve, "GridID", ['SUM_Join_C','SUM_Vol_in'])
    # arcpy.CalculateField_management(sn, "df_up", "!SUM_Join_C!","PYTHON","")
    # arcpy.CalculateField_management(sn, "df_vol_up", "!SUM_Vol_in!","PYTHON","")
    # arcpy.DeleteField_management(sn,'SUM_Join_C;SUM_Vol_in')
    
      ## Wood Volume
    #with arcpy.da.UpdateCursor(in_table=disc_poly, field_names=['FID','ORIG_FID']) as updater:
        #for row in updater:
            #row[1]=row[0]
            #updater.updateRow(row)
    # if wood[i]=="Empty":
    #     pass
    # else:
    #     arcpy.DeleteField_management(sn, ["lw_density","Lw_volume"])
    #     disc_poly3=out_tmp_id+"disc_poly3.shp"
    #     arcpy.analysis.SpatialJoin(lw_file, disc_poly, disc_poly3, match_option="CLOSEST_GEODESIC", search_radius="100 Meters")
    #     stats=out_tmp_id+"stats.dbf"
    #     arcpy.analysis.Statistics(disc_poly3, stats, [["Volume","SUM"]], "GridID")
    #     arcpy.management.JoinField(disc_poly, "GridID", stats, "GridID", ['SUM_Volume'])
    #     arcpy.DeleteField_management(sn, "wood_den")
    #     arcpy.management.AddField(disc_poly, "wood_den", "FLOAT")
    #     arcpy.management.AddGeometryAttributes(disc_poly, "AREA_GEODESIC")
    #     arcpy.CalculateField_management(disc_poly, "wood_den", "!SUM_Volume!/(!AREA_GEO!/10000)", "PYTHON_9.3")
    #     arcpy.management.JoinField(sn, "GridID", disc_poly, "GridID", ['SUM_Volume','wood_den'])
    #     arcpy.CalculateField_management(sn, "lw_volume", "!SUM_Volume!","PYTHON","")
    #     arcpy.CalculateField_management(sn, "lw_density", "!wood_den!","PYTHON","")
    #     arcpy.DeleteField_management(sn,'SUM_Volume;wood_den')
    
    #   ## Sediment Balance
    # arcpy.DeleteField_management(sn, ["Sed_Trans"])
    # arcpy.management.AddField(sn, "Sed_Trans", "FLOAT")
    # with arcpy.da.UpdateCursor(in_table=sn, field_names=['srb_vol_up','df_vol_up','Sed_Trans']) as updater:
    #     for row in updater:
    #         row[2]=row[1]-row[0]
    #         updater.updateRow(row)
    #%%     ## Finishing Touches
    #   ##Add Fire Name
    # arcpy.DeleteField_management(sn, ["Fire_name"])
    # arcpy.management.AddField(sn, "Fire_name", "TEXT")
    # updater=arcpy.da.UpdateCursor(sn, field_names=["Fire_name"])
    # for row in updater:
    #     row[0]=fire[i]
    #     updater.updateRow(row)
        
    #     ##Add Watershed Name
    # arcpy.DeleteField_management(sn, ["Wtshd"])
    # arcpy.management.AddField(sn, "Wtshd", "TEXT")
    # updater=arcpy.da.UpdateCursor(sn, field_names=["Wtshd"])
    # for row in updater:
    #     row[0]=outid[i]
    #     updater.updateRow(row)
    
    
    
    
    