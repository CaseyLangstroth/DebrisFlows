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
from scipy import stats
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

# Code snippet to make a new directory
def chk_mk_dir(input_dir):
    try:
        os.mkdir(input_dir)
    except:
        pass
 
# Load in excel input data for each watershed
df = pd.read_excel(inData)
c = df.columns
for col in c:
   globals()[col]=np.array(df[col]) 

# Set watershed and area of interest shapefiles
indem=r'D:\GIS\Utah_DEM_utm_tunnel.tif' # Set DEM file
pp=r'D:\Alec\WEST\pp.shp' # Set pour point file
aoi=r'D:\Alec\WEST\wtshd_buff.shp' # Set watershed file

# Create output directory for each watershed and run the USUAL metrics using this watershed abbreviation (start of loop)
for i in np.arange(len(ppfid)):
    print(outid[i])
    outpath=os.path.join(out_folder,outid[i])
    chk_mk_dir(outpath)
    
    #%% Run the USUAL Tools: These are tools developed by David et al., 2022 to prepare stream networks for sediment routing. We use these tools to create a discretized stream network with the attributes of slope, reach length, elevation and valley bottom width.
    # Delineate Watershed and Stream Network
    arcpy.USUAL.WatershedRiverDelineation(indem, 2.5, pp, "FID", str(ppfid[i]), 50, aoi, "FID", str(wtshdfid[i]), outpath, outid[i], False, '', False)
    print ("Watershed Delineation Done")
    indemf=fr"{outpath}\{outid[i]}demf.tif"    # output filled dem file
    fdr=fr"{outpath}\{outid[i]}fdr.tif"        # output flow direction file
    fac=fr"{outpath}\{outid[i]}fac.tif"        # output flow accumulation file
    stream_5k_ras=fr"{outpath}\{outid[i]}_stream_2km.tif"      # output stream network raster file
    stream_5k_shp=fr"{outpath}\{outid[i]}_stream_2km.shp"      # output stream network shapefile
    wtshd=fr"{outpath}\{outid[i]}_watershed.shp"               # output watershed from area of interest
    # Delineate SubCatchments and Interfluves
    arcpy.USUAL.SubCatchmentInterfluveDelineation(indemf, fdr, fac, stream_5k_shp, wtshd, 0.1, outpath, outid[i],'',True,True,False,False)
    print ("SubCatchment Delineation Done")
    # Discretize Stream Network
    arcpy.USUAL.NetworkDiscretizationAttribution(stream_5k_shp, fac, indemf, 500, 0.001, outpath, outid[i], False)
    print ("Network Discretization Done")
    sn=fr"{outpath}\{outid[i]}_network.shp"     # output discretized stream network
    vb= fr'D:\Alec\Input\VBs\{outid[i]}.shp' #Valley Bottom     #valley bottom delineated using the fluvial corridor tool
    # Delineate Valley Bottom using the Fluvial Corridor Tool
    arcpy.FluvialCorridor.ValleyBottom(False, sn, indemf, 100, 100, 10, 100, -4, 4, None, 100, 20, 1000, 100, vb, False)   # run tool manually/fails in arcpy
    print ("Valley Bottom Delineation Done")         
    # Delineate centerline of valley bottom
    centerline=fr'D:\Alec\WEST\Output\centerline.gdb\{outid[i]}'  
    arcpy.topographic.PolygonToCenterline(vb, centerline, None)  
    # Create Transects across Valley Bottom
    arcpy.USUAL.FluvialPolygonTransects(vb, centerline, 10, outpath, outid[i], None, 30, False)
    transects=fr"{outpath}\{outid[i]}transects.shp"
    # Calculate Width Using Transects
    arcpy.USUAL.FluvialLinkAverageWidth(sn, "VB_width", vb, transects, 10, outpath, outid[i], 10, False)
    print ("Width Calculations Done")
    
    #%% Metrics ####
   
   # input variables 
    srb_points=r'D:\Alec\WEST\Input\SRB_model_points_2022.shp' #Sediment Bottleneck points
    dnbr=r'D:\Alec\Input\mtbs8_18_reclassify.tif'#Burn Severity raster from MTBS
    perimeters='D:\Alec\Input\Final_fire_perimeters.shp' # Fire perimeters
    modhigh=r'D:\Alec\Input\mod_high.tif' # Moderate High Burn Severity raster 
    cover=r'D:\Alec\Input\us110\US_110evt_Utah_reclass_new.tif' # Reclassified existing vegetation type raster from Landfire
    conifer=r'D:\Alec\Input\us110\Conifer.tif' # Reclassified conifer raster from Landfire
    lw_file=fr'D:\Alec\Input\Wood\{wood[i]}.shp'  # Corresponding Wood Dataset
    debrisflow=fr'D:\Alec\Input\DF_volume_all2.shp' # Utah Debris Flow Dataset
    tpi=r'D:\Alec\Input\tpi_Ut.tif' # Topographic Position Index Raster for Utah
    slope=r'D:\Alec\Input\slope_utm.tif' # Slope raster for Utah
    herbup=r'D:\Alec\Input\nlcd_92_herbaceous_upland.tif' # Landfire existing vegetation type reclassified as herbaceous upland (for USGS 2 yr flow regression)
    
    #input values ####
    theta=concavity[i] # Concavity measured using Topotoolbox
    referencetheta=0.45
    a=1.505 #channel width power function variables
    b=0.273
    
    out_temp_dir= os.path.join(outpath,'temp','metrics') # make directory name
    chk_mk_dir(out_temp_dir) # make directory
    print(out_temp_dir)
   
    #.......Set up directory with base file id for outputs.......
    out_fid=outid[i]
    out_tmp_id=out_temp_dir+"\\"+out_fid
    out_id=outpath+"\\"+out_fid
    
    disc_poly=outpath+"\\"+outid[i]+"disc_poly.shp"
    wtshds=outpath+"\\"+outid[i]+"wtshds.shp"
    sn=fr"{outpath}\{outid[i]}_network.shp"
    vb= fr'D:\Alec\Input\VBs\{outid[i]}.shp' #Valley Bottom
    
    # %% Geomorphic Metrics ####
      #Normalized Steepness
      #Source- The resilience of logjams to floods- Wohl 2021
      #Calculation- a=Floodplain width b=
      #Importance
    arcpy.DeleteField_management(sn, ["NormSteep"])
    arcpy.management.AddField(sn, "NormSteep", "FLOAT")
    with arcpy.da.UpdateCursor(in_table=sn, field_names=['Slope','usarea_km2','NormSteep']) as updater:
        for row in updater:
            row[2]=row[0]/math.pow(abs(row[1])*1000000,-theta)
            updater.updateRow(row)
    
    arcpy.DeleteField_management(sn, ["NormSteep2"])
    arcpy.management.AddField(sn, "NormSteep2", "FLOAT")
    ksn_list=[]
    with arcpy.da.SearchCursor(in_table=sn, field_names=['Slope','usarea_km2','NormSteep2']) as searcher:
        for row in searcher:
            ksn=row[0]/math.pow(abs(row[1])*1000000,-referencetheta)
            ksn_list.append(ksn)
    ksn_percentile=(iter(stats.rankdata(ksn_list, "average")/len(ksn_list)))
    #ksn_iter=iter(ksn_percentile)
    with arcpy.da.UpdateCursor(in_table=sn, field_names=['NormSteep2']) as updater:
        for row in updater:
            item=next(ksn_percentile)
            row[0]=item
            updater.updateRow(row)
           
                    
      ##Stream Power
      #Source- The resilience of logjams to floods- Wohl 2021
      #Calculation- a=Floodplain width b=
      #Importance
    arcpy.DeleteField_management(sn, ["StreamPow"])
    arcpy.management.AddField(sn, "StreamPow", "FLOAT")
    with arcpy.da.UpdateCursor(in_table=sn, field_names=['Slope','usarea_km2','StreamPow']) as updater:
        for row in updater:
            row[2]=row[0]*row[1]
            updater.updateRow(row)      
    
    
      ## Finite Differencing
      # Source- The resilience of logjams to floods- Wohl 2021
      # Calculation- a=Floodplain width b=
      # Importance
    arcpy.DeleteField_management(sn, ["S_Change","VB_Change","S_Ratio"])
    arcpy.management.AddField(sn, "S_Change", "FLOAT")
    arcpy.management.AddField(sn, "VB_Change", "FLOAT")
    arcpy.management.AddField(sn, "S_Ratio", "FLOAT")
    updater=arcpy.da.UpdateCursor(sn, field_names=['Slope', 'VB_width','S_Change', 'VB_Change','GRIDID','S_Ratio'])
    for row in updater:
        drainsfrom=row[4]
        us_slope=[]
        us_width=[]
        searcher= arcpy.da.SearchCursor(sn, field_names=['ToLink','Slope', 'VB_width'],where_clause=f'ToLink={drainsfrom}')
        for bow in searcher:
            us_slope.append(bow[1])
            us_width.append(bow[2])
        if len(us_slope)==0:
            row[2]=0
        else:
            avg_slope=min(us_slope)
            ds_slope=row[0]
            row[2]=(ds_slope-avg_slope)/ds_slope
            row[5]=avg_slope/ds_slope
        if len(us_width)==0:
            row[3]=0
        else:
            avg_width=max(us_width)
            ds_width=row[1]
            if ds_width==0:
                row[3]=0
            else:
                row[3]=(ds_width-avg_width)/ds_width
        updater.updateRow(row)
    
      ##Floodplain width vs Channel Width
    #Source- The resilience of logjams to floods- Wohl 2021
    #Calculation- a=Floodplain width b=
    #Importance  
    arcpy.DeleteField_management(sn, ["FpChanRat","ChanWidth"])
    arcpy.management.AddField(sn, "FpChanRat", "FLOAT")
    arcpy.management.AddField(sn, "ChanWidth", "FLOAT")
    with arcpy.da.UpdateCursor(in_table=sn, field_names=['usarea_km2','VB_width','FpChanRat','ChanWidth']) as updater:
        for row in updater:
            row[2]=row[1]/(a*math.pow(abs(row[0]),b))
            row[3]=a*math.pow(abs(row[0]),b)
            updater.updateRow(row)
    
      ##Confinement
      #Source- The resilience of logjams to floods- Wohl 2021
      #Calculation- a=Floodplain width b=
      #Importance
    arcpy.DeleteField_management(sn, ["Confine"])
    arcpy.management.AddField(sn, "Confine", "FLOAT")
    arcpy.management.AddField(sn, "chan_buf", "FLOAT")
    arcpy.CalculateField_management(sn, "chan_buf", "!ChanWidth!+15", "PYTHON_9.3")
    sn_buffer=out_temp_dir+"\sn_buffer_"+outid[i]+".shp"
    arcpy.analysis.Buffer(sn, sn_buffer, "chan_buf", "FULL", "FLAT", "NONE", None, "PLANAR")
    boundary_file=out_tmp_id+"boundary_file.shp"
    arcpy.PolygonToLine_management (vb, boundary_file) # turn the valley bottom polygon into a polyline boundary
    sn_intersect=out_temp_dir+"\sn_intersect_"+outid[i]+".shp"
    arcpy.analysis.Intersect(f"{boundary_file};{sn_buffer}", sn_intersect, "ONLY_FID", None, "INPUT")
    arcpy.management.AddGeometryAttributes(sn_intersect, "LENGTH_GEODESIC")
    confine_stats=out_tmp_id+"confine_stats.dbf"
    arcpy.analysis.Statistics(sn_intersect, confine_stats, [["LENGTH_GEO","SUM"]], "FID_sn_buf")
    arcpy.management.JoinField(sn, "FID", confine_stats, "FID_sn_buf", ['SUM_LENGTH'])
    arcpy.CalculateField_management(sn, "Confine", "!SUM_LENGTH!/!Length_m!","PYTHON","")
    arcpy.DeleteField_management(sn, "SUM_LENGTH;chan_buf")

      ##Sinuosity
      #Source- The resilience of logjams to floods- Wohl 2021
      #Calculation- a=Floodplain width b=
      #Importance
    arcpy.DeleteField_management(sn, ["Sinuosity"])
    arcpy.management.AddField(sn, "Sinuosity", "FLOAT")
    sin_line=out_tmp_id+"sin_line.shp"
    sin_points=out_tmp_id+"sin_points.shp"
    updater=arcpy.da.UpdateCursor(sn, field_names=["Length_m",'shape@','Sinuosity'])
    for row in updater:
        sin_line=row[1]
        arcpy.management.FeatureVerticesToPoints(sin_line,sin_points,  "BOTH_ENDS")
        searcher= arcpy.da.SearchCursor(in_table=sin_points, field_names=['shape@'])
        for bow in searcher:
            item=bow[0]
            item2=next(searcher)
            item3=item2[0]
            new_distance=item.distanceTo(item3)
        row[2]=row[0]/new_distance
        updater.updateRow(row)
    
     # %% Local Landcover Variables 
    
    ## Develop Discretized Valley Bottom
    arcpy.CalculateField_management(sn, "vb_buf", "!VB_width!+25", "PYTHON_9.3")
    sn_vb_buffer=out_tmp_id+"sn_vb_buffer.shp"
    arcpy.analysis.Buffer(sn, sn_vb_buffer, "vb_buf", "FULL", "FLAT", "NONE", None, "PLANAR")
    disc_poly=outpath+"\\"+outid[i]+"disc_poly.shp"
    arcpy.analysis.Clip(sn_vb_buffer, vb, disc_poly, None)
    arcpy.DeleteField_management(sn, "vb_buf")
     
    ##Local Dominant Cover Type
    #Source- The resilience of logjams to floods- Wohl 2021
    #Calculation- a=Floodplain width b=
    #Importance
    arcpy.DeleteField_management(sn, ["Cover"])
    arcpy.management.AddField(sn, "Cover", "FLOAT")
    covertable=out_tmp_id+"covertable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(disc_poly, "FID", cover, covertable, "DATA", "MAJORITY", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.management.JoinField(disc_poly, "FID", covertable, "FID_", "MAJORITY")
    arcpy.management.JoinField(sn, "GridID", disc_poly, "GridID","MAJORITY")
    arcpy.DeleteField_management(disc_poly, ["MAJORITY"])
    arcpy.CalculateField_management(sn, "Cover", "!MAJORITY!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MAJORITY"])
    
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
    arcpy.sa.ZonalStatisticsAsTable(disc_poly, "FID", dnbr, burnsevtable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.DeleteField_management(disc_poly, ["MEAN"])
    arcpy.management.JoinField(disc_poly, "FID", burnsevtable, "FID_", "MEAN")
    #arcpy.SpatialJoin_analysis(rivnetmidp, disc_poly, joindata)
    arcpy.management.JoinField(sn, "GridID", disc_poly, "GridID","MEAN")
    arcpy.CalculateField_management(sn, "Burn", "!MEAN!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MEAN"])

    
     # %% Upstream Landcover Variables
    wtshds=outpath+"\\"+outid[i]+"wtshds.shp"
        ##Upstream 
    srs=arcpy.SpatialReference(26912)
    arcpy.CreateFeatureclass_management(out_path=outpath, out_name=f"{outid[i]}wtshds.shp", geometry_type='POLYGON',spatial_reference=srs)
    arcpy.management.AddField(wtshds, "GridID", "FLOAT")
    # Collapse stream network to points
    rivnetmidp=out_tmp_id+"networkmidpoints.shp"
    arcpy.FeatureVerticesToPoints_management(sn,rivnetmidp,"MID")
    searcher= arcpy.da.SearchCursor(in_table=rivnetmidp, field_names=['GridID','shape@'])
    for row in searcher:
          id=row[0]
          pourpoint=row[1]
          flowdir=out_id+"Fdr.tif"
          in_raster=arcpy.gp.Watershed_sa(flowdir, pourpoint)
          wtshd=out_tmp_id+"wtshd.shp"
          polygon=arcpy.conversion.RasterToPolygon(in_raster, wtshd)
          wtshd_dissolve=out_tmp_id+"wtshd_dissolve.shp"
          arcpy.management.Dissolve(wtshd,wtshd_dissolve , "","","SINGLE_PART", "DISSOLVE_LINES")
          searcher2= arcpy.da.SearchCursor(in_table=wtshd_dissolve, field_names=['shape@','SHAPE@AREA'])
          for row in searcher2:
              polygon=row[0]
              area=row[1]
          if area==0:
              pass
          else:
              icur=arcpy.da.InsertCursor(in_table=wtshds,field_names=['GridID','shape@'])
              icur.insertRow([id,polygon])
              arcpy.CalculateField_management(wtshds, "Id", "!GridID!", "PYTHON_9.3")
    print('Done with Upstream!')
    
      ##Cover Type    
      #Source- The resilience of logjams to floods- Wohl 2021
      #Calculation- a=Floodplain width b=
      #Importance
    arcpy.DeleteField_management(sn, ["Cover_shd"])
    arcpy.management.AddField(sn, "Cover_shd", "FLOAT")
    covershdtable=out_tmp_id+"covershdtable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", cover, covershdtable, "DATA", "MAJORITY", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.management.JoinField(sn, "GridID", covershdtable, "Id", "MAJORITY")
    arcpy.CalculateField_management(sn, "Cover_Shd", "!MAJORITY!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MAJORITY"])
    
      ##Conifer
      #Source- The resilience of logjams to floods- Wohl 2021
      #Calculation- a=Floodplain width b=
      #Importance
    arcpy.DeleteField_management(sn, ["Conifer"])
    arcpy.management.AddField(sn, "Conifer", "FLOAT")
    conifertable=out_tmp_id+"conifertable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", conifer, conifertable, "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.management.AddField(conifertable, "conifer_up", "FLOAT")
    arcpy.CalculateField_management(conifertable, "conifer_up", "!SUM!*900/!AREA!*100", "PYTHON_9.3")
    arcpy.management.JoinField(sn, "GridID",conifertable , "Id", "conifer_up")
    arcpy.CalculateField_management(sn, "Conifer", "!conifer_up!","PYTHON","")
    arcpy.DeleteField_management(sn, ["conifer_up"])
    
      ## TPI Upstream
      #Source- Weiss (2000) Topographic Position and Landforms Analysis
      #Calculation- Average topographic position index in upstream draining area
      #Importance- The value represents how inset the stream network upstream of the reach, potentially indicating higher sediment transport
    arcpy.DeleteField_management(sn, ["tpi_up"])
    arcpy.management.AddField(sn, "tpi_up", "FLOAT")
    wtshds_tpi=out_tmp_id+"wtshds_tpi.shp"
    arcpy.analysis.Intersect(f"{wtshds};{sn}", wtshds_tpi, "ALL", None, "INPUT")
    #arcpy.SpatialJoin_analysis(wtshds, sn, wtshds_tpi, "JOIN_ONE_TO_MANY")
    tpi_dissolve=out_tmp_id+"tpi_dissolve"
    arcpy.management.Dissolve(wtshds_tpi,tpi_dissolve , "GridID", "TPI MEAN", "MULTI_PART", "DISSOLVE_LINES")
    arcpy.management.JoinField(sn, "GridID",tpi_dissolve , "GridID", "MEAN_TPI")
    arcpy.CalculateField_management(sn, "tpi_up", "!MEAN_TPI!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MEAN_TPI"])
    # #arcpy.analysis.SpatialJoin("sf_watershed", "sf_network", r"C:\Users\BelmontLab\Documents\ArcGIS\Projects\MyProject7\MyProject7.gdb\sf_watershed_SpatialJoin", "JOIN_ONE_TO_MANY", "KEEP_ALL", 'Id "Id" true true false 10 Long 0 10,First,#,sf_watershed,Id,-1,-1;gridcode "gridcode" true true false 10 Long 0 10,First,#,sf_watershed,gridcode,-1,-1;Id_1 "Id" true true false 6 Long 0 6,First,#,sf_network,Id,-1,-1;ORIG_FID "ORIG_FID" true true false 10 Long 0 10,First,#,sf_network,ORIG_FID,-1,-1;ORIG_SEQ "ORIG_SEQ" true true false 10 Long 0 10,First,#,sf_network,ORIG_SEQ,-1,-1;GridID "GridID" true true false 19 Double 0 0,First,#,sf_network,GridID,-1,-1;Length_m "Length_m" true true false 19 Double 0 0,First,#,sf_network,Length_m,-1,-1;ToLink "ToLink" true true false 19 Double 0 0,First,#,sf_network,ToLink,-1,-1;usarea_km2 "usarea_km2" true true false 19 Double 0 0,First,#,sf_network,usarea_km2,-1,-1;uselev_m "uselev_m" true true false 19 Double 0 0,First,#,sf_network,uselev_m,-1,-1;dselev_m "dselev_m" true true false 19 Double 0 0,First,#,sf_network,dselev_m,-1,-1;Slope "Slope" true true false 19 Double 0 0,First,#,sf_network,Slope,-1,-1;VB_width "VB_width" true true false 19 Double 0 0,First,#,sf_network,VB_width,-1,-1;StreamPow "StreamPow" true true false 13 Float 0 0,First,#,sf_network,StreamPow,-1,-1;S_Change "S_Change" true true false 13 Float 0 0,First,#,sf_network,S_Change,-1,-1;VB_Change "VB_Change" true true false 13 Float 0 0,First,#,sf_network,VB_Change,-1,-1;FpChanRat "FpChanRat" true true false 13 Float 0 0,First,#,sf_network,FpChanRat,-1,-1;ChanWidth "ChanWidth" true true false 13 Float 0 0,First,#,sf_network,ChanWidth,-1,-1;Confine "Confine" true true false 13 Float 0 0,First,#,sf_network,Confine,-1,-1;Sinuosity "Sinuosity" true true false 13 Float 0 0,First,#,sf_network,Sinuosity,-1,-1;Cover "Cover" true true false 13 Float 0 0,First,#,sf_network,Cover,-1,-1;TPI "TPI" true true false 13 Float 0 0,First,#,sf_network,TPI,-1,-1;Burn "Burn" true true false 13 Float 0 0,First,#,sf_network,Burn,-1,-1;Cover_shd "Cover_shd" true true false 13 Float 0 0,First,#,sf_network,Cover_shd,-1,-1;Conifer "Conifer" true true false 13 Float 0 0,First,#,sf_network,Conifer,-1,-1;Mod_High "Mod_High" true true false 13 Float 0 0,First,#,sf_network,Mod_High,-1,-1;mh_sqm "mh_sqm" true true false 13 Float 0 0,First,#,sf_network,mh_sqm,-1,-1;per_burn "per_burn" true true false 13 Float 0 0,First,#,sf_network,per_burn,-1,-1;srb_up "srb_up" true true false 13 Float 0 0,First,#,sf_network,srb_up,-1,-1;srb_vol_up "srb_vol_up" true true false 13 Float 0 0,First,#,sf_network,srb_vol_up,-1,-1;df_vol "df_vol" true true false 13 Float 0 0,First,#,sf_network,df_vol,-1,-1;df_vol_up "df_vol_up" true true false 13 Float 0 0,First,#,sf_network,df_vol_up,-1,-1;df_up "df_up" true true false 254 Text 0 0,First,#,sf_network,df_up,0,254;Sed_Trans "Sed_Trans" true true false 13 Float 0 0,First,#,sf_network,Sed_Trans,-1,-1;Fire_name "Fire_name" true true false 254 Text 0 0,First,#,sf_network,Fire_name,0,254;Wtshd "Wtshd" true true false 254 Text 0 0,First,#,sf_network,Wtshd,0,254;tpi_up "tpi_up" true true false 13 Float 0 0,First,#,sf_network,tpi_up,-1,-1;NormSteep "NormSteep" true true false 13 Float 0 0,First,#,sf_network,NormSteep,-1,-1;srb_vol "srb_vol" true true false 13 Float 0 0,First,#,sf_network,srb_vol,-1,-1;srb "srb" true true false 13 Float 0 0,First,#,sf_network,srb,-1,-1;VBM_SF "VBM_SF" true true false 254 Text 0 0,First,#,sf_network,VBM_SF,0,254;Cause "Cause" true true false 254 Text 0 0,First,#,sf_network,Cause,0,254', "INTERSECT", None, '')
    # #arcpy.analysis.SpatialJoin("sfwtshds", r"Networks\sf_network", r"D:\Alec\WEST\WEST.gdb\sfwtshds_SpatialJoin", "JOIN_ONE_TO_MANY", "KEEP_ALL", r'Id "Id" true true false 6 Long 0 6,First,#,sfwtshds,Id,-1,-1;GridID "GridID" true true false 13 Float 0 0,First,#,sfwtshds,GridID,-1,-1;AREA_GEO "AREA_GEO" true true false 19 Double 0 0,First,#,sfwtshds,AREA_GEO,-1,-1;perburn "perburn" true true false 13 Float 0 0,First,#,sfwtshds,perburn,-1,-1;Id_1 "Id" true true false 6 Long 0 6,First,#,Networks\sf_network,Id,-1,-1;ORIG_FID "ORIG_FID" true true false 10 Long 0 10,First,#,Networks\sf_network,ORIG_FID,-1,-1;ORIG_SEQ "ORIG_SEQ" true true false 10 Long 0 10,First,#,Networks\sf_network,ORIG_SEQ,-1,-1;GridID_1 "GridID" true true false 19 Double 0 0,First,#,Networks\sf_network,GridID,-1,-1;Length_m "Length_m" true true false 19 Double 0 0,First,#,Networks\sf_network,Length_m,-1,-1;ToLink "ToLink" true true false 19 Double 0 0,First,#,Networks\sf_network,ToLink,-1,-1;usarea_km2 "usarea_km2" true true false 19 Double 0 0,First,#,Networks\sf_network,usarea_km2,-1,-1;uselev_m "uselev_m" true true false 19 Double 0 0,First,#,Networks\sf_network,uselev_m,-1,-1;dselev_m "dselev_m" true true false 19 Double 0 0,First,#,Networks\sf_network,dselev_m,-1,-1;Slope "Slope" true true false 19 Double 0 0,First,#,Networks\sf_network,Slope,-1,-1;VB_width "VB_width" true true false 19 Double 0 0,First,#,Networks\sf_network,VB_width,-1,-1;StreamPow "StreamPow" true true false 13 Float 0 0,First,#,Networks\sf_network,StreamPow,-1,-1;S_Change "S_Change" true true false 13 Float 0 0,First,#,Networks\sf_network,S_Change,-1,-1;VB_Change "VB_Change" true true false 13 Float 0 0,First,#,Networks\sf_network,VB_Change,-1,-1;FpChanRat "FpChanRat" true true false 13 Float 0 0,First,#,Networks\sf_network,FpChanRat,-1,-1;ChanWidth "ChanWidth" true true false 13 Float 0 0,First,#,Networks\sf_network,ChanWidth,-1,-1;Confine "Confine" true true false 13 Float 0 0,First,#,Networks\sf_network,Confine,-1,-1;Sinuosity "Sinuosity" true true false 13 Float 0 0,First,#,Networks\sf_network,Sinuosity,-1,-1;Cover "Cover" true true false 13 Float 0 0,First,#,Networks\sf_network,Cover,-1,-1;TPI "TPI" true true false 13 Float 0 0,First,#,Networks\sf_network,TPI,-1,-1;Cover_shd "Cover_shd" true true false 13 Float 0 0,First,#,Networks\sf_network,Cover_shd,-1,-1;Conifer "Conifer" true true false 13 Float 0 0,First,#,Networks\sf_network,Conifer,-1,-1;Mod_High "Mod_High" true true false 13 Float 0 0,First,#,Networks\sf_network,Mod_High,-1,-1;mh_sqm "mh_sqm" true true false 13 Float 0 0,First,#,Networks\sf_network,mh_sqm,-1,-1;per_burn "per_burn" true true false 13 Float 0 0,First,#,Networks\sf_network,per_burn,-1,-1;srb_up "srb_up" true true false 13 Float 0 0,First,#,Networks\sf_network,srb_up,-1,-1;srb_vol_up "srb_vol_up" true true false 13 Float 0 0,First,#,Networks\sf_network,srb_vol_up,-1,-1;df_vol "df_vol" true true false 13 Float 0 0,First,#,Networks\sf_network,df_vol,-1,-1;df_vol_up "df_vol_up" true true false 13 Float 0 0,First,#,Networks\sf_network,df_vol_up,-1,-1;df_up "df_up" true true false 254 Text 0 0,First,#,Networks\sf_network,df_up,0,254;Sed_Trans "Sed_Trans" true true false 13 Float 0 0,First,#,Networks\sf_network,Sed_Trans,-1,-1;Fire_name "Fire_name" true true false 254 Text 0 0,First,#,Networks\sf_network,Fire_name,0,254;Wtshd "Wtshd" true true false 254 Text 0 0,First,#,Networks\sf_network,Wtshd,0,254;srb_vol "srb_vol" true true false 13 Float 0 0,First,#,Networks\sf_network,srb_vol,-1,-1;srb "srb" true true false 13 Float 0 0,First,#,Networks\sf_network,srb,-1,-1;VBM_SF "VBM_SF" true true false 254 Text 0 0,First,#,Networks\sf_network,VBM_SF,0,254;Cause "Cause" true true false 254 Text 0 0,First,#,Networks\sf_network,Cause,0,254;NormSteep "NormSteep" true true false 13 Float 0 0,First,#,Networks\sf_network,NormSteep,-1,-1;Burn "Burn" true true false 13 Float 0 0,First,#,Networks\sf_network,Burn,-1,-1;tpi_up "tpi_up" true true false 13 Float 0 0,First,#,Networks\sf_network,tpi_up,-1,-1', "INTERSECT", None, '')
        
        ##Mod High
        #Source- The resilience of logjams to floods- Wohl 2021
        #Calculation- a=Floodplain width b=
        #Importance
    arcpy.DeleteField_management(sn, ["Mod_High"])
    arcpy.management.AddField(sn, "Mod_High", "FLOAT")
    modhightable=out_tmp_id+"modhightable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", modhigh, modhightable, "DATA", "SUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.management.AddField(modhightable, "modhigh", "FLOAT")
    arcpy.CalculateField_management(modhightable, "modhigh", "!SUM!/!AREA!*100", "PYTHON_9.3")
    arcpy.management.JoinField(sn, "GridID", modhightable, "Id", "modhigh")
    arcpy.CalculateField_management(sn, "Mod_High", "!modhigh!","PYTHON","")
    arcpy.DeleteField_management(sn, ["modhigh"])
    
        ##Mod High Square Meters
        #Source- The resilience of logjams to floods- Wohl 2021
        #Calculation- a=Floodplain width b=
        #Importance
    arcpy.DeleteField_management(sn, ["mh_sqm"])
    arcpy.management.AddField(sn, "mh_sqm", "FLOAT")
    arcpy.CalculateField_management(modhightable, "mhsqm", "!SUM!", "PYTHON_9.3")
    arcpy.management.JoinField(sn, "GridID", modhightable, "Id", "mhsqm")
    arcpy.CalculateField_management(sn, "mh_sqm", "!mhsqm!","PYTHON","")
    arcpy.DeleteField_management(sn, ["mhsqm"])
   
        ##PRISM Data
        #Source-Beechie and Imacki (2014) Predicting natural channel patterns based on landscape and geomorphic controls in the Columbia River basin, USA: Predicting Channel Patterns in the Columbia Basin
        #Calculation- Power function to calculate channel width
        #Importance- Alternative way to calculate channel width
    arcpy.sa.ZonalStatisticsAsTable(wtshds, "GridID", prism, prism_shd, "DATA", "MAXIMUM", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.management.AddField(prism_shd, "prism_shd", "FLOAT")
    arcpy.CalculateField_management(prism_shd, "prism_shd", "!MAX!", "PYTHON_9.3")
    arcpy.management.JoinField(sn, "GridID", prism_shd, "GridID", "prism_shd")
    print("Help")
    
        #Total Burn Percentage
        #Source-
        #Calculation- Intersects the burn perimeter with the upstream draining area, and divides that intersected area by the upstream draining area
        #Importance- Helpful for checking for site criteria
    arcpy.DeleteField_management(sn, ["per_burn"])
    arcpy.management.AddField(sn, "per_burn", "FLOAT")
    arcpy.management.AddGeometryAttributes(wtshds, "AREA_GEODESIC")
    arcpy.management.AddField(wtshds, "perburn", "FLOAT")
    perintersect=out_tmp_id+"perintersect.shp"
    water_feature_name=out_fid+"water_feature"
    water_feature=out_tmp_id+"water_feature"
    updater=arcpy.da.UpdateCursor(wtshds, field_names=['AREA_GEO', 'perburn','SHAPE@'])
    for row in updater:
        waterarea=row[0]
        water=row[2]
        arcpy.conversion.FeatureClassToFeatureClass(water,out_temp_dir, water_feature_name)
        arcpy.analysis.Intersect(f"{water_feature};{perimeters}", perintersect)
        arcpy.management.AddGeometryAttributes(perintersect, "AREA_GEODESIC")
        searcher= arcpy.da.SearchCursor(perintersect, field_names=['AREA_GEO'])
        for bow in searcher:
            if bow[0] is None:
                burnarea=0
            else:
                burnarea=bow[0]
        row[1]=burnarea/waterarea*100
        #print(row[1])
        updater.updateRow(row)
    arcpy.management.JoinField(sn, "GridID", wtshds, "GridID", ['perburn'])
    arcpy.CalculateField_management(sn, "per_burn", "!perburn!","PYTHON","")
    arcpy.DeleteField_management(sn, "perburn")
    
    #%% Bottlenecks, Debris Flows and Wood
    
      ##Bottlenecks Present 
      #Source-
      #Calculation- Associates the attributes of sediment bottlenecks with the reach they are located in, including Presence, Cause, Volume 
      #Importance- Sediment Bottleneck Presence and Volume Response Variables
    arcpy.DeleteField_management(sn, ["srb_vol","srb","VBM_SF","Cause","srb_vbm_sf","srb_cause"])
    arcpy.management.AddField(sn, "srb_vol", "FLOAT")
    arcpy.management.AddField(sn, "srb", "FLOAT")
    arcpy.management.AddField(sn, "VBM_SF", "TEXT")
    arcpy.management.AddField(sn, "Cause", "TEXT")
    disc_poly=outpath+"\\"+outid[i]+"disc_poly.shp"
    disc_poly2=out_tmp_id+"disc_poly2"
    arcpy.SpatialJoin_analysis(disc_poly, srb_points, disc_poly2, "JOIN_ONE_TO_MANY")
    disc_poly2_dissolve=out_tmp_id+"disc_poly2_dissolve"
    arcpy.management.Dissolve(disc_poly2,disc_poly2_dissolve , "GridID", "Join_Count SUM;Volume SUM;VBM_SF FIRST;Cause FIRST", "MULTI_PART", "DISSOLVE_LINES")
    arcpy.management.JoinField(sn, "GridID", disc_poly2_dissolve, "GridID",['SUM_Join_C','FIRST_VBM_','FIRST_Caus','SUM_Volume'])
    arcpy.CalculateField_management(sn, "srb_vol", "!SUM_Volume!","PYTHON","")
    arcpy.CalculateField_management(sn, "srb", "!SUM_Join_C!","PYTHON","")
    arcpy.CalculateField_management(sn, "VBM_SF", "!FIRST_VBM_!","PYTHON","")
    arcpy.CalculateField_management(sn, "Cause", "!FIRST_Caus!","PYTHON","")
    arcpy.DeleteField_management(sn,'SUM_Join_C;SUM_Volume;FIRST_VBM_;FIRST_Caus')
    
      ##Bottleneck Upstream
     #Source-
     #Calculation- The sum of sediment bottleneck volume delivered in the area draining to the reach
     #Importance- Sediment bottlenecks upstream might mean less sediment is delivered downstream
    arcpy.DeleteField_management(sn, ["srb_up","srb_vol_up"])
    arcpy.management.AddField(sn, "srb_up", "FLOAT")
    arcpy.management.AddField(sn, "srb_vol_up", "FLOAT")
    srb_up=out_tmp_id+"srb_up"
    arcpy.SpatialJoin_analysis(wtshds, srb_points, srb_up, "JOIN_ONE_TO_MANY")
    srb_up_dissolve=out_tmp_id+"srb_up_dissolve"
    arcpy.management.Dissolve(srb_up,srb_up_dissolve , "GriDID", "Join_Count SUM;Volume SUM", "MULTI_PART", "DISSOLVE_LINES")
    arcpy.management.JoinField(sn, "GridID", srb_up_dissolve, "GridID", ['SUM_Join_C','SUM_Volume'])
    arcpy.CalculateField_management(sn, "srb_up", "!SUM_Join_C!","PYTHON","")
    arcpy.CalculateField_management(sn, "srb_vol_up", "!SUM_Volume!","PYTHON","")
    arcpy.DeleteField_management(sn,'SUM_Join_C;SUM_Volume')
    
      ##Debris Flow Delivery
      #Source-Murphy et al., 2019 Post‐wildfire sediment cascades: A modeling framework linking debris flow generation and network‐scale sediment routing
      #Calculation- The proportion of the debris flow area intersecting the modeled channel, multiplied by the modeled volume of the debris flow (Wall et al., 2022)
      #Importance- Helps calculate the next metric
    arcpy.DeleteField_management(sn, ["df_vol"])
    arcpy.management.AddField(sn, "df_vol", "FLOAT")
    sn_buffer2=out_tmp_id+"sn_buffer2.shp"
    arcpy.analysis.Buffer(sn, sn_buffer2, "Chanwidth", "FULL", "FLAT", "NONE", None, "PLANAR")
    sn_intersect2=out_tmp_id+"sn_intersect2.shp"
    arcpy.analysis.PairwiseIntersect(f"{debrisflow};{sn_buffer2}", sn_intersect2, "ALL", None, "INPUT")
    arcpy.management.AddField(sn_intersect2, "Vol_input", "FLOAT")
    arcpy.management.AddGeometryAttributes(sn_intersect2, "AREA_GEODESIC")
    arcpy.CalculateField_management(sn_intersect2, "Vol_input", "!AREA_GEO!/!m2!*!m3!", "PYTHON_9.3")
    arcpy.management.JoinField(sn, "FID", sn_intersect2, "ORIG_FID", ['Vol_input'])
    arcpy.CalculateField_management(sn, "df_vol", "!Vol_input!","PYTHON","")
    arcpy.DeleteField_management(sn,'Vol_input')
    
      ##Debris Flow Delivery Upstream
      #Source-Murphy et al., 2019 Post‐wildfire sediment cascades: A modeling framework linking debris flow generation and network‐scale sediment routing
      #Calculation- The sum of debris flows volume delivered in the area draining to the reach
      #Importance- Debris flows are the primary sources of sediment following wildfire
    arcpy.DeleteField_management(sn, ["df_up","df_vol_up"])
    arcpy.management.AddField(sn, "df_vol_up", "FLOAT")
    df_up=out_tmp_id+"df_up.shp"
    arcpy.SpatialJoin_analysis(wtshds, sn_intersect2 , df_up, "JOIN_ONE_TO_MANY")
    df_up_dissolve=out_tmp_id+"df_up_dissolve.shp"
    arcpy.management.Dissolve(df_up,df_up_dissolve , "GridID", "Join_Count SUM;Vol_input SUM", "MULTI_PART", "DISSOLVE_LINES")
    arcpy.management.JoinField(sn, "GridID", df_up_dissolve, "GridID", ['SUM_Join_C','SUM_Vol_in'])
    arcpy.CalculateField_management(sn, "df_up", "!SUM_Join_C!","PYTHON","")
    arcpy.CalculateField_management(sn, "df_vol_up", "!SUM_Vol_in!","PYTHON","")
    arcpy.DeleteField_management(sn,'SUM_Join_C;SUM_Vol_in')
    
      # Wood Volume
      #Source- See Wood Volume Calculation Code
      #Calculation- Assigns a volume of wood to a reach and then sums the total volume of wood for that reach, divide the volume of wood by the area of the reach to get the density
      #Importance- Response Variable
    with arcpy.da.UpdateCursor(in_table=disc_poly, field_names=['FID','ORIG_FID']) as updater:
        for row in updater:
            row[1]=row[0]
            updater.updateRow(row)
    if wood[i]=="Empty":
        pass
    else:
        arcpy.DeleteField_management(sn, ["Lw_density","lw_volume","lw_density"])
        disc_poly3=out_tmp_id+"disc_poly3.shp"
        arcpy.analysis.SpatialJoin(lw_file, disc_poly, disc_poly3, match_option="CLOSEST_GEODESIC", search_radius="100 Meters")
        stats=out_tmp_id+"stats.dbf"
        arcpy.analysis.Statistics(disc_poly3, stats, [["Volume","SUM"]], "GridID")
        arcpy.management.JoinField(disc_poly, "GridID", stats, "GridID", ['SUM_Volume'])
        arcpy.DeleteField_management(sn, "wood_den")
        arcpy.management.AddField(disc_poly, "wood_den", "FLOAT")
        arcpy.management.AddGeometryAttributes(disc_poly, "AREA_GEODESIC")
        arcpy.CalculateField_management(disc_poly, "wood_den", "!SUM_Volume!/(!AREA_GEO!/10000)", "PYTHON_9.3")
        arcpy.management.JoinField(sn, "GridID", disc_poly, "GridID", ['SUM_Volume','wood_den'])
        arcpy.CalculateField_management(sn, "lw_volume", "!SUM_Volume!","PYTHON","")
        arcpy.CalculateField_management(sn, "lw_density", "!wood_den!","PYTHON","")
        arcpy.DeleteField_management(sn,'SUM_Volume;wood_den')
    
        ## Sediment Balance
        #Source-
        #Calculation- The amount of debris flows upstream minus the amount of sediment bottlenecks upstream-a surrogate for how much sediment is being transported through the reach
        #Importance- Potentially an interesting response variable, looking at where sediment is transported vs not
    arcpy.DeleteField_management(sn, ["Sed_Trans"])
    arcpy.management.AddField(sn, "Sed_Trans", "FLOAT")
    with arcpy.da.UpdateCursor(in_table=sn, field_names=['srb_vol_up','df_vol_up','Sed_Trans']) as updater:
        for row in updater:
            row[2]=row[1]-row[0]
            updater.updateRow(row)
   
      #%% Bankfull Flow and Grain Size
    
        ## Calculate 2 year/Bankfull Discharge
        #Source- Kenney et al. (2007) Methods for estimating magnitude and frequency of peak flows for natural streams in Utah
        #Calculation- Regression equations developed for each region using raster variables calculated below
        #Importance- 2 yr bankfull flows represent the average maximum flow of the year which does the most to shape the channel
    arcpy.DeleteField_management(sn, ["q2"])
    arcpy.management.AddField(sn, "q2", "FLOAT")
    ## Average Basin Elevation
    arcpy.DeleteField_management(sn, ["mean_elev"])
    arcpy.management.AddField(sn, "mean_elev", "FLOAT")
    meanelevtable=out_tmp_id+"meanelevtable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", indem, meanelevtable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.CalculateField_management(meanelevtable, "meanelevup", "!MEAN!*3.28", "PYTHON_9.3")
    arcpy.management.JoinField(sn, "GridID",meanelevtable , "Id", "MEAN")
    arcpy.CalculateField_management(sn, "mean_elev", "!MEAN!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MEAN"])
    
    ## Average Upstream Slope
    arcpy.DeleteField_management(sn, ["mean_slope"])
    arcpy.management.AddField(sn, "mean_slope", "FLOAT")
    meanslopetable=out_tmp_id+"meanslopetable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", slope, meanslopetable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.management.JoinField(sn, "GridID",meanslopetable , "Id", "MEAN")
    arcpy.CalculateField_management(sn, "mean_slope", "!MEAN!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MEAN"])
   
    ## Average Upstream Herbaceous Area
    arcpy.DeleteField_management(sn, ["mean_herb"])
    arcpy.management.AddField(sn, "mean_herb", "FLOAT")
    meanherbtable=out_tmp_id+"meaneherbtable.dbf"
    arcpy.sa.ZonalStatisticsAsTable(wtshds, "Id", herbup, meanherbtable, "DATA", "MEAN", "CURRENT_SLICE", 90, "AUTO_DETECT")
    arcpy.management.JoinField(sn, "GridID",meanherbtable , "Id", "MEAN")
    arcpy.CalculateField_management(sn, "mean_herb", "!MEAN!","PYTHON","")
    arcpy.DeleteField_management(sn, ["MEAN"])
    
    ## Regional Regression Calculations
    if region[i] == 1:
        
        updater=arcpy.da.UpdateCursor(sn, field_names=['q2', 'usarea_km2','mean_elev'])
        for row in updater:
                row[0]=1.52*math.pow(row[1],0.677)*math.pow(1.39,row[2]*3.28/1000)
                updater.updateRow(row)
    elif region[i] == 3:
        updater=arcpy.da.UpdateCursor(sn, field_names=['q2', 'usarea_km2'])
        for row in updater:
                row[0]=14.5*math.pow(row[1],0.328)
                updater.updateRow(row)
    elif region[i] == 4:
        updater=arcpy.da.UpdateCursor(sn, field_names=['q2', 'usarea_km2','mean_elev','mean_slope'])
        for row in updater:
                row[0]=0.083*math.pow(row[1],0.822)*math.pow(2.72,0.656*(row[2]*3.28/1000)-0.039*row[3])
                updater.updateRow(row)
    elif region[i] == 5:
        updater=arcpy.da.UpdateCursor(sn, field_names=['q2', 'usarea_km2','mean_herb'])
        for row in updater:
                row[0]=4.32*math.pow(row[1],0.623)*math.pow(row[2]*100+1,0.503)
                updater.updateRow(row)
    elif region[i] == 6:
        updater=arcpy.da.UpdateCursor(sn, field_names=['q2', 'usarea_km2','mean_elev'])
        for row in updater:
                row[0]=4150*math.pow(row[1],0.553)*math.pow(row[2]*3.28/1000,-2.45)
                updater.updateRow(row)
    elif region[i] == 7:
        updater=arcpy.da.UpdateCursor(sn, field_names=['q2', 'usarea_km2'])
        for row in updater:
                row[0]=18.4*math.pow(row[1],0.630)
                updater.updateRow(row)
    else:
        pass
    
      ##Stream Power
      #Source- Common Geomorphic Metric
      #Calculation- density of water*gravity*2 yr flow*CFS to CMS conversion factor*slope
      #Importance- Energy/Ability of stream to transport sediment
    arcpy.DeleteField_management(sn, ["StreamPow"])
    arcpy.management.AddField(sn, "StreamPow", "FLOAT")
    with arcpy.da.UpdateCursor(in_table=sn, field_names=['Slope','q2','StreamPow']) as updater:
        for row in updater:
            row[2]=row[0]*row[1]*0.028*1000*9.8
            updater.updateRow(row)   
    
    ## Pre-fire Grain Size Predictions
    #Source- Snyder et al. (2013) Predicting grain size in gravel-bedded rivers using digital elevation models: Application to three Maine watersheds
    #Calculation- predict grain size based on 2 yr flow, drainage area, slope and channel width (Equation 2)
    #Importance- Surrogate for pre-fire grain size, compare against measured D50
    arcpy.DeleteField_management(sn, ["D50"])
    arcpy.management.AddField(sn, "D50", "FLOAT")
    updater=arcpy.da.UpdateCursor(sn, field_names=['D50','q2', 'usarea_km2','Slope','ChanWidth'])
    for row in updater:
        row[0]=(1000*9.8*math.pow(0.04,3/5)*math.pow(row[2]*0.028,3/5)*math.pow(row[4],-3/5)*math.pow(row[3],7/10))/((2650-1000)*9.8*0.04)*1000
        updater.updateRow(row)
    
    #%%     ## Finishing Touches
      ##Add Fire Name to Field
    arcpy.DeleteField_management(sn, ["Fire_name"])
    arcpy.management.AddField(sn, "Fire_name", "TEXT")
    updater=arcpy.da.UpdateCursor(sn, field_names=["Fire_name"])
    for row in updater:
        row[0]=fire[i]
        updater.updateRow(row)
        
        ##Add Watershed Name to Field
    arcpy.DeleteField_management(sn, ["Wtshd"])
    arcpy.management.AddField(sn, "Wtshd", "TEXT")
    updater=arcpy.da.UpdateCursor(sn, field_names=["Wtshd"])
    for row in updater:
        row[0]=outid[i]
        updater.updateRow(row)
    
    
    
    
    