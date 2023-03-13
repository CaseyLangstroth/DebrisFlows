# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 17:44:51 2022

@author: clang
"""
import arcpy
import os


arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')


def CatchStats(wtsd, concov, contrue, appendyr, in_dem, slopetrue1, slopetrue2, dbnr, subcatch, subcatchID, workspace, outdir, basename):

    arcpy.env.workspace = workspace
    path = os.path.join(outdir, 'temp\\')
    os.mkdir(path)
    #### Clip Rasters to watershed
    
    arcpy.management.Clip(in_raster = concov,
                        out_raster = path + 'covwsclip',
                        in_template_dataset = wtsd,
                        clipping_geometry = True)
    arcpy.management.Clip(in_raster = in_dem,                
                        out_raster = path + 'demwsclip',
                        in_template_dataset = wtsd,
                        clipping_geometry = True)
  
    ####Conditional Rasters
    conwhere = 'VALUE = ' + contrue

    CovCON = arcpy.sa.Con(in_conditional_raster = path +'covwsclip',
                            in_true_raster_or_constant = 1,
                            where_clause = conwhere)
    outcon = 'CovCon' + basename + '.tif'
    
    CovCON.save(path + outcon) #cover conditional raster

    
    Slope_Ras = arcpy.sa.Slope(in_raster = path +'demwsclip', 
                           output_measurement = 'DEGREE')
    outslope = 'Slope' + basename + '.tif'
    
    Slope_Ras.save(path+ outslope) #basic slope raster
    
    slopewhere1 = 'VALUE >=' + slopetrue1
    
    SlopeCon = arcpy.sa.Con(in_conditional_raster = path +outslope,
                            in_true_raster_or_constant = 1,
                            where_clause = slopewhere1)  
    outslopecon = 'SlopeCon' + slopetrue1 + basename + '.tif'
    
    SlopeCon.save(path +outslopecon) #slope con raster
    
    if len(slopetrue2) > 0:
        slopewhere2 = 'VALUE >=' + slopetrue2
        SlopeCon2 = arcpy.sa.Con(in_conditional_raster = path + outslope,
                                    in_true_raster_or_constant = 1,
                                    where_clause = slopewhere2)
        outslopecon2 = 'SlopeCon' + slopetrue2 + basename + '.tif'
        
        SlopeCon2.save(path +outslopecon2) #second slope con raster
        
 
    bmhCon = arcpy.sa.Con(in_conditional_raster = dbnr,
                          in_true_raster_or_constant = 1,
                          where_clause = 'VALUE = 3 or VALUE = 4')
    outbmhcon = 'bmh' + basename + '.tif'
    bmhCon.save(path + outbmhcon) #bmh con raster
    
   ####Raster to polygon
      
    arcpy.conversion.RasterToPolygon(in_raster = path + outcon,
                                                 out_polygon_features = path + 'CoverPoly') #cover poly
    outsp1 = path +'SlopePoly' +slopetrue1
    arcpy.conversion.RasterToPolygon(in_raster = path + outslopecon,
                                                 out_polygon_features = outsp1) #slope poly
    if len(slopetrue2) > 0:
        outsp2 = path + 'SlopePoly' +slopetrue2
        arcpy.conversion.RasterToPolygon(in_raster = path + outslopecon2,
                                                     out_polygon_features = outsp2) #slope2 poly
    
        arcpy.conversion.RasterToPolygon(in_raster = path + outbmhcon,
                                               out_polygon_features = path +'BmhPoly') #bmh poly
    
    
    ####Dissolve
    
    arcpy.management.Dissolve(in_features = path + 'CoverPoly',
                              out_feature_class = path + 'CoverPolyDiss')
    sp1diss = path + 'SlopeDiss' + slopetrue1
    arcpy.management.Dissolve(in_features = outsp1,
                              out_feature_class = sp1diss)
    if len(slopetrue2) > 0:
        sp2diss = path + 'SlopeDiss' +slopetrue2
        arcpy.management.Dissolve(in_features = outsp2,
                                  out_feature_class = sp2diss)
    arcpy.management.Dissolve(in_features = path + 'BmhPoly',
                              out_feature_class = path + 'BmhPolyDiss')
    
    ####Clip dissolved features
    
    arcpy.analysis.Clip(in_features= subcatch,
                        clip_features = path + 'CoverPolyDiss',
                        out_feature_class = path + 'CovSubCatchClip')
    
    arcpy.analysis.Clip(in_features= subcatch,
                        clip_features = path + 'BmhPolyDiss',
                        out_feature_class = path + 'BmhSubCatchClip')
    
    sc1 = 'SlopeSubCatchClip' +slopetrue1
    arcpy.analysis.Clip(in_features= subcatch,
                        clip_features = sp1diss,
                        out_feature_class = path + sc1)
    if len(slopetrue2) > 0:
        sc2 = 'SlopeSubCatchClip' +slopetrue2
        arcpy.analysis.Clip(in_features= subcatch,
                            clip_features = sp2diss,
                            out_feature_class = path + sc2)
    
    #### Get area data
    arcpy.management.AddGeometryAttributes(Input_Features = path +'CovSubCatchClip.shp',
                                           Geometry_Properties = 'AREA_GEODESIC',
                                           Area_Unit = 'SQUARE_KILOMETERS')
    arcpy.management.AddGeometryAttributes(Input_Features=  path + sc1 + '.shp',
                                           Geometry_Properties = 'AREA_GEODESIC',
                                           Area_Unit = 'SQUARE_KILOMETERS')
    if len(slopetrue2) > 0:
        arcpy.management.AddGeometryAttributes(Input_Features = path + sc2 +'.shp',
                                       Geometry_Properties = 'AREA_GEODESIC',
                                       Area_Unit = 'SQUARE_KILOMETERS')   
    arcpy.management.AddGeometryAttributes(Input_Features = path + 'BmhSubCatchClip.shp',
                                           Geometry_Properties = 'AREA_GEODESIC',
                                           Area_Unit = 'SQUARE_KILOMETERS')
    arcpy.management.AddGeometryAttributes(Input_Features = subcatch,
                                           Geometry_Properties = 'AREA_GEODESIC',
                                           Area_Unit = 'SQUARE_KILOMETERS')
    ConCovArea = []
    with arcpy.da.SearchCursor(in_table = path + 'CovSubCatchClip.shp', field_names = 'AREA_GEO') as searcher:
        for row in searcher:
            ConCovArea.append(row[0])

    SlopeArea1 = []
    with arcpy.da.SearchCursor(in_table = path + sc1 + '.shp', field_names = 'AREA_GEO') as searcher:
        for row in searcher:
            SlopeArea1.append(row[0])
    if len(slopetrue2) > 0:
        SlopeArea2 = []
        with arcpy.da.SearchCursor(in_table = path +sc2 + '.shp', field_names = 'AREA_GEO') as searcher:
            for row in searcher:
                SlopeArea2.append(row[0])

    BmhArea = []        
    with arcpy.da.SearchCursor(in_table = path+ 'BmhSubCatchClip.shp', field_names = 'AREA_GEO') as searcher:
        for row in searcher:
            BmhArea.append(row[0])
  
    CatchArea = []
    with arcpy.da.SearchCursor(in_table = subcatch, field_names = 'AREA_GEO') as searcher:
        for row in searcher:
            CatchArea.append(row[0])
   
    arcpy.sa.ZonalStatisticsAsTable(in_zone_data = subcatch, 
                                    zone_field = subcatchID, 
                                    in_value_raster = path + 'demwsclip', 
                                    out_table = path + 'DEMStats')
    Elevation = []
    with arcpy.da.SearchCursor(in_table = path + 'DEMStats', field_names = 'MEAN') as searcher:
        for row in searcher:
            Elevation.append(row[0])
    
    relief = []
    with arcpy.da.SearchCursor(in_table = path + 'DEMStats', field_names = 'RANGE') as searcher:
        for row in searcher:
            relief.append(row[0])
    
    arcpy.sa.ZonalStatisticsAsTable(in_zone_data = subcatch,
                                    zone_field = subcatchID,
                                    in_value_raster = path + outslope,
                                    out_table = path + 'SlopeStats')
    Avg_slope = []
    with arcpy.da.SearchCursor(in_table = path + 'SlopeStats', field_names = 'MEAN') as searcher:
        for row in searcher:
            Avg_slope.append(row[0])
     
    #Slope Percent
    slope_perc1 = []
    for i in range(len(CatchArea)):
        slope_perc1.append((SlopeArea1[i]/CatchArea[i])*100)
    
    if len(slopetrue2)>0:
        slope_perc2 = []
        for i in range(len(CatchArea)):
            slope_perc2.append((SlopeArea2[i]/CatchArea[i])*100)
    
                
    #Conifer Cover percent
    Cp = []
    for i in range(len(CatchArea)):
        Cp.append((ConCovArea[i]/CatchArea[i]) *100)
    
    ###slope field names
    slopeperc1 = 'Sp' + slopetrue1 
    if len(slopetrue2)>0:
        slopeperc2 = 'Sp' + slopetrue2 
    Sl1 = 'S' + slopetrue1 + 'Area'
    if len(slopetrue2)>0:
        Sl2 = 'S' + slopetrue2 + 'Area'
    
    
    #### Create table
    out_table = basename + '_CatchmentStats'
    arcpy.management.CreateTable(out_path = workspace,
                                  out_name = out_table)
    
    fields = ['Catch_Area','Mean_Elev', 'Mean_Slope', 'Relief', slopeperc1, 'Cp', 'BmhArea', Sl1]
    
    if len(slopetrue2)> 0:
        fields.insert(5, slopeperc2)
        fields.append(Sl2)
    
    arcpy.management.AddField(in_table = out_table,
                              field_name = 'Site',
                              field_type = 'TEXT')
    for i in fields:
        arcpy.management.AddField(in_table = out_table,
                                  field_name = i,
                                  field_type = 'FLOAT')
    
    #### get site names
    site_names = []
    with arcpy.da.SearchCursor(in_table = subcatch, field_names = subcatchID) as searcher:
        for row in searcher:
            site_names.append(row[0])
    ##input data
    if len(slopetrue2) > 0:
        all_fields = ['Site','Catch_Area','Mean_Elev', 'Mean_Slope', 'Relief', slopeperc1, slopeperc2, 'Cp', 'BmhArea', Sl1, Sl2]
        data = []
        for i in range(len(site_names)):    
            data = []
            data.append(site_names[i])
            data.append(round(CatchArea[i],2))
            data.append(round(Elevation[i],2))
            data.append(round(Avg_slope[i],2))
            data.append(round(relief[i],2))
            data.append(round(slope_perc1[i],2))
            data.append(round(slope_perc2[i],2))
            data.append(round(Cp[i],2))
            data.append(round(BmhArea[i],2))
            data.append(round(SlopeArea1[i],2))
            data.append(round(SlopeArea2[i],2))
            with arcpy.da.InsertCursor(in_table = out_table, field_names = all_fields) as inserter:
                inserter.insertRow(data)
    else:
        all_fields = ['Site','Catch_Area','Mean_Elev', 'Mean_Slope', 'Relief', slopeperc1, 'Cp', 'BmhArea', Sl1]
        for i in range(len(site_names)): 
            data = []
            data.append(site_names[i])
            data.append(round(CatchArea[i],2))
            data.append(round(Elevation[i],2))
            data.append(round(Avg_slope[i],2))
            data.append(round(relief[i],2))
            data.append(round(slope_perc1[i],2))
            data.append(round(Cp[i],2))
            data.append(round(BmhArea[i],2))
            data.append(round(SlopeArea1[i],2))
            with arcpy.da.InsertCursor(in_table = out_table, field_names = all_fields) as inserter:
                inserter.insertRow(data)


    if len(appendyr) > 0:
        new_name = 'Cp' + appendyr
        arcpy.management.AlterField(in_table = out_table, 
                                    field = 'Cp',
                                    new_field_name = new_name,
                                    new_field_alias = new_name)

            
    arcpy.management.DeleteField(in_table = subcatch,
                                 drop_field = 'AREA_GEO')

    #Watershed Stats
    wtsd_table = basename + '_WatershedStats'
    arcpy.management.CreateTable(out_path = workspace,
                                  out_name = wtsd_table)
    
    wtsd_fields = ['Watershed_Area','Mean_Elev', 'Mean_Slope', 'Relief']
    
    for i in wtsd_fields:
        arcpy.management.AddField(in_table = wtsd_table,
                                  field_name = i,
                                  field_type = 'FLOAT')
    arcpy.sa.ZonalStatisticsAsTable(in_zone_data = wtsd,
                                    zone_field = 'Id',
                                    in_value_raster = path + outslope,
                                    out_table = path + 'SlopeStatsWtsd')
    arcpy.sa.ZonalStatisticsAsTable(in_zone_data = wtsd,
                                    zone_field = 'Id',
                                    in_value_raster = path + 'demwsclip',
                                    out_table = path + 'DEMStatsWtsd')
    wtsdslope = []
    with arcpy.da.SearchCursor(in_table = path + 'SlopeStatsWtsd', field_names = 'MEAN') as searcher:
        for row in searcher:
            wtsdslope.append(row[0])
    wtsdelev = []
    with arcpy.da.SearchCursor(in_table = path + 'DEMStatsWtsd', field_names = 'MEAN') as searcher:
        for row in searcher:
            wtsdelev.append(row[0])
    wtsdrelief = []
    with arcpy.da.SearchCursor(in_table = path + 'DEMStatsWtsd', field_names = 'RANGE') as searcher:
        for row in searcher:
            wtsdrelief.append(row[0])
    
    arcpy.management.AddGeometryAttributes(Input_Features = wtsd,
                                   Geometry_Properties = 'AREA_GEODESIC',
                                   Area_Unit = 'SQUARE_KILOMETERS')
    wtsdarea = []
    with arcpy.da.SearchCursor(in_table = wtsd, field_names = 'AREA_GEO') as searcher:
        for row in searcher:
            wtsdarea.append(row[0])
    arcpy.management.DeleteField(in_table = wtsd,
                                 drop_field = 'AREA_GEO')
    
    for i in range(len(wtsdslope)):
        wtsddata = []
        wtsddata.append(round(wtsdarea[i],2))
        wtsddata.append(round(wtsdslope[i],2))
        wtsddata.append(round(wtsdelev[i],2))
        wtsddata.append(round(wtsdrelief[i],2))
        with arcpy.da.InsertCursor(in_table = wtsd_table, field_names = wtsd_fields) as inserter:
            inserter.insertRow(wtsddata)

    return

# This is used to execute code if the file was run but not imported
if __name__ == '__main__':

    # Tool parameter accessed with GetParameter or GetParameterAsText
    wtsd = arcpy.GetParameterAsText(0)
    concov = arcpy.GetParameterAsText(1)
    contrue =arcpy.GetParameterAsText(2)
    appendyr = arcpy.GetParameterAsText(3)
    in_dem =arcpy.GetParameterAsText(4)
    slopetrue1 =arcpy.GetParameterAsText(5)
    slopetrue2 =arcpy.GetParameterAsText(6)
    dbnr =arcpy.GetParameterAsText(7)
    subcatch =arcpy.GetParameterAsText(8)
    subcatchID =arcpy.GetParameterAsText(9)
    workspace =arcpy.GetParameterAsText(10)
    outdir = arcpy.GetParameterAsText(11)
    basename = arcpy.GetParameterAsText(12)
    
    CatchStats(wtsd, concov, contrue, appendyr, in_dem, slopetrue1, slopetrue2, dbnr, subcatch, subcatchID, workspace, outdir, basename)
    
