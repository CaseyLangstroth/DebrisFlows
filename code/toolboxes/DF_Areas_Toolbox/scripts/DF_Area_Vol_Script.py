#Import Packages

import arcpy
from arcpy import env
import math


aprx = arcpy.mp.ArcGISProject('CURRENT')
arcpy.env.overwriteOutput = True #allow overwrite



def AreaVol(param0, param1, param2, param3, param4):
    
### ---------------------------------------------------Calculate geometry and add to table (indv input)-----------------------------------------------------------------------------------------------------###
    if len(param1) == 0:  #param1 not chosen (indv input)
        for x in param0: #multiple inputs allowed, do them all            
        ##--------------------------------------------------------------------------------Area-------------------------------------------------------------------------##
            if len(param3) == 0: #if no update area field chosen             
                    arcpy.management.AddGeometryAttributes(Input_Features = x,
                                                           Geometry_Properties = 'AREA_GEODESIC',
                                                           Area_Unit = param2) #Add geometry attributes to table
                   
                    arcpy.management.AlterField(in_table = x,
                                                field = 'AREA_GEO',
                                                new_field_name = 'Area',
                                                new_field_alias = 'Area') #change field name to Area from AREA_GEO
                    
            if len(param3) > 0: #if update area field chosen
                    arcpy.management.DeleteField(in_table = x, drop_field = param3) #delete the field since I don't know data type

                    arcpy.management.AddGeometryAttributes(Input_Features = x,
                                                           Geometry_Properties = 'AREA_GEODESIC',
                                                           Area_Unit = param2)
                    arcpy.management.AlterField(in_table = x,
                                                field = 'AREA_GEO',
                                                new_field_name = param3,
                                                new_field_alias = param3)
                            
            ##-----------------------------------------------------------------------Volume------------------------------------------------------------------------------##
        for x in param0:
            if len(param4) == 0: #if no update volume field chosen
                    arcpy.management.AddField(in_table = x, #make new volume field
                                              field_name = 'EstVol',
                                              field_type = 'FLOAT')
                    
                    searcher = arcpy.da.SearchCursor(in_table = x, field_names = '*') #search if area field is == Area or param3
                    field_names = searcher.fields #get field names
                    for field in field_names:
                        if field == 'Area':
                            arcpy.management.CalculateField(in_table = x,
                                                            field = 'EstVol',
                                                            expression = '!Area!**0.9',
                                                            field_type = 'FLOAT')                   
                        if field == param3: #if area field was updated to set name of param3
                            arcpy.management.AlterField(in_table = x, #change this for the expression
                                                        field = param3,
                                                        new_field_name = 'Area',
                                                        new_field_alias = 'Area')
                            
                            arcpy.management.CalculateField(in_table = x, #calculate volume based on area 
                                                            field = 'EstVol',
                                                            expression = '!Area!**0.9',
                                                            field_type = 'FLOAT')
                            
                           
                            arcpy.management.AlterField(in_table = x, #change area field back to input param3 field name
                                                        field = 'Area',
                                                        new_field_name = param3,
                                                        new_field_alias = param3)
                   
            if len(param4) > 0: #if replacement volume field chosen
                    searcher = arcpy.da.SearchCursor(in_table = x, field_names = '*') #search if area field is == Area or param3
                    field_names = searcher.fields #get field names
                    for field in field_names:
                        if field == 'Area':
                            arcpy.management.CalculateField(in_table = x,
                                                            field = param4,
                                                            expression = '!Area!**0.9',
                                                            field_type = 'FLOAT')                   
                        if field == param3: #if area field was updated to set name of param3
                            arcpy.management.AlterField(in_table = x, #change this for the expression
                                                        field = param3,
                                                        new_field_name = 'Area',
                                                        new_field_alias = 'Area')
                            
                            arcpy.management.CalculateField(in_table = x, #calculate volume based on area 
                                                            field = param4,
                                                            expression = '!Area!**0.9',
                                                            field_type = 'FLOAT')
                            
                           
                            arcpy.management.AlterField(in_table = x, #change area field back to input param3 field name
                                                        field = 'Area',
                                                        new_field_name = param3,
                                                        new_field_alias = param3)

# ---------------------------------------------------Calculate geometry and add to table (group layer input)-----------------------------------------------------------------------------------------------------###
    
    if len(param0) == 0: #param0 not chosen, group layer input
        for x in param1:#multiple inputs allowed, do them all
            input_layers = x.listLayers() #group layer as input, so I need the layers in the group layer       
            ##--------------------------------------------------------------------------------Area-------------------------------------------------------------------------##

            if len(param3) == 0: #update area field not chosen                   
                for a in input_layers:
                    if a.isFeatureLayer== True:
                        if (arcpy.Describe(a).shapeType) == 'Polygon': #this avoids trying to calculate geometry if there are points in the group layer (ex. pour points for other analysis)
                                
                             arcpy.management.AddGeometryAttributes(Input_Features = a,
                                                                   Geometry_Properties = 'AREA_GEODESIC',
                                                                   Area_Unit = param2) #Add geometry attributes to table

                             arcpy.management.AlterField(in_table = a,
                                                        field = 'AREA_GEO',
                                                        new_field_name = 'Area',
                                                        new_field_alias = 'Area') #change field name to Area from AREA_GEO
                    

            if len(param3) > 0: #update area field chosen
                for a in input_layers:
                    if a.isFeatureLayer== True:
                        if (arcpy.Describe(a).shapeType) == 'Polygon': #this avoids trying to calculate geometry if there are points in the group layer (ex. pour points for other analysis)
                            
                            arcpy.management.DeleteField(in_table = a, drop_field = param3) #delete the field since I don't know data type

                            arcpy.management.AddGeometryAttributes(Input_Features = a,
                                                                   Geometry_Properties = 'AREA_GEODESIC',
                                                                   Area_Unit = param2)
                            arcpy.management.AlterField(in_table = a,
                                                        field = 'AREA_GEO',
                                                        new_field_name = param3,
                                                        new_field_alias = param3)
                 ##-----------------------------------------------------------------------Volume------------------------------------------------------------------------------##               
            if len(param4) == 0: #no replacement vol field chosen
                for a in input_layers:
                    if a.isFeatureLayer== True:
                        if (arcpy.Describe(a).shapeType) == 'Polygon': #this avoids trying to calculate geometry if there are points in the group layer (ex. pour points for other analysis)
                            arcpy.management.AddField(in_table = a, #make new volume field
                                                      field_name = 'EstVol',
                                                      field_type = 'FLOAT')
                
                            searcher = arcpy.da.SearchCursor(in_table = a, field_names = '*') #search if area field is == Area or param3
                            field_names = searcher.fields #get field names
                            for field in field_names:
                                if field == 'Area':
                                    arcpy.management.CalculateField(in_table = a,
                                                                    field = 'EstVol',
                                                                    expression = '!Area!**0.9')                   
                                if field == param3: #if area field was updated to set name of param3
                                    arcpy.management.AlterField(in_table = a, #change this for the expression
                                                                field = param3,
                                                                new_field_name = 'Area',
                                                                new_field_alias = 'Area')
                                    
                                    arcpy.management.CalculateField(in_table = a, #calculate volume based on area 
                                                                    field = 'EstVol',
                                                                    expression = '!Area!**0.9',
                                                                    field_type = 'FLOAT')
                                    
                                   
                                    arcpy.management.AlterField(in_table = a, #change area field back to input param3 field name
                                                                field = 'Area',
                                                                new_field_name = param3,
                                                                new_field_alias = param3)
                
            if len(param4) > 0: #if param 3 chosen and replacement volume field chosen
                  for a in input_layers:
                    if a.isFeatureLayer== True:
                        if (arcpy.Describe(a).shapeType) == 'Polygon':#this avoids trying to calculate geometry if there are points in the group layer (ex. pour points for other analysis)
                            searcher = arcpy.da.SearchCursor(in_table = a, field_names = '*') #search if area field is == Area or param3
                            field_names = searcher.fields #get field names
                            for field in field_names:
                                if field == 'Area':
                                    arcpy.management.CalculateField(in_table = a,
                                                                    field = param4,
                                                                    expression = '!Area!**0.9',
                                                                    field_type = 'FLOAT')                   

                                if field == param3: #if area field was updated to set name of param3
                                    arcpy.management.AlterField(in_table = a, #change this for the expression
                                                                field = param3,
                                                                new_field_name = 'Area',
                                                                new_field_alias = 'Area')
                                    
                                    arcpy.management.CalculateField(in_table = a, #calculate volume based on area 
                                                                    field = param4,
                                                                    expression = '!Area!**0.9',
                                                                    field_type = 'FLOAT')
                                    
                                   
                                    arcpy.management.AlterField(in_table = a, #change area field back to input param3 field name
                                                                field = 'Area',
                                                                new_field_name = param3,
                                                                new_field_alias = param3)
    return

# This is used to execute code if the file was run but not imported
if __name__ == '__main__':
    
    # Tool parameter accessed with GetParameter or GetParameterAsText
    param0 = arcpy.GetParameter(0)
    param1 = arcpy.GetParameter(1)
    param2 =arcpy.GetParameterAsText(2)
    param3 =arcpy.GetParameterAsText(3)
    param4= arcpy.GetParameterAsText(4)
    
    
    AreaVol(param0, param1, param2, param3, param4)
