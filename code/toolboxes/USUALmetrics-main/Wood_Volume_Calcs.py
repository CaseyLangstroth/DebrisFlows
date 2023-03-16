# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 11:31:08 2022

@author: BelmontLab
"""

# LW Volume

# import things

import arcpy
import math
arcpy.env.overwriteOutput=True
arcpy.CheckOutExtension("Spatial")


# set input variables
arcpy.env.workspace = r'C:\Users\BelmontLab\Documents\Alec\Research\Valley Bottom Setting\Valley Bottom Metrics'
lw_file=r'C:\Users\BelmontLab\Documents\Alec\Research\LW\LW Data Collection\LW Data\LW_after_fire.gdb\Brianhead\Brian_Head_2021HexagonFixed' # wood file

# Add measurement fields/calculate area and perimeter
arcpy.management.AddGeometryAttributes(lw_file, "PERIMETER_LENGTH_GEODESIC", "METERS")
arcpy.management.AddGeometryAttributes(lw_file, "AREA_GEODESIC", "METERS")
arcpy.management.AddField(lw_file, "Length", "FLOAT")
arcpy.management.AddField(lw_file, "Width", "FLOAT")
arcpy.management.AddField(lw_file, "Volume", "FLOAT")
arcpy.management.AddField(lw_file, "Orient", "FLOAT")

#Individuals
# Create list for holding diameters
diam_list=[]
# Iterate through each piece of wood (>20 cm in diameter, 1 meter in length)
updater=arcpy.da.UpdateCursor(lw_file, field_names=["Length",'Width','AREA_GEO','PERIM_GEO',"Volume"], where_clause='quantity=1')
for row in updater:
    Area=row[2]
    Perim=row[3]
    row[0]=Perim/4+1/4*(Perim**2-16*Area)**(1/2) # Calculate length
    Length=row[0]
    row[1]=Area/Length #Calculate diameter
    Diameter=row[1]
    diam_list.append(Diameter) # Add diamters to list
    row[4]=math.pi*(Diameter/2)**2*Length # Calculate volume, approximating as a hexagonal prism
    updater.updateRow(row)

# Calculate the average diameter for pieces of wood within this fire    
avg_diam=sum(diam_list)/len(diam_list)

#Jams
#Iterate through each wood jam (3 or more pieces of wood)
updater=arcpy.da.UpdateCursor(lw_file, field_names=['AREA_GEO','Covered_dcml',"Volume"], where_clause='quantity>2')
for row in updater:
    Area=row[0]    
    Porosity=row[1]
    row[2]=avg_diam*2*Area*Porosity #Calculate volume of jam using the area, porosity and estimate depth to be the 2 average pieces of wood
    updater.updateRow(row) 
    
#Errors
#Check to make sure jams have more than 2 pieces of wood
searcher=arcpy.da.SearchCursor(lw_file, field_names=["Length",'Width','AREA_GEO','PERIM_GEO',"Volume"], where_clause='Quantity=2')
for row in searcher:
    print("QUANTITY VALUE OF 2")