# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 19:04:05 2022

@author: clang
"""
import arcpy
from arcpy import env
import math
import os
import csv

aprx = arcpy.mp.ArcGISProject('CURRENT')
arcpy.env.overwriteOutput = True #allow overwrite


def AreaVolTbl(group_layer, fire_name, fire_year, outpath, param4):
    
    data_layer = []
    ##Get layers from group layer for site names
    DF_layers = aprx.listLayers()
        
    for i in DF_layers:
        if i.isFeatureLayer == True and i.name == group_layer:
            if arcpy.Describe(i).shapeType == 'Polygon':
                start = i.name.find("20")
                end = start - 1
                for y in i.name:
                    if i.name.endswith('Assumed_Original'):
                        pass
                    else:
                        x = i.name[0:end]
                        data_layer.append(x)
    #make a list of unique values
    data_layers = set(data_layer)
     #get name of fire          
    if fire_name.find(' ') == True:
        fire = fire_name.replace(' ', '_') #replace spaces with _
    else:
        fire = fire_name
    
    fire_year = int(fire_year)
    
    var1 = 'EstVol'
    var2 = 'Area'
    
    data_table_name1_vol = fire + '_'+ 'DF_Eroded_' + var1
    data_table_name2_vol = fire + '_'+ 'DF_Exposed_' + var1
    data_table_name3_vol = fire + '_'+ 'DF_Gross_' + var1
    data_table_name4_vol = fire + '_'+ 'DF_Main_Fan_Gross_' + var1
    data_table_name5_vol = fire + '_'+ 'DF_Main_Fan_Eroded_' + var1
    data_table_name6_vol = fire + '_'+ 'DF_Main_Fan_Exposed_' + var1
    data_table_name7_vol = fire + '_'+ 'DF_Channel_' + var1
    data_table_name8_vol = fire + '_'+ 'DF_Small_Fan_Eroded_' + var1
    data_table_name9_vol = fire + '_'+ 'DF_Small_Fan_Exposed_' + var1
    data_table_name10_vol = fire + '_'+ 'DF_Small_Fan_Gross_' + var1

    groups = []
    years = []

    for i in DF_layers:
        if i.name.startswith(fire):
            DF = i.listLayers()
    for i in DF:
        if i.isFeatureLayer == False:
            name = i.name
            groups.append(name)

    for i in groups:
        year_start = i.find("20")
        if year_start > 0:
            year_si = year_start
            year_end = year_si + 4 #index year in group name
            year = i[year_si:year_end]
            years.append(year)
            
    years.append('Assumed_Original')
    years.reverse() #list it in reverse
    year_len = len(years)
    fields = ['Max']
    t_str = ['0'] #string of time values
    for i in years[1:year_len]: #skip the assumed original
        i_int = int(i)
        t = i_int - fire_year
        t_str.append(str(t))
    for i in t_str:
        i = 't_' + i
        fields.append(i)
     
    #Create volume tables
    data_table1_vol = arcpy.management.CreateTable(outpath, data_table_name1_vol)
    data_table2_vol = arcpy.management.CreateTable(outpath, data_table_name2_vol)
    data_table3_vol = arcpy.management.CreateTable(outpath, data_table_name3_vol)
    data_table4_vol = arcpy.management.CreateTable(outpath, data_table_name4_vol)
    data_table5_vol = arcpy.management.CreateTable(outpath, data_table_name5_vol)
    data_table6_vol = arcpy.management.CreateTable(outpath, data_table_name6_vol)
    data_table7_vol = arcpy.management.CreateTable(outpath, data_table_name7_vol)
    data_table8_vol = arcpy.management.CreateTable(outpath, data_table_name8_vol)
    data_table9_vol = arcpy.management.CreateTable(outpath, data_table_name9_vol)
    data_table10_vol = arcpy.management.CreateTable(outpath, data_table_name10_vol)
    arcpy.management.AddField(data_table1_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table2_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table3_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table4_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table5_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table6_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table7_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table8_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table9_vol, 'Site', 'TEXT')
    arcpy.management.AddField(data_table10_vol, 'Site', 'TEXT')
    for i in fields:
        arcpy.management.AddField(in_table = data_table1_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table2_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table3_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table4_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table5_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table6_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table7_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table8_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table9_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table10_vol, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )

    data_table_name1_area = fire + '_'+ 'DF_Eroded_' + var2
    data_table_name2_area = fire + '_'+ '_DF_Exposed_' + var2
    data_table_name3_area = fire + '_'+ '_DF_Gross_' + var2
    data_table_name4_area = fire + '_'+ '_DF_Main_Fan_Gross_' + var2
    data_table_name5_area = fire + '_'+ '_DF_Main_Fan_Eroded_' + var2
    data_table_name6_area = fire + '_'+ '_DF_Main_Fan_Exposed_' + var2
    data_table_name7_area = fire + '_'+ '_DF_Channel_' + var2
    data_table_name8_area = fire + '_'+ '_DF_Small_Fan_Eroded_' + var2
    data_table_name9_area = fire + '_'+ '_DF_Small_Fan_Exposed_' + var2
    data_table_name10_area = fire + '_'+ '_DF_Small_Fan_Gross_' + var2
    data_table1_area = arcpy.management.CreateTable(outpath, data_table_name1_area)
    data_table2_area = arcpy.management.CreateTable(outpath, data_table_name2_area)
    data_table3_area = arcpy.management.CreateTable(outpath, data_table_name3_area)
    data_table4_area = arcpy.management.CreateTable(outpath, data_table_name4_area)
    data_table5_area = arcpy.management.CreateTable(outpath, data_table_name5_area)
    data_table6_area = arcpy.management.CreateTable(outpath, data_table_name6_area)
    data_table7_area = arcpy.management.CreateTable(outpath, data_table_name7_area)
    data_table8_area = arcpy.management.CreateTable(outpath, data_table_name8_area)
    data_table9_area = arcpy.management.CreateTable(outpath, data_table_name9_area)
    data_table10_area = arcpy.management.CreateTable(outpath, data_table_name10_area)
    arcpy.management.AddField(data_table1_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table2_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table3_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table4_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table5_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table6_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table7_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table8_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table9_area, 'Site', 'TEXT')
    arcpy.management.AddField(data_table10_area, 'Site', 'TEXT')
    for i in fields:
        arcpy.management.AddField(in_table = data_table1_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table2_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table3_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table4_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table5_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table6_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table7_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table8_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table9_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table10_area, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        
    data_layer = []
    for i in DF_layers:
        if i.isFeatureLayer == True:
            if arcpy.Describe(i).shapeType == 'Polygon':
                start = i.name.find("20")
                end = start - 1
                for y in i.name:
                    if i.name.endswith('Assumed_Original'):
                        pass
                    else:
                        x = i.name[0:end]
                        data_layer.append(x)
    data_layers = set(data_layer)

    for i in data_layers:
        row_name = i.replace('_',' ')
        channel_name = row_name + ' ' 'Channel'
        small_fan_name = row_name + '.1'
        start = len(i) - 1
        end = len(i)
        exposed_name = i[start:end] + ')'
        exposed_name_A = i[start:end] + 'A' +')'

    searcher_fields1 = ['Name', 'Area']
    searcher_fields2 = ['Name', 'EstVol']
    inserter_fields = ['Site']
    for i in fields:
        inserter_fields.append(i)
        

    DR = []
    gross_sum = []
    gross_sum_r = []
    eroded_sum = []
    eroded_sum_r = []
    exposed_sum =[]
    exposed_sum_r = []
    main_fan_gross_sum = []
    main_fan_gross_sum_r = []
    main_fan_eroded_sum = []
    main_fan_eroded_sum_r = []
    main_fan_exposed_sum = []
    main_fan_exposed_sum_r =[]
    small_fan_eroded_sum = []
    small_fan_eroded_sum_r =[]
    small_fan_exposed_sum = []
    small_fan_exposed_sum_r = []
    small_fan_gross_sum =[]
    small_fan_gross_sum_r = []
    channel_sum = []
    channel_sum_r = []

    for i in years:
        DR_year = data_layer + '_' + i
        for i in DF_layers:
            if i.isFeatureLayer == True:
                if i.name == DR_year:
                    DR.append(i)

    for i in DR: 
        gross = []
        gross_int = []
        eroded_gross = []
        eroded_gross_int = []
        exposed_gross =[]
        exposed_gross_int = []
        main_fan_gross = []
        main_fan_gross_int = []
        main_fan_eroded = []
        main_fan_eroded_int = []
        main_fan_exposed = []
        main_fan_exposed_int = []
        small_fan_gross = []
        small_fan_gross_int = []
        small_fan_eroded = []
        small_fan_eroded_int = []
        small_fan_exposed = []
        small_fan_exposed_int = []
        channel = []
        channel_int = []
        with arcpy.da.SearchCursor(i, searcher_fields1) as searcher:
           for row in searcher:
              #exposed gross 
               if row[0].startswith('exposed') == True:
                   exposed_gross.append(row)
                   gross.append(row)
                   #print('exposed_gross', exposed_gross)
               #small fan 
               if row[0].endswith('.1') == True:
                   small_fan_gross.append(row)
                   small_fan_eroded.append(row)
                   gross.append(row)
               if row[0].endswith('.1A') == True:
                   small_fan_gross.append(row)
                   small_fan_eroded.append(row)
                   gross.append(row)
                   #print('small_fan_eroded',small_fan_eroded)
               #small fan gross
               if row[0].endswith('.1)') == True:
                   small_fan_gross.append(row)
                   small_fan_exposed.append(row)
               if row[0].endswith('.1A') == True:
                   main_fan_exposed.append(row)
                   #print('small_fan_gross', small_fan_gross)
               #main fan exposed
               if row[0].endswith(exposed_name) == True:
                   main_fan_exposed.append(row) 
               if row[0].endswith(exposed_name_A) == True:
                   main_fan_exposed.append(row)
               #channel
               if row[0].endswith('Slice') == True:
                   channel.append(row)
                   gross.append(row)
                   #print('channel',channel)
               #eroded 
               if row[0].startswith(row_name) == True:
                   eroded_gross.append(row)
                   #print(eroded_gross)
               #main fan eroded
               if row[0].endswith(row_name) == True:
                   main_fan_eroded.append(row)
                   main_fan_gross.append(row)
                   gross.append(row)
                   #print('main_fan_eroded',main_fan_eroded)
               #main fan gross
               if row[0].endswith(exposed_name) == True:
                   main_fan_gross.append(row)
                   gross.append(row)
               if row[0].endswith(exposed_name_A) == True:
                   main_fan_gross.append(row)
                   gross.append(row)
                   #print('main_fan_gross',main_fan_gross)               
        with arcpy.da.SearchCursor(i, searcher_fields2) as searcher:
            for row in searcher:
               #exposed gross 
                if row[0].startswith('exposed') == True:
                    exposed_gross.append(row)
                    gross.append(row)
                    #print('exposed_gross', exposed_gross)
                #small fan 
                if row[0].endswith('.1') == True:
                    small_fan_gross.append(row)
                    small_fan_eroded.append(row)
                    gross.append(row)
                if row[0].endswith('.1A') == True:
                    small_fan_gross.append(row)
                    small_fan_eroded.append(row)
                    gross.append(row)
                    #print('small_fan_eroded',small_fan_eroded)
                #small fan gross
                if row[0].endswith('.1A)') == True:
                    small_fan_gross.append(row)
                    small_fan_exposed.append(row)
                if row[0].endswith('.1)') == True:
                    small_fan_gross.append(row)
                    small_fan_exposed.append(row)
                    #print('small_fan_gross', small_fan_gross)
                #main fan exposed
                if row[0].endswith(exposed_name) == True:
                    main_fan_exposed.append(row)
                if row[0].endswith(exposed_name_A) == True:
                    main_fan_gross.append(row)
                    gross.append(row)
                #channel
                if row[0].endswith('Slice') == True:
                    channel.append(row)
                    gross.append(row)
                    #print('channel',channel)
                #eroded 
                if row[0].startswith(row_name) == True:
                    eroded_gross.append(row)
                    #print(eroded_gross)
                #main fan eroded
                if row[0].endswith(row_name) == True:
                    main_fan_eroded.append(row)
                    main_fan_gross.append(row)
                    gross.append(row)
                    #print('main_fan_eroded',main_fan_eroded)
                #main fan gross
                if row[0].endswith(exposed_name) == True:
                    main_fan_gross.append(row)
                    gross.append(row)
                if row[0].endswith(exposed_name_A) == True:
                    main_fan_gross.append(row)
                    gross.append(row)
                    #print('main_fan_gross',main_fan_gross)
        for x in exposed_gross:
            exposed_gross_int.append(x[1])
        exposed_sum.append(sum(exposed_gross_int))
        for x in eroded_gross:
            eroded_gross_int.append(x[1])
        eroded_sum.append(sum(eroded_gross_int))
        for x in channel:
            channel_int.append(x[1])
        channel_sum.append(sum(channel_int))
        for x in small_fan_gross:
            small_fan_gross_int.append(x[1])
        small_fan_gross_sum.append(sum(small_fan_gross_int))
        for x in small_fan_exposed:
            small_fan_exposed_int.append(x[1])
        small_fan_exposed_sum.append(sum(small_fan_exposed_int))
        for x in main_fan_gross:
            main_fan_gross_int.append(x[1])
        main_fan_gross_sum.append(sum(main_fan_gross_int))
        for x in main_fan_exposed:
            main_fan_exposed_int.append(x[1])
        main_fan_exposed_sum.append(sum(main_fan_exposed_int))
        for x in gross:
            gross_int.append(x[1])
        gross_sum.append(sum(gross_int))

        for x in main_fan_eroded:
            main_fan_eroded_int.append(x[1])
        main_fan_eroded_sum.append(sum(main_fan_eroded_int))
        for x in small_fan_eroded:
            small_fan_eroded_int.append(x[1])
        small_fan_eroded_sum.append(sum(small_fan_eroded_int))

    for x in exposed_sum:
        exposed_sum_r.append(round(x, 1))
    for x in eroded_sum:
        eroded_sum_r.append(round(x,1))
    for x in channel_sum:
        channel_sum_r.append(round(x,1))
    for x in small_fan_gross_sum:
        small_fan_gross_sum_r.append(round(x,1))
    for x in main_fan_gross_sum:
        main_fan_gross_sum_r.append(round(x,1))
    for x in gross_sum:
        gross_sum_r.append(round(x,1))
    for x in small_fan_eroded_sum:
        small_fan_eroded_sum_r.append(round(x,1))
    for x in main_fan_eroded_sum:
        main_fan_eroded_sum_r.append(round(x,1))
    for x in main_fan_exposed_sum:
        main_fan_exposed_sum_r.append(round(x,1))
    for x in small_fan_exposed_sum:
        small_fan_exposed_sum_r.append(round(x,1))

    exposed_dat = list(exposed_sum_r)
    exposed_dat.insert(0,row_name)
    eroded_dat = list(eroded_sum_r)
    eroded_dat.insert(0,row_name)
    channel_dat =list(channel_sum_r)
    channel_dat.insert(0,channel_name)
    small_fan_gross_dat = list(small_fan_gross_sum_r)
    small_fan_gross_dat.insert(0,small_fan_name)
    main_fan_gross_dat = list(main_fan_gross_sum_r)
    main_fan_gross_dat.insert(0, row_name)
    main_fan_eroded_dat = list(main_fan_eroded_sum_r)
    main_fan_eroded_dat.insert(0,row_name)
    small_fan_eroded_dat = list(small_fan_eroded_sum_r)
    small_fan_eroded_dat.insert(0,small_fan_name)
    small_fan_exposed_dat = list(small_fan_exposed_sum_r)
    small_fan_exposed_dat.insert(0,small_fan_name)
    main_fan_exposed_dat = list(main_fan_exposed_sum_r)
    main_fan_exposed_dat.insert(0, (row_name + ' '+ 'exposed'))
    gross_dat = list(gross_sum_r)
    gross_dat.insert(0, row_name)


    with arcpy.da.InsertCursor(data_table1_vol, inserter_fields) as inserter:
        inserter.insertRow(eroded_dat)
    with arcpy.da.InsertCursor(data_table2_vol, inserter_fields) as inserter:
        inserter.insertRow(exposed_dat)
    with arcpy.da.InsertCursor(data_table3_vol, inserter_fields) as inserter:
        inserter.insertRow(gross_dat)
    with arcpy.da.InsertCursor(data_table4_vol, inserter_fields) as inserter:
        inserter.insertRow(main_fan_gross_dat)
    with arcpy.da.InsertCursor(data_table5_vol, inserter_fields) as inserter:
        inserter.insertRow(main_fan_eroded_dat)
    with arcpy.da.InsertCursor(data_table6_vol, inserter_fields) as inserter:
        inserter.insertRow(main_fan_exposed_dat)
    with arcpy.da.InsertCursor(data_table7_vol, inserter_fields) as inserter:
        inserter.insertRow(channel_dat)
    with arcpy.da.InsertCursor(data_table8_vol, inserter_fields) as inserter:
        inserter.insertRow(small_fan_eroded_dat)
    with arcpy.da.InsertCursor(data_table9_vol, inserter_fields) as inserter:
        inserter.insertRow(small_fan_exposed_dat)
    with arcpy.da.InsertCursor(data_table10_vol, inserter_fields) as inserter:
        inserter.insertRow(small_fan_gross_dat)

    with arcpy.da.InsertCursor(data_table1_area, inserter_fields) as inserter:
        inserter.insertRow(eroded_dat)
    with arcpy.da.InsertCursor(data_table2_area, inserter_fields) as inserter:
        inserter.insertRow(exposed_dat)
    with arcpy.da.InsertCursor(data_table3_area, inserter_fields) as inserter:
        inserter.insertRow(gross_dat)
    with arcpy.da.InsertCursor(data_table4_area, inserter_fields) as inserter:
        inserter.insertRow(main_fan_gross_dat)
    with arcpy.da.InsertCursor(data_table5_area, inserter_fields) as inserter:
        inserter.insertRow(main_fan_eroded_dat)
    with arcpy.da.InsertCursor(data_table6_area, inserter_fields) as inserter:
        inserter.insertRow(main_fan_exposed_dat)
    with arcpy.da.InsertCursor(data_table7_area, inserter_fields) as inserter:
        inserter.insertRow(channel_dat)
    with arcpy.da.InsertCursor(data_table8_area, inserter_fields) as inserter:
        inserter.insertRow(small_fan_eroded_dat)
    with arcpy.da.InsertCursor(data_table9_area, inserter_fields) as inserter:
        inserter.insertRow(small_fan_exposed_dat)
    with arcpy.da.InsertCursor(data_table10_area, inserter_fields) as inserter:
        inserter.insertRow(small_fan_gross_dat)
        
    fields_write = ['ID']
    for i in fields:
        fields_write.append(i)

    table_list = arcpy.ListTables()
    for i in table_list:
        data = []
        y = arcpy.Describe(i).name
        filename = y +'.csv'
        completename = os.path.join(param4, filename)
        with arcpy.da.SearchCursor(i, '*') as searcher:
            for row in searcher:
                data.append(row)
        with open(completename, 'wt', newline ='') as out_file:
            csv_writer = csv.writer(out_file)
            csv_writer.writerow(fields_write)
            csv_writer.writerows(data)
            
            
    return

     # This is used to execute code if the file was run but not imported
if __name__ == '__main__':
    
    # Tool parameter accessed with GetParameter or GetParameterAsText
    group_layer = arcpy.GetParameter(0)
    fire_name = arcpy.GetParameterAsText(1)
    fire_year = arcpy.GetParameterAsText(2)
    outpath = arcpy.GetParameterAsText(3)
    param4 = arcpy.GetParameterAsText(4)
    
    
    
    AreaVolTbl(group_layer, fire_name, fire_year, outpath, param4)
