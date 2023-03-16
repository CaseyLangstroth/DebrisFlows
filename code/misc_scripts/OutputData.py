# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 16:52:37 2023

@author: clang
"""

import arcpy
import os
import csv
arcpy.env.overwriteOutput = True

def OutputData(param0, param1, param2, param3, param4, param5):
    
    ############################################### Get the Layers #########################################################
    DF_layers_all = []
    DF_layers = []
    aprx = arcpy.mp.ArcGISProject('CURRENT')
    
    m1 = aprx.listMaps()
    for i in m1:
        if i.name == 'DF_Areas':
            l1 = i.listLayers()
        else:
            l1 = i.listLayers()
    for i in l1:
        if param0 in i.longName:
            DF_layers_all.append(i)
    for i in DF_layers_all:
        if i.isGroupLayer == False:
            DF_layers.append(i) #actual layers with data
    
    DF_layers_groups = [] # group layers for years
    
    for i in DF_layers_all:
        if i.isGroupLayer == True:
            if i.name == param0: #don't want the top layer
                pass
            else:
                DF_layers_groups.append(i) #groups to get years
                
    ##################################################### Get the site names and years ###########################################       
    sites = []

    for i in DF_layers:
        end = i.name.find('_20')
        for y in i.name:
            if i.name.endswith('Assumed_Original'):
                pass
            else:
                x = i.name[0:end]
                sites.append(x)
    
    sitesset = list(set(sites)) #sites unique
    
    fire_year = int(param2)
    
    years1 = []
    
    for i in DF_layers_groups:
        year_start = i.name.find("20")
        year_si = year_start
        year_end = year_si + 4
        year = i.name[year_si:year_end]
        years1.append(year)
    
    yearsint = []
    for i in years1:
        yearsint.append(int(i))
    yearsint.sort()
    
    years = []
    for i in yearsint:
        years.append(str(i))
    years.insert(0,'Assumed_Original')

    ######################################################### Get field and row/site names #########################################
    fields = ['Max']
    
    t_str = []
    t_str2 = []
    for i in yearsint:
        t = i- fire_year
        t_str.append(str(t))

    for i in t_str:
        y = 't_' + i
        t_str2.append(y)
        
    for i in t_str2:
        fields.append(i)
            
    row_names_A = []        
    row_names_all = []
    row_names_ = []
    row_names_A_ = []
    
    for i in sitesset:
        row_names_all.append(i.replace('_', ' '))
        if i.endswith('A') == False:
            row_names_.append(i + '_')
        else:
            row_names_A.append(i.replace('_', ' '))
            row_names_A_.append(i + '_')
            
    channel_names = []
    for i in row_names_all:
        channel_names.append(i + ' Channel')
    
    small_fan_names = []
    for i in row_names_all:
        name_len = len(i)
        if i.endswith('A') == True:
            new_name = i[:name_len-1] + '.1' + 'A'
            small_fan_names.append(new_name)
        else:
            small_fan_names.append(i + '.1') 
            
    exposed_names = []
    exposed_names_A = []
    
    for i in row_names_all:
        name_len = len(i)
        if i.endswith('A'):
            exposed_names_A.append(i[name_len-2:] + ')')
        else:
            exposed_names.append(i[name_len-1:] + ')')
    
    inserter_fields = ['Site']
    
    for i in fields:
        inserter_fields.append(i)

    ############################################################ Create Tables ########################################################
    fire = param1

    var1 = 'EstVol'
    
    outpath = param3
    arcpy.AddMessage('creating volume tables')
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
    
    data_table1 = arcpy.management.CreateTable(outpath, data_table_name1_vol)
    data_table2 = arcpy.management.CreateTable(outpath, data_table_name2_vol)
    data_table3 = arcpy.management.CreateTable(outpath, data_table_name3_vol)
    data_table4 = arcpy.management.CreateTable(outpath, data_table_name4_vol)
    data_table5 = arcpy.management.CreateTable(outpath, data_table_name5_vol)
    data_table6 = arcpy.management.CreateTable(outpath, data_table_name6_vol)
    data_table7 = arcpy.management.CreateTable(outpath, data_table_name7_vol)
    data_table8 = arcpy.management.CreateTable(outpath, data_table_name8_vol)
    data_table9 = arcpy.management.CreateTable(outpath, data_table_name9_vol)
    data_table10 = arcpy.management.CreateTable(outpath, data_table_name10_vol)
    arcpy.management.AddField(data_table1, 'Site', 'TEXT')
    arcpy.management.AddField(data_table2, 'Site', 'TEXT')
    arcpy.management.AddField(data_table3, 'Site', 'TEXT')
    arcpy.management.AddField(data_table4, 'Site', 'TEXT')
    arcpy.management.AddField(data_table5, 'Site', 'TEXT')
    arcpy.management.AddField(data_table6, 'Site', 'TEXT')
    arcpy.management.AddField(data_table7, 'Site', 'TEXT')
    arcpy.management.AddField(data_table8, 'Site', 'TEXT')
    arcpy.management.AddField(data_table9, 'Site', 'TEXT')
    arcpy.management.AddField(data_table10, 'Site', 'TEXT')
    
    arcpy.AddMessage('adding fields:')
    for i in fields:
        arcpy.AddMessage(i)
        arcpy.management.AddField(in_table = data_table1, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table2, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table3, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table4, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table5, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table6, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table7, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table8, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table9, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        arcpy.management.AddField(in_table = data_table10, 
                                  field_name = i, 
                                  field_type = 'FLOAT'  
                                 )
        ############################################## Prep layers for search cursor, seperate by site #####################################
    DR = []
    DR_years = []
    for i in years:
        for y in sitesset:
            DR_years.append(y + '_' + i)
    for i in DF_layers:
        for y in DR_years:
            if i.name == y:
                 DR.append(i)
    lA = []
    l2 = []    

    for i in row_names_A:
        name_len = len(i)
        if i.endswith('A'):
            lA.append(i[name_len-2:])

    for i in row_names_:
        name_len = len(i)
        l2.append(i[name_len-2:])

    
    lA_AO = []
    for i in lA:
        lA_AO.append(i)
        lA_AO.append(i + '_Assumed_Original')
    lA_AO_set = list(set(lA_AO))
    l2_AO = []
    for i in l2:
        l2_AO.append(i)
        l2_AO.append(i + 'Assumed_Original')
    l2_AO_set = list(set(l2_AO))
    
    
    DF_A = [] #seperate lists of layers for A (added) and original sites
    DF = []
    for y in DR:
        end = len(y.name)
        eA = end - 7
        e1 = end - 6
        for x in lA_AO_set:
            if y.name[eA:eA+2] == x: #get the site number
                DF_A.append(y)
            if y.name.endswith(x) == True: #get the assumed original
                DF_A.append(y)
        for a in l2_AO_set:
            if y.name[e1:e1+2] == a: #get the site number
                DF.append(y)
            if y.name.endswith(a) == True: #get the assumed original
                DF.append(y)

    #################################################### Get and insert data (A sites)####################################################
    if len(lA) > 0:
        for i in row_names_A_:
            arcpy.AddMessage('getting data from:')
            data_layer = []
            gross_sum = []
            exposed_gross_sum = []
            small_fan_gross_sum = []
            small_fan_eroded_sum = []
            small_fan_exposed_sum = []
            main_fan_exposed_sum = []
            main_fan_gross_sum = []
            channel_sum = []
            eroded_sum = []
            main_fan_eroded_sum = []
            for y in DF_A:
                if i in y.name:
                    data_layer.append(y)
            for x in data_layer:
                arcpy.AddMessage(x.name)
                gross = []
                gross_flt = []
                exposed_gross = []
                exposed_gross_flt = []
                small_fan_gross = []
                small_fan_gross_flt = []
                small_fan_eroded = []
                small_fan_eroded_flt = []
                small_fan_exposed = []
                small_fan_exposed_flt = []
                main_fan_exposed = []
                main_fan_exposed_flt = []
                main_fan_gross = []
                main_fan_gross_flt = []
                channel = []
                channel_flt = []
                eroded_gross = []
                eroded_gross_flt = []
                main_fan_eroded = []
                main_fan_eroded_flt = []
                with arcpy.da.SearchCursor(in_table =x, field_names = ['Name', 'EstVol']) as searcher:
                    for row in searcher:
                       #exposed gross 
                        if row[0].startswith('exposed') == True:
                            exposed_gross.append(row[1])
                            gross.append(row[1])
                        #small fan 
                        if row[0].endswith('.1A') == True:
                            small_fan_gross.append(row[1])
                            small_fan_eroded.append(row[1])
                            gross.append(row[1])
                        #small fan gross
                        if row[0].endswith('.1A)') == True:
                            small_fan_gross.append(row[1])
                            small_fan_exposed.append(row[1])
                        #main fan exposed
                        for y in exposed_names:
                            if row[0].endswith(y) == True:
                                main_fan_exposed.append(row[1])
                        for y in exposed_names_A:
                            if row[0].endswith(y) == True:
                                main_fan_gross.append(row[1])
                                gross.append(row[1])
                        #channel
                        if row[0].endswith('Slice') == True:
                            channel.append(row[1])
                            gross.append(row[1])
                        if row[0].endswith('(channel)') == True:
                            channel.append(row[1])
                        #eroded
                        for y in row_names_all:
                            if row[0].startswith(y) == True:
                                eroded_gross.append(row[1])
                        #main fan eroded
                            if row[0].endswith(y) == True:
                                main_fan_eroded.append(row[1])
                                main_fan_gross.append(row[1])
                                gross.append(row[1])
                        #main fan gross
                        for y in exposed_names_A:
                            if row[0].endswith(y) == True:
                                main_fan_gross.append(row[1])
                
                for i in gross:
                    gross_flt.append(float(i))
                gross_sum.append(round(sum(gross_flt),1))
                for i in exposed_gross:
                    exposed_gross_flt.append(float(i))
                exposed_gross_sum.append(round(sum(exposed_gross_flt),1))
                for i in small_fan_gross:
                    small_fan_gross_flt.append(float(i))
                small_fan_gross_sum.append(round(sum(small_fan_gross_flt),1))
                for i in small_fan_eroded:
                    small_fan_eroded_flt.append(float(i))
                small_fan_eroded_sum.append(round(sum(small_fan_eroded_flt),1))
                for i in small_fan_exposed:
                    small_fan_exposed_flt.append(float(i))
                small_fan_exposed_sum.append(round(sum(small_fan_exposed_flt),1))
                for i in main_fan_exposed:
                    main_fan_exposed_flt.append(float(i))
                main_fan_exposed_sum.append(round(sum(main_fan_exposed_flt),1))
                for i in main_fan_gross:
                    main_fan_gross_flt.append(float(i))
                main_fan_gross_sum.append(round(sum(main_fan_gross_flt),1))
                for i in channel:
                    channel_flt.append(float(i))
                channel_sum.append(round(sum(channel_flt),1))
                for i in eroded_gross:
                    eroded_gross_flt.append(float(i))
                eroded_sum.append(round(sum(eroded_gross_flt),1))
                for i in main_fan_eroded:
                    main_fan_eroded_flt.append(float(i))
                main_fan_eroded_sum.append(round(sum(main_fan_eroded_flt),1))
             
            for i in row_names_A_:
                 if i in x.name:
                     name_len = len(i)
                     y = i[:name_len-1]
                     c = y.replace('_', ' ')
                     gross_sum.append(c)
                     exposed_gross_sum.append(c)
                     main_fan_eroded_sum.append(c)
                     main_fan_gross_sum.append(c)
                     eroded_sum.append(c)
                     main_fan_exposed_sum.append(c)
                     small_fan_gross_sum.append(c)
                     small_fan_exposed_sum.append(c)
                     small_fan_eroded_sum.append(c)
                     channel_sum.append(c)
        
            gross_sum.reverse()
            exposed_gross_sum.reverse()
            main_fan_eroded_sum.reverse()
            main_fan_gross_sum.reverse()
            eroded_sum.reverse()
            main_fan_exposed_sum.reverse()
            small_fan_gross_sum.reverse()
            small_fan_exposed_sum.reverse()
            small_fan_eroded_sum.reverse()
            channel_sum.reverse()
            
            ## Correct for A's not having 2020 data
            
            yr2020 = 't_' + str(2020 - fire_year)
            i2020 = fields.index(yr2020) + 1
            
            field_len = len(inserter_fields)
            
            if len(gross_sum) == field_len - 1: #should only be 1 short because 2020 is missing
                gross_sum.insert(i2020, -9999) # -9999 = NULL
            if len(exposed_gross_sum) == field_len - 1:
                exposed_gross_sum.insert(i2020, -9999)
            if len(exposed_gross_sum) == 1: # = 1 not 0 because name is already appended
                for i in range(len(fields)):
                    exposed_gross_sum.append(0) #if there's no exposed data, its 0, not NULL
            if len(main_fan_eroded_sum) == field_len - 1:
                main_fan_eroded_sum.insert(i2020, -9999)
            if len( main_fan_gross_sum) == field_len - 1:
                main_fan_gross_sum.insert(i2020, -9999)
            if len(eroded_sum) == field_len - 1:
                eroded_sum.insert(i2020, -9999)
            if len(main_fan_exposed_sum) == field_len - 1:
                main_fan_exposed_sum.insert(i2020, -9999)
            if len(main_fan_exposed_sum) == 1:
                for i in range(len(fields)):
                    main_fan_exposed_sum.append(0)    
            if len(small_fan_gross_sum) == field_len - 1:
                small_fan_gross_sum.insert(i2020, -9999)
            if len(small_fan_gross_sum) == 1:
                for i in range(len(fields)):
                    small_fan_gross_sum.append(0)
            if len(small_fan_exposed_sum) == field_len - 1:
                small_fan_exposed_sum.insert(i2020, -9999)
            if len(small_fan_exposed_sum) == 1:
                for i in range(len(fields)):
                    small_fan_exposed_sum.append(0)
            if len(small_fan_eroded_sum) == field_len - 1:
                small_fan_eroded_sum.insert(i2020, -9999)
            if len(small_fan_eroded_sum) == 1:
                for i in range(len(fields)):
                    small_fan_eroded_sum.append(0)    
            if len(channel_sum) == field_len - 1:
                channel_sum.insert(i2020, -9999)
            if len(channel_sum) == 1:
                for i in range(len(fields)):
                    channel_sum.append(0)
            arcpy.AddMessage('inserting data to table')
            with arcpy.da.InsertCursor(data_table1, inserter_fields) as inserter:
                inserter.insertRow(eroded_sum)
            with arcpy.da.InsertCursor(data_table2, inserter_fields) as inserter:
                inserter.insertRow(exposed_gross_sum)
            with arcpy.da.InsertCursor(data_table3, inserter_fields) as inserter:
                inserter.insertRow(gross_sum)
            with arcpy.da.InsertCursor(data_table4, inserter_fields) as inserter:
                inserter.insertRow(main_fan_gross_sum)
            with arcpy.da.InsertCursor(data_table5, inserter_fields) as inserter:
                inserter.insertRow(main_fan_eroded_sum)
            with arcpy.da.InsertCursor(data_table6, inserter_fields) as inserter:
                inserter.insertRow(main_fan_exposed_sum)
            with arcpy.da.InsertCursor(data_table7, inserter_fields) as inserter:
                inserter.insertRow(channel_sum)
            with arcpy.da.InsertCursor(data_table8, inserter_fields) as inserter:
                inserter.insertRow(small_fan_eroded_sum)
            with arcpy.da.InsertCursor(data_table9, inserter_fields) as inserter:
                inserter.insertRow(small_fan_exposed_sum)
            with arcpy.da.InsertCursor(data_table10, inserter_fields) as inserter:
                inserter.insertRow(small_fan_gross_sum)
                
            
    ################################################## Get and insert data for original sites ##########################################                   
    
    for i in row_names_:
        arcpy.AddMessage('getting data from:')
        data_layer = []
        gross_sum = []
        exposed_gross_sum = []
        small_fan_gross_sum = []
        small_fan_eroded_sum = []
        small_fan_exposed_sum = []
        main_fan_exposed_sum = []
        main_fan_gross_sum = []
        channel_sum = []
        eroded_sum = []
        main_fan_eroded_sum = []
        for y in DF:
            if i in y.name:
                data_layer.append(y)
        for x in data_layer:
            arcpy.AddMessage(x.name)
            gross = []
            gross_flt = []
            exposed_gross = []
            exposed_gross_flt = []
            small_fan_gross = []
            small_fan_gross_flt = []
            small_fan_eroded = []
            small_fan_eroded_flt = []
            small_fan_exposed = []
            small_fan_exposed_flt = []
            main_fan_exposed = []
            main_fan_exposed_flt = []
            main_fan_gross = []
            main_fan_gross_flt = []
            channel = []
            channel_flt = []
            eroded_gross = []
            eroded_gross_flt = []
            main_fan_eroded = []
            main_fan_eroded_flt = []
            with arcpy.da.SearchCursor(in_table =x, field_names = ['Name', 'EstVol']) as searcher:
                for row in searcher:
                   #exposed gross 
                    if row[0].startswith('exposed') == True:
                        exposed_gross.append(row[1])
                        gross.append(row[1])
                    #small fan 
                    if row[0].endswith('.1') == True:
                        small_fan_gross.append(row[1])
                        small_fan_eroded.append(row[1])
                        gross.append(row[1])
                    #small fan gross
                    if row[0].endswith('.1)') == True:
                        small_fan_gross.append(row[1])
                        small_fan_exposed.append(row[1])
                    #main fan exposed
                    for y in exposed_names:
                        if row[0].endswith(y) == True:
                            main_fan_exposed.append(row[1])
                    #channel
                    if row[0].endswith('Slice') == True:
                        channel.append(row[1])
                        gross.append(row[1])
                    #eroded
                    for y in row_names_all:
                        if row[0].startswith(y) == True:
                            eroded_gross.append(row[1])
                    #main fan eroded
                        if row[0].endswith(y) == True:
                            main_fan_eroded.append(row[1])
                            main_fan_gross.append(row[1])
                            gross.append(row[1])
                    #main fan gross
                    for y in exposed_names:
                        if row[0].endswith(y) == True:
                            main_fan_gross.append(row[1])
                    
    
            for i in gross:
                gross_flt.append(float(i))
            gross_sum.append(round(sum(gross_flt),1))
            for i in exposed_gross:
                exposed_gross_flt.append(float(i))
            exposed_gross_sum.append(round(sum(exposed_gross_flt),1))
            for i in small_fan_gross:
                small_fan_gross_flt.append(float(i))
            small_fan_gross_sum.append(round(sum(small_fan_gross_flt),1))
            for i in small_fan_eroded:
                small_fan_eroded_flt.append(float(i))
            small_fan_eroded_sum.append(round(sum(small_fan_eroded_flt),1))
            for i in small_fan_exposed:
                small_fan_exposed_flt.append(float(i))
            small_fan_exposed_sum.append(round(sum(small_fan_exposed_flt), 1))
            for i in main_fan_exposed:
                main_fan_exposed_flt.append(float(i))
            main_fan_exposed_sum.append(round(sum(main_fan_exposed_flt),1))
            for i in main_fan_gross:
                main_fan_gross_flt.append(float(i))
            main_fan_gross_sum.append(round(sum(main_fan_gross_flt),1))
            for i in channel:
                channel_flt.append(float(i))
            channel_sum.append(round(sum(channel_flt),1))
            for i in eroded_gross:
                eroded_gross_flt.append(float(i))
            eroded_sum.append(round(sum(eroded_gross_flt),1))
            for i in main_fan_eroded:
                main_fan_eroded_flt.append(float(i))
            main_fan_eroded_sum.append(round(sum(main_fan_eroded_flt),1))
    
        for i in row_names_:
            if i in x.name:
                name_len = len(i)
                y = i[:name_len-1]
                c = y.replace('_', ' ')
                gross_sum.append(c)
                exposed_gross_sum.append(c)
                main_fan_eroded_sum.append(c)
                main_fan_gross_sum.append(c)
                eroded_sum.append(c)
                main_fan_exposed_sum.append(c)
                small_fan_gross_sum.append(c)
                small_fan_exposed_sum.append(c)
                small_fan_eroded_sum.append(c)
                channel_sum.append(c)
    
        
    
        gross_sum.reverse()
        exposed_gross_sum.reverse()
        main_fan_eroded_sum.reverse()
        main_fan_gross_sum.reverse()
        eroded_sum.reverse()
        main_fan_exposed_sum.reverse()
        small_fan_gross_sum.reverse()
        small_fan_exposed_sum.reverse()
        small_fan_eroded_sum.reverse()
        channel_sum.reverse()
        
        arcpy.AddMessage(gross_sum)
        with arcpy.da.InsertCursor(data_table1, inserter_fields) as inserter:
            inserter.insertRow(eroded_sum)
        with arcpy.da.InsertCursor(data_table2, inserter_fields) as inserter:
            inserter.insertRow(exposed_gross_sum)
        with arcpy.da.InsertCursor(data_table3, inserter_fields) as inserter:
            inserter.insertRow(gross_sum)
        with arcpy.da.InsertCursor(data_table4, inserter_fields) as inserter:
            inserter.insertRow(main_fan_gross_sum)
        with arcpy.da.InsertCursor(data_table5, inserter_fields) as inserter:
            inserter.insertRow(main_fan_eroded_sum)
        with arcpy.da.InsertCursor(data_table6, inserter_fields) as inserter:
            inserter.insertRow(main_fan_exposed_sum)
        with arcpy.da.InsertCursor(data_table7, inserter_fields) as inserter:
            inserter.insertRow(channel_sum)
        with arcpy.da.InsertCursor(data_table8, inserter_fields) as inserter:
            inserter.insertRow(small_fan_eroded_sum)
        with arcpy.da.InsertCursor(data_table9, inserter_fields) as inserter:
            inserter.insertRow(small_fan_exposed_sum)
        with arcpy.da.InsertCursor(data_table10, inserter_fields) as inserter:
            inserter.insertRow(small_fan_gross_sum)
                
    ############################################################ write csv ######################################################
    arcpy.AddMessage('writing to csv')
    filepath = param4
    writer_fields = ['ID']
    for i in inserter_fields:
        writer_fields.append(i) 
    arcpy.env.workspace = param3
    table_list = arcpy.ListTables()
    for i in table_list:
        data = []
        y = arcpy.Describe(i).name
        if y.startswith(param1) == True:
            filename = y +'.csv'
            completename = (os.path.join(filepath, filename))
            with arcpy.da.SearchCursor(i, '*') as searcher:
                for row in searcher:
                        data.append(row)
            with open(completename, 'wt', newline ='') as out_file:
                csv_writer = csv.writer(out_file)
                csv_writer.writerow(writer_fields)
                csv_writer.writerows(data)
        else:
            pass
                
    ######################################## Do the same for area data if chosen to be done #############################################
    
    if param5 > 0:
        arcpy.AddMessage('creating area tables')
        var2 = 'Area'
        data_table_name1_area= fire + '_'+ 'DF_Eroded_' + var2
        data_table_name2_area = fire + '_'+ 'DF_Exposed_' + var2
        data_table_name3_area = fire + '_'+ 'DF_Gross_' + var2
        data_table_name4_area = fire + '_'+ 'DF_Main_Fan_Gross_' + var2
        data_table_name5_area = fire + '_'+ 'DF_Main_Fan_Eroded_' + var2
        data_table_name6_area = fire + '_'+ 'DF_Main_Fan_Exposed_' + var2
        data_table_name7_area = fire + '_'+ 'DF_Channel_' + var2
        data_table_name8_area = fire + '_'+ 'DF_Small_Fan_Eroded_' + var2
        data_table_name9_area = fire + '_'+ 'DF_Small_Fan_Exposed_' + var2
        data_table_name10_area = fire + '_'+ 'DF_Small_Fan_Gross_' + var2
        
        data_table1 = arcpy.management.CreateTable(outpath, data_table_name1_area)
        data_table2 = arcpy.management.CreateTable(outpath, data_table_name2_area)
        data_table3 = arcpy.management.CreateTable(outpath, data_table_name3_area)
        data_table4 = arcpy.management.CreateTable(outpath, data_table_name4_area)
        data_table5 = arcpy.management.CreateTable(outpath, data_table_name5_area)
        data_table6 = arcpy.management.CreateTable(outpath, data_table_name6_area)
        data_table7 = arcpy.management.CreateTable(outpath, data_table_name7_area)
        data_table8 = arcpy.management.CreateTable(outpath, data_table_name8_area)
        data_table9 = arcpy.management.CreateTable(outpath, data_table_name9_area)
        data_table10 = arcpy.management.CreateTable(outpath, data_table_name10_area)
        arcpy.management.AddField(data_table1, 'Site', 'TEXT')
        arcpy.management.AddField(data_table2, 'Site', 'TEXT')
        arcpy.management.AddField(data_table3, 'Site', 'TEXT')
        arcpy.management.AddField(data_table4, 'Site', 'TEXT')
        arcpy.management.AddField(data_table5, 'Site', 'TEXT')
        arcpy.management.AddField(data_table6, 'Site', 'TEXT')
        arcpy.management.AddField(data_table7, 'Site', 'TEXT')
        arcpy.management.AddField(data_table8, 'Site', 'TEXT')
        arcpy.management.AddField(data_table9, 'Site', 'TEXT')
        arcpy.management.AddField(data_table10, 'Site', 'TEXT')
        arcpy.AddMessage('adding fields:')
        for i in fields:
            arcpy.AddMessage(i)
            arcpy.management.AddField(in_table = data_table1, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table2, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table3, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table4, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table5, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table6, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table7, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table8, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table9, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
            arcpy.management.AddField(in_table = data_table10, 
                                      field_name = i, 
                                      field_type = 'FLOAT'  
                                     )
    
        for i in row_names_A_:
            arcpy.AddMessage('getting data from:')
            data_layer = []
            gross_sum = []
            exposed_gross_sum = []
            small_fan_gross_sum = []
            small_fan_eroded_sum = []
            small_fan_exposed_sum = []
            main_fan_exposed_sum = []
            main_fan_gross_sum = []
            channel_sum = []
            eroded_sum = []
            main_fan_eroded_sum = []
            for y in DF_A:
                if i in y.name:
                    data_layer.append(y)
            for x in data_layer:
                arcpy.AddMessage(x.name)
                gross = []
                gross_flt = []
                exposed_gross = []
                exposed_gross_flt = []
                small_fan_gross = []
                small_fan_gross_flt = []
                small_fan_eroded = []
                small_fan_eroded_flt = []
                small_fan_exposed = []
                small_fan_exposed_flt = []
                main_fan_exposed = []
                main_fan_exposed_flt = []
                main_fan_gross = []
                main_fan_gross_flt = []
                channel = []
                channel_flt = []
                eroded_gross = []
                eroded_gross_flt = []
                main_fan_eroded = []
                main_fan_eroded_flt = []
                with arcpy.da.SearchCursor(in_table =x, field_names = ['Name', 'Area']) as searcher:
                    for row in searcher:
                       #exposed gross 
                        if row[0].startswith('exposed') == True:
                            exposed_gross.append(row[1])
                            gross.append(row[1])
                        #small fan 
                        if row[0].endswith('.1A') == True:
                            small_fan_gross.append(row[1])
                            small_fan_eroded.append(row[1])
                            gross.append(row[1])
                        #small fan gross
                        if row[0].endswith('.1A)') == True:
                            small_fan_gross.append(row[1])
                            small_fan_exposed.append(row[1])
                        #main fan exposed
                        for y in exposed_names:
                            if row[0].endswith(y) == True:
                                main_fan_exposed.append(row[1])
                        for y in exposed_names_A:
                            if row[0].endswith(y) == True:
                                main_fan_gross.append(row[1])
                                gross.append(row[1])
                        #channel
                        if row[0].endswith('Slice') == True:
                            channel.append(row[1])
                            gross.append(row[1])
                        if row[0].endswith('(channel)') == True:
                            channel.append(row[1])
                        #eroded
                        for y in row_names_all:
                            if row[0].startswith(y) == True:
                                eroded_gross.append(row[1])
                        #main fan eroded
                            if row[0].endswith(y) == True:
                                main_fan_eroded.append(row[1])
                                main_fan_gross.append(row[1])
                                gross.append(row[1])
                        #main fan gross
                        for y in exposed_names_A:
                            if row[0].endswith(y) == True:
                                main_fan_gross.append(row[1])
                
                for i in gross:
                    gross_flt.append(float(i))
                gross_sum.append(round(sum(gross_flt),1))
                for i in exposed_gross:
                    exposed_gross_flt.append(float(i))
                exposed_gross_sum.append(round(sum(exposed_gross_flt),1))
                for i in small_fan_gross:
                    small_fan_gross_flt.append(float(i))
                small_fan_gross_sum.append(round(sum(small_fan_gross_flt),1))
                for i in small_fan_eroded:
                    small_fan_eroded_flt.append(float(i))
                small_fan_eroded_sum.append(round(sum(small_fan_eroded_flt),1))
                for i in small_fan_exposed:
                    small_fan_exposed_flt.append(float(i))
                small_fan_exposed_sum.append(round(sum(small_fan_exposed_flt),1))
                for i in main_fan_exposed:
                    main_fan_exposed_flt.append(float(i))
                main_fan_exposed_sum.append(round(sum(main_fan_exposed_flt),1))
                for i in main_fan_gross:
                    main_fan_gross_flt.append(float(i))
                main_fan_gross_sum.append(round(sum(main_fan_gross_flt),1))
                for i in channel:
                    channel_flt.append(float(i))
                channel_sum.append(round(sum(channel_flt),1))
                for i in eroded_gross:
                    eroded_gross_flt.append(float(i))
                eroded_sum.append(round(sum(eroded_gross_flt),1))
                for i in main_fan_eroded:
                    main_fan_eroded_flt.append(float(i))
                main_fan_eroded_sum.append(round(sum(main_fan_eroded_flt),1))
            
            for i in row_names_A_:
                if i in x.name:
                     name_len = len(i)
                     y = i[:name_len-1]
                     c = y.replace('_', ' ')
                     gross_sum.append(c)
                     exposed_gross_sum.append(c)
                     main_fan_eroded_sum.append(c)
                     main_fan_gross_sum.append(c)
                     eroded_sum.append(c)
                     main_fan_exposed_sum.append(c)
                     small_fan_gross_sum.append(c)
                     small_fan_exposed_sum.append(c)
                     small_fan_eroded_sum.append(c)
                     channel_sum.append(c)
            gross_sum.reverse()
            exposed_gross_sum.reverse()
            main_fan_eroded_sum.reverse()
            main_fan_gross_sum.reverse()
            eroded_sum.reverse()
            main_fan_exposed_sum.reverse()
            small_fan_gross_sum.reverse()
            small_fan_exposed_sum.reverse()
            small_fan_eroded_sum.reverse()
            channel_sum.reverse()
            
            ## Correct for A's not having 2020 data
            
            yr2020 = 't_' + str(2020 - fire_year)
            i2020 = fields.index(yr2020) + 1
            
            field_len = len(inserter_fields)
            
            if len(gross_sum) == field_len - 1: #should only be 1 short because 2020 is missing
                gross_sum.insert(i2020, -9999) # -9999 = NULL
            if len(exposed_gross_sum) == field_len - 1:
                exposed_gross_sum.insert(i2020, -9999)
            if len(exposed_gross_sum) == 1: # = 1 not 0 because name is already appended
                for i in range(len(fields)):
                    exposed_gross_sum.append(0) #if there's no exposed data, its 0, not NULL
            if len(main_fan_eroded_sum) == field_len - 1:
                main_fan_eroded_sum.insert(i2020, -9999)
            if len( main_fan_gross_sum) == field_len - 1:
                main_fan_gross_sum.insert(i2020, -9999)
            if len(eroded_sum) == field_len - 1:
                eroded_sum.insert(i2020, -9999)
            if len(main_fan_exposed_sum) == field_len - 1:
                main_fan_exposed_sum.insert(i2020, -9999)
            if len(main_fan_exposed_sum) == 1:
                for i in range(len(fields)):
                    main_fan_exposed_sum.append(0)    
            if len(small_fan_gross_sum) == field_len - 1:
                small_fan_gross_sum.insert(i2020, -9999)
            if len(small_fan_gross_sum) == 1:
                for i in range(len(fields)):
                    small_fan_gross_sum.append(0)
            if len(small_fan_exposed_sum) == field_len - 1:
                small_fan_exposed_sum.insert(i2020, -9999)
            if len(small_fan_exposed_sum) == 1:
                for i in range(len(fields)):
                    small_fan_exposed_sum.append(0)
            if len(small_fan_eroded_sum) == field_len - 1:
                small_fan_eroded_sum.insert(i2020, -9999)
            if len(small_fan_eroded_sum) == 1:
                for i in range(len(fields)):
                    small_fan_eroded_sum.append(0)    
            if len(channel_sum) == field_len - 1:
                channel_sum.insert(i2020, -9999)
            if len(channel_sum) == 1:
                for i in range(len(fields)):
                    channel_sum.append(0)
           
            with arcpy.da.InsertCursor(data_table1, inserter_fields) as inserter:
                inserter.insertRow(eroded_sum)
            with arcpy.da.InsertCursor(data_table2, inserter_fields) as inserter:
                inserter.insertRow(exposed_gross_sum)
            with arcpy.da.InsertCursor(data_table3, inserter_fields) as inserter:
                inserter.insertRow(gross_sum)
            with arcpy.da.InsertCursor(data_table4, inserter_fields) as inserter:
                inserter.insertRow(main_fan_gross_sum)
            with arcpy.da.InsertCursor(data_table5, inserter_fields) as inserter:
                inserter.insertRow(main_fan_eroded_sum)
            with arcpy.da.InsertCursor(data_table6, inserter_fields) as inserter:
                inserter.insertRow(main_fan_exposed_sum)
            with arcpy.da.InsertCursor(data_table7, inserter_fields) as inserter:
                inserter.insertRow(channel_sum)
            with arcpy.da.InsertCursor(data_table8, inserter_fields) as inserter:
                inserter.insertRow(small_fan_eroded_sum)
            with arcpy.da.InsertCursor(data_table9, inserter_fields) as inserter:
                inserter.insertRow(small_fan_exposed_sum)
            with arcpy.da.InsertCursor(data_table10, inserter_fields) as inserter:
                inserter.insertRow(small_fan_gross_sum)
                           
        for i in row_names_:
            arcpy.AddMessage('getting data from')
            data_layer = []
            gross_sum = []
            exposed_gross_sum = []
            small_fan_gross_sum = []
            small_fan_eroded_sum = []
            small_fan_exposed_sum = []
            main_fan_exposed_sum = []
            main_fan_gross_sum = []
            channel_sum = []
            eroded_sum = []
            main_fan_eroded_sum = []
            for y in DF:
                if i in y.name:
                    data_layer.append(y)
            for x in data_layer:
                arcpy.AddMessage(x.name)
                gross = []
                gross_flt = []
                exposed_gross = []
                exposed_gross_flt = []
                small_fan_gross = []
                small_fan_gross_flt = []
                small_fan_eroded = []
                small_fan_eroded_flt = []
                small_fan_exposed = []
                small_fan_exposed_flt = []
                main_fan_exposed = []
                main_fan_exposed_flt = []
                main_fan_gross = []
                main_fan_gross_flt = []
                channel = []
                channel_flt = []
                eroded_gross = []
                eroded_gross_flt = []
                main_fan_eroded = []
                main_fan_eroded_flt = []
                with arcpy.da.SearchCursor(in_table =x, field_names = ['Name', 'Area']) as searcher:
                    for row in searcher:
                       #exposed gross 
                        if row[0].startswith('exposed') == True:
                            exposed_gross.append(row[1])
                            gross.append(row[1])
                        #small fan 
                        if row[0].endswith('.1') == True:
                            small_fan_gross.append(row[1])
                            small_fan_eroded.append(row[1])
                            gross.append(row[1])
                        #small fan gross
                        if row[0].endswith('.1)') == True:
                            small_fan_gross.append(row[1])
                            small_fan_exposed.append(row[1])
                        #main fan exposed
                        for y in exposed_names:
                            if row[0].endswith(y) == True:
                                main_fan_exposed.append(row[1])
                        #channel
                        if row[0].endswith('Slice') == True:
                            channel.append(row[1])
                            gross.append(row[1])
                        #eroded
                        for y in row_names_all:
                            if row[0].startswith(y) == True:
                                eroded_gross.append(row[1])
                        #main fan eroded
                            if row[0].endswith(y) == True:
                                main_fan_eroded.append(row[1])
                                main_fan_gross.append(row[1])
                                gross.append(row[1])
                        #main fan gross
                        for y in exposed_names:
                            if row[0].endswith(y) == True:
                                main_fan_gross.append(row[1])
                        
        
                for i in gross:
                    gross_flt.append(float(i))
                gross_sum.append(round(sum(gross_flt),1))
                for i in exposed_gross:
                    exposed_gross_flt.append(float(i))
                exposed_gross_sum.append(round(sum(exposed_gross_flt),1))
                for i in small_fan_gross:
                    small_fan_gross_flt.append(float(i))
                small_fan_gross_sum.append(round(sum(small_fan_gross_flt),1))
                for i in small_fan_eroded:
                    small_fan_eroded_flt.append(float(i))
                small_fan_eroded_sum.append(round(sum(small_fan_eroded_flt),1))
                for i in small_fan_exposed:
                    small_fan_exposed_flt.append(float(i))
                small_fan_exposed_sum.append(round(sum(small_fan_exposed_flt), 1))
                for i in main_fan_exposed:
                    main_fan_exposed_flt.append(float(i))
                main_fan_exposed_sum.append(round(sum(main_fan_exposed_flt),1))
                for i in main_fan_gross:
                    main_fan_gross_flt.append(float(i))
                main_fan_gross_sum.append(round(sum(main_fan_gross_flt),1))
                for i in channel:
                    channel_flt.append(float(i))
                channel_sum.append(round(sum(channel_flt),1))
                for i in eroded_gross:
                    eroded_gross_flt.append(float(i))
                eroded_sum.append(round(sum(eroded_gross_flt),1))
                for i in main_fan_eroded:
                    main_fan_eroded_flt.append(float(i))
                main_fan_eroded_sum.append(round(sum(main_fan_eroded_flt),1))
           
            for i in row_names_:
                if i in x.name:
                    name_len = len(i)
                    y = i[:name_len-1]
                    c = y.replace('_', ' ')
                    gross_sum.append(c)
                    exposed_gross_sum.append(c)
                    main_fan_eroded_sum.append(c)
                    main_fan_gross_sum.append(c)
                    eroded_sum.append(c)
                    main_fan_exposed_sum.append(c)
                    small_fan_gross_sum.append(c)
                    small_fan_exposed_sum.append(c)
                    small_fan_eroded_sum.append(c)
                    channel_sum.append(c)
           
        
            gross_sum.reverse()
            exposed_gross_sum.reverse()
            main_fan_eroded_sum.reverse()
            main_fan_gross_sum.reverse()
            eroded_sum.reverse()
            main_fan_exposed_sum.reverse()
            small_fan_gross_sum.reverse()
            small_fan_exposed_sum.reverse()
            small_fan_eroded_sum.reverse()
            channel_sum.reverse()
            arcpy.AddMessage('inserting data to table')
            with arcpy.da.InsertCursor(data_table1, inserter_fields) as inserter:
                inserter.insertRow(eroded_sum)
            with arcpy.da.InsertCursor(data_table2, inserter_fields) as inserter:
                inserter.insertRow(exposed_gross_sum)
            with arcpy.da.InsertCursor(data_table3, inserter_fields) as inserter:
                inserter.insertRow(gross_sum)
            with arcpy.da.InsertCursor(data_table4, inserter_fields) as inserter:
                inserter.insertRow(main_fan_gross_sum)
            with arcpy.da.InsertCursor(data_table5, inserter_fields) as inserter:
                inserter.insertRow(main_fan_eroded_sum)
            with arcpy.da.InsertCursor(data_table6, inserter_fields) as inserter:
                inserter.insertRow(main_fan_exposed_sum)
            with arcpy.da.InsertCursor(data_table7, inserter_fields) as inserter:
                inserter.insertRow(channel_sum)
            with arcpy.da.InsertCursor(data_table8, inserter_fields) as inserter:
                inserter.insertRow(small_fan_eroded_sum)
            with arcpy.da.InsertCursor(data_table9, inserter_fields) as inserter:
                inserter.insertRow(small_fan_exposed_sum)
            with arcpy.da.InsertCursor(data_table10, inserter_fields) as inserter:
                inserter.insertRow(small_fan_gross_sum)
        arcpy.AddMessage('writing csv')
        filepath = param4
        writer_fields = ['ID']
        for i in inserter_fields:
            writer_fields.append(i)
        arcpy.env.workspace = param3
        table_list = arcpy.ListTables()
        for i in table_list:
            data = []
            y = arcpy.Describe(i).name
            if y.startswith(param1) == True:
                filename = y +'.csv'
                completename = os.path.join(filepath, filename)
                with arcpy.da.SearchCursor(i, '*') as searcher:
                    for row in searcher:
                            data.append(row)
                with open(completename, 'wt', newline ='') as out_file:
                        csv_writer = csv.writer(out_file)
                        csv_writer.writerow(writer_fields)
                        csv_writer.writerows(data)
            else:
                pass
    return

# This is used to execute code if the file was run but not imported
if __name__ == '__main__':

    # Tool parameter accessed with GetParameter or GetParameterAsText
    param0 = arcpy.GetParameterAsText(0)
    param1 = arcpy.GetParameterAsText(1)
    param2 = arcpy.GetParameterAsText(2)
    param3 = arcpy.GetParameterAsText(3)
    param4 = arcpy.GetParameterAsText(4)
    param5 = arcpy.GetParameter(5)
    
    OutputData(param0, param1, param2, param3, param4, param5)
    
    # Update derived parameter values using arcpy.SetParameter() or arcpy.SetParameterAsText()
