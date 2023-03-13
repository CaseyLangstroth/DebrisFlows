# import libraries 
from urllib.request import urlopen
from datetime import date
import os
import arcpy

# toogle on and off using GUI interface
usegui=False


#%% user inputs
if not usegui:
    inshp=r'C:\Users\Scott\Downloads\Scott\Scott\utah_gages_test.shp'
    gaugeidfield="gageid"
    startdatefield="Begin_date"#[] if not using
    enddatefield=[]#"End_date"# [] if not using
    outdir=r'D:\test\gagedata'

## gui based inputs
#%% user inputs
if  usegui:
    inshp=arcpy.GetParameterAsText(0)
    gaugeidfield=arcpy.GetParameterAsText(1)
    
    try:
        startdatefield=arcpy.GetParameterAsText(2)
    except:
        startdatefield=[]
    try:
        enddatefield=arcpy.GetParameterAsText(3)
    except:
        enddatefield=[]
    
    outdir=arcpy.GetParameterAsText(4)



#%% Functions
def getusgsdischarges(sid,sdate,edate,outdir):
    print('Downloading Site ' +sid)
    url="https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&" \
        "site_no="+sid+"&referred_module=sw&period=&begin_date="+sdate+"&end_date="+edate
    htmlfile=urlopen(url)
    htmltext=htmlfile.read()
    fid=outdir+'\\'+sid+'.txt'
    file=open(fid,"wb")
    file.write(htmltext)
    file.close
    return;
    
# Check if a directory exists if not make it    
def chk_mk_dir(input_dir):
    try:
        os.mkdir(input_dir)
    except:
        pass


def extractshpfield(shp,field):
    try:
        arr = arcpy.da.FeatureClassToNumPyArray(shp, field ,skip_nulls=True)
    except:
        errmsg='Check all attributs are in the attributed pour point shapefile'
        raise SyntaxError(errmsg)        
    return arr[field]

#%% Execution
# Get the gage names from the shapefile
chk_mk_dir(outdir)

#---------------setup fields----------------------------
gageid=extractshpfield(inshp,gaugeidfield)



if startdatefield:
    sdate=extractshpfield(inshp,startdatefield)
else:
    sdate="1900-01-01"

if enddatefield:
    edate=extractshpfield(inshp,enddatefield)
else:
    today = date.today()
    edate = today.strftime("%Y-%m-%d")
    
#---------Extract and Download the data to a txt file----------
# but what if we have a start date or end date and not the other
# how should we edit this ....
if enddatefield and startdatefield:
    for i in range(0,len(gageid)):
        getusgsdischarges(gageid[i],sdate[i],edate[i],outdir)

elif enddatefield and not startdatefield:
    for i in range(0,len(gageid)):
        getusgsdischarges(gageid[i],sdate,edate[i],outdir)
        
elif startdatefield  and not enddatefield:
    for i in range(0,len(gageid)):
        getusgsdischarges(gageid[i],sdate[i],edate,outdir)    
else:    
    for i in range(0,len(gageid)):
        getusgsdischarges(gageid[i],sdate,edate,outdir)
        
        
        
""" how to bypass arcgis
import shapefile as shp
def extractshpfield(inshp,fieldid):
    sf=shp.Reader(inshp)#open shapefile
    fields=sf.fields# extract all field information
    field_names = [field[0] for field in fields]# get just fielnd names
    res = any(ele in fieldid for ele in field_names) # check if field name exists
    
    if res==1:
        loc=field_names.index(fieldid)# find field of interested location 
    else:
        raise SyntaxError('No field exists with that name')
    
    records=sf.records()# get all records in the shapefile
    num_rec=len(records)
    
    outdata=[records[i][loc] for i in range(0,num_rec)] # extract data of interest
    return outdata;
"""