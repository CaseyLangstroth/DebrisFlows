# import libraries 
from urllib.request import urlopen
from datetime import date
import os
import arcpy

# toogle on and off using GUI interface
usegui=False


#%% user inputs
if not usegui:
    gageid='0000'
    sdate=[]# if not using if using '2010-02-27'
    edate=[]#"End_date"# [] if not using
    outdir=r'D:\test\gagedata'

## gui based inputs
#%% user inputs
if  usegui:
    inshp=arcpy.GetParameterAsText(0)
    gaugeidfield=arcpy.GetParameterAsText(1)
    
    try:
        startdate=arcpy.GetParameterAsText(2)
    except:
        startdate=[]
    try:
        enddate=arcpy.GetParameterAsText(3)
    except:
        enddate=[]
    
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

#%% Execution
# Get the gage names from the shapefile
chk_mk_dir(outdir)


if not sdate:
    sdate="1900-01-01"

if not edate:
    today = date.today()
    edate = today.strftime("%Y-%m-%d")
    
#---------Extract and Download the data to a txt file----------
getusgsdischarges(gageid,sdate,edate,outdir)


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