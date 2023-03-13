# USDA Forest Service - RMRS - Boise Aquatic Sciences Lab
# Author: David Nagel
# 322 E. Front St., Boise, ID  83702
# dnagel [at] fs.fed.us
# Last Modified: 1-25-13
## ---------------------------------------------------------------------------
## ValleyConfinementNHDPlus30m_Ver10-02.py
##
## See the following website for a full description of this program:
## http://www.fs.fed.us/rm/boise/AWAE/projects/valley_confinement.shtml
##
## ------------------------------------------------------------------------------
##
## Input data:
##
## DEM - 30 m DEM from NHDPlus dataset or other source
## Stream lines - NHDPlus flow lines shapefile
## Waterbodies - NHDPlus waterbodes shapefile
##
## Parameters:
##
## Workspace - The input working directory
## DEM - The input DEM. Should be 30 m DEM in rectangular coordinate system
## Streams - Input stream lines. Should be NHDPlus ver. 1 stream lines 
## Use waterbody - Parameter to record if waterbodies are being used
## Input waterbody - NHDPlus waterbody shapefile
## Output type -  1 = valley bottoms only  2 = valley bottom and distance
## Use ground slope - Use ground slope = 1. Do not use ground slope = 0
## Ground slope threshold  - Slope value for maximum slope included in valley bottom
## Use flood factor - Use flood depth = 1. Do not use flood depth = 0
## Flood factor - A multiplier of bankfull depth to raise flood level above channel
## Average annual precip - Precip in cm for highest point in watershed
## Max valley width - Threshold for valley width
## Drainage area - A threshold that requires valley bottom polygons to have a minimum size
## Stream length - Threshold requiring a minimum amount of stream lenght in a valley bottom
## Poly area - Threshold for minimum size of a valley bottom polygon 
## outputfile - Output file name
##
## Output:
##
## A polygon shapefile that 1) outlines flat valley bottoms only or 2) computes
## distance along the NHD stream channel from valley bottoms. Output classes:
## 0 - Valley bottom polygon
## 1-30 - Distance from valley bottom (km)
## 31 - Distance greater than 30 km
## 50 - Lakes and reservoirs
## ---------------------------------------------------------------------------

# Import modules
import sys, arcpy, os, arcgisscripting
from arcpy import env
from arcpy.sa import *
import arcpy.cartography as CA

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out any necessary licenses
arcpy.CheckOutExtension("Spatial")

# Set overwrite
arcpy.env.overwriteOutput = True

# env.workspace = "C:\\AllWorkspace\\Spatial1\\ValleyBottomPythonDevelopment\\SmallSFBoise"
# inputDEM = "a_elev"
# inputStreams = "a_NHDFlowline.shp"
# useWaterbody = 1   # 1 = waterbody shapefile exists 0 = no waterbody shapefile
# inputWaterbody = "a_Waterbody.shp"
# valleyType = 1    # 1 = valley bottoms only  2 = valley bottom and distance
# useGroundSlope = 1  # Use ground slope = 1. Do not use ground slope = 0
# groundSlope = 9    # 10 should be default - seems about right
# useFloodFactor = 1   # Use flood depth = 1. Do not use flood depth = 0
# floodFactor = 6      # 6 should be default - 1/2 contour interval
# aveAnnualPrecip = 100
# maxValleyWidth = 1000
# drainageArea = 1
# streamLength = 1
# polyArea = 10000
# outputfile = "ValleyBottom.shp"

env.workspace = sys.argv[1]
inputDEM = sys.argv[2]
inputStreams = sys.argv[3]
useWaterbody = sys.argv[4]
inputWaterbody = sys.argv[5]
valleyType = sys.argv[6]
useGroundSlope = sys.argv[7]
groundSlope = sys.argv[8]
useFloodFactor = sys.argv[9]
floodFactor = sys.argv[10]
aveAnnualPrecip = sys.argv[11]
maxValleyWidth = sys.argv[12]
drainageArea = sys.argv[13]
streamLength = sys.argv[14]
polyArea = sys.argv[15]
outputfile = sys.argv[16]

Workspace = env.workspace

# convert all types to float so arcpy works correctly
useWaterbody = float(useWaterbody)
valleyType = float(valleyType)
groundSlope = float(groundSlope)
useGroundSlope = float(useGroundSlope)
floodFactor = float(floodFactor)
useFloodFactor = float(useFloodFactor)
maxValleyWidth = float(maxValleyWidth)
drainageArea = float(drainageArea)
streamLength = float(streamLength)
polyArea = float(polyArea)
aveAnnualPrecip = float(aveAnnualPrecip)


# Check to see if temporary folder exists already, if so delete
exists = arcpy.Exists("Temp")
if exists == True:
    arcpy.Delete_management("Temp")

# Set up the temporary workspace
workspaceTemp = env.workspace + "\\Temp"

# Now create the new folder \Temp in the original working directory. This is for
# temporary data. Create the folder using the Python os module.
os.chdir(Workspace)
os.mkdir("Temp")

arcpy.AddMessage("Temporary workspace created")

# Get the DEM cell size
desc = arcpy.Describe(inputDEM)
cellSize = desc.MeanCellHeight

# Set the output extent
arcpy.env.extent = inputDEM

arcpy.AddMessage("Environments set")

try:

    x = 1
    if x == 1:

        # Copy the NHD+ input streams and waterbodies to a new file so nothing is overwritten
        # Copy into a temporary workspace.

        # copy streams to \Temp
        # from the full path name, find the shapefile name after the last backslash
        # by finding the offset, then add 1 to the offset
        offsetbackend = inputStreams.rfind("\\")
        offsetbackendplus = offsetbackend + 1
        shapefilename = inputStreams[offsetbackendplus:]
        # recreate the filename with \Temp inserted
        pathnameTemp = workspaceTemp
        pathnameTempStreams = pathnameTemp + "\\" + shapefilename
        # print
        arcpy.AddMessage("Input streams: " + inputStreams)
        arcpy.AddMessage("Output streams: " + pathnameTempStreams)
        # Do the copy to \Temp
        arcpy.Copy_management(inputStreams, pathnameTempStreams)
        arcpy.AddMessage("Copy streams done")
        arcpy.AddMessage("")

        # copy waterbody to \Temp - see comments above under streams
        # only do so if there actually is a waterbody file. If there is not one,
        # see "if Waterbody == 0" below after the workspace is changed to Temp.
        if useWaterbody == 1:
            offsetbackend = inputWaterbody.rfind("\\")
            offsetbackendplus = offsetbackend + 1
            shapefilename = inputWaterbody[offsetbackendplus:]

            pathnameTemp = workspaceTemp
            pathnameTempWaterbody = pathnameTemp + "\\" + shapefilename

            arcpy.AddMessage("Input waterbody: " + inputWaterbody)
            arcpy.AddMessage("Output waterbody: " + pathnameTempWaterbody)
            
            arcpy.Copy_management(inputWaterbody, pathnameTempWaterbody)

            arcpy.AddMessage("Copy waterbodies done")
            arcpy.AddMessage("")

        # copy DEM to \Temp - See comments above under streams
        offsetbackend = inputDEM.rfind("\\")
        offsetbackendplus = offsetbackend + 1
        dem = inputDEM[offsetbackendplus:]

        pathnameTemp = workspaceTemp
        pathnameTempDEM = pathnameTemp + "\\" + dem

        arcpy.AddMessage("Input dem: " + inputDEM)
        arcpy.AddMessage("Output dem: " + pathnameTempDEM)
        
        arcpy.Copy_management(inputDEM, pathnameTempDEM)

        arcpy.AddMessage("Copy DEM done")
        arcpy.AddMessage("")

        # after copies to \Temp are complete, change workspace to \Temp
        env.workspace = workspaceTemp

        # If no waterbody file is present, create a dummy shapefile. This will create one
        # point along a polyline, buffered out from the streams. The point is buffered to
        # create a polygon. This is the fake waterbody.
        if useWaterbody == 0:
            arcpy.Buffer_analysis(inputStreams, "tempbuffer.shp", "90", "FULL", "ROUND", "ALL")
            arcpy.FeatureToLine_management("tempbuffer.shp", "templines.shp")
            arcpy.CreateRandomPoints_management(workspaceTemp, "temppoint.shp", "templines.shp", "", "1")
            arcpy.Buffer_analysis("temppoint.shp", "fakewaterbody.shp", "30")
            inputWaterbody = "fakewaterbody.shp"
        
        arcpy.AddMessage("Starting Drainage Area Threshold Section")    

        # Copy input streams to the hard coded name mfkb_demnet, this is a legacy name.
        arcpy.Copy_management(inputStreams, "mfkb_demnet")

        # First be sure that the field CUMDRAINAG exists
        existCUMDRAINAG = 0
        fieldList = arcpy.ListFields("mfkb_demnet.shp")
        for field in fieldList:
            if field.name == "CUMDRAINAG":
                existCUMDRAINAG = 1
        if existCUMDRAINAG == 0:
                arcpy.AddMessage("The field CUMDRAINAG does not exist in the stream line shapefile, exiting program")
                sys.exit()

        # Make a new shapefile from the NHD streams that have contrib. area > X HA.
        # The shapefile is named StreamsCAThresh.shp.  This layer is used later
        # in the algorithm to select polygons that meet the threshold criteria.
        # First create the layer file
        arcpy.MakeFeatureLayer_management("mfkb_demnet.shp", "Layer", "", "")
        # Select based on drainage area threshold
        arcpy.SelectLayerByAttribute_management("Layer", "NEW_SELECTION", "CUMDRAINAG >= "+ str(drainageArea))

        # Copy the drainage area threshold streams to a new layer    
        arcpy.CopyFeatures_management("Layer", "StreamsCAThresh.shp", "", "0", "0", "0")
        arcpy.management.Delete("Layer", "FeatureLayer")

        arcpy.AddMessage("Finished Drainage Area Threshold Section")

        # ----------------- Slope, Range, and Width Processing Section -----------------------

    
        arcpy.AddMessage("Starting Slope and Width Processing Section")

        
        # Compute slope in percent using the input DEM
        outGrid = Slope(inputDEM, "PERCENT_RISE", "1")
        outGrid.save("slope")

        # Reclassify the slope file so that flat areas are == 1 and all else == 2 using
        # the user defined ground slope.
        # Create a string variable for the Reclassify command
        slopeInputs = "0 "+ str(groundSlope) +" 1; "+ str(groundSlope) +" 100000 2"
        # Reclassify the slope grid
        outGrid = Reclassify("slope", "Value", str(slopeInputs), "DATA")
        outGrid.save("slope10")
        
        # Use vector stream network to create GRID of same network. This will be used as the
        # basis for Euclidean Distance.  "slope" is used as the base grid template.
        outGrid = ExtractByMask("slope", "mfkb_demnet.shp")
        outGrid.save("strlinegrd")

        # Reclass strlinegrd to all 1s for use with Euclidean Distance
        outGrid = Reclassify("strlinegrd", "Value", "0 10000000 1", "DATA")
        outGrid.save("strline1")

        # Run Euclidean distance on strline1
        outGrid = EucDistance("strline1", "", cellSize)
        outGrid.save("eucdist")

##        # Multiply Euclidean Distance by Slope
##        # outGrid = Times("eucdist", "slope")
##        # outGrid.save("eucxslope")
##
##        # Euclidean * slope
##        # Recode Euclidean distance x slope to create a slope-distance restriction
##        # This is a course level restriction that will be further refined by the user
##        # inputs.
##        # outGrid = Reclassify("eucxslope", "Value", "0 4000 1; 4000 1000000 2", "DATA")
##        # outGrid.save("eucxslp850")

        # Replace Euclidean distance * slope procedures above with the Path Distance function
        # outPathDist = PathDistance("strline1", "slope", inputDEM)
        # outPathDist.save("pathdist")

        # Replace Euclidean distance * slope procedures above with the Path Distance function
        outCostDist = CostDistance("strline1", "slope")
        outCostDist.save("costdist")

        # Threshold the path distance grid at 2500
        outGrid = Reclassify("costdist", "Value", "0 2500 1; 2500 10000000 2", "DATA")
        outGrid.save("eucxslp850")

        # Incorporate the maximum width restriction
        maxValleyWidth = float(float(maxValleyWidth)/2)
        widthExpression = "0 "+ str(maxValleyWidth) +" 1; "+ str(maxValleyWidth) +" 100000 2"
        outGrid = Reclassify("eucdist", "Value", str(widthExpression), "DATA")
        outGrid.save("eucdist500")

        arcpy.AddMessage("Finished Slope and Width Processing Section")

        # ----------------- Flood Processing Section -----------------------

        # This section uses precipitation and contributing area to compute bankfull width
        # and then bankfull depth. The flood factor is multiplied by bankfull depth to
        # compute the flood depth above the channel.

        arcpy.AddMessage("Starting Flood Processing Section")

        # Compute the precipitation variable
        precipX = aveAnnualPrecip**0.355

        # Add new fields
        arcpy.AddField_management("StreamsCAThresh.shp", "BANKFWTMP", "FLOAT", "8", "2")
        arcpy.AddField_management("StreamsCAThresh.shp", "PRECIPX", "FLOAT", "8", "2")
        arcpy.AddField_management("StreamsCAThresh.shp", "BANKFW", "FLOAT", "8", "2")
        arcpy.AddField_management("StreamsCAThresh.shp", "BANKFD", "FLOAT", "8", "2")
        arcpy.AddField_management("StreamsCAThresh.shp", "FLOOD_DEP", "FLOAT", "8", "2")

        # Calculate new fields
        arcpy.CalculateField_management("StreamsCAThresh.shp", "BANKFWTMP", "float(0.196*!CUMDRAINAG!**0.280)", "PYTHON", "")
        arcpy.CalculateField_management("StreamsCAThresh.shp", "PRECIPX", precipX, "PYTHON", "")
        arcpy.CalculateField_management("StreamsCAThresh.shp", "BANKFW", "float(!BANKFWTMP!*!PRECIPX!)", "PYTHON", "")
        arcpy.CalculateField_management("StreamsCAThresh.shp", "BANKFD", "float(!BANKFW!**0.607*0.145)", "PYTHON", "")

        # Calculate depth above DEM to flood. This is the computed flood depth to flood above the bankfull
        # elevation. It is stored in the FLOOD_DEP field in the shapefile StreamsCAThresh.shp
        calc = "float(!BANKFD!*" + str(floodFactor) + ")"
        arcpy.CalculateField_management("StreamsCAThresh.shp", "FLOOD_DEP", calc, "PYTHON", "")

        # These next steps are preparetory to get ready for the nibble procedure.

        # First, convert the original DEM to integer format to get ready for nibble        
        # Convert float to integer
        mathStatement = "int(float(" + inputDEM + ") * 1000)"
        gp.SingleOutputMapAlgebra_sa(mathStatement, "demint1000")

        # Next covert the flood depth value from StreamCAThresh.shp to a raster line
        arcpy.env.snapRaster = inputDEM
        arcpy.env.extent = inputDEM
        arcpy.PolylineToRaster_conversion("StreamsCAThresh.shp", "FLOOD_DEP", "flood_line", "MAXIMUM_LENGTH", "NONE", inputDEM)

        # Now multiply by 1000 and convert the flood depth line to integer
        gp.SingleOutputMapAlgebra_sa("int(float(flood_line * 1000))", "fline1000")

        
        # Create the file that has the flood depth that will be nibbled out
        # fdep_nibme = flood depth nibble me
        # This is the flood depth line (integer) with the DEM (integer) underlying it.
        # This is a confusing grid but Arc requires that the grid being nibbled has some type
        # of data in the space that will be nibbled over. The important point is to have
        # the correct flood depth along the stream line. 
        mathStatement = "con(isnull(fline1000),demint1000,fline1000 + demint1000)"
        gp.SingleOutputMapAlgebra_sa(mathStatement, "fdep_nibme")


        # Here is the nibble. strline1 is the grid mask that represents the location of
        # the pixels that will be nibbled out. fdep_nibme respresents the pixels values
        # that will be nibbled out. So, fdep_nibme is the flood values to nibble and strline1
        # is the location mask to nibble.
        # This step takes the flooded thalweg data and nibbles out the thalweg elevations to the
        # edges of the study area.  This will be used to find pixels that are too high
        # above the stream channel. The output is nibstrelev.
        outGrid = Nibble("fdep_nibme", "strline1")
        outGrid.save("nibstrelev")
        # Convert back to float
        gp.SingleOutputMapAlgebra_sa("float(float(nibstrelev) / 1000)", "nibstrflt")
        # Smooth the data a bit
        outGrid = FocalStatistics("nibstrflt", "Rectangle 3 3 CELL", "MEAN", "DATA")
        outGrid.save("str_nib3x3")

        # At this point we have nibbled out the flooded thalweg elevation and smoothed
        # it a bit (str_nib3x3). Now wherever this nibbled layer intersects the original DEM is
        # the extent of the flooded valley.

        # Subtract nibbled stream channel elevation from the input DEM
        mathStatement = "(" + inputDEM + " - str_nib3x3)"
        gp.SingleOutputMapAlgebra_sa(mathStatement, "dem-nib")

        # Preserve areas that are <= the flood depth and write to floodArea
        conStatement = "con(dem-nib <= 0, 1, 0)"
        gp.SingleOutputMapAlgebra_sa(conStatement, "floodarea")

        ## Note: The entire flood processing section including computing the flood depth
        ## with precipitation and output of the final floodarea raster has been checked
        ## and double checked and is computing properly. DN 3-5-2013


        arcpy.AddMessage("Finished Flood Processing Section")        

        # Combine slope10 (slope threshold), eucdist850 (eucdist * slope), eucdist500 (valley width), and floodArea (flood depth grid)
        # This is the preliminary valley bottom grid without any clean up
        if useGroundSlope == 1 and useFloodFactor == 1:
            gp.SingleOutputMapAlgebra_sa("con(slope10 == 1 AND eucdist500 == 1 AND eucxslp850 == 1 AND floodarea == 1, 1, 2)", "rawflat", "")  
        if useGroundSlope == 0 and useFloodFactor == 0:
            gp.SingleOutputMapAlgebra_sa("con(eucdist500 == 1 AND eucxslp850 == 1, 1, 2)", "rawflat", "")
        if useGroundSlope == 1 and useFloodFactor == 0:
            gp.SingleOutputMapAlgebra_sa("con(slope10 == 1 AND eucdist500 == 1 AND eucxslp850 == 1, 1, 2)", "rawflat", "")
        if useGroundSlope == 0 and useFloodFactor == 1:
            gp.SingleOutputMapAlgebra_sa("con(eucdist500 == 1 AND eucxslp850 == 1 AND floodarea == 1, 1, 2)", "rawflat", "")

        arcpy.RefreshCatalog(Workspace)

        # ---------------- Start cleaning and filtering section -------------------------
        
       
        arcpy.AddMessage("Start Cleaning and Filtering Section")

        # Do 3x3 majority filter to generalize the rawflat layer
        outGrid = FocalStatistics("rawflat", "Rectangle 3 3 CELL", "MAJORITY", "DATA")
        outGrid.save("slp3x3_tmp")

        # Convert filtered slope layer to shapefile polygon format to eliminate polygons < certain size
        arcpy.RasterToPolygon_conversion("slp3x3_tmp", "SlopePoly_tmp.shp", "NO_SIMPLIFY", "Value")
        #Calculate polygon area for the slope shapefile.  This will be used to eliminate small polygons.
        arcpy.AddField_management("SlopePoly_tmp.shp", "SHAPE_AREA", "DOUBLE")
        arcpy.CalculateField_management("SlopePoly_tmp.shp", "SHAPE_AREA", "float(!SHAPE.AREA!)", "PYTHON")
        arcpy.MakeFeatureLayer_management("SlopePoly_tmp.shp", "Layer", "", "")

        # Select by attribute to create a selection set for the area threshold
        arcpy.SelectLayerByAttribute_management("Layer", "NEW_SELECTION", "SHAPE_AREA < "+ str(polyArea))
        # Eliminate polygons selected above
        arcpy.Eliminate_management("Layer", "Eliminate.shp", "LENGTH")
        # Make a Feature Layer from the eliminate file.  These are the larger polygons
        arcpy.MakeFeatureLayer_management("Eliminate.shp", "Layer", "", "")
        # Process: Select Layer By Attribute to select all polygons where GRIDCODE = 1
        arcpy.SelectLayerByAttribute_management("Layer", "NEW_SELECTION", "GRIDCODE = 1")

        ## Do the stream size threshold process

        # Copy selection to shapefile
        arcpy.CopyFeatures_management("Layer", "GridCode1.shp", "", "0", "0", "0")
        # Make GridCode1.shp a feature layer
        arcpy.MakeFeatureLayer_management("GridCode1.shp", "Layer", "", "")
        # Make StreamsCAThresh.shp a feature layer - these are streams with the
        # minumum magnitude threshold
        arcpy.MakeFeatureLayer_management("StreamsCAThresh.shp", "net_Layer", "", "")
        # Select polygons from GridCode1.shp (flat areas) that intersect StreamsCAThresh.shp
        # but only if the NHD stream is > the contributing area threshold
        arcpy.SelectLayerByLocation_management("Layer", "INTERSECT", "net_Layer", "", "NEW_SELECTION")

        # Copy output - These are the flat areas that intersect the NHD stream layer with
        # contributing area threshold and minimum size threshold
        arcpy.CopyFeatures_management("Layer", "FlatAreas.shp", "", "0", "0", "0")

        # Compute shape area for the flat polygons
        gp.CalculateField_management("FlatAreas.shp", "SHAPE_AREA", "float(!SHAPE.AREA!)", "PYTHON")

        # Get rid of small polyons less than .1 HA. This cleans up the flat area file and eliminates
        # polygons that are too small. This is an unnecessary legacy step that was retained.
        arcpy.MakeFeatureLayer_management("FlatAreas.shp", "Layer", "", "")
        arcpy.SelectLayerByAttribute_management("Layer", "NEW_SELECTION", "SHAPE_AREA > 1")
        arcpy.CopyFeatures_management("Layer", "FlatAreasElim.shp", "", "0", "0", "0")

        # Buffer the polygons 30 m (1 pixel).  This will widen the polygons a bit so that the NHD
        # stream lines will fall within the polygons in the many places were they run up on the hillslopes
        # rather than staying in the VB like they should.  This will allow us to more accurately compute
        # the total length of stream lines within a polygon.
        arcpy.AddMessage("Buffering with dissolve")
        distance = str(cellSize) + " Meters"
        arcpy.Buffer_analysis("FlatAreasElim.shp", "FlatAreasElimBuf.shp", distance, "FULL", "ROUND", "ALL", "")
        arcpy.MultipartToSinglepart_management("FlatAreasElimBuf.shp","FlatAreasElimBufsp.shp")
   
        
        arcpy.AddMessage("Finished Cleaning and Filtering Section")

        arcpy.RefreshCatalog(Workspace)

        # ----------------- Stream Length per Polygon Section-----------------------

        # This section deletes valley bottom polygons that don't have enough stream
        # length based on the minimum stream length threshold

        
        arcpy.AddMessage("Starting Stream Length Computation Section")

        # Add a field to the FlatAreas shapefile called POLY_ID.  This will give each polygon a
        # unique ID that can be used later.
        arcpy.AddField_management("FlatAreasElimBufsp.shp", "POLY_ID", "LONG")
        arcpy.CalculateField_management("FlatAreasElimBufsp.shp", "POLY_ID", "int(!FID! + 1)", "PYTHON", "")

        # Clip the NHD stream lines using the buffer polygons.  This is done in order to
        # compute the total length of stream line within each polygon, to be used to
        # eliminate polyons with stream length that is too small to support a fish population.

        # Erase lakes using the waterbody shapefile.        
        arcpy.Erase_analysis("mfkb_demnet.shp", inputWaterbody, "EraseWaterbody.shp")

        # Clip and calculate shape length.        
        arcpy.Clip_analysis("EraseWaterbody.shp", "FlatAreasElimBufsp.shp", "streamsclip.shp")
        arcpy.MultipartToSinglepart_management("streamsclip.shp","streamsclipsp.shp")
        arcpy.AddField_management("streamsclipsp.shp", "SHAPE_LEN", "DOUBLE")
        arcpy.CalculateField_management("streamsclipsp.shp", "SHAPE_LEN", "float(!SHAPE.LENGTH!)", "PYTHON")

        # Overlay the clipped streams with the polygons in order to compute total stream length for each
        # polygon by doing an identity.  Then compute statistics for each polygon to yield length.
        arcpy.Identity_analysis("streamsclipsp.shp", "FlatAreasElimBufsp.shp", "polyidentity.shp")
        arcpy.Statistics_analysis("polyidentity.shp", "PolyLengthSummary.dbf", "SHAPE_LEN SUM", "POLY_ID")

        # Join the length statistics to the polygon layer using the unique ID generated above.
        arcpy.MakeFeatureLayer_management("FlatAreasElimBufsp.shp", "Layer", "", "")
        arcpy.AddJoin_management("Layer", "POLY_ID", "PolyLengthSummary.dbf", "POLY_ID", "KEEP_ALL", )
        arcpy.CopyFeatures_management("Layer", "PolyJoin.shp")
        arcpy.RemoveJoin_management("Layer", "PolyLengthSummary")

        arcpy.RefreshCatalog(Workspace)

        # Delete unnecessary fields and create new field names that make sense.  Calculate attributes
        # into new fields.
       
        arcpy.AddMessage("Compute stream length")
        arcpy.AddField_management("PolyJoin.shp", "STRM_LEN", "DOUBLE")
        arcpy.CalculateField_management("PolyJoin.shp", "STRM_LEN", "float(!PolyLeng_3!)", "PYTHON", "")
        
        gp.RefreshCatalog(Workspace)

        # Eliminate the polygons with stream length that is too small.  To do so, actually select
        # the polygons that have enough stream length and copy to new shapefile.
        arcpy.MakeFeatureLayer_management("PolyJoin.shp", "Layer", "", "")
        string = "\"STRM_LEN\" > " + str(streamLength)
        arcpy.SelectLayerByAttribute_management("Layer", "NEW_SELECTION", string)
        arcpy.CopyFeatures_management("Layer", "STRMLENGT1500.shp")

        
        arcpy.AddMessage("Finished Stream Length Computation Section")

        # ---------------------- Start Flat Area Distance Section ---------------------------

        # This section creates the output where distance from valley bottom is computed along
        # with the valley bottom polygon.

        
        arcpy.AddMessage("Starting Flat Area Distance Section")

        arcpy.Copy_management("STRMLENGT1500.shp", "LargeFlat.shp")

        # Create a field called GRIDCODE and calculate to 1.  These are sources for
        # distance analysis
        arcpy.AddField_management("LargeFlat.shp", "GRIDVALUE", "SHORT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management("LargeFlat.shp", "GRIDVALUE", "1")

        # Make a shapefile that can be used as the cost grid
        arcpy.Copy_management("mfkb_demnet.shp", "StrLines.shp")
        arcpy.AddField_management("StrLines.shp", "GRIDVALUE", "SHORT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management("StrLines.shp", "GRIDVALUE", "1")

        #Also create field for waterbodies and calculate to 1
        arcpy.Copy_management(inputWaterbody, "WaterbodyCopy.shp")
        distance = str(cellSize * 3) + " Meters"
        arcpy.Buffer_analysis("WaterbodyCopy.shp", "WaterbodyBuf.shp", distance, "FULL", "ROUND", "ALL")
        arcpy.AddField_management("WaterbodyBuf.shp", "GRIDVALUE", "SHORT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management("WaterbodyBuf.shp", "GRIDVALUE", "1")

        # Set environment variables for raster processing
        arcpy.CellSize = cellSize
        arcpy.Extent = inputDEM

        # Delete the DEM to avoid schema lock problems
        # killDEM = workspaceTemp + "\\" + inputDEM
        # arcpy.Delete_management(killDEM)

        # Convert shapefiles to raster
        # This is the source grid
        arcpy.FeatureToRaster_conversion("LargeFlat.shp", "GRIDVALUE", "lrgflatlines", cellSize)
        # This is the cost grid
        arcpy.FeatureToRaster_conversion("StrLines.shp", "GRIDVALUE", "strlines", cellSize)
        # This is the waterbody raster
        arcpy.FeatureToRaster_conversion("WaterbodyBuf.shp", "GRIDVALUE", "waterbody", cellSize)

        # Do cost distance and set distance threshold using reclassify
        outGrid = CostDistance("lrgflatlines", "strlines", "", "")
        outGrid.save("distance")

        # This is the newest reclass scheme
        outGrid = Reclassify("distance", "Value", "0 50 0;50 1000 1;1000 2000 2;2000 3000 3;3000 4000 4;4000 5000 5;5000 6000 6;6000 7000 7;7000 8000 8;8000 9000 9;9000 10000 10;10000 11000 11;11000 12000 12;12000 13000 13;13000 14000 14;14000 15000 15;15000 16000 16;16000 17000 17;17000 18000 18;18000 19000 19;19000 20000 20;20000 21000 21;21000 22000 22;22000 23000 23;23000 24000 24;24000 25000 25;25000 26000 26;26000 27000 27;27000 28000 28;28000 29000 29;29000 30000 30;30000 10000000 31", "NODATA")
        outGrid.save("distreclas")

        # Do conditional statement to convert waterbody to 50 
        distreclas = "distreclas"
        waterbody = "waterbody"
        dist_water = "dist_water"
        # Process: Single Output Map Algebra...
        gp.SingleOutputMapAlgebra_sa("con(isnull(waterbody), distreclas, con(waterbody == 1 and distreclas >= 0, 50))", dist_water, "distreclas;waterbody")
        
 
        arcpy.AddMessage("Finished Flat Area Distance Section")

        ## New section for valleyType = 1
        arcpy.AddMessage("Starting valley type mapping")
        valley0 = "valley0"
        smoothTolerance = cellSize * 3
        if valleyType == 1:
            conStatement = "con(dist_water == 0, 1, setnull(dist_water))"
            gp.SingleOutputMapAlgebra_sa(conStatement, valley0)
            # Shrink the valley bottom 1 cell because it was buffered 1 cell earlier.
            outGrid = Shrink("valley0", "1", "1")
            outGrid.save("valley0sh")
            arcpy.RasterToPolygon_conversion("valley0sh", "distreclaspoly_a.shp", "NO_SIMPLIFY", "VALUE")
            # Smooth the output polygon.
            CA.SmoothPolygon("distreclaspoly_a.shp", "distreclaspoly_b.shp", "PAEK", smoothTolerance, "", "")
            #Calculate polygon area for the slope shapefile.  This will be used to eliminate small polygons.
            arcpy.AddField_management("distreclaspoly_b.shp", "SHAPE_AREA", "DOUBLE")
            arcpy.CalculateField_management("distreclaspoly_b.shp", "SHAPE_AREA", "float(!SHAPE.AREA!)", "PYTHON")
            arcpy.MakeFeatureLayer_management("distreclaspoly_b.shp", "Layer", "", "")
            # Select by attribute to create a selection set for the area threshold
            arcpy.SelectLayerByAttribute_management("Layer", "NEW_SELECTION", "SHAPE_AREA >= "+ str(polyArea))
            arcpy.CopyFeatures_management("Layer", "distreclaspoly.shp")
            arcpy.management.Delete("Layer", "FeatureLayer")

        if valleyType == 2:
            # Convert distance raster with waterbody to polygon shapefile
            arcpy.RasterToPolygon_conversion("dist_water", "distreclaspoly.shp", "NO_SIMPLIFY", "Value")

        arcpy.AddMessage("Finished valley type mapping")

        # Fix the fields        
        arcpy.AddField_management("distreclaspoly.shp", "VB_CLASS", "LONG", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management("distreclaspoly.shp", "VB_CLASS", "int(!GRIDCODE!)", "PYTHON", "")
        arcpy.DeleteField_management("distreclaspoly.shp", "GRIDCODE")

        # ----------------- Final Cleanup Section ------------------------


        arcpy.AddMessage("Starting Final Cleanup Section")

        # Set workspace
        env.workspace = workspaceTemp
        os.chdir(workspaceTemp)

        # **** Copy final output to working folder - For Toolbox version only
        string1 = workspaceTemp + "\\" + "distreclaspoly.shp"
        string2 = outputfile
        arcpy.Copy_management(string1, string2, "Shapefile")

        # **** Copy final output to working folder - For PythonWin version only
        # string1 = workspaceTemp + "\\" + "distreclaspoly.shp"
        # string2 = Workspace + "\\" + outputfile
        # arcpy.Copy_management(string1, string2, "Shapefile")   

        # Delete all raster datasets
        try:
            dsList = arcpy.ListDatasets("", "All")
            for dataset in dsList:
                arcpy.Delete_management(dataset)
        except:
            pass

        # Delete shapefiles
        try:
            fcList = arcpy.ListFeatureClasses("", "All")
            for fc in fcList:
                arcpy.Delete_management(fc)
        except:
            pass

        arcpy.Delete_management("PolyLengthSummary.dbf")

        ## Note: The \Temp folder cannot be removed due to a schema lock that cannot be removed.
        ## This is a system level problem

except:
    arcpy.AddError("An error was encountered")
    arcpy.AddMessage(arcpy.GetMessages())
 

arcpy.AddMessage("Program End")
arcpy.AddMessage("")
