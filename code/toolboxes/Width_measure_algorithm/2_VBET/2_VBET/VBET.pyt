import os
import arcpy
import VBET
import NHDNetworkBuilder
from importlib import reload



class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Riparian Area Condition Assessments"
        self.alias = "Riparian Area Condition Assessments"

        # List of tool classes associated with this toolbox
        self.tools = [VBETtool, NHDNetworkBuildertool]


class VBETtool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Valley Bottom Extraction Tool"
        self.description = "Uses a DEM and stream network to extract a valley bottom polygon"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            displayName="Input DEM",
            name="inDEM",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input")

        param1 = arcpy.Parameter(
            displayName="Input Stream Network",
            name="inNetwork",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input")
        param1.filter.list = ["Polyline"]

        param2 = arcpy.Parameter(
            displayName="Input Drainage Area Raster",
            name="inDA",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input")

        param3 = arcpy.Parameter(
            displayName="Valley Bottom Output",
            name="outValleyBottom",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        param4 = arcpy.Parameter(
            displayName="High Drainage Area Threshold",
            name="high_da_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param4.value = 250

        param5 = arcpy.Parameter(
            displayName="Low Drainage Area Threshold",
            name="low_da_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param5.value = 25

        param6 = arcpy.Parameter(
            displayName="Large Buffer Size",
            name="lg_buf_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        param7 = arcpy.Parameter(
            displayName="Medium Buffer Size",
            name="med_buf_size",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        param8 = arcpy.Parameter(
            displayName="Small Buffer Size",
            name="sm_buf_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        param9 = arcpy.Parameter(
            displayName="Minimum Buffer Size",
            name="min_buf_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        param10 = arcpy.Parameter(
            displayName="Large Slope Threshold",
            name="lg_slope_thresh",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        param10.value = 5

        param11 = arcpy.Parameter(
            displayName="Medium Slope Threshold",
            name="med_slope_thresh",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        param11.value = 7

        param12 = arcpy.Parameter(
            displayName="Small Slope Threshold",
            name="sm_slope_thresh",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param12.value = 12

        param13 = arcpy.Parameter(
            displayName="Scratch Workspace",
            name="scratchWS",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        param13.filter.list = ["Local Database"]
        param13.value = arcpy.env.scratchWorkspace

        param14 = arcpy.Parameter(
            displayName="Aggregation Distance",
            name="ag_distance",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param14.value = 100

        param15 = arcpy.Parameter(
            displayName="Minimum Polygon Area to Keep in Output",
            name="min_area",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param15.value = 30000

        param16 = arcpy.Parameter(
            displayName="Minimum Hole Area to Keep in Output",
            name="min_hole",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param16.value = 50000

        return [param0, param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12, param13, param14, param15, param16]

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, p, messages):
        """The source code of the tool."""
        reload(VBET)
        VBET.main(p[0].valueAsText,
                  p[1].valueAsText,
                  p[2].valueAsText,
                  p[3].valueAsText,
                  p[4].valueAsText,
                  p[5].valueAsText,
                  p[6].valueAsText,
                  p[7].valueAsText,
                  p[8].valueAsText,
                  p[9].valueAsText,
                  p[10].valueAsText,
                  p[11].valueAsText,
                  p[12].valueAsText,
                  p[13].valueAsText,
                  p[14].valueAsText,
                  p[15].valueAsText,
                  p[16].valueAsText)
        return



class NHDNetworkBuildertool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "NHD Network Builder"
        self.description = "Creates a user specified stream network using attributes from NHD hydrography data"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            displayName="Input NHD Flowline",
            name="inFlowline",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input")
        param0.filter.list = ["Polyline"]

        param1 = arcpy.Parameter(
            displayName="Input NHD Waterbody",
            name="inWaterbody",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input")
        param1.filter.list = ["Polygon"]

        param2 = arcpy.Parameter(
            displayName="Input NHD Area",
            name="inArea",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input")
        param2.filter.list = ["Polygon"]

        param3 = arcpy.Parameter(
            displayName="Check to subset artifical paths",
            name="ap_fix",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param4 = arcpy.Parameter(
            displayName="Waterbody threshold size (sq km)",
            name="subsize",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")
        param4.value = 0.001

        param5 = arcpy.Parameter(
            displayName="Remove Artifical Paths",
            name="boolArtPath",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param6 = arcpy.Parameter(
            displayName="Remove Canals",
            name="boolCanals",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param7 = arcpy.Parameter(
            displayName="Remove Aqueducts",
            name="boolAqueducts",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param8 = arcpy.Parameter(
            displayName="Remove Stormwater",
            name="boolStormwater",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param9 = arcpy.Parameter(
            displayName="Remove Connectors",
            name="boolConnectors",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param10 = arcpy.Parameter(
            displayName="Remove General Streams",
            name="boolStreams",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param11 = arcpy.Parameter(
            displayName="Remove Intermittent Streams",
            name="boolIntermittent",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param12 = arcpy.Parameter(
            displayName="Remove Perennial Streams",
            name="boolPerennial",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param13 = arcpy.Parameter(
            displayName="Remove Ephemeral Streams",
            name="boolEphemeral",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")

        param14 = arcpy.Parameter(
            displayName="Output Stream Network",
            name="outNetwork",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        param15 = arcpy.Parameter(
            displayName="Select a Projection",
            name="proj",
            datatype="GPCoordinateSystem",
            parameterType="Required",
            direction="Input")

        param16 = arcpy.Parameter(
            displayName="Scratch Workspace",
            name="scratchWS",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
        param16.filter.list = ["File System"]

        return [param0, param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12, param13, param14, param15, param16]

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, p, messages):
        """The source code of the tool."""
        reload(NHDNetworkBuilder)
        NHDNetworkBuilder.main(p[0].valueAsText,
                               p[1].valueAsText,
                               p[2].valueAsText,
                               p[3].valueAsText,
                               p[4].valueAsText,
                               p[5].valueAsText,
                               p[6].valueAsText,
                               p[7].valueAsText,
                               p[8].valueAsText,
                               p[9].valueAsText,
                               p[10].valueAsText,
                               p[11].valueAsText,
                               p[12].valueAsText,
                               p[13].valueAsText,
                               p[14].valueAsText,
                               p[15].valueAsText,
                               p[16].valueAsText)
        return

