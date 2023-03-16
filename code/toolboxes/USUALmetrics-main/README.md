# USUALmetrics
Metrics Developed to Measure Reached-Based Geomorphic, Local and Upstream variables across Fluvial Networks after Wildfire

## Preprocessing

### Fluvial Network

This set of metrics can be run on fluvial networks delineated using the USUAL Watershed Tools.
https://github.com/WatershedsWildfireResearchCollaborative/USUAL (David et al., 2022)
Watersheds and pour point shapefile FID are entered into the code using an excel datasheet and are ran through the USUAL tools in line 54-72 to delineate the fluvial network and measure the valley bottom width.
https://github.com/WatershedsWildfireResearchCollaborative/USUAL
### Valley Bottom/Floodplain Delineation

The Fluvial Corridor Tool is run to delineate the valley bottom using the delineated stream network prior to calculating valley bottom width.
https://github.com/EVS-GIS/Fluvial-Corridor-Toolbox-ArcGIS 
Roux C, Alber A, Bertrand M, Vaudor L, Piégay H. 2015. “FluvialCorridor”: a ArcGIS toolbox 	package for multiscaleriverscape exploration. Geomorphology 242:29–37. 	https://doi.org/10.1016/j.geomorph.2014.04.018.
Because parameters were optimized for each stream network, the Fluvial Corridor was run manually for each individual site and valley bottoms for each site were manually set as variables in each USUALmetrics run.  

## Inputs

### Raster Layers for Calculation of Local and Upstream Variables

-DEM: The elevation raster used in this study is a 10 m DEM sourced from USGS https://gis.utah.gov/data/elevation-and-terrain/
-RDNBR: This is derived from the MTBS database for the selected wildfires of interest. https://www.mtbs.gov/
-Moderate/High: This is derived from the RDNBR layer above. Cells as classified as being burned at Moderate/High Severity (1) or not (0)
-Cover: This raster sources from the Landfire Existing Vegetation Type Layer (USEVT130). https://landfire.gov/ It was reclassified into the following categories:

-Conifer: This raster is a reclassified version of the USEVT130 layer above, where cells are designated as either conifer (1) or not (0)
-Topographic Position Index: This raster is calculated from the DEM. http://www.jennessent.com/downloads/tpi-poster-tnc_18x22.pdf

### Shapefiles

Valley Bottom: This is the preprocessed layer created using the Fluvial Corridor Tool.
Sediment Bottlenecks: This is the layer that contains our field-measured sediment bottlenecks, attributed with Cause, Type and Volume
Debris Flows: This shapefile is a compiled dataset of debris flows mapped across Utah. Volume is calculated using Wall et al., 2022 area volume relationship.
Large Woody Debris: This shapefile is a mapped dataset of large wood within the valley bottom at each study site. 

### Output  

Background, calculation and importance of each variable are commented out in the code.

### Required Software

- ArcGIS Pro 3.0 w/ Spatial Analyst and 3D Analyst license

### Accessing Toolkit

This code can be cloned or downloaded and run in any Python interface. While this tools uses ArcGIS tools, a GUI has not been developed at this time for ArcGIS Pro.

### Recommended Citation
 Arditti A., Murphy B.P., & Belmont P. (2023). Controls on Sediment Connectivity in Fluvial Networks Impacted by Wildfire Across Utah. Earth Surface Processes and Landforms

### Known Issues
ArcGIS does not always allows for tables to be overwritten. Temporary solution: In this case delete the temp file and rerun. 

### Updates/Fixes

# MOVE TO WIKI AFTER MAKING PUBLIC
# USUAL User Manual

## Inputs and Outputs For Each Tool


### Overview of Output Data Structure


#### Watershed and River Delineation

##### Inputs

| Input | Data Type | Required/Optional|
| ----------- | ----------- |----------- |
| DEM | Raster | Required |
| Pour Point| Point Shapefile|Required
