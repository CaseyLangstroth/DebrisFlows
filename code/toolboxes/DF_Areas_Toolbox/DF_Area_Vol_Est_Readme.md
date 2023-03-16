### Tool Overview:

The purpose of this tool is for producing area and estimated volume data for debris flow polygons. After delineating a debris flow, the polygon(s) area will be calculated with ArcGIS Add Geometry tool box. Area will be calculated as geodesic and added to the attribute table. The estimated volume tool is a simple log10 math equation derived from Wall et al. 2022, which estimates the volume based on the plan area of the debris flow. There is an option to output a .csv file for the new attribute table, if this data is desired outside of ArcGIS for further analysis. 

---

### Inputs:

Inputs for this tool are feature layers. 

##### **Working outside a group layer:**

User has the option to input feature layers outside of a group layer. This is useful if the user has updated the shape of a polygon and wants to simply rerun the area and volume estimates on the edited polygon. The tool allows for multiple feature layer inputs to be run at once.

##### **Working inside a group layer:**

If the user has their contents pane organized by group layers, the tool allows for group layers to be inputted. A typical content organization may look like:

                      "Fire
                        Year
                          Location
                            Debris Flow Sites"
In any case, when choosing to input a group layer, the layer selected from the contents must be the lowest tier before the Debris Flow site polygons. In the example above, the inputted group layer would end at the Location group layer ("Fire/Year/Location"). The tool allows for multiple group layer inputs to be run at once.

While the required input data type is Feature Layer, the tool will only search for layers in which the shape type is Polygon. Therefore, if there are other shape types within the layer (ex. Pour Points), the tool will not return an error and instead will skip these data types when calculating an area.



##### **Area Units:**

The user must specify what units to calculate area via the drop down arrow. The same units will be used for the volume calculations.

##### **Area and Volume fields to be populated:**

The user has the ability to input Area and Volume fields which they would like to be populated. The tool allows for many options with these fields.

  The defaults for the fields which will be added to the attribute table are "Area" and "EstVol". If these fields are left blank when assigning parameters in the GUI, the default fields will be added to the attribute tables. However, if there are already fields titled "Area" or "EstVol", the tool will return an error stating these fields already exist.
    
  If there are fields already in the attribute table, the user can specify for these fields to be updated with the new data that will output by the tool. To do so, the user will type the fields into each of these boxes. The tool will then update the existing fields with the area and volume data calculated in the tool. It is best practice to have all fields with the same label in each input, as this parameter assumes the fields inputted exist in each attribute table.
    
  If the user would like to add an Area and EstVol field to their attribute tables, but prefers different field names, the preferred names can also be input in this parameter, and the tool will create the fields accordingly.

The area field can be blank while the volume field can be populated and visa versa. 

---

### Outputs:

The tool allows the user to chose to output csv files for further analysis. The output fields will be "Name", "Area", "EstVol". In cases where the user has specified the names of the Area or Volume fields to be different, the field names in the csv will be these specified names instead of the default "Area" or "EstVol" fields.

##### **Single Output file:**
If the user is running a feature layer outside of a group, and is only running one feature layer, the output file option can be chosen. This allows the user to chose the name of the file as well as the output location. This is useful in cases where the user has simply edited the polygon shape and would like to update an existing file with the new shape area/volume. The tool has overwrite capabilities if the user choses to overwrite existing files.

##### **Multiple Output files:**
If the user is inputing multiple layers to run the tool for either group layers or feature layers outside of groups, the user must specify a directory to save these files. The file names will be the same as the feature layer name (ex. SiteName_SiteNo_Year.csv) and will be saved in the specified directory.

---

##### Citations:
 Wall, S., Murphy, B.P., Belmont, P. & Yocom, L. (2022) Predicting post-fire debris flow grain sizes and depositional volumes in the Intermountain West, United States. Earth Surface Processes and Landforms, 1– 19. Available from: https://doi.org/10.1002/esp.5480



