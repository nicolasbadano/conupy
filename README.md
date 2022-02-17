![alt text](https://github.com/esteldunedain/conupy/blob/master/logo-sm.png?raw=true)

**conuPy** is a python program designed to create complete SWMM schematizations of complex urban catchments from easily available geographic data.

### Core features
- Convert streams, conduits and street data into SWMM models.
- Support spatially distributed rainfall.
- Post-process flood maps, both from event envelopes and individual time steps.
- Supports QGIS 2.x and ArcGIS 10.1 for geographic backends

# Installing
Currently using QGIS 2.X is recommended, given it's the development platform used at the moment. **conuPy** can be loaded as a QGIS plugin by extracting all the source files under the following folder:

On Windows:

```
C:\Users\MyUser\.qgis2\python\plugins\conuPy
```

On Linux:

```
/home/MyUser/.qgis2/python/plugins/conuPy
```

# Basic usage

## Required spatial datasets
The software requires the following geographic datasets. All geographic data should use the same **projected** coordinate system. The coordinate system should have **meters** as the unit. More flexibility may be included in the future.

### Drainage network dataset
This is a polyline shape detailing any stream or conduit present in the catchment. The data can be provided either on un-prepared or prepared conditions. The un-prepared data may be loaded once and converted to prepared condition data using conuPy. The prepared dataset may then be loaded in any GIS system and edited to improve the representation if so desired.

#### Un-prepared Drainage network
The un-prepared dataset shape only requires fields detailing the following information:
- Type (field called `Tipo`): indicates what type of drainage is the current polyline, either a `"conduit"` or a `"channel"`.
- Height (field called `Alto`): in meters.
- Width (field called `Ancho`): in meters.

Additionally, if available, the following data may be provided:
- Initial and final depths of the bottom of the stream represented by each polyline, relative to the terrain (fields called `depthIni` and `depthFin`).
- Initial and final elevation of the bottom of the stream represented by of each polyline (fields called `levelIni` and `levelFin`). If this data is provided, the depth data is ignored.

When this dataset is prepared, **conuPy** samples the terrain elevation on each node and computes the bottom level on each node from the best available data:
- If depth and elevation were not provided, the bottom level is computed by substracting the height of the stream from the terrain elevation. For conduits, a coverage depth is added (0.3m by default)
- If depth was provided, the bottom level is computed by substracting it from the terrain elevation.
- If elevation was provided, no additional processing is done for that node.

The prepared data is exported in shape format to be used in the subsequent stages of the process.

#### Prepared Drainage network
The prepared drainage network dataset contains:
- Type (field called `type`): indicates what type of drainage is the current polyline, either a `"conduit"` or a `"channel"`.
- Height (field called `h`): in meters.
- Width (field called `w`): in meters.
- Initial and final elevation of the bottom of the stream represented by of each polyline (fields called `levelIni` and `levelFin`). If this data is provided, the depth data is ignored.
- Slope of the stream (field called `slope`). This value is reported as positive when the levels decrease in the direction of the polyline.

### Street dataset
This dataset represents the street network, as a polyline shape. It has to be edited to be topologically correct: Each street block needs to connect as best as possible on corners (node) with the rest of the incoming streets; However, a certain tolerance is allowed (5m). The street dataset may be edited to remove small features that may not be important enough to be included in the hydrologic model. Each street needs the following data:
- Width (field called `Ancho`): in meters.

### Catchment dataset
This dataset represents the entire catchment to simulate, as a polygon shape with a single feature. All the area contained by the polygon will be included in the model divided in subcatchments.

### Outfall nodes dataset
This dataset represents the exit point/s of the catchment, where water is allowed to leave the model; it should be a point shape. Each point should be located in close proximity (< 50m) to the end node of the stream or street that connects to them.

### Digital elevation model
The elevation model should be provided as a raster in GRID format, covering the entirety of the catchment.

### General Terrain slopes
The general terrain slope should be provided as a raster in GRID format, covering the entirety of the catchment. This data represents an average value of the subcatchments slope, that is sampled by **conuPy** on each subcatchment center. It should therefore not be computed as a pixel slope.

### Rainfall coefficient
This dataset should be provided as a raster in GRID format. It represents a map of coefficients (typically betwen 0 and 1) that are multiplied by the main rainfall time series to simulate the spatial distribution of the storm.

### Impermeability
This dataset should be provided as a raster in GRID format. It represents the coefficient of impermeability on each zone; it is sampled by **conuPy** on each subcatchment center. The raster should be in the 0 to 100 range, with 100 representng full impermeability.