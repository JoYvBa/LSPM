# Land Surface Process Modelling personal project: Land degradation modelling

Model script and supporting files for the LSPM course personal project.

- `ErosionModel.py`     --> PCRaster Python script with the model code and resulting graphs.
- `canopy_tbl.txt`      --> ASCII text table used for giving the proper values to the parameter for canopy cover. Values based on Morgan (2004).
- `classes.map`         --> Raster map with the land cover classification of the studied region. 0 is bare soil, 1 is other vegetation (grasses and shrubbery), 2 is forest. Based on Shug et al. (2020).
- `clone.map`           --> Raster map defining the catchment.
- `dem.map`             --> Elevation map of the study site, for source, see `readme.txt`.
- `dis/`                --> Folder for discharge output from the model.
- `ground_tbl.txt`      --> ASCII text table used for giving the proper values to the parameter for ground cover. Values based on Morgan (2004).
- `intercept_tbl.txt`   --> ASCII text table used for giving the proper values to the parameter for rainfall interception. Values based on Morgan (2004).
- `ldd.map`             --> Drainage map of the area, derived from the elevation map.
- `ldis/`               --> Folder for logarithmic discharge output from the model.
- `pltheight_tbl`       --> ASCII text table used for giving the proper values to the parameter for plant height. Values based on Gelabert et al. (2020).
- `precipitation.txt`   --> Rainfall data for the catchment area, for more information and source, see `readme.txt`.
- `readme.txt`          --> text file containing additional information about model input and data source.
- `sample_location.map` --> Location of catchment outflow point.
- `snow/`               --> Folder for snow deck output from the model.
- `temperature.txt`     --> Temperature data for the catchment area, for more information and source, see `readme.txt`.
- `totd/`               --> Folder for total soil detachment per cell model output from the model.

