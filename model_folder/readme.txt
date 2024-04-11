This is a data set and model to simulate the hydrology of the Dorferbach, Austria.


Observed streamflow data
--------------------------

streamflow.txt
Streamflow from Dorferbach. Downloaded from the Global Runoff Data
Centre, https://www.bafg.de/GRDC/
It represents the time span from 1 jan 1990 up to 31 dec 1993

Meteorological data
---------------------

temperature.txt and precipitation.txt
Meteorological data for the catchment areas, reanalysis data from
the NCEP, downloaded from https://globalweather.tamu.edu
It represents the time span from 1 jan 1990 up to 31 dec 1993

Map data
-----------

dem.map
Resampled from the data source:
DGM_Tirol_10m_epsg31254
Downloaded from https://www.data.gv.at/katalog/dataset/land-tirol_tirolgelnde
It is the gelande model.

ldd.map
derived from dem

sample_location.map
location of measurement for streamflow (see above)

Model
-------------------

runoff.py
PCRaster Python model
