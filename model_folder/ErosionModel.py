from pcraster import *
from pcraster.framework import *
import numpy as np
import matplotlib.pyplot as plt

################################################################
###  Coupled precipitation and erosion model for predicting  ###
###     erosion values in a catchment on a yearly basis      ###
################################################################

#%% Model function code

# helper function 1 to read values from a map
def getCellValue(Map, Row, Column):
    Value, Valid=cellvalue(Map, Row, Column)
    if Valid:
      return Value
    else:
      raise RuntimeError('missing value in input of getCellValue')

# helper function 2 to read values from a map
def getCellValueAtBooleanLocation(location,map):
    # map can be any type, return value always float
    valueMap=mapmaximum(ifthen(location,scalar(map)))
    value=getCellValue(valueMap,1,1)
    return value


class MyFirstModel(DynamicModel):
    def __init__(self, snow_toggle, cover_toggle):
        DynamicModel.__init__(self)
        setclone('clone.map')
        
        # assign 'external' input to the model variable
        self.snow_toggle  = snow_toggle
        self.cover_toggle = cover_toggle
        

    def initial(self):
        self.clone = self.readmap('clone')
        
        # load in the elevation map
        dem = self.readmap('dem')

        # elevation (m) of the observed meteorology, this is taken from the
        # reanalysis input data set
        elevationMeteoStation = 1180.0
        elevationAboveMeteoStation = dem - elevationMeteoStation
        temperatureLapseRate = 0.005
        
        # Correct temperture of the area for the difference in elevation 
        # compared to the meteorology station
        self.temperatureCorrection = elevationAboveMeteoStation * temperatureLapseRate
        
        # set the melt rate parameter
        self.meltRateParameter = 0.04
        
        # potential loss of water to the atmosphere (m/day)
        self.atmosphericLoss = 0.002 

        # infiltration capacity, m/day
        self.infiltrationCapacity = scalar(0.0018 * 24.0)

        # proportion of subsurface water that seeps out to surface water per day
        self.seepageProportion =  0.06 

        # amount of water in the subsurface water (m), initial value
        self.subsurfaceWater = 0.0

        # amount of upward seepage from the subsurface water (m/day), initial value
        self.upwardSeepage = 0.0
 
        # snow thickness (m), initial value
        self.snow = 0.0

        # flow network
        self.ldd = self.readmap('ldd')

        # location where streamflow is measured (and reported by this model)
        self.sampleLocation = self.readmap("sample_location")

        # initialize streamflow timeseries for directly writing to disk
        self.runoffTss = TimeoutputTimeseries("streamflow_modelled", self, self.sampleLocation, noHeader=True)
        
        # initialize streamflow timeseries as numpy array for directly writing to disk
        self.simulation = numpy.zeros(self.nrTimeSteps())
        
        # initialize total rainfall
        self.totalrain = scalar(0)
        
        # initialize total discharge
        self.totaldis = scalar(0)
        
        #-------------------------------#
        # soil erosion model parameters #
        #-------------------------------#
        
        # spatially constant variables
        self.intensity      = 10        # Erosivity of rainfall, set to temperate 10 (mm/h)
        self.detach         = 0.3       # Soil detachability, set to sand; 1.2 (g / J)
        self.cohesion       = 2         # Soil cohesion, set to sand; 2 (kPa)

        # spatially dependant variables
        
        # read land cover map and set every location outside the drainage basin to -1
        self.areacover      = self.readmap("classes")
        self.areacover      = ifthenelse(self.clone == 1, self.areacover, -1)
        
        # set model parameters to those for specific land covers using ascii text tables when running the model with land cover
        if self.cover_toggle:
            self.intercept  = lookupscalar('intercept_tbl.txt', self.areacover)
            self.canopy     = lookupscalar('canopy_tbl.txt', self.areacover)
            self.pltheight  = lookupscalar('pltheight_tbl.txt', self.areacover)
            self.ground     = lookupscalar('ground_tbl.txt', self.areacover)
            
        # otherwise, set those parameters to 0
        else:
            self.intercept  = 0         # Fraction of interception by crops or vegetation (0-1)
            self.canopy     = 0         # Fraction of canopy cover (0-1)
            self.pltheight  = 0         # Height of plant from which rain falls on ground from plants (m)
            self.ground     = 0         # Fraction of ground cover (0-1)
        
        self.slope          = atan(slope(dem)) # Slope steepness (º)
        
        # determine amount of cells in the catchment
        self.totalcells_r   = int(maptotal(ifthenelse(self.clone == 1, cellarea(), 0)) / cellarea())
        
        # initiallize output arrays
        
        self.nrOfYears     = self.nrTimeSteps() // 365 
        
        self.detach_np     = np.zeros(self.nrOfYears)
        self.bare_np       = np.zeros(self.nrOfYears)
        self.veg_np        = np.zeros(self.nrOfYears)
        self.forrest_np    = np.zeros(self.nrOfYears)
        
        self.raind_np      = np.zeros(self.nrOfYears)
        self.rund_np       = np.zeros(self.nrOfYears)
        
    def dynamic(self):
        # load in precipitation data, which is giving in mm, is converted to m
        precipitation = timeinputscalar('precipitation.txt',self.clone)/1000.0
        
        # load in temperature data
        temperatureObserved = timeinputscalar('temperature.txt',self.clone)
        
        # correct temperature for altitude of cell
        temperature = temperatureObserved - self.temperatureCorrection

        # determine if a particular cell has a temperature below 0ºC
        freezing = temperature < 0.0
        
        # if temperature below 0ºC, precipitation falls as snow, otherwise, it falls as rain
        snowFall = ifthenelse(freezing,precipitation,0.0)
        rainFall = ifthenelse(pcrnot(freezing),precipitation,0.0)
        
        # if snow forming is enabled, the snow fall is added to the total amount of snow on a cell [m]
        if snow_toggle:
            self.snow = self.snow+snowFall
            
        # otherwise, the amount of snow on a cell is always 0 and the snow fall is added to the amount of rain that has fallen, equivalent to the snow instantaneously melting once it has hit the ground
        else:
            self.snow = 0
            rainFall = rainFall + snowFall
        
        # determine the amount of melting taking place when a cell is not freezing
        potentialMelt = ifthenelse(pcrnot(freezing),temperature * self.meltRateParameter, 0)
        actualMelt = min(self.snow, potentialMelt)
        
        # update total amoung of snow with the melting
        self.snow = self.snow - actualMelt

        # sublimate first from atmospheric loss
        self.sublimation = min(self.snow,self.atmosphericLoss)
        self.snow = self.snow - self.sublimation
   
        # potential evapotranspiration from subsurface water (m/day)
        self.potential_evapotranspiration = max(self.atmosphericLoss - self.sublimation,0.0)

        # actual evapotranspiration from subsurface water (m/day)
        self.evapotranspiration = min(self.subsurfaceWater, self.potential_evapotranspiration)

        # subtract actual evapotranspiration from subsurface water
        self.subsurfaceWater = max(self.subsurfaceWater - self.evapotranspiration, 0)

        # available water on surface (m/day) and infiltration
        availableWater = actualMelt + rainFall
        infiltration = min(self.infiltrationCapacity,availableWater)
        self.runoffGenerated = availableWater - infiltration

        # streamflow in m water depth per day
        discharge = accuflux(self.ldd,self.runoffGenerated + self.upwardSeepage)

        # upward seepage (m/day) from subsurface water
        self.upwardSeepage = self.seepageProportion * self.subsurfaceWater 

        # update subsurface water
        self.subsurfaceWater = max(self.subsurfaceWater + infiltration - self.upwardSeepage, 0)

        # convert streamflow from m/day to mm/day
        self.dischargeMilimPerDay = discharge * 1000
        
        # if there is any snow cover on a cell, rainfall does not cause soil erosion and is thus set to zero
        rainFall = ifthenelse(self.snow > 0, 0, rainFall)
        
        # as snowfall does not cause splash detachment, substract the snowfall amount previously added to rainFall.
        if not snow_toggle:
            rainFall = rainFall - snowFall

        # convert rainfall to mm/day and add to total rainfall 
        self.totalrain = self.totalrain + rainFall * 1000
        
        # add current discharge to total discharge
        self.totaldis = self.totaldis + self.dischargeMilimPerDay
        
        # report every 10 days
        if self.currentTimeStep() % 10 == 0:

            self.report(self.snow,'snow/snow')
        
        # calculate erosion every year
        if self.currentTimeStep() % 365 == 0:
            
            logdis = log10(self.totaldis)
            self.report(self.totaldis, 'dis/dis')
            self.report(logdis, 'ldis/ldis')
            
            # only during the first year, calculate amount of cells for each land cover type
            if self.currentTimeStep() == 365:
                # high discharges are removed for more realstic representation of the runoff erosion
                self.not_river = ifthenelse(self.totaldis < 100000, boolean(1), boolean(0))

                # calculate the amount of cells for each cover type 
                self.barecells     = int(maptotal(ifthenelse(pcrand(self.areacover == 0, self.not_river == 1), cellarea(), 0)) / cellarea())
                self.vegcells      = int(maptotal(ifthenelse(pcrand(self.areacover == 1, self.not_river == 1), cellarea(), 0)) / cellarea())
                self.forrestcells  = int(maptotal(ifthenelse(pcrand(self.areacover == 2, self.not_river == 1), cellarea(), 0)) / cellarea())
                self.totalcells    = int(maptotal(ifthenelse(self.not_river == 1, cellarea(), 0)) / cellarea())
            
            # take only the discharge on locations not determined to be a river
            self.totaldis = ifthenelse(self.not_river == 1, self.totaldis, 0)
                
            # calculate effective rainfall (mm)
            effective_rain        = self.totalrain  * (1 - self.intercept)
            
            # calculate leaf drainage by multiplying with the relative canopy cover (mm)
            leaf_drainage         = effective_rain * self.canopy
            
            # calculate the direct throughfall (mm) 
            direct_throughfall    = effective_rain - leaf_drainage
            
            # calculate the kinetic energy of raindrops that fall through the canopy (J / m2)
            kinetic_throughfall   = direct_throughfall * (11.9 + 8.7 * log10(self.intensity))
            
            # calculate the kinetic energy of raindrops that fall on the canopy (J / m2)
            kinetic_leaf_drainage = leaf_drainage * (15.8 - (self.pltheight ** 0.5) - 5.87) 
            
            # calculate the total kinetic energy of raindrops hitting the ground (J / m2)
            total_kinetic         = kinetic_throughfall + kinetic_leaf_drainage
            
            # annual rate of soil particle detachment by raindrop impact (kg / m2)
            rain_detachment       = self.detach * total_kinetic / 1000 
            
            # factor based on the cohesion of the soil (1 / kPa)
            Z_factor              = 1 / (0.5 * self.cohesion)
            
            # annual rate of soil particle detachment by runoff (kg / m2)
            runoff_detachment     = Z_factor * (self.totaldis ** 1.5) * sin(self.slope) * (1 - self.ground) / 1000
            
            # annual rate of total soil particle detachment (kg m-2)
            total_detachment      = rain_detachment + runoff_detachment
            
            # separate detachment for different land covers (kg m-2)
            bare_detachment       = ifthenelse(self.areacover == 0, total_detachment, 0)
            veg_detachment        = ifthenelse(self.areacover == 1, total_detachment, 0)
            forrest_detachment    = ifthenelse(self.areacover == 2, total_detachment, 0)
            
            # total amount of annual detachment of soil over the whole catchment (kg m-2)
            map_detachment        = maptotal(total_detachment) / self.totalcells
            
            # convert to numpy 
            self.detach_np[int(self.currentTimeStep() / 365 - 1)] = pcraster.numpy_operations.pcr2numpy(map_detachment, 0)[0][0]
            
            # total amoung of annual detachment of soil for different land covers (kg m-2)
            baremap_detachment    = maptotal(bare_detachment) / self.barecells
            vegmap_detachment     = maptotal(veg_detachment) / self.vegcells
            forrestmap_detachment = maptotal(forrest_detachment) / self.forrestcells
            
            # convert to numpy
            self.bare_np[int(self.currentTimeStep() / 365 - 1)]    = pcraster.numpy_operations.pcr2numpy(baremap_detachment, 0)[0][0]
            self.veg_np[int(self.currentTimeStep() / 365 - 1)]     = pcraster.numpy_operations.pcr2numpy(vegmap_detachment, 0)[0][0]
            self.forrest_np[int(self.currentTimeStep() / 365 - 1)] = pcraster.numpy_operations.pcr2numpy(forrestmap_detachment, 0)[0][0]
            
            # calculate the rain and runoff detachment over the whole catchment (kg m-2)
            map_raindetachment = maptotal(rain_detachment) / self.totalcells_r
            map_runoffdetachment = maptotal(runoff_detachment) / self.totalcells
            
            # convert to numpy
            self.raind_np[int(self.currentTimeStep() / 365 - 1)] = pcraster.numpy_operations.pcr2numpy(map_raindetachment, 0)[0][0]
            self.rund_np[int(self.currentTimeStep() / 365 - 1)] = pcraster.numpy_operations.pcr2numpy(map_runoffdetachment, 0)[0][0]
            
            # save the detachment to disk
            self.report(total_detachment, 'totd/totd')

            # reset the total rainfall and discharge to 0 for next year
            self.totalrain = scalar(0)
            self.totaldis = scalar(0)

# Amount of time steps in the model            
nrOfTimeSteps=1461


#%% Model scenarios for report results: with and without landcover

# initialize variables and output arrays
toggle_list     = [True, False]
snow_toggle     = True
cover_toggle    = True
bare_erosion    = np.zeros(2)
veg_erosion     = np.zeros(2)
forrest_erosion = np.zeros(2)

# run model twice, for both the scenario with and without snowdeck formation

for i, toggle in enumerate(toggle_list):

    cover_toggle = toggle
    myModel      = MyFirstModel(snow_toggle, cover_toggle)
    
    dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
    dynamicModel.setQuiet()
    dynamicModel.run()
    
    bare_erosion[i]    = np.mean(myModel.bare_np)
    veg_erosion[i]     = np.mean(myModel.veg_np)
    forrest_erosion[i] = np.mean(myModel.forrest_np)

#%% Model scenarios for report results: with and without snowdeck

# initialize variables and output arrays
toggle_list     = [True, False]
snow_toggle     = True
cover_toggle    = True
erosion         = [0,0]

# run model twice, for both the scenario with and without concideration of landcover
for i, toggle in enumerate(toggle_list):
    
    snow_toggle  = toggle
    myModel      = MyFirstModel(snow_toggle, cover_toggle)
    
    dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
    dynamicModel.setQuiet()
    dynamicModel.run()
    
    erosion[i]   = myModel.detach_np
    
#%% Model values for presentation discussion: total amount of rainsplash and runoff erosion

snow_toggle     = True
cover_toggle    = True

# run the model with both snow cover and landcover
myModel      = MyFirstModel(snow_toggle, cover_toggle)
   
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.setQuiet()
dynamicModel.run()

rainero   = myModel.raind_np
runoffero = myModel.rund_np

#%% plotting of the result figure for the model with and without snow deck

# set resolution for plot
plt.figure(dpi = 600)

# make list with the different model scenarios
names = ['with landcover', 'without landcover']

# setting the width of the bars 
bar_width = 0.2

# determining bar positions the three land cover categories
bar_positions1 = np.arange(len(names))
bar_positions2 = bar_positions1 + bar_width
bar_positions3 = bar_positions1 + 2*bar_width

# plotting the barplot
plt.bar(bar_positions1, bare_erosion, width=bar_width, label='Bare soil', color='grey')
plt.bar(bar_positions2, veg_erosion, width=bar_width, label='Low vegetation', color='lawngreen')
plt.bar(bar_positions3, forrest_erosion, width=bar_width, label='Forest', color='darkgreen')

# setting plot layout parameters
plt.xticks(bar_positions2, names)
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize = 14)
plt.rc('legend', fontsize = 11)

# axis labels
plt.xlabel('Model run')
plt.ylabel('Average yearly erosion $\mathregular{kg\;m^{-2}}$')

# including the legend and showing the plot
plt.legend()
plt.show()

#%% plotting of the result figure for the model with and without snow deck

# set resolution for plot
plt.figure(dpi = 600)

# make list with the years for which the model is run
names = ['1990', '1991', '1992', '1993']

# setting the width of the bars 
bar_width = 0.35

# determining bar positions the scenarios with and without snowdeck
bar_positions1 = np.arange(len(names))
bar_positions2 = bar_positions1 + bar_width
text_position =  bar_positions1 + 0.5 * bar_width

# plotting the barplot
plt.bar(bar_positions1, erosion[0], width=bar_width, label='With snowdeck', color='deepskyblue')
plt.bar(bar_positions2, erosion[1], width=bar_width, label='Without snowdeck', color='sandybrown')

# setting plot layout parameters
plt.xticks(text_position, names)
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize = 16)
plt.rc('legend', fontsize = 13)

# axis labels
plt.xlabel('Year')
plt.ylabel('Yearly soil erosion $\mathregular{kg\;m^{-2}}$')

# including the legend and showing the plot
plt.legend()
plt.show()

#%% plotting of the presentation discussion figure of splash and runoff erosion

# set resolution for plot
plt.figure(dpi = 600)

# make list with the years for which the model is run
names = ['1990', '1991', '1992', '1993']

# setting the width of the bars 
bar_width = 0.35

# determinining bar location for splash and runoff erosion
bar_positions1 = np.arange(len(names))
bar_positions2 = bar_positions1 + bar_width
text_positions = bar_positions1 + 0.5 * bar_width

# plotting the barplot
plt.bar(bar_positions1, rainero, width=bar_width, label='Rain erosion', color='deepskyblue')
plt.bar(bar_positions2, runoffero, width=bar_width, label='Runoff erosion', color='sandybrown')

# setting plot layout parameters
plt.xticks(text_positions, names)
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize = 16)
plt.rc('legend', fontsize = 13)

# axis labels
plt.xlabel('Year')
plt.ylabel('Yearly soil erosion $\mathregular{kg\;m^{-2}}$')

# including the legend and showing the plot
plt.legend()
plt.show()
