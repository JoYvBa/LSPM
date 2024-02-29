from pcraster import *
from pcraster.framework import *

import os
os.chdir(r"/Users/jorrit/DynamicModelling/probab/fire")

class Fire(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone('clone.map')

  def initial(self):
    self.dem = readmap('dem')
    self.fire = readmap('start')
    self.veg = readmap('veg')
    self.ldd = lddcreate(self.dem, 1e31, 1e31, 1e31, 1e31)
    aguila(self.ldd, self.veg, self.dem)
    self.report(self.veg, 'veggi')

  def dynamic(self):
    scale_fire = scalar(self.fire)
    burning_neighbours = window4total(scale_fire)
    neighbour_burns = burning_neighbours > 0
    catch_fire = pcrand(neighbour_burns, pcrnot(self.fire))
    downstream_fire = downstream(self.ldd, self.fire)
    realiz = uniform(1)
    realiz = ifthenelse(downstream_fire, realiz < 0.8, realiz < 0.1)

    new_fire = pcrand(catch_fire, realiz)
    self.fire = pcror(new_fire, self.fire)
    
    self.report(neighbour_burns, 'brn')
    self.report(catch_fire, 'ctch')
    self.report(new_fire, 'nwf')
    self.report(self.fire, 'fir')

nrOfTimeSteps=200
myModel = Fire()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()

