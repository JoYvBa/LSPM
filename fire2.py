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
    #aguila(self.dem, self.fire)

  def dynamic(self):
    scale_fire = scalar(self.fire)
    burning_neighbours = window4total(scale_fire)
    neighbour_burns = burning_neighbours > 0
    catch_fire = pcrand(neighbour_burns, pcrnot(self.fire))
    realiz = uniform(1) < 0.1
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

