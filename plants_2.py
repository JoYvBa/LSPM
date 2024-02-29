from pcraster import *
from pcraster.framework import *

import os
os.chdir(r"/Users/jorrit/DynamicModelling/probab/plants")

class Plants(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone('clone.map')

  def initial(self):
    self.vegg = readmap("veg")
    self.d_herbs = self.vegg == 5
    self.biotope = pcror(self.vegg == 5, self.vegg == 3)
    self.r = 40
    aguila(self.d_herbs)
    
    aguila(self.biotope, self.vegg)

  def dynamic(self):
    herb_distance = spread(self.d_herbs, 0, 1)
    probs = 0.1**(herb_distance / self.r)
    uni = uniform(1)
    seeded = uni < probs
    growable = seeded & self.biotope
    self.d_herbs = growable | self.d_herbs
    #self.report(herb_distance, 'hdis')
    self.report(probs, 'prb')
    self.report(seeded, 'sdd')
    self.report(self.d_herbs, 'hrb')

    
nrOfTimeSteps=1000
myModel = Plants()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()

