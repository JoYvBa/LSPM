from pcraster import *
from pcraster.framework import *

import os
os.chdir(r"/Users/jorrit/DynamicModelling/probab/fire")

class RandomModel(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone('clone.map')

  def initial(self):
    pass

  def dynamic(self):
    #normal_d = normal(1)
    #uniform_d = uniform(1)
    #mapuniform_d = mapuniform()
    #gaus_10_5 = normal(1) * 5 + 10
    unif = uniform(1)
    #distr = ifthenelse(unif < 0.8, boolean(1), boolean(0))
    classes = lookupnominal('threes.txt', unif)
    
    self.report(classes, 'rcl')
    #self.report(gaus_10_5, 'gaus')
    #self.report(normal_d, 'nrm')
    #self.report(uniform_d, 'uni')
    #self.report(mapuniform_d, 'muni')
    #self.report(distr, 'dist')
    
nrOfTimeSteps=10
myModel = RandomModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()

''' make ASCII file containing the following:
[0,0.01> 1
[0.01,0.91> 2
[0.91,1.0] 3
'''

