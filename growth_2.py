from pcraster import *
from pcraster.framework import *
import os
os.chdir(r"/Users/jorrit/DynamicModelling/neighbourhood/growth")

class Growth(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone('clone.map')

  def initial(self):
    self.x = spatial(scalar(8.5))
    self.r = 0.08
    self.k = 10
    self.c = 0.1
    self.cinc = 0.00006
    self.random = normal(1) / 10.0
    self.d = 0.01
    self.numberOfNeighbours = window4total(spatial(scalar(1)))
    #aguila(self.random)

  def dynamic(self):
      
     diffusion = self.d*(window4total(self.x)-self.numberOfNeighbours*self.x)
     dBdt = self.r * self.x * (1 - (self.x / self.k)) - self.c * (self.x**2 / (self.x**2 + 1))
     self.x = dBdt + self.x + self.random + self.d * diffusion
     self.x = ifthenelse(self.x < 0, 0, self.x)
     self.c = self.c + self.cinc
     B_mean = maptotal(self.x) / (40 * 40)
     B_var = maptotal((self.x - B_mean)**2) / (40 * 40)
     
     self.report(self.x, 'b')
     self.report(B_mean, 'mean')
     self.report(B_var, 'var')

nrOfTimeSteps=2500
myModel = Growth()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()

  




