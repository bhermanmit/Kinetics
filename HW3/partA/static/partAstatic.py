#! /usr/bin/python

# module imports
import math
import matplotlib.pyplot as plt
import numpy as np
import h5py

# begin module classes

################################################################################
# CLASS RESULT 
################################################################################

class result: 

  power_ref = 0.0

  def __init__(self,filename,factor):

    # open hdf5 file
    f = h5py.File(filename,'r')

    # read in keff
    self.keff = f['keff'].value[0]

#-------------------------------------------------------------------------------

def process_all():

  rodin = result('./RODIN/output.h5',1)
  rodout = result('./RODOUT/output.h5',1)

  beta = 0.006648

  static = (rodout.keff - rodin.keff)/rodout.keff/beta
  print "UNRODDED"+str(rodout.keff)
  print "RODDED"+str(rodin.keff)
  print static

################################################################################

if __name__ == "__main__":

  r = process_all() 
