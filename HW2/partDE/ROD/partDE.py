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

  def __init__(self,filename,factor):

    # open hdf5 file
    f = h5py.File(filename,'r')

    # read in iter 
    self.iters = f['iter'].value[0]

    # grid factor
    self.factor = factor

################################################################################

def gnuplot_file(unrod,rod):

  datastr = ""
  for item in unrod:
    datastr += "{factor} {iters} \n".format(factor=1.0/item.factor,       \
                                          iters=item.iters)

  datastr += "e\n"
  for item in rod:
    datastr += "{factor} {iters} \n".format(factor=1.0/item.factor,       \
                                          iters=item.iters)

  datastr += "e\n"

  outstr = """
set terminal wxt persist
set key horizontal
set key bmargin center
set log x
plot '-' using 1:2 with linespoint axes x1y1 title "UN-RODDED", '-' using 1:2 with linespoint axes x1y2 title "RODDED"
  """

  outstr += datastr

  with open("gnuplot.dat",'w') as f:
    f.write(outstr)

#-------------------------------------------------------------------------------

def process_all():

  # write filenames manually to process
  unrod = []
  rod = []
  temp = result('./UNROD/tol1/output.h5',1.e-1)
  unrod.append(temp)
  temp = result('./UNROD/tol2/output.h5',1.e-2)
  unrod.append(temp)
  temp = result('./UNROD/tol3/output.h5',1.e-3)
  unrod.append(temp)
  temp = result('./UNROD/tol4/output.h5',1.e-4)
  unrod.append(temp)
  temp = result('./UNROD/tol5/output.h5',1.e-5)
  unrod.append(temp)
  temp = result('./ROD/tol1/output.h5',1.e-1)
  rod.append(temp)
  temp = result('./ROD/tol2/output.h5',1.e-2)
  rod.append(temp)
  temp = result('./ROD/tol3/output.h5',1.e-3)
  rod.append(temp)
  temp = result('./ROD/tol4/output.h5',1.e-4)
  rod.append(temp)
  temp = result('./ROD/tol5/output.h5',1.e-5)
  rod.append(temp)

  # generate gnuplot file
  gnuplot_file(results)

  return results

################################################################################

if __name__ == "__main__":

  r = process_all() 
