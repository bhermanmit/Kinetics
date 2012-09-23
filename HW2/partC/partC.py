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
    self.dr = f['dr'].value[0]

    # grid factor
    self.factor = factor

################################################################################

def gnuplot_file(unrod,rod):

  datastr = ""
  for item in unrod:
    datastr += "{factor} {dr} \n".format(factor=1.0/item.factor,       \
                                          dr=item.dr)

  datastr += "e\n"
  for item in rod:
    datastr += "{factor} {dr} \n".format(factor=1.0/item.factor,       \
                                          dr=item.dr)

  datastr += "e\n"

  outstr = """
set terminal pdf
set output "dr.pdf"
set key horizontal
set key bmargin center
set xlabel "Coarse Mesh to Fine Mesh Ratio [-]"
set ylabel "Dominance Ratio [-]
set log x
plot '-' using 1:2 with linespoint title "UN-RODDED", '-' using 1:2 with linespoint title "RODDED"
  """

  outstr += datastr

  with open("gnuplot.dat",'w') as f:
    f.write(outstr)

#-------------------------------------------------------------------------------

def process_all():

  # write filenames manually to process
  unrod = []
  rod = []
  temp = result('./mesh10/output.h5',10.0)
  unrod.append(temp)
  temp = result('./mesh5/output.h5',5.0)
  unrod.append(temp)
  temp = result('./mesh1/output.h5',1.0)
  unrod.append(temp)
  temp = result('./mesh05/output.h5',0.5)
  unrod.append(temp)
  temp = result('./mesh02/output.h5',0.2)
  unrod.append(temp)
  temp = result('./mesh013/output.h5',0.133333333333333333)
  unrod.append(temp)
# temp = result('./mesh011/output.h5',0.111111111111111111)
# unrod.append(temp)
# temp = result('./mesh01/output.h5',0.1)
# unrod.append(temp)
  temp = result('./ROD/mesh10/output.h5',10.0)
  rod.append(temp)
  temp = result('./ROD/mesh5/output.h5',5.0)
  rod.append(temp)
  temp = result('./ROD/mesh1/output.h5',1.0)
  rod.append(temp)
  temp = result('./ROD/mesh05/output.h5',0.5)
  rod.append(temp)
  temp = result('./ROD/mesh02/output.h5',0.2)
  rod.append(temp)
  temp = result('./ROD/mesh013/output.h5',0.133333333333333333)
  rod.append(temp)
# temp = result('./ROD/mesh011/output.h5',0.111111111111111111)
# rod.append(temp)
# temp = result('./ROD/mesh01/output.h5',0.1)
# rod.append(temp)

  # generate gnuplot file
  gnuplot_file(unrod,rod)

  return

################################################################################

if __name__ == "__main__":

  r = process_all()
