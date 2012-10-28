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
    self.iters = f['iterations'].value[0]

    # grid factor
    self.factor = factor

################################################################################

def gnuplot_file(gs_unrod,pj_unrod,gs_rod,pj_rod):

  datastr = ""
  for item in gs_unrod:
    datastr += "{factor} {iters} \n".format(factor=item.factor,       \
                                          iters=item.iters)

  datastr += "e\n"
  for item in pj_unrod:
    datastr += "{factor} {iters} \n".format(factor=item.factor,       \
                                          iters=item.iters)

  datastr += "e\n"
  for item in gs_rod:
    datastr += "{factor} {iters} \n".format(factor=item.factor,       \
                                          iters=item.iters)

  datastr += "e\n"
  for item in pj_rod:
    datastr += "{factor} {iters} \n".format(factor=item.factor,       \
                                          iters=item.iters)

  datastr += "e\n"



  outstr = """
set terminal pdf
set output "inner.pdf"
set key horizontal
set key bmargin center
set log x
set log y
set xlabel "Inner Iteration Tolerance [-]"
set ylabel "Number of Fission Source Iterations [-]"
plot '-' using 1:2 with linespoint title "G-S UNRODDED", '-' using 1:2 with linespoint title "P-J UNRODDED", '-' using 1:2 with linespoint title "G-S RODDED", '-' using 1:2 with linespoint title "P-J RODDED"
  """

  outstr += datastr

  with open("gnuplot.dat",'w') as f:
    f.write(outstr)

#-------------------------------------------------------------------------------

def process_all():

  # write filenames manually to process
  gs_unrod = []
  pj_unrod = []
  gs_rod = []
  pj_rod = []
  temp = result('./GS/tol1/output.h5',1.e-1)
  gs_unrod.append(temp)
  temp = result('./GS/tol2/output.h5',1.e-2)
  gs_unrod.append(temp)
  temp = result('./GS/tol3/output.h5',1.e-3)
  gs_unrod.append(temp)
  temp = result('./GS/tol4/output.h5',1.e-4)
  gs_unrod.append(temp)
  temp = result('./GS/tol5/output.h5',1.e-5)
  gs_unrod.append(temp)
  temp = result('./GS/tol6/output.h5',1.e-6)
  gs_unrod.append(temp)
  temp = result('./GS/tol7/output.h5',1.e-7)
  gs_unrod.append(temp)
  temp = result('./GS/tol8/output.h5',1.e-8)
  gs_unrod.append(temp)
  temp = result('./GS/tol9/output.h5',1.e-9)
  gs_unrod.append(temp)
  temp = result('./GS/tol10/output.h5',1.e-10)
  gs_unrod.append(temp)


  temp = result('./PJ/tol1/output.h5',1.e-1)
  pj_unrod.append(temp)
  temp = result('./PJ/tol2/output.h5',1.e-2)
  pj_unrod.append(temp)
  temp = result('./PJ/tol3/output.h5',1.e-3)
  pj_unrod.append(temp)
  temp = result('./PJ/tol4/output.h5',1.e-4)
  pj_unrod.append(temp)
  temp = result('./PJ/tol5/output.h5',1.e-5)
  pj_unrod.append(temp)
  temp = result('./PJ/tol6/output.h5',1.e-6)
  pj_unrod.append(temp)
  temp = result('./PJ/tol7/output.h5',1.e-7)
  pj_unrod.append(temp)
  temp = result('./PJ/tol8/output.h5',1.e-8)
  pj_unrod.append(temp)
  temp = result('./PJ/tol9/output.h5',1.e-9)
  pj_unrod.append(temp)
  temp = result('./PJ/tol10/output.h5',1.e-10)
  pj_unrod.append(temp)


  temp = result('./ROD/GS/tol1/output.h5',1.e-1)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol2/output.h5',1.e-2)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol3/output.h5',1.e-3)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol4/output.h5',1.e-4)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol5/output.h5',1.e-5)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol6/output.h5',1.e-6)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol7/output.h5',1.e-7)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol8/output.h5',1.e-8)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol9/output.h5',1.e-9)
  gs_rod.append(temp)
  temp = result('./ROD/GS/tol10/output.h5',1.e-10)
  gs_rod.append(temp)


  temp = result('./ROD/PJ/tol1/output.h5',1.e-1)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol2/output.h5',1.e-2)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol3/output.h5',1.e-3)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol4/output.h5',1.e-4)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol5/output.h5',1.e-5)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol6/output.h5',1.e-6)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol7/output.h5',1.e-7)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol8/output.h5',1.e-8)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol9/output.h5',1.e-9)
  pj_rod.append(temp)
  temp = result('./ROD/PJ/tol10/output.h5',1.e-10)
  pj_rod.append(temp)




  # generate gnuplot file
  gnuplot_file(gs_unrod,pj_unrod,gs_rod,pj_rod)

  return

################################################################################

if __name__ == "__main__":

  r = process_all() 
