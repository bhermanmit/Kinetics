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

    # read in maps
    self.nfx = f['nfx'].value[0]
    self.nfy = f['nfy'].value[0]
    self.nfz = f['nfz'].value[0]
    self.nfg = f['nfg'].value[0]

    # read in flux
    self.phi = f['phi'].value

    # reshape all
    self.phi = np.reshape(self.phi,(self.nfg,self.nfx,self.nfy,self.nfz),order='F')

################################################################################

def gnuplot_file(for_unrod,adj_unrod,for_rod,adj_rod):

  xx = np.arange(0.0,420.0,0.2)
  datastr = ""
  i = 0
  while i < for_unrod.nfx:
    datastr += "{x} {y}\n".format(x=xx[i],y=for_unrod.phi[0,i,0,0])
    i += 1

  datastr += "e\n"

  i = 0
  while i < for_unrod.nfx:
    datastr += "{x} {y}\n".format(x=xx[i],y=for_unrod.phi[1,i,0,0])
    i += 1

  datastr += "e\n"

  i = 0
  while i < adj_unrod.nfx:
    datastr += "{x} {y}\n".format(x=xx[i],y=adj_unrod.phi[0,i,0,0])
    i += 1
  
  datastr += "e\n"

  i = 0
  while i < adj_unrod.nfx:
    datastr += "{x} {y}\n".format(x=xx[i],y=adj_unrod.phi[1,i,0,0])
    i += 1

  datastr += "e\n"

  outstr = """
set terminal pdf dashed
set output "adjoint_unrod.pdf"
set key horizontal
set key bmargin center
set xlabel "Core Length [cm]"
set ylabel "Forward Flux [-]"
set y2label "Adjoint Flux [-]"
set xrange[0:420]
set title "Unrodded Case"
set ytics nomirror
set y2tics
plot '-' using 1:2 with lines axes x1y1 linewidth 3 linetype 1 linecolor rgb "red" title "Group 1 Forward", '-' using 1:2 with lines axes x1y1 linetype 1 linewidth 3 linecolor rgb "blue" title "Group 2 Forward", '-' using 1:2 with lines axes x1y2 linewidth 3 linetype 2 linecolor rgb "red" title "Group 1 Adjoint", '-' using 1:2 with lines axes x1y2 linewidth 3 linetype 2 linecolor rgb "blue" title "Group 2 Adjoint"
  """

  outstr += datastr

  with open("gnuplot_unrod.dat",'w') as f:
    f.write(outstr)

  datastr = ""
  i = 0
  while i < for_unrod.nfx:
    datastr += "{x} {y}\n".format(x=xx[i],y=for_rod.phi[0,i,0,0])
    i += 1

  datastr += "e\n"

  i = 0
  while i < for_unrod.nfx:
    datastr += "{x} {y}\n".format(x=xx[i],y=for_rod.phi[1,i,0,0])
    i += 1

  datastr += "e\n"

  i = 0
  while i < adj_unrod.nfx:
    datastr += "{x} {y}\n".format(x=xx[i],y=adj_rod.phi[0,i,0,0])
    i += 1

  datastr += "e\n"

  i = 0
  while i < adj_unrod.nfx:
    datastr += "{x} {y}\n".format(x=xx[i],y=adj_rod.phi[1,i,0,0])
    i += 1

  datastr += "e\n"

  outstr = """
set terminal pdf dashed
set output "adjoint_rod.pdf"
set key horizontal
set key bmargin center
set xlabel "Core Length [cm]"
set ylabel "Forward Flux [-]"
set y2label "Adjoint Flux [-]"
set xrange[0:420]
set title "Rodded Case"
set ytics nomirror
set y2tics
plot '-' using 1:2 with lines axes x1y1 linewidth 3 linetype 1 linecolor rgb "red" title "Group 1 Forward", '-' using 1:2 with lines axes x1y1 linetype 1 linewidth 3 linecolor rgb "blue" title "Group 2 Forward", '-' using 1:2 with lines axes x1y2 linewidth 3 linetype 2 linecolor rgb "red" title "Group 1 Adjoint", '-' using 1:2 with lines axes x1y2 linewidth 3 linetype 2 linecolor rgb "blue" title "Group 2 Adjoint"
  """

  outstr += datastr

  with open("gnuplot_rod.dat",'w') as f:
    f.write(outstr)


#-------------------------------------------------------------------------------

def process_all():

  # write filenames manually to process
  for_unrod = result('./UNROD/FOR/output.h5',1.)
  adj_unrod = result('./UNROD/ADJ/output.h5',1.)
  for_rod = result('./ROD/FOR/output.h5',1.)
  adj_rod = result('./ROD/ADJ/output.h5',1.)

  # generate gnuplot file
  gnuplot_file(for_unrod,adj_unrod,for_rod,adj_rod)

  return

################################################################################

if __name__ == "__main__":

  r = process_all() 
