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

    # read in power
    self.power = f['power'].value
    self.power = self.power/np.sum(self.power)*np.sum(self.power > 1.e-11)

    # read in maps
    self.nfx = f['nfx'].value[0]
    self.nfy = f['nfy'].value[0]
    self.nfz = f['nfz'].value[0]
    self.nfg = f['nfg'].value[0]
    self.mat = f['mat'].value
    self.reg = f['reg'].value

    # read in flux
    self.phi = f['phi'].value

    # reshape all
    self.phi = np.reshape(self.phi,(self.nfg,self.nfx,self.nfy,self.nfz),order='F')
    self.mat = np.reshape(self.mat,(self.nfx,self.nfy,self.nfz,self.nfg),order='F')
    self.reg = np.reshape(self.reg,(self.nfx,self.nfy,self.nfz,self.nfg),order='F')

    # grid factor
    self.factor = factor

################################################################################

def gnuplot_file(results):

  datastr = ""
  j = 0
  while(j < results.nfy):

    i = 0
    while(i < results.nfx):

      datastr += "{val} ".format(val = results.phi[0,i,j,0])

      i += 1

    datastr += "\n"
    j += 1

  datastr += "e\n"

  outstr = """
set terminal pdf
set output "LRA_flux1.pdf" 
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
set view map
set size ratio -1
set lmargin at screen 0.10
set rmargin at screen 0.90
set bmargin at screen 0.15
set tmargin at screen 0.90
set xrange [0.0:165.0]
set yrange [0.0:165.0]
set xlabel "Core x dimeimsion [-]"
set ylabel "Core y dimension [-]"
set title "Group 1 Flux Distribution"
splot '-' matrix with image 
  """

  outstr += datastr

  with open("gnuplot_flux1.dat",'w') as f:
    f.write(outstr)

  datastr = ""
  j = 0
  while(j < results.nfy):

    i = 0
    while(i < results.nfx):

      datastr += "{val} ".format(val = results.phi[1,i,j,0])

      i += 1

    datastr += "\n"
    j += 1

  datastr += "e\n"

  outstr = """
set terminal pdf
set output "LRA_flux2.pdf" 
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
set view map
set size ratio -1
set lmargin at screen 0.10
set rmargin at screen 0.90
set bmargin at screen 0.15
set tmargin at screen 0.90
set xrange [0.0:165.0]
set yrange [0.0:165.0]
set xlabel "Core x dimeimsion [-]"
set ylabel "Core y dimension [-]"
set title "Group 2 Flux Distribution"
splot '-' matrix with image 
  """

  outstr += datastr

  with open("gnuplot_flux2.dat",'w') as f:
    f.write(outstr)

  datastr = ""
  j = 0
  while(j < results.nfy):

    i = 0
    while(i < results.nfx):

      datastr += "{val} ".format(val = results.mat[i,j,0,0])

      i += 1

    datastr += "\n"
    j += 1

  datastr += "e\n"

  outstr = """
set terminal pdf 
set output "LRA_mat.pdf"
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
set view map
set size ratio -1
set lmargin at screen 0.10
set rmargin at screen 0.90
set bmargin at screen 0.15
set tmargin at screen 0.90
set xrange [0.0:166.0]
set yrange [0.0:166.0]
set xlabel "Core x dimeimsion [-]"
set ylabel "Core y dimension [-]"
set title "Material Map"
splot '-' matrix with image 
  """

  outstr += datastr

  with open("gnuplot_mat.dat",'w') as f:
    f.write(outstr)

  datastr = ""
  j = 0
  while(j < results.nfy):

    i = 0
    while(i < results.nfx):

      datastr += "{val} ".format(val = results.reg[i,j,0,0])

      i += 1

    datastr += "\n"
    j += 1

  datastr += "e\n"

  outstr = """
set terminal pdf 
set output "LRA_reg.pdf"
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
set view map
set size ratio -1
set lmargin at screen 0.10
set rmargin at screen 0.90
set bmargin at screen 0.15
set tmargin at screen 0.90
set xrange [0.0:166.0]
set yrange [0.0:166.0]
set xlabel "Core x dimeimsion [-]"
set ylabel "Core y dimension [-]"
set title "Region Map"
splot '-' matrix with image 
  """

  outstr += datastr

  with open("gnuplot_reg.dat",'w') as f:
    f.write(outstr)


#-------------------------------------------------------------------------------

def process_all():

  # write filenames manually to process
  results = result('./output.h5',1.0)

  # generate gnuplot file
  gnuplot_file(results)

  return results

################################################################################

if __name__ == "__main__":

  r = process_all() 
