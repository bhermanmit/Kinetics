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
    self.power = f['core_power'].value
    self.time = f['time'].value

    # grid factor
    self.factor = factor

#-------------------------------------------------------------------------------

  def process_max_power(self):

    self.max_power = abs(max(self.power) - self.power_ref) / self.power_ref

#-------------------------------------------------------------------------------

  def process_end_power(self):

    self.end_power = abs(self.power[len(self.power)-1] - self.power_ref) / self.power_ref

################################################################################

def gnuplot_file(results):

  datastr = ""
  i = 0
  while i < len(results):
    item = results[i] 
    datastr += "{factor} {maxpower}\n".format(factor=item.factor,       \
                                          maxpower=item.max_power)
    i += 1

  datastr += "e\n"

  i = 0
  while i < len(results):
    item = results[i]
    datastr += "{factor} {endpower}\n".format(factor=item.factor,       \
                                          endpower=item.end_power)
    i += 1

  datastr += "e\n"

  outstr = """
set terminal pdf 
set output "partA_timeerror.pdf"
set key left top
set xlabel "Time step multiplier [-]"
set ylabel "Error in Power from Reference [-]"
set title "Fraction Error in Core Power vs. Time step (Ref. dt = 0.005s)"
plot '-' using 1:2 with linespoint linewidth 3.0 title "Max Core Power", '-' using 1:2 with linespoint linewidth 3.0 title "End Core Power"
  """

  outstr += datastr

  with open("gnuplot_time.dat",'w') as f:
    f.write(outstr)

#-------------------------------------------------------------------------------

def process_all():

  # write filenames manually to process
  results = []
  temp = result('./time_1/output.h5',1.0)
  results.append(temp)
  temp = result('./time_2/output.h5',2.0)
  results.append(temp)
  temp = result('./time_4/output.h5',4.0)
  results.append(temp)
  temp = result('./time_8/output.h5',8.0)
  results.append(temp)
  temp = result('./time_16/output.h5',16.0)
  results.append(temp)
  temp = result('./time_32/output.h5',32.0)
  results.append(temp)

  # set reference
  result.power_ref = max(results[0].power)

  # begin loop around results
  i = 0
  while i < len(results):
    results[i].process_max_power()
    i += 1

  # set reference
  result.power_ref = results[0].power[len(results[0].power)-1]

  # begin loop around results
  i = 0
  while i < len(results):
    results[i].process_end_power()
    i += 1

  # generate gnuplot file
  gnuplot_file(results)

  return results

################################################################################

if __name__ == "__main__":

  r = process_all() 
