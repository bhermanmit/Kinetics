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

    self.max_power = abs(max(self.power) - self.power_ref) / max(self.power) * 100.0

#-------------------------------------------------------------------------------

  def process_end_power(self):

    self.end_power = abs(self.power[len(self.power)-1] - self.power_ref) / self.power[len(self.power)-1] * 100.0

################################################################################

def gnuplot_file(results):

  datastr = ""
  i = 0
  while i < len(results)-1:
    item = results[i] 
    datastr += "{factor} {maxpower}\n".format(factor=item.factor,       \
                                          maxpower=item.max_power)
    i += 1

  datastr += "e\n"

  i = 0
  while i < len(results)-1:
    item = results[i]
    datastr += "{factor} {endpower}\n".format(factor=item.factor,       \
                                          endpower=item.end_power)
    i += 1

  datastr += "e\n"


  outstr = """
set terminal pdf 
set output "partB_timestep.pdf"
set key horizontal
set key bmargin center
set log x
set log y
set grid
set xlabel "Inverse of time step [1/s]"
set ylabel "Power Deviation [%]"
set title "Power Deviation from Grid Refinement (ref dt = 2e-6 s)"
plot '-' using 1:2 with linespoint linewidth 3.0 title "Max Power", '-' using 1:2 with linespoint linewidth 3.0 title "Final Power"
  """

  outstr += datastr

  with open("gnuplot_timestep.dat",'w') as f:
    f.write(outstr)

  datastr = ""
  j = 0
  while j < len(results):

    i = 0
    while i < len(results[j].power):
      datastr += "{time} {power}\n".format(time=results[j].time[i],power=results[j].power[i])
      i += 1

    datastr += "e\n"
    j += 1

  outstr = """
set terminal pdf
set output "partB_powers.pdf"
set key horizontal
set key bmargin
set grid
set xlabel "Time [s]"
set ylabel "Relative Core Power [-]"
set log y
set title "Core Power Evolution During Transient"
set log y
plot '-' using 1:2 with lines linewidth 2.0 title "dt = 1e-2 s", '-' using 1:2 with lines linewidth 2.0 title "dt = 2e-2 s", '-' using 1:2 with lines linewidth 2.0 title "dt = 1e-3 s", '-' using 1:2 with lines linewidth 2.0 title "dt = 2e-3 s", '-' using 1:2 with lines linewidth 2.0 title "dt = 1e-4 s", '-' using 1:2 with lines linewidth 2.0 title "dt = 2e-4 s", '-' using 1:2 with lines linewidth 2.0 title "dt = 1e-5 s"
  """
  outstr += datastr

  with open("gnuplot_power.dat",'w') as f:
    f.write(outstr)

  datastr = ""
  i = 0
  while i < len(results[9].power):
    datastr += "{time} {power}\n".format(time=results[9].time[i],power=results[9].power[i])
    i += 1

  datastr += "e\n"

  outstr = """
set terminal pdf
set output "partB_powerconv.pdf"
set xlabel "Time [s]"
set ylabel "Relative Core Power [-]"
set log y
set grid
set title "Core Power Evolution During Transient"
plot '-' using 1:2 with lines linewidth 2.0 title "dt = 1e-5 s"
  """
  outstr += datastr

  with open("gnuplot_powerconv.dat",'w') as f:
    f.write(outstr)



#-------------------------------------------------------------------------------

def process_all():

  # write filenames manually to process
  results = []
  temp = result('./time100/output.h5',100.0)
  results.append(temp)
  temp = result('./time200/output.h5',200.0)
  results.append(temp)
  temp = result('./time500/output.h5',500.0)
  results.append(temp)
  temp = result('./time1000/output.h5',1000.0)
  results.append(temp)
  temp = result('./time2000/output.h5',2000.0)
  results.append(temp)
  temp = result('./time5000/output.h5',5000.0)
  results.append(temp)
  temp = result('./time10000/output.h5',10000.0)
  results.append(temp)
  temp = result('./time20000/output.h5',20000.0)
  results.append(temp)
  temp = result('./time50000/output.h5',50000.0)
  results.append(temp)
  temp = result('./time100000/output.h5',100000.0)
  results.append(temp)
  temp = result('./time200000/output.h5',200000.0)
  results.append(temp)
  temp = result('./time500000/output.h5',500000.0)
  results.append(temp)

  # begin loop around results
  i = 0
  while i < len(results)-1:
#   result.power_ref = max(results[i+1].power)
    result.power_ref = max(results[11].power)
    results[i].process_max_power()
#   result.power_ref = results[i+1].power[len(results[i+1].power)-1]
    result.power_ref = results[11].power[len(results[11].power)-1]
    results[i].process_end_power()
    i += 1

  # generate gnuplot file
  gnuplot_file(results)

  return results

################################################################################

if __name__ == "__main__":

  r = process_all() 
