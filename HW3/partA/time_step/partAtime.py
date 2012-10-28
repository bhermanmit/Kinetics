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
set output "partA_timestep.pdf"
set key horizontal
set key bmargin center
set log x
set log y
set grid
set xlabel "Inverse of time step [1/s]"
set ylabel "Power Deviation [%]"
set title "Power Deviation from Grid Refinement (ref dt = 1e-4 s)"
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
set output "partA_powers.pdf"
set xlabel "Time [s]"
set ylabel "Relative Core Power [-]"
set title "Core Power Evolution During Transient"
set grid
plot '-' using 1:2 with lines linewidth 2.0 title "dt = 1 s", '-' using 1:2 with lines linewidth 2.0 title "dt = 0.5 s", '-' using 1:2 with lines linewith 2.0 title "dt = 0.1 s",'-' using 1:2 with lines linewidth 2.0 title "dt = 0.005 s", '-' using 1:2 with lines linewidth 2.0 title "dt = 0.001s"
  """
  outstr += datastr

  with open("gnuplot_power.dat",'w') as f:
    f.write(outstr)

  datastr = ""
  i = 0
  while i < len(results[3].power):
    datastr += "{time} {power}\n".format(time=results[3].time[i],power=results[3].power[i])
    i += 1

  datastr += "e\n"

  outstr = """
set terminal pdf
set output "partA_powerconv.pdf"
set xlabel "Time [s]"
set ylabel "Relative Core Power [-]"
set title "Core Power Evolution During Transient"
set grid
plot '-' using 1:2 with lines linewidth 2.0 title "dt = 0.1 s"
  """
  outstr += datastr

  with open("gnuplot_powerconv.dat",'w') as f:
    f.write(outstr)



#-------------------------------------------------------------------------------

def process_all():

  # write filenames manually to process
  results = []
  temp = result('./time1/output.h5',1.0)
  results.append(temp)
  temp = result('./time2/output.h5',2.0)
  results.append(temp)
  temp = result('./time5/output.h5',5.0)
  results.append(temp)
  temp = result('./time10/output.h5',10.0)
  results.append(temp)
  temp = result('./time20/output.h5',20.0)
  results.append(temp)
  temp = result('./time50/output.h5',50.0)
  results.append(temp)
  temp = result('./time100/output.h5',100.0)
  results.append(temp)
  temp = result('./time10000/output.h5',10000.0)
  results.append(temp)

  # begin loop around results
  i = 0
  while i < len(results)-1:
#   result.power_ref = max(results[i+1].power)
    result.power_ref = max(results[7].power)
    results[i].process_max_power()
#   result.power_ref = results[i+1].power[len(results[i+1].power)-1]
    result.power_ref = results[7].power[len(results[7].power)-1]
    results[i].process_end_power()
    i += 1

  # generate gnuplot file
  gnuplot_file(results)

  return results

################################################################################

if __name__ == "__main__":

  r = process_all() 
