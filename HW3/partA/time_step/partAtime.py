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

    self.max_power = abs(max(self.power) - self.power_ref) / self.power_ref * 100.0

#-------------------------------------------------------------------------------

  def process_end_power(self):

    self.end_power = abs(self.power[len(self.power)-1] - self.power_ref) / self.power_ref * 100.0

################################################################################

def gnuplot_file(results):

  datastr = ""
  i = 1
  while i < len(results):
    item = results[i] 
    datastr += "{factor} {maxpower}\n".format(factor=item.factor,       \
                                          maxpower=item.max_power)
    i += 1

  datastr += "e\n"

  outstr = """
set terminal pdf 
set output "partA_max.pdf"
set key horizontal
set key bmargin center
set log x
set xlabel "Coarse Mesh to Fine Mesh Ratio [-]"
set ylabel "Max Power Deviation [%]"
set title "Max Power Deviation from Grid Refinement"
plot '-' using 1:2 with linespoint linewidth 3.0
  """

  outstr += datastr

  with open("gnuplot_max.dat",'w') as f:
    f.write(outstr)

  datastr = ""
  i = 1
  while i < len(results):
    item = results[i]
    datastr += "{factor} {endpower}\n".format(factor=item.factor,       \
                                          endpower=item.end_power)
    i += 1

  datastr += "e\n"

  outstr = """
set terminal pdf 
set output "partA_end.pdf"
set key horizontal
set key bmargin center
set log x
set xlabel "Coarse Mesh to Fine Mesh Ratio [-]"
set ylabel "End Power Deviation [%]"
set title "End Power Deviation from Grid Refinement"
plot '-' using 1:2 with linespoint linewidth 3.0
  """

  outstr += datastr

  with open("gnuplot_end.dat",'w') as f:
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
plot '-' using 1:2 with lines linewidth 2.0 title "dt = 0.005 s"
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
  temp = result('./time05/output.h5',5.0)
  results.append(temp)
  temp = result('./time01/output.h5',10.0)
  results.append(temp)
  temp = result('./time005/output.h5',50.0)
  results.append(temp)
  temp = result('./time001/output.h5',100.0)
  results.append(temp)

  # begin loop around results
  i = 1
  while i < len(results):
    result.power_ref = max(results[i-1].power)
    results[i].process_max_power()
    result.power_ref = results[i-1].power[len(results[i-1].power)-1]
    results[i].process_end_power()
    i += 1

  # generate gnuplot file
  gnuplot_file(results)

  return results

################################################################################

if __name__ == "__main__":

  r = process_all() 
