#!/usr/bin/env python

from corefig import CoreFig

import sys
import os
import h5py
import numpy as np

CATAWBA=True
filepath = sys.argv[1]
num = sys.argv[2]

################################################################################
################################################################################
################################################################################

if CATAWBA:

################################################################################


  caption = "Layout of 2-D LRA Benchmark"
  label = "LRA_layout"
  out = 'LRA_layout' 

  fig = CoreFig(caption,label,scale=1.8,scalebox=0.9)

  fig.set_legend()

  fig.write_fig(out+'.tex')

###################################################################################


  LETS = ['A','B','C','D','E','F','G','H','J','K','L']
  NUMS = ['1','2','3','4','5','6','7','8','9','10','11']
  f = h5py.File(filepath,'r')
  assy_pow = f['time'+num]['assy_pow']
  assy_temp = f['time'+num]['assy_temp']
  time = f['time'+num]['time']
  pdens = f['time'+num]['pow']
  temp = f['time'+num]['avg_T']

  caption = "Normalized Assembly Power Distribution @ time {thetime:.4f}s, fuel avg. power density {powden:.3e}W/cc and fuel avg. temperature {temper:.1f}K".format(thetime=time[0], powden = pdens[0], temper=temp[0])
  altcap = 'Normalized Assembly Power Distribution'
  label = "power"+num
  out = "Power" + num

  fig = CoreFig(caption,label,altcap=altcap,scale=1.8,scalebox=0.9,colorbyval=True)

  i = 0
  while i < 11:

    j = 0
    while j < 11:

      assy = LETS[j]+NUMS[i]
      idx = j + 11*i
      text = "{0:.3f}".format(assy_pow[idx])
      v = assy_pow[idx]
      fig.set_pos(assy,text,val=v)
      j += 1
    i += 1
  fig.write_fig(out+'.tex')


  caption = "Assembly Temperature Distribution @ time {thetime:.4f}s, fuel avg. power density {powden:.3e}W/cc and fuel avg. temperature {temper:.1f}K".format(thetime=time[0], powden = pdens[0], temper=temp[0])
  altcap = 'Assembly Temperature Distribution'
  label = "temp"+num
  out = "Temperature" + num

  fig = CoreFig(caption,label,altcap=altcap,scale=1.8,scalebox=0.9,colorbyval=True)

  i = 0
  while i < 11:

    j = 0
    while j < 11:

      assy = LETS[j]+NUMS[i]
      idx = j + 11*i
      text = "{0:.1f}".format(assy_temp[idx])
      v = assy_temp[idx]
      fig.set_pos(assy,text,val=v)
      j += 1
    i += 1
  fig.write_fig(out+'.tex')
  f.close()
