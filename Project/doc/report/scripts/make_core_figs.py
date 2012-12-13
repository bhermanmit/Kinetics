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

  fig = CoreFig(caption,label)

  fig.set_legend()

  fig.write_fig(out+'.tex')

###################################################################################

  caption = "test"
  altcap = ''
  label = "fig_test"
  out = "test"

  fig = CoreFig(caption,label,altcap=altcap,scale=1.1,scalebox=0.9,colorbyval=True)

  text="1.0"
  v = 1.0
  fig.set_pos("A1",text,val=v)
  text="0.0"
  v = 0.1
  fig.set_pos("D10",text,val=v)
  text="0.5"
  v = 0.5
  fig.set_pos("K7",text,val=v)

  fig.write_fig(out+'.tex')

###################################################################################

  caption = "Time " + num
  altcap = ''
  label = "fig_test"
  out = "Time" + num

  fig = CoreFig(caption,label,altcap=altcap,scale=1.1,scalebox=0.9,colorbyval=True)

  LETS = ['A','B','C','D','E','F','G','H','J','K','L']
  NUMS = ['1','2','3','4','5','6','7','8','9','10','11']
  f = h5py.File(filepath,'r')
  assy_pow = f['time'+num]['assy_pow']
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
