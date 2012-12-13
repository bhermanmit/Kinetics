#!/usr/bin/env python

from corefig import CoreFig

import sys
import os

try:
  base = sys.argv[1]
except:
  base = ".."

CATAWBA=True

################################################################################
################################################################################
################################################################################

if CATAWBA:

################################################################################


  caption = "Layout of fuel assemblies showing enrichment loading pattern and burnable absorber positions. Source: \\ref{num:sheet_CA}"
  altcap = "Core enrichment zones and burnable absorber positions"
  label = "fig_enr_ba_pos"
  out = 'enr' 

  fig = CoreFig(caption,label,altcap=altcap)

  fig.set_legend()

  fig.write_fig(out+'.tex')

  caption = "test"
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
