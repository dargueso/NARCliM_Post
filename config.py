#!/usr/bin/env python
"""config.py script
	Configuration file for NARCliM postprocessing (postprocess_NARCliM.py)
	
   Authors: Daniel Argueso (d.argueso@unsw.edu.au), Alejandro Di Luca (a.diluca@unsw.edu.au)
   Institution: CoECSS, UNSW. CCRC, UNSW.
   Created: 19 September 2013
"""
import datetime as dt

#Defines wich calendar uses each GCM
GCM_calendar={}
GCM_calendar['MIROC3.2']='standard'
GCM_calendar['ECHAM5']='standard'
GCM_calendar['MK3.0']='standard'
GCM_calendar['CCCMA3.1']='no_leap'


#Reference date for output files
ref_date=dt.datetime(1949,12,1,00,00,00)