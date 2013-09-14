#!/usr/bin/env python
"""postprocess_NARCliM.py script
   Script to postprocess WRF outputs from NARCliM project
   It reads an input file (NARClIM_post.input), where the input arguments are provided
	
   Authors: Daniel Argueso (d.argueso@unsw.edu.au), Alejandro Di Luca (a.diluca@unsw.edu.au)
   Institution: CoECSS, UNSW. CCRC, UNSW.
   Created: 13 September 2013
"""

import netCDF4 as nc
import numpy as np
import sys
import os 
import datetime as dt


#### READING INPUT FILE ######

##############################

pathin=input['pathin']
pathout=input['pathout']
GCM=input['GCM']
RCM=input['RCM']
syear=input['start_year']
eyear=input['end_year']
out_variables=input['out_variables']


#CREATE OUTPUT DIR IF IT DOESN'T EXIST
if not os.path.exists(pathout):
	os.makedirs(pathout)

#CREATE A TEMPORAL DIR WITHIN THE OUTPUT DIR IF IT DOESN'T EXIST
if not os.path.exist("%s/temp/" %(pathout)):
	os.makedirs("%s/temp/" %(pathout))


	
# Loop over all output variables
for file in filetype:
    for var in out_variables:
    
    
    
    
    
    
    
    
    
    
    
    
    