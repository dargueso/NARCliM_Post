#!/usr/bin/env python
"""postprocess_NARCliM.py script
   Script to postprocess WRF outputs from NARCliM project
   It reads an input file (NARClIM_post.input), where the input arguments are provided
	
   Authors: Daniel Argueso (d.argueso@unsw.edu.au), Alejandro Di Luca (a.diluca@unsw.edu.au)
   Institution: CoECSS, UNSW. CCRC, UNSW.
   Created: 13 September 2013
   Modified: 14 Septemvber 2013
"""

import netCDF4 as nc
import numpy as np
import sys
import os 
import datetime as dt
import glob
from postprocess_modules import *
import calendar

# Check initial time
ctime_i=checkpoint(0)
ctime=checkpoint(0)

#### READING INPUT FILE ######

##############################

#pathin=input['pathin']
#pathout=input['pathout']
#GCM=input['GCM']
#RCM=input['RCM']
#syear=input['start_year']
#eyear=input['end_year']
#out_variables=input['out_variables']

pathin='/home/z3393020/WRFouts/NARCliM/MK30/R1/1990-2010/'
pathout='/srv/ccrc/data28/z3444417/NARCliM/new_post/'
GCM='MK3.0'
RCM='R1'
syear=1990
eyear=2009
out_variables=['tas']

file_type=['wrfhrly', 'wrfout', 'wrfxtrm', 'wrfdly']
domain='d02'
out_variables=['T2']

if GCM=='CCMA':
	calendar=='noleap'
else:
	calendar='standard'

#CREATE OUTPUT DIR IF IT DOESN'T EXIST
if not os.path.exists(pathout):
	os.makedirs(pathout)

#CREATE A TEMPORAL DIR WITHIN THE OUTPUT DIR IF IT DOESN'T EXIST
if not os.path.exists("%s/temp/" %(pathout)):
	os.makedirs("%s/temp/" %(pathout))


	
#***********************************************
# Loop over all types of WRF output files (i.e., wrfhrly, wrfout, etc) 
for file in file_type:

	if file=='wrfhrly':
		print ' Processing ', file, ' WRF outputs'
		

		#***********************************************
		# LOOP over years
		for year in np.arange(syear,eyear+1):
			print ' -> Processing year: ', year
			
			loadfiles = pathin+'%s_%s_%s*' % (file,domain,year) # Specify path
			files_in=sorted(glob.glob(loadfiles))

			print '   Number of files to read:', len(files_in)

			# CHECKING: Check if the number of files is right
			if len(files_in)!=12:
				print 'ERROR: the number of ',file, ' files in year ', year,' is INCORRECT'
				print ' ---- SOME FILES ARE MISSING ---'
				print 'SCRIPT stops running '
				sys.exit(0)
			
			# READ FILES
			fin=nc.MFDataset(files_in) #Read all files
			time = fin.variables['Times'][:] # Get time variable
			lat=fin.variables['XLAT']
			lon=fin.variables['XLONG']
			print '   READ LATITUDE, LONGITUDE AND TIMES'
			
			# CHECKING: Check if the of time steps is right
			n_timesteps = dt.datetime(year+1,01,01,00)-dt.datetime(year,01,01,00)*24
			if calendar=='noleap' and calendar.isleap(year)==True:
				n_timesteps = (dt.datetime(year+1,01,01,00)-dt.datetime(year,01,01,00)-1)*24
			
			if len(files_in)!=12:
				print 'ERROR: the number of ',file, ' files in year ', year,' is INCORRECT'
				print ' ---- SOME FILES ARE MISSING ---'
				print 'SCRIPT stops running '
				sys.exit(0)

			# ***********************************************
			# LOOP over variables
			for var in out_variables:
				print '   READ VARIABLE ', var
				varval=fin.variables['T2']
				varatt=fin.variables[var].ncattrs()

				create_netcdf(file_out, var, lat, lon, times, rotpole, varatt, overwrite=None)
				
				
				
    

    
    
    
    
    
    
