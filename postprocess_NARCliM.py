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
from optparse import OptionParser
from postprocess_modules import *
import calendar as cal

# Check initial time
ctime_i=checkpoint(0)
ctime=checkpoint(0)

#### READING INPUT FILE ######
### Options 

parser = OptionParser()

parser.add_option("-i", "--infile", dest="infile",
help="file with the input arguments", metavar="INPUTFILE")
(opts, args) = parser.parse_args()

###

inputinf=read_input(opts.infile)

##############################

pathin=inputinf['pathin']
pathout=inputinf['pathout']
GCM=inputinf['GCM']
RCM=inputinf['RCM']
syear=inputinf['start_year']
eyear=inputinf['end_year']
out_variables=inputinf['out_variables']
domain=inputinf['domain']
outfile_patt=inputinf['outfile_patt']

#CREATE OUTPUT DIR IF IT DOESN'T EXIST
fullpathout='%s/%s/%s/%s-%s/%s' %(pathout,GCM,RCM,syear,eyear,domain,)
if not os.path.exists(fullpathout):
	os.makedirs(fullpathout)

#CREATE A TEMPORAL DIR WITHIN THE OUTPUT DIR IF IT DOESN'T EXIST
if not os.path.exists("%s/temp/" %(fullpathout)):
	os.makedirs("%s/temp/" %(fullpathout))



#### Reading variable info file ######

varinfo=read_varinfo("./info_files/variables.inf")
file_type=varinfo.keys()


######################################


out_variables=['tas']
out_variables=['T2']

if GCM=='CCCMA3.1':
	calendar=='noleap'
else:
	calendar='standard'




	
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

			# -------------------
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
			
			dates_day = dt.datetime(year+1,01,01,00)-dt.datetime(year,01,01,00)
			dates = [datetime(int(year_i),01,01,00)+ timedelta(hours=x) for x in xrange(0,int(slp.shape[0])*ts,ts)]
			
			if calendar=='noleap' and calendar.isleap(year)==True:
				dates = (dt.datetime(year+1,01,01,00)-dt.datetime(year,01,01,00)-1)
				months_all=np.asarray([dates[i].month for i in xrange(len(dates))]) 
				days_all=np.asarray([dates[i].day for i in xrange(len(dates))])
				dates=dates[((months_all==2) & (days_all==29))==False]
			
			n_timesteps=n_timesteps_days*24

			months_all=np.asarray([dates[i].month for i in xrange(len(dates))]) 
			days_all=np.asarray([dates[i].day for i in xrange(len(dates))])
			dates=dates[((months_all==2) & (days_all==29))==False]

			# -------------------
			# CHECKING: Check if the of time steps is right
			if n_timesteps!=time.shape[0]:
				print 'ERROR: the number of timesteps in year ', year,' is INCORRECT'
				print 'There should be: ', n_timesteps
				print 'There are: ', time.shape[0]
				print 'SCRIPT stops running '
				sys.exit(0)


			# ***********************************************
			# LOOP over variables
			for var in out_variables:
				print '   READ VARIABLE ', var
				varval=fin.variables['T2']
				varatt=fin.variables[var].ncattrs()

				create_netcdf(file_out, var, lat, lon, times, rotpole, varatt, overwrite=None)
				
				
				
    

    
    
    
    
    
    
