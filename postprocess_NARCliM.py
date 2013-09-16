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
import postprocess_modules as pm
import compute_vars as cv
import calendar as cal

# Check initial time
ctime_i=pm.checkpoint(0)
ctime=pm.checkpoint(0)

#### READING INPUT FILE ######
### Options 

parser = OptionParser()

parser.add_option("-i", "--infile", dest="infile",
help="file with the input arguments", metavar="INPUTFILE")
(opts, args) = parser.parse_args()

###


inputinf,out_variables=pm.read_input(opts.infile)



##############################

pathin=inputinf['pathin']
pathout=inputinf['pathout']
GCM=inputinf['GCM']
RCM=inputinf['RCM']
syear=int(inputinf['syear'])
eyear=int(inputinf['eyear'])
domain=inputinf['domain']
outfile_patt=inputinf['outfile_patt']
overwrite=True

#CREATE OUTPUT DIR IF IT DOESN'T EXIST
fullpathout='%s/%s/%s/%s-%s/%s' %(pathout,GCM,RCM,syear,eyear,domain)
if not os.path.exists(fullpathout):
	os.makedirs(fullpathout)

#CREATE A TEMPORAL DIR WITHIN THE OUTPUT DIR IF IT DOESN'T EXIST
if not os.path.exists("%s/temp/" %(fullpathout)):
	os.makedirs("%s/temp/" %(fullpathout))


#CREATE A LOG IFLE TO PUT OUTPUT FROM THE MAIN SCRIPT
datenow=dt.datetime.now().strftime("%Y-%m-%d_%H:%M")
logfile = '%s/postprocess_%s_%s_%s-%s_%s_%s.log' %(fullpathout,GCM,RCM,syear,eyear,domain,datenow)
print 'The output messages are written to %s' %(logfile)
#sys.stdout = open('%s' %(logfile), "w") 


#### Reading variable info file ######
varinfo=pm.read_varinfo("./info_files/variables.inf")
file_type=varinfo.keys()

######################################
if GCM=='CCCMA3.1':
	calendar=='noleap'
else:
	calendar='standard'


#***********************************************
# Loop over all types of WRF output files (i.e., wrfhrly, wrfout, etc) 
for filet in file_type:

	if filet=='wrfhrly':
		print '\n', ' PROCESSING ', filet, ' OUTPUTS'
	

		#***********************************************
		# LOOP over years
		for year in np.arange(syear,eyear+1):
			print '\n', ' -> PROCESSING YEAR: ', year
			time_bounds=False

			loadfiles = pathin+'%s_%s_%s*' % (filet,domain,year) # Specify path
			files_in=sorted(glob.glob(loadfiles))

			print '  -->  Number of files to read:', len(files_in)

			# -------------------
			# CHECKING: Check if the number of files is right
			if len(files_in)!=12:
				print '\n', 'ERROR: the number of ',filet, ' files in year ', year,' is INCORRECT'
				print ' ---- SOME FILES ARE MISSING ---'
				print 'SCRIPT stops running ','\n' 
				sys.exit(0)
			
			# READ FILES
			fin=nc.MFDataset(files_in) #Read all files
			time_old = fin.variables['Times'][:] # Get time variable
					
	
			# DEFINE DATES TAKING INTO ACCOUNT IF LEAP-NOLEAP YEAR
			dates_day = dt.datetime(year+1,01,01,00)-dt.datetime(year,01,01,00)
			n_timesteps=dates_day.days*24
			dates = [dt.datetime(year,01,01,00)+ dt.timedelta(hours=x) for x in xrange(0,n_timesteps,1)]

			if calendar=='noleap' and cal.isleap(year)==True:
				months_all=np.asarray([dates[i].month for i in xrange(len(dates))]) 
				days_all=np.asarray([dates[i].day for i in xrange(len(dates))])	
				dates=dates[((months_all==2) & (days_all==29))==False]

			time=nc.date2num(dates[:],units="hours since 1949-12-01 00:00:00")
			

			# -------------------
			# CHECKING: Check if the of time steps is right
			if n_timesteps!=time.shape[0]:
				print '\n', 'ERROR: the number of timesteps in year ', year,' is INCORRECT'
				print 'There should be: ', n_timesteps
				print 'There are: ', time.shape[0]
				print 'SCRIPT stops running ','\n'
				sys.exit(0)


			# ***********************************************
			# LOOP over variables
			for var in out_variables:
				print '  -->  READING VARIABLE ', var
				wrfvar=(pm.getwrfname(var)[0]).split('-')
		

				# ***********************************************************
				# BEFORE READING AND PROCESSING THE VARIABLE OF INTEREST CHECK 
				# IF THE FILE ALREADY EXISTS
				# If it does then go to the next one...
				file_out=pathout+'%s%s_%s-%s_%s.nc' % (outfile_patt,'01H',year,year,var) # Specify output file
				a=os.path.exists(file_out)
				print '  --> OUTPUT FILE:'
				print '                 ', file_out
				if a==True and overwrite==False:
					print '                  +++ FILE ALREADY EXISTS +++'
				else:
					if  a==True and overwrite==True:
						print '                   +++ FILE EXISTS AND WILL BE OVERWRITE +++'
					else:
						print '                   +++ FILE DOES NOT EXISTS YET +++'
				# ***********************************************************


					# ***********************************************
					# LOOP over wrf variables need to compute the variable of interest
					count_v=0
					for wrfv in wrfvar:
						if count_v==0:
							varval=np.array(fin.variables[wrfv][:], dtype='d')
						if count_v==1:
							varval1=np.array(fin.variables[wrfv][:], dtype='d')
						if count_v==2:
							varval2=np.array(fin.variables[wrfv][:], dtype='d')
					fin.close()
					count_v=count_v+1

					ctime=pm.checkpoint(ctime)

					# GET ATTRIBUTES AND MODIFY OUTPUT VARIABLE IF NEEDED
					compute=getattr(cv,'compute_'+var)
					sys.exit(0)
					# CALL COMPUTE_ MODULE
					if var!='pracc':
						result, varatt=compute(varval,dates)
					else:
						# PRECIPITATION NEEDS ONE TIME STEP MORE
						last_file = sorted(glob.glob(pathin+'%s_%s_%s-12-01_*' % (filet,domain,year-1)))
						fin=nc.Dataset(last_file,mode='r')
						count_v=0
						for wrfv in wrfvar:
							if count_v==0:
								varval=np.array(fin.variables[wrfv][:], dtype='d')
							if count_v==1:
								varval1=np.array(fin.variables[wrfv][:], dtype='d')

							result, varatt=compute(varval,last_val,dates)

					# INFO NEEDED TO WRITE THE OUTPUT NETCDF
					nctcdf_info=[file_out, var, varatt, calendar, domain, files_in[0],\
							     GCM, RCM, time_bounds]
					
					# CREATE NETCDF FILE
					aa=pm.create_netcdf(netcdf_info, result, time)

					print aa,'\n'
					ctime=pm.checkpoint(ctime)
					sys.exit(0)
				
		# ***********************************************
		# LOOP over variables
		for var in out_variables:
		
			# ***********************************************
			# LOOP over FREQUENCIES
			freqlist=(varinfo[filet][var]).split(',')
			for freq in freqlist:
				print freq

    
