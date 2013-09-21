#!/usr/bin/env python
"""postprocess_NARCliM.py script
   Script to postprocess WRF outputs from NARCliM project
   It reads an input file (NARClIM_post.input), where the input arguments are provided
	
   Authors: Daniel Argueso (d.argueso@unsw.edu.au), Alejandro Di Luca (a.diluca@unsw.edu.au)
   Institution: CoECSS, UNSW. CCRC, UNSW.
   Created: 13 September 2013
   Modified: 17 September 2013
"""

import netCDF4 as nc
import numpy as np
import sys
import os 
import datetime as dt
import glob
from optparse import OptionParser
import postprocess_modules as pm
import compute_vars as comv
import calendar as cal
import compute_stats as coms
from dateutil.relativedelta import relativedelta

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
overwrite=False
gvars=pm.gvar(inputinf)

fullpathout=pm.create_outdir(gvars)
 
#CREATE A LOG IFLE TO PUT OUTPUT FROM THE MAIN SCRIPT
datenow=dt.datetime.now().strftime("%Y-%m-%d_%H:%M")
logfile = '%s/postprocess_%s_%s_%s-%s_%s_%s.log' %(fullpathout,gvars.GCM,gvars.RCM,gvars.syear,gvars.eyear,gvars.domain,datenow)
print 'The output messages are written to %s' %(logfile)
#sys.stdout = open('%s' %(logfile), "w") 


#### Reading variable info file ######
varinfo=pm.read_varinfo("./info_files/variables.inf")
file_type=varinfo.keys()
#***********************************************
# LOOP OVER ALL TYPES OF WRF FILE OUTPUTS (i.e., wrfhrly, wrfout, etc) 
for filet in file_type:
	ctime_filet=pm.checkpoint(0)
	
	print '\n','\n', '*************************************'
	print '  PROCESSING ', filet, ' FILE OUTPUTS'
	print '*************************************'
	
	perstep=1

	if filet=='wrfhrly':
		n_files=12	
		time_step=1 #hours between two time steps
		file_freq='01H'
		tbounds=False

	if filet=='wrfout':
		time_step=3 #hours between two time steps
		file_freq='03H'
		tbounds=False

	if filet=='wrfxtrm':
		n_files=12	
		time_step=24 #hours between two time steps
		file_freq='DAY'
		tbounds=True

	if filet=='wrfdly':
		n_files=12	
		time_step=24 #hours between two time steps
		file_freq='DAY'
		tbounds=True

	#***********************************************
	# LOOP OVER YEARS
	for year in np.arange(gvars.syear,gvars.eyear+1):
		ctime_year=pm.checkpoint(0)
		if filet=='wrfout':
			n_files=365
			if cal.isleap(year):
				n_files=perstep*366

				
		# SELECTING FILES TO READ
		print '\n', ' -> PROCESSING PERIOD: ', str(year)+' - '+str(year)
		loadfiles = '%s%s_%s_%s*' % (gvars.pathin,filet,gvars.domain,year) # Specify path
		files_in=np.array(sorted(glob.glob(loadfiles)))
		print '  -->  Number of files to read:', files_in.shape[0]
		files_list=list(files_in)

		# -------------------
		# CHECKING: Check if the number of files is right
		if len(files_in)!=n_files:
			print '\n', 'ERROR: the number of ',filet, ' files in period ', year,' is INCORRECT'
			print ' ---- SOME FILES ARE MISSING ---'
			print 'SCRIPT stops running ','\n' 
			sys.exit(0)

		# ***********************************************
		# LOOP OVER VARIABLES IN THE GIVEN KIND OF FILE
		for var in varinfo[filet].keys():
			if var in out_variables:
				ctime_var=pm.checkpoint(0)

				# ***********************************************************
				# BEFORE READING AND PROCESSING THE VARIABLE OF INTEREST CHECK 
				# IF THE FILE ALREADY EXISTS
				# If it does then go to the next one...
				file_out='%s/%s%s_%s-%s_%s.nc' % (fullpathout,gvars.outfile_patt,file_freq,year,year,var) # Specify output file
				filewrite=pm.checkfile(file_out,overwrite)
				if filewrite==True:
					
					
					
					# # READ FILES FROM THE CORRESPONDING PERIOD
					# print '    -->  READING FILES '
					# fin=nc.MFDataset(files_list) # Read all files
					# print '    -->  EXTRACTING VARIABLE Time'
					# time_old = fin.variables['Times'][:] # Get time variable

					# # FIRST YEAR, MONTH, DAY AND HOUR OF ALL READ FILES
					# year_i, month_i, day_i, hour_i = pm.get_wrfdate(time_old[0,:])
					# mins=0

					# # LAST YEAR, MONTH, DAY AND HOUR OF ALL READ FILES
					# year_f, month_f, day_f, hour_f = pm.get_wrfdate(time_old[-1,:])

					# # DEFINE TIME BOUNDS VARIABLES 
					# time_bounds=tbounds
					# time_bnds=pm.const.missingval

					# # DEFINE DATES USING STANDARD CALENDAR
					# n_days = dt.datetime(year+1,month_i,day_i,hour_i)-dt.datetime(year,month_i,day_i,hour_i)
					# n_timesteps=n_days.days*int(24./time_step)
					# date = pm.get_dates(year,month_i,day_i,hour_i,mins,time_step,n_timesteps)
					# time=pm.date2hours(date,gvars.ref_date)

					# # -------------------
					# # CHECKING: Check if the of time steps is right
					# if n_timesteps!=time_old.shape[0]:
					# 	print '\n', 'ERROR: the number of timesteps in period ', year,' is INCORRECT'
					# 	print 'There should be: ', n_timesteps
					# 	print 'There are: ', time_old.shape[0]
					# 	print 'SCRIPT stops running ','\n'
					# 	sys.exit(0)


					# ***********************************************
					# LOOP over wrf variables need to compute the variable of interest
					print '    -->  EXTRACTING VARIABLE ', var
					wrfvar=(pm.getwrfname(var)[0]).split('-')
					
					varvals,time_old=pm.get_wrfvars(wrfvar,files_list)
					# FIRST YEAR, MONTH, DAY AND HOUR OF ALL READ FILES
					year_i, month_i, day_i, hour_i = pm.get_wrfdate(time_old[0,:])
					mins=0
					# LAST YEAR, MONTH, DAY AND HOUR OF ALL READ FILES
					year_f, month_f, day_f, hour_f = pm.get_wrfdate(time_old[-1,:])
					# DEFINE TIME BOUNDS VARIABLES 
					time_bounds=tbounds
					time_bnds=pm.const.missingval
					# DEFINE DATES USING STANDARD CALENDAR
					n_days = dt.datetime(year+1,month_i,day_i,hour_i)-dt.datetime(year,month_i,day_i,hour_i)
					n_timesteps=n_days.days*int(24./time_step)
					date = pm.get_dates(year,month_i,day_i,hour_i,mins,time_step,n_timesteps)
					time=pm.date2hours(date,gvars.ref_date)
					# -------------------
					# CHECKING: Check if the of time steps is right
					if n_timesteps!=time_old.shape[0]:
						print '\n', 'ERROR: the number of timesteps in period ', year,' is INCORRECT'
						print 'There should be: ', n_timesteps
						print 'There are: ', time_old.shape[0]
						print 'SCRIPT stops running ','\n'
						sys.exit(0)
					
					if gvars.GCM_calendar=='noleap' and cal.isleap(year)==True:
						varvals=add_leapdays(varvals)

					
					
					# ***********************************************
					# ACCUMULATED VARIABLES NEED ONE TIME STEP MORE TO COMPUTE DIFFERENCES
					if var=='pracc' or var=='potevp' or var=='evspsbl':

						# DEFINE TIME BOUNDS FOR ACCUMULATED VARIABLES
						time_bounds=True
						if filet=='wrfhrly' or filet=='wrfout':
							time=pm.create_outtime(date,gvars)
							time_bnds=pm.create_timebnds(time)
							varvals=pm.add_timestep_acc(wrfvar,varvals,year,gvars,filet)

					# ***********************************************
					# DEFINE TIME BOUNDS FOR XTRM AND DAILY VARIABLES
					if filet=='wrfxtrm' or filet=='wrfdly':
						time=pm.date2hours(date,gvars.ref_date)
						time_bnds=pm.create_timebnds(time)


					# CALL COMPUTE_VAR MODULE
					compute=getattr(comv,'compute_'+var) # FROM STRING TO ATTRIBUTE
					varval, varatt=compute(varvals,date)

					# INFO NEEDED TO WRITE THE OUTPUT NETCDF
					netcdf_info=[file_out, var, varatt, time_bounds]

					# CREATE NETCDF FILE
					pm.create_netcdf(netcdf_info, gvars, varval, time, time_bnds)
					ctime=pm.checkpoint(ctime_var)
					print '=====================================================', '\n', '\n', '\n'
	print ' =======================  YEAR:',year, ' FINISHED ==============', '\n', '\n',
	ctime=pm.checkpoint(ctime_year)
print ' =======================  FILE TYPE :',filet, ' FINISHED ==============', '\n', '\n',
ctime=pm.checkpoint(ctime_filet)


#***********************************************
# DAILY STATISTICS
# Loop over all types of WRF output files (i.e., wrfhrly, wrfout, etc) 
for filet in file_type:
	if (filet!='wrfxtrm') and (filet!='wrfdly'): # These files are already daily
		for varname in varinfo[filet].keys():
			if varname in out_variables:
				stat_all=varinfo[filet][varname].split(',')
				pm.create_dailyfiles(gvars,varname,stat_all)
					
					
#***********************************************
# MONTHLY STATISTICS
for filet in file_type:
	for varname in varinfo[filet].keys():
		if varname in out_variables:
			stat_all=varinfo[filet][varname].split(',')
			pm.create_monthlyfiles(gvars,varname,stat_all)    
