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


#### Reading input info file ######
inputinf,out_variables=pm.read_input(opts.infile)

#### Reading variable info file ######
varinfo=pm.read_varinfo("./variables.inf")
file_type=varinfo.keys()

#### Creating global variables ####
gvars=pm.gvar(inputinf)
fullpathout=pm.create_outdir(gvars)
 
#CREATE A LOG IFLE TO PUT OUTPUT FROM THE MAIN SCRIPT
datenow=dt.datetime.now().strftime("%Y-%m-%d_%H:%M")
logfile = '%spostprocess_%s_%s_%s-%s_%s_%s.log' %(fullpathout,gvars.GCM,gvars.RCM,gvars.syear,gvars.eyear,gvars.domain,datenow)
print 'The output messages are written to %s' %(logfile)
#sys.stdout = open('%s' %(logfile), "w") 

#***********************************************
# LOOP OVER ALL TYPES OF WRF FILE OUTPUTS (i.e., wrfhrly, wrfout, etc) 
for filet in file_type:
	ctime_filet=pm.checkpoint(0)
	print '\n','\n', '*************************************'
	print '  PROCESSING ', filet, ' FILE OUTPUTS'
	print '*************************************'
	sper=gvars.syear
	eper=gvars.eyear
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
		perstep=5

	if filet=='wrfxtrm':
		n_files=60      
		time_step=24 #hours between two time steps
		file_freq='DAY'
		tbounds=True
		perstep=5

	if filet=='wrfdly':
		n_files=60      
		time_step=24 #hours between two time steps
		file_freq='DAY'
		tbounds=True
		perstep=5

        # LOOP OVER PERIODS
	for per in np.arange(sper,eper+1,perstep):
		ctime_year=pm.checkpoint(0)
		per_f=per+perstep-1

		n_leap=0
		for pp in np.arange(per,per_f+1):
			if cal.isleap(pp):
				n_leap=n_leap+1	
			if filet=='wrfout':
				if gvars.GCM_calendar!='no_leap':
					n_files=perstep*365+n_leap
				else:
					n_files=perstep*365
				

		# SELECTING FILES TO READ
		files_list=pm.file_list(gvars, per, per_f, filet, n_files)

		# LOOP OVER VARIABLES IN THE GIVEN KIND OF FILE
		for var in pm.intersect(varinfo[filet].keys(),out_variables):
			ctime_var=pm.checkpoint(0)
			
			# CHECK IF THE FILE ALREADY EXISTS
			file_out='%s%s%s_%s-%s_%s.nc' % (fullpathout,gvars.outfile_patt,file_freq,per,per_f,var) # Specify output file
			filewrite=pm.checkfile(file_out,gvars.overwrite)
			if filewrite==True:
				
				# READ FILES FROM THE CORRESPONDING PERIOD
				wrfvar=(pm.getwrfname(var)[0]).split('-')
				time_old, varvals=pm.read_list(files_list, var)	

				# FIRST/LAST YEAR, MONTH, DAY AND HOUR OF ALL READ FILES
				year_i, month_i, day_i, hour_i = pm.get_wrfdate(time_old[0,:])
				year_f, month_f, day_f, hour_f = pm.get_wrfdate(time_old[-1,:])

				# DEFINE TIME BOUNDS VARIABLES 
				time_bounds=tbounds
				time_bnds=pm.const.missingval

				# DEFINE DATES USING STANDARD CALENDAR
				n_days = dt.datetime(per_f+1,month_i,day_i,hour_i)-dt.datetime(per,month_i,day_i,hour_i)
				n_timesteps=n_days.days*int(24./time_step)
				date = pm.get_dates(year_i,month_i,day_i,hour_i,0,time_step,n_timesteps)
				time=pm.date2hours(date,gvars.ref_date)

				# ADD LEAP DAY FOR MODELS WITHOU IT 
				if gvars.GCM_calendar=='no_leap' and n_leap>=1:
					varvals=pm.add_leapdays(varvals,wrfvar,date,time_step)

				# ***********************************************
				# ACCUMULATED VARIABLES NEED ONE TIME STEP MORE TO COMPUTE DIFFERENCES
				if var=='pracc' or var=='potevp' or var=='evspsbl':

					# DEFINE TIME BOUNDS FOR ACCUMULATED VARIABLES
					time_bounds=True
					if filet=='wrfhrly' or filet=='wrfout':
						time=pm.create_outtime(date,gvars)
						time_bnds=pm.create_timebnds(time)
						varvals=pm.add_timestep_acc(wrfvar,varvals,per_f,gvars,filet)

				# ***********************************************
				# DEFINE TIME BOUNDS FOR XTRM AND DAILY VARIABLES
				if filet=='wrfxtrm' or filet=='wrfdly':
					time=pm.date2hours(date,gvars.ref_date)
					time=[time[i]+time_step/2 for i in xrange(len(time))]
					time_bnds=pm.create_timebnds(time)
					varvals=pm.mv_timestep(wrfvar,varvals,per_f,gvars,filet)

				# CALL COMPUTE_VAR MODULE
				compute=getattr(comv,'compute_'+var) # FROM STRING TO ATTRIBUTE
				varval, varatt=compute(varvals,date,gvars)

				# INFO NEEDED TO WRITE THE OUTPUT NETCDF
				netcdf_info=[file_out, var, varatt, time_bounds]

				# CREATE NETCDF FILE
				pm.create_netcdf(netcdf_info, gvars, varval, time, time_bnds)
				ctime=pm.checkpoint(ctime_var)
				print '=====================================================', '\n', '\n', '\n'
			print ' =======================  PERIOD: ',per, ' - ', per_f, ' FINISHED ==============', '\n', '\n',
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
