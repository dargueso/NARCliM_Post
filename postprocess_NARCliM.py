#!/usr/bin/env python
"""postprocess_NARCliM.py script
   Script to postprocess WRF outputs from NARCliM project
   It reads an input file (NARClIM_post.input), where the input arguments are provided
	
   Authors: Daniel Argueso (d.argueso@unsw.edu.au), Alejandro Di Luca (a.diluca@unsw.edu.au)
   Institution: CoECSS, UNSW. CCRC, UNSW.
   Created: 13 September 2013
   Modified: 17 Septemvber 2013
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



##############################

pathin=inputinf['pathin']
pathout=inputinf['pathout']
GCM=inputinf['GCM']
RCM=inputinf['RCM']
syear=int(inputinf['syear'])
eyear=int(inputinf['eyear'])
domain=inputinf['domain']
outfile_patt=inputinf['outfile_patt']
overwrite=False


### Reference file (attributes and other) #########

fileref_att=pathin+'/wrfhrly_%s_%s-01-01_00:00:00' %(domain,syear)
sch_info=pm.read_schemes(fileref_att) #Information about the physics schemes used in the simulation

### Time units (reference time for netcdf time units) ###########

time_units="hours since 1949-12-01 00:00:00"

#CREATE OUTPUT DIR IF IT DOESN'T EXIST
fullpathout='%s%s/%s/%s-%s/%s/' %(pathout,GCM,RCM,syear,eyear,domain)
if not os.path.exists(fullpathout):
	os.makedirs(fullpathout)

#CREATE A TEMPORAL DIR WITHIN THE OUTPUT DIR IF IT DOESN'T EXIST
if not os.path.exists("%stemp/" %(fullpathout)):
	os.makedirs("%stemp/" %(fullpathout))


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
# LOOP OVER ALL TYPES OF WRF FILE OUTPUTS (i.e., wrfhrly, wrfout, etc) 
for filet in file_type:
	ctime_filet=pm.checkpoint(0)
	filet='wrfout'
	filet='wrfhrly'
	
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
	for year in np.arange(syear,eyear+1):
		ctime_year=pm.checkpoint(0)
		if filet=='wrfout':
			n_files=365
			if cal.isleap(year):
				n_files=perstep*366

				
		# SELECTING FILES TO READ
		print '\n', ' -> PROCESSING PERIOD: ', str(year)+' - '+str(year)
		loadfiles = pathin+'%s_%s_%s*' % (filet,domain,year) # Specify path
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
				file_out=fullpathout+'/%s%s_%s-%s_%s.nc' % (outfile_patt,file_freq,year,year,var) # Specify output file
				filewrite=pm.checkfile(file_out,overwrite)
				if filewrite==True:

					# READ FILES FROM THE CORRESPONDING PERIOD
					print '    -->  READING FILES '
					fin=nc.MFDataset(files_list) # Read all files
					print '    -->  EXTRACTING VARIABLE Time'
					time_old = fin.variables['Times'][:] # Get time variable

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
					time=nc.date2num(date[:],units=time_units)

					# -------------------
					# CHECKING: Check if the of time steps is right
					if n_timesteps!=time_old.shape[0]:
						print '\n', 'ERROR: the number of timesteps in period ', year,' is INCORRECT'
						print 'There should be: ', n_timesteps
						print 'There are: ', time_old.shape[0]
						print 'SCRIPT stops running ','\n'
						sys.exit(0)


					# ***********************************************
					# LOOP over wrf variables need to compute the variable of interest
					print '    -->  EXTRACTING VARIABLE ', var
					wrfvar=(pm.getwrfname(var)[0]).split('-')
					for cv,wrfv in enumerate(wrfvar):
						if cv==0:
							varval=np.array(fin.variables[wrfv][:], dtype='d')
						if cv==1:
							varval1=np.array(fin.variables[wrfv][:], dtype='d')
						if cv==2:
							varval2=np.array(fin.variables[wrfv][:], dtype='d')


					# ***********************************************
					# FOR LEAP YEARS AND MODELS WITH NO LEAP YEARS REPLACE THE 29TH FEBRUARY BY 
					# MISSING VALUES
					if calendar=='noleap' and cal.isleap(year)==True:
						months_all=np.asarray([date[i].month for i in xrange(len(date))]) 
						days_all=np.asarray([date[i].day for i in xrange(len(date))])
						index=np.where(np.logical_and(months_all==2,days_all==29))

						for cv,wrfv in enumerate(wrfvar):
							if cv==0:
								varval = pm.add_leap(varval,index)
							if cv==1:
								varval1 = pm.add_leap(varval1,index)
							if cv==2:
								varval2 = pm.add_leap(varval2,index)
					fin.close() # CLOSE FILES

					
					# ***********************************************
					# ACCUMULATED VARIABLES NEED ONE TIME STEP MORE TO COMPUTE DIFFERENCES
					if var=='pracc' or var=='potevp' or var=='evspsbl':

						# DEFINE TIME BOUNDS FOR ACCUMULATED VARIABLES
						time_bounds=True
						if filet=='wrfhrly' or filet=='wrfout':
							date_inf=date
							time_inf=nc.date2num(date_inf[:],units=time_units)
							date_sup=pm.get_dates(year,month_i,day_i,hour_i+time_step,mins,time_step,n_timesteps)
							time_sup=nc.date2num(date_sup[:],units=time_units)
							time_bnds=np.reshape(np.concatenate([time_inf,time_sup]), (time.shape[0],2))

							# Time variable is defined in the middle of time bounds
							date=pm.get_dates(year,month_i,day_i,hour_i,30,time_step,n_timesteps)
							time=nc.date2num(date[:],units=time_units)

						# READ ONE MORE TIME STEP FOR ACCUMULATED COMPUTATIONS
						if year<2009:
							last_file = sorted(glob.glob(pathin+'%s_%s_%s-01-01_*' % (filet,domain,year+1)))
							fin2=nc.Dataset(last_file[0],mode='r')
							for cv,wrfv in enumerate(wrfvar):
								if cv==0:
									last_varval=np.reshape(np.array(fin2.variables[wrfv][0,:,:], \
														dtype='d'),(1,varval.shape[1],varval.shape[2]))
									varval=np.concatenate((varval,last_varval))
								if cv==1:
									last_varval=np.reshape(np.array(fin2.variables[wrfv][0,:,:], \
														dtype='d'),(1,varval.shape[1],varval.shape[2]))
									varval1=np.concatenate((varval1,last_varval))
							fin2.close()
						else:
							last_varval=np.zeros((1,varval.shape[1],varval.shape[2]))
							last_varval[:]=pm.const.missingval
							for cv,wrfv in enumerate(wrfvar):
								if cv==0:
									varval=np.concatenate((varval,last_varval))
								if cv==1:
									varval1=np.concatenate((varval1,last_varval))
						

					# ***********************************************
					# DEFINE TIME BOUNDS FOR XTRM AND DAILY VARIABLES
					if filet=='wrfxtrm' or filet=='wrfdly':
						datei=date[0]-dt.timedelta(hours=int(float(time_step)/2.))
						date_inf=[datei+dt.timedelta(hours=x) for x in xrange(0,n_timesteps*time_step,time_step)]
						time_inf=nc.date2num(date_inf[:],units=time_units)
						datei=date[0]+dt.timedelta(hours=int(float(time_step)/2.))
						date_sup=[datei+dt.timedelta(hours=x) for x in xrange(0,n_timesteps*time_step,time_step)]
						time_sup=nc.date2num(date_sup[:],units=time_units)
						time_bnds=np.reshape(np.concatenate([time_inf,time_sup]), (time.shape[0],2))


					# CALL COMPUTE_VAR MODULE
					compute=getattr(comv,'compute_'+var) # FROM STRING TO ATTRIBUTE
					if len(wrfvar)==1:
						varval, varatt=compute(varval,date)
					if len(wrfvar)==2:
						varval, varatt=compute(varval,varval1,date)
					if len(wrfvar)==3:
						varval, varatt=compute(varval,varval1,varval2,date)

					# INFO NEEDED TO WRITE THE OUTPUT NETCDF
					netcdf_info=[file_out, var, varatt, 'standard', domain, GCM, RCM, time_bounds]

					# CREATE NETCDF FILE
					pm.create_netcdf(netcdf_info, varval, time, time_bnds, sch_info,time_units)
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
				print '%s/%s0?H_*_%s.nc' %(fullpathout,outfile_patt,varname)
				fileall=sorted(glob.glob('%s/%s0?H_*_%s.nc' %(fullpathout,outfile_patt,varname)))
				eyfile=np.asarray([int(fileall[i].split('_%s.nc' %(varname))[0][-4:]) for i in xrange(len(fileall))])
				syfile=np.asarray([int(fileall[i].split('_%s.nc' %(varname))[0][-9:-5]) for i in xrange(len(fileall))])

				fileref=nc.Dataset(fileall[0],'r')

				syp=syear
				while syp<eyear:
					ctime_var=pm.checkpoint(0)
					eyp=((int(syp)/5)+1)*5
					loadfile=False
					for stat in stat_all:
						if varname=='pracc':
							varstat=varname
						else:
							varstat=varname+stat
						file_out=fullpathout+'/%sDAY_%s-%s_%s.nc' % (outfile_patt,syp,eyp-1,varstat) # Specify output file
						filewrite=pm.checkfile(file_out,overwrite)
						if filewrite==True:
							loadfile=True
					if loadfile==True:

						print 'PROCESSING PERIOD %s-%s for variable %s' %(syp,eyp,varname)
						sel_files=[fileall[i] for i in xrange(len(syfile)) if ((syfile[i]>=syp) & (eyfile[i]<eyp))]
						files=nc.MFDataset(sel_files)
						time=nc.num2date(files.variables['time'][:],units=files.variables['time'].units)
						var=files.variables[varname][:]
						ctime=pm.checkpoint(ctime_var)

						for stat in stat_all:
							ctime_var=pm.checkpoint(0)
							if varname=='pracc':
								varstat=varname
							else:
								varstat=varname+stat
							file_out=fullpathout+'/%sDAY_%s-%s_%s.nc' % (outfile_patt,syp,eyp-1,varstat) # Specify output file
							filewrite=pm.checkfile(file_out,overwrite)
							if filewrite==True:
								dvar,dtime=coms.compute_daily(var,time,stat)
								dtime_nc=nc.date2num(dtime,units=time_units)
								time_bnds_inf=dtime_nc-12
								time_bnds_sup=dtime_nc+12
								time_bnds=np.reshape(np.concatenate([time_bnds_inf,time_bnds_sup]), (len(dtime),2))
								varatt={}
								for att in fileref.variables[varname].ncattrs():
									varatt[att]=getattr(fileref.variables[varname],att)

								# INFO NEEDED TO WRITE THE OUTPUT NETCDF
								netcdf_info=[file_out, varstat, varatt, calendar, domain, GCM, RCM, True]

								# CREATE NETCDF FILE
								pm.create_netcdf(netcdf_info, dvar, dtime_nc, time_bnds, sch_info,time_units)
								ctime=pm.checkpoint(ctime_var)
								print '=====================================================', '\n', '\n', '\n'
					syp=eyp
					
					
#***********************************************
# MONTHLY STATISTICS
for filet in file_type:
	for varname in varinfo[filet].keys():
		if varname in out_variables:
			stat_all=varinfo[filet][varname].split(',')
			for stat in stat_all:
				if varname=='pracc':
					varstat=varname
				else:
					varstat=varname+stat
				
				fileall=sorted(glob.glob('%s/%sDAY_*_%s.nc' %(fullpathout,outfile_patt,varstat)))
				print np.asarray([(fileall[i].split('_%s.nc' %(varstat))[0][-4:]) for i in xrange(len(fileall))])
				eyfile=np.asarray([int(fileall[i].split('_%s.nc' %(varstat))[0][-4:]) for i in xrange(len(fileall))])
				syfile=np.asarray([int(fileall[i].split('_%s.nc' %(varstat))[0][-9:-5]) for i in xrange(len(fileall))])

				fileref=nc.Dataset(fileall[0],'r')

				syp=syear
				while syp<eyear:
					print varstat
					ctime_var=pm.checkpoint(0)
					eyp=((int(syp)/10)+1)*10
					print 'PROCESSING PERIOD %s-%s for variable %s' %(syp,eyp,varname)
					sel_files=[fileall[i] for i in xrange(len(syfile)) if ((syfile[i]>=syp) & (eyfile[i]<eyp))]
					files=nc.MFDataset(sel_files)
					time=nc.num2date(files.variables['time'][:],units=files.variables['time'].units)
					var=files.variables[varstat][:]
					ctime=pm.checkpoint(ctime_var)
					file_out=fullpathout+'/%sMON_%s-%s_%s.nc' % (outfile_patt,syp,eyp-1,varstat) # Specify output file
					filewrite=pm.checkfile(file_out,overwrite)
					if filewrite==True:
						mvar,mtime=coms.compute_monthly(var,time,stat)
						mtime_nc=nc.date2num(mtime,units=time_units)
						mtime_bnds_inf=nc.date2num([dt.datetime(mtime[i].year,mtime[i].month,1,0,0,0) for i in xrange(len(mtime))],units=time_units)
						mtime_bnds_sup=nc.date2num([dt.datetime(mtime[i].year,mtime[i].month,1,0,0,0)+relativedelta(months=1) for i in xrange(len(mtime))],units=time_units)
						mtime_bnds=np.reshape(np.concatenate([mtime_bnds_inf,mtime_bnds_sup],axis=0), (len(mtime),2))
						print mtime_bnds
						varatt={}
		
						for att in fileref.variables[varstat].ncattrs():
							varatt[att]=getattr(fileref.variables[varstat],att)
						
							# INFO NEEDED TO WRITE THE OUTPUT NETCDF
							netcdf_info=[file_out, varstat, varatt, calendar, domain, GCM, RCM, True]

							# CREATE NETCDF FILE
							pm.create_netcdf(netcdf_info, mvar, mtime_nc, mtime_bnds, sch_info,time_units)
							ctime=pm.checkpoint(ctime_var)
							print '=====================================================', '\n', '\n', '\n'
					syp=eyp    
