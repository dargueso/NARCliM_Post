#!/usr/bin/env python

"""compute_stats.py
   Methods to compute daily and monthly statistics from high-frequency variables
"""
import netCDF4 as nc
import numpy as np
import datetime as dt
import sys
import postprocess_modules as pm
def compute_daily(var,time,stat):
	"""Method to compute daily statistics
	   var: variable that will be processed
	   time: list of times corresponding to var 1st dimension
	   stat: requested stat (now 'acc', 'mean', 'min' or 'max')
	   ---
	   dvar: daily statistic of the variable
	   dtime: list of times corresponding to dailyvar 1st dimension
	"""
	if stat not in ['acc', 'mean','min','max']:
		sys.exit("ERROR the requested daily statistic %s does not exist. Please choose between 'acc', 'mean', 'min' or 'max'" %(stat))
	if len(time)!=var.shape[0]:
	    sys.exit('ERROR in compute_daily: The lenght of time variable does not correspond to var first dimension')
	#Remove hours from time list, keep only dates
	ndays=(time[-1]-time[0]).days+1
	sdate=time[0].date()
	stime=dt.datetime(sdate.year,sdate.month,sdate.day,12,0,0)
	dtime=[stime+dt.timedelta(days=x) for x in xrange(ndays)]
	nsteps=len(time)/ndays
	#Remove degenerated dimensions
	var=np.squeeze(var)
	#In the reshape, keep right order for the last two dimensions (lat,lon)
	var_r=np.reshape(var,[nsteps,ndays,var.shape[-2],var.shape[-1]],order='F')
	if stat == 'acc':
	    dvar=np.ma.sum(var_r,axis=0)
	elif stat == 'mean':
	    dvar=np.ma.mean(var_r,axis=0)
	elif stat == 'max':
	    dvar=np.ma.max(var_r,axis=0)
	elif stat == 'min':
	    dvar=np.ma.min(var_r,axis=0)

	return dvar,dtime
    
def compute_monthly(var,time,stat):
	"""Method to compute monthly statistics
	   var: variable that will be processed
	   time: list of times corresponding to var 1st dimension
	   stat: requested stat (now 'acc', 'mean', 'min' or 'max')
	   ---
	   mvar: monthly statistic of the variable
	   mtime: list of times corresponding to mvar 1st dimension
	"""
	if stat not in ['acc', 'mean','min','max']:
		sys.exit("ERROR the requested monthly statistic %s does not exist. Please choose between 'acc', 'mean', 'min' or 'max'" %(stat))
	if len(time)!=var.shape[0]:
		sys.exit('ERROR in compute_monthly: The lenght of time variable does not correspond to var first dimension')
	#Remove hours from time list, keep only dates
	years=np.asarray([time[i].year for i in xrange(len(time))])-time[0].year
	climmonths=np.asarray([time[i].month for i in xrange(len(time))])
	months=years*12+climmonths
	

	var=np.squeeze(var)
	mvar=np.ma.ones((max(months),)+var.shape[1:],dtype=np.float64)*pm.const.missingval

	#Calculating the middle of each month. Data is provided in the mid point between the time_bounds
	mtime=[0]*max(months)
	for mo in xrange(max(months)):
		time_m=time[months==mo+1]
		tdiference=(time_m[-1]-time_m[0]).total_seconds()/2
		mtime[mo]=time_m[0]+dt.timedelta(seconds=tdiference)


	if stat == 'acc':
	    for mo in xrange(max(months)):
	        mvar[mo,:,:]=np.ma.sum(var[months==mo+1,:,:],axis=0)
	if stat == 'mean':
	    for mo in xrange(max(months)):
	        mvar[mo,:,:]=np.ma.mean(var[months==mo+1,:,:],axis=0)
	if stat == 'max':
	    for mo in xrange(max(months)):
	        mvar[mo,:,:]=np.ma.max(var[months==mo+1,:,:],axis=0)
	if stat == 'min':
	    for mo in xrange(max(months)):
	        mvar[mo,:,:]=np.ma.min(var[months==mo+1,:,:],axis=0)


	return mvar,mtime


            
        