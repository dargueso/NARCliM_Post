#!/usr/bin/env python
## 
# python script to launch the ncl script calculate_mslp_arguments.ncl that extract 
# mean sea level pressure from wrfout files using command line
# arguments.
# This script goes through all NARCliM simulations. If the mslp field already 
# exists then it goes to the next one.
#
# The way the ncl script is called is better explained here: 
#    https://wiki.c2sm.ethz.ch/Wiki/VisNCLBasicExamples
#
# Author: Alejandro Di Luca
# Created: 30/04/2014

from netCDF4 import Dataset as ncdf
from netCDF4 import num2date, date2num
import numpy as np
import glob
import os
import copy
import sys		
import ccrc_utils as cu
import postprocess_modules as pm

# Loop parameters
domain='d02'
GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5','CSIRO-MK3.0']
GCM_out={'MIROC3.2':'MIROC3.2','CCCMA3.1':'WRF_CCCMA','ECHAM5':'WRF_ECHAM5','CSIRO-MK3.0':'WRF_MK30'}
RCM_names=['R1','R2','R3']
Periods=['1990-2010','2060-2080']
runcommand='./runncl.sh calculate_mslp_arguments.ncl '

for gind,gname in enumerate(GCM_names):
    for rind,rname in enumerate(RCM_names):
        for pers,pername in enumerate(Periods):
            index_rg=gind*3+rind
            
            indir=cu.get_raw_location(gname,rname,pername)[0]
            print ' -->  INPUT DIRECTORY: ',indir
            
            outdir = "/srv/ccrc/data28/z3444417/Data/WRF/"+gname+"/"+rname+"/"+pername+"/psl/raw/"+domain+"/"
            os.system('mkdir -p '+outdir)
            print '  --> OUTPUT DIRECTORY: ',outdir,'\n'
            for year in np.arange(int(pername.split('-')[0]),int(pername.split('-')[1])):
            	for month in np.arange(1,13):
			month_name=str(month)
	    		if month<10:
				month_name='0'+str(month)
			file_out=outdir+'WRF_mslp_'+rname+'_'+domain+'_'+str(year)+'-'+month_name+'.nc'

            		filewrite=pm.checkfile(file_out,False)
            		if filewrite==True:
                		os.system(runcommand+indir+' '+outdir+' '+str(year)+' '+str(month))
	sys.exit(0)
