"""add_global_att_postprocess_files.py

This script adds new global attributes to ALL postprocessed files.
Four new global attributes are added:
 - CF convention version
 - Comment on conventions
 - version number
 - version date in the repo
 

Author: Alejandro Di Luca @ CCRC, UNSW. Sydney (Australia)
email: a.diluca@unsw.edu.au
Created: 17/03/2014
Modified: 13/06/2014
    - I added the version number and I added variable for the global attributes and their values.

"""
import netCDF4 as nc
import numpy as np
import sys
import os
import glob
import ccrc_utils as cu


GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5','CSIRO-MK30']
RCM_names=['R1','R2','R3']
Period_names=['1990-2010','2020-2040','2060-2080']
Domain_names=['d01','d02']

# Names and values of global att to be added
gatt_name=['Conventions','Comments_on_Conventions','version']
gatt_value={'Conventions':'CF-1.6',\
              'Comments_on_Conventions':'Some variables do not strictly follow CF-1.6 conventions (e.g., precipitation is given as an accumulated quantity rather than as a flux quantity).',\
'version':'2.0'}


for g in xrange(len(GCM_names)):
  for r in xrange(len(RCM_names)):
    for p in xrange(len(Period_names)):
      for d in xrange(len(Domain_names)):

        print "PROCESSING SIMULATION: %s - %s - %s - %s" \
            %(GCM_names[g],RCM_names[r],Period_names[p],Domain_names[d])
        
        # Retrieving the location where postprocessed data are stored
        pathin=cu.get_postproc_location(GCM_names[g],RCM_names[r],Period_names[p])[0]

        # Output files are going to the same directory...
        pathout=pathin

        # Check if we have permissions to write in output folder
        if os.access(pathin, os.W_OK)==True:

          # List files to modify
          files_in=sorted(glob.glob('%s/%s/*.nc' % (pathin,Domain_names[d])))

          for file in files_in:
              print ' --->> Adding attributes in: ', file

              # Open and read file
              fout = nc.Dataset(file,mode='a')

              # Read global atts in the file
              globatt={}
              for attname in fout.ncattrs():
                globatt[attname]=getattr(fout, attname)

              for gatt in gatt_name:

                # Check if the attribute already exists
                if np.any(np.asarray(globatt.keys())==gatt)==False:
                  setattr(fout,gatt,gatt_value[gatt])   
              
              fout.close()

