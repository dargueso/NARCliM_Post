# ======================================================================
# PURPOSE
# =========
# To calculate tasmean_bc (mean temperature) from tasminmean_bc
# (mean minimum temperature) and tasmaxmean_bc (mean maximum
# temperature) files for the NARCliM project and output the result
# into a similarly formatted NetCDF file. Author and email of the
# resulting file are specified as inputs.

# STATUS
# ========
# Working

# CAUTION! Only expected to be tested on the monthly files.

# NOTE! We change the data type to float64. Regardless of the initial
# _FillValue, the output data is filled with 1E+20. 

# INPUT
# =======
# minfile => NARCliM NetCDF file containing tasminmean_bc. String. 
# maxfile => NARCliM NetCDF file containing tasmaxmean_bc. String. 
# outfile => NetCDF file where to output tasmean_bc. String ending with "*.nc".
# author  => File's author (for global attributes). String.
# email   => Author's email (for global attributes). String. 

# CALLS
# =======
# read_narclim_var and create_netcdf_flexi from netcdf_utils

# OUTPUT
# ========
# outfile => Mean temperature file. NARCliM NetCDF4 file. Minor
# modifications are made to global and variable attributes, which are
# generally copied from the minimum temperature file. 

# REQUIRES
# =========
# Python 2.7.5SL, numpy 1.7.1, netCDF 1.0.4,  netcdf/4.2.1-intel library

# HISTORY
# ========
# Jan 29 2014 => Written by Roman Olson, UNSW.
# Feb 02 2015 => Improved variable and global comments fields.
# ======================================================================
def make_tasmean_bc(minfile, maxfile, outfile, author, email):


   # SOURCE MODULES
   from netcdf_utils import create_netcdf_flexi, read_narclim_var
   import netCDF4 as ncdf
   from datetime import datetime 
   import numpy as np
   import numpy.ma as ma
   import pdb

   
   # OPEN MINFILE
   (nclat1, nclon1, nctime1, tmin) = read_narclim_var(minfile, "tasminmean_bc")
   (nclat2, nclon2, nctime2, tmax) = read_narclim_var(maxfile, "tasmaxmean_bc")

   # CHECKS
   if (np.any(nclat1 != nclat2)):
      raise ValueError, "Latitudes don't match"
   if (np.any(nclon1 != nclon2)):
      raise ValueError, "Longitudes don't match"
   if (np.any(nctime1 != nctime2)):
      raise ValueError, "Times don't match" 

   # OBTAIN MEAN TEMPERATURE AND CONVERT TO MASKED ARRAY
   tasmean_reg = (tmin + tmax)/2.0 
   tasmean_pt = ma.array(tasmean_reg, dtype=np.float64, 
      mask=np.isnan(tasmean_reg), fill_value=1.0E20)

   # CREATE ATTRIBUTES DICTIONARY
   # Globals
   minfileobj              = ncdf.Dataset(minfile)
   minfile_glob            = minfileobj.__dict__.copy()
   minfile_glob["author"]  = author
   minfile_glob["contact"] = email 
   minfile_glob["history"] = "\n On " + datetime.utcnow().isoformat(' ') + \
   " UTC created mean temperature from tasminmean_bc and tasmaxmean_bc." 
   minfile_glob.update(comments="tasmean_bc is not an original NARCliM\
 variable. It is calculated as an average of tasminmean_bc and tasmaxmean_bc\
 variables")

   # Tasmean
   tasmean_attrs = dict(minfileobj.variables["tasminmean_bc"].__dict__.copy())
   del tasmean_attrs["_FillValue"]
   tasmean_attrs.update(history="Created by averaging tasminmean_bc (monthly\
 mean minimum temperature) and tasmaxmean_bc (monthly mean maximum temperature)\
.")
   tasmean_attrs.update(comment="tasmean_bc is not an original NARCliM\
 variable. It is calculated as an average of tasminmean_bc and tasmaxmean_bc\
 variables")
   attrs = {"tasmean_bc":tasmean_attrs, "global":minfile_glob} 

   # WRITE OUTPUT FILE
   create_netcdf_flexi(outfile, tasmean_pt, ["lon", "lat", "time"], 
      "Projection run", varname="tasmean_bc", time=nctime1, lon=nclon1, 
      lat=nclat1, attr_file=minfile, attrs=attrs)

   
