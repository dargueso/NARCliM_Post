#!/usr/bin/env python

"""compute_vars.py
   Methods to compute the output variables from input variables
"""
import netCDF4 as nc
import numpy as np
import attributes as at
import datetime as dt

def compute_tas(var,time):
    """Method to compute tas
    var: T2 from wrf files
    time: list of times corresponding to var 1st dimension
    ---
    varout: output variable
    attributes: attributes of the output variable to be used in the output netcdf 
    """
    if len(time)!=var.shape[0]:
        sys.exit('ERROR in compute_tas: The lenght of time variable does not correspond to var first dimension')     
    
    attribs={}
    attribs[tseconds]=(time[-1]-time[0]).seconds/len(time)
    attribs[standard_name] = "air_temperature"
    attribs[long_name] = "Surface air temperature"
    attribs[units] = "K"
    attribs[coordinates] = "lon lat"
    attribs[cell_method] = "time: point values %s second" %(tseconds)
    attribs[grid_mapping] = "Rotated_pole"
        
    return attribs