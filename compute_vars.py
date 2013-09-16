#!/usr/bin/env python

"""compute_vars.py
   Methods to compute the output variables from input variables
"""
import netCDF4 as nc
import numpy as np
import datetime as dt
import sys
import postprocess_modules as pm

def compute_tas(t2,time):
    """Method to compute 2-m temperature
    t2: T2 from wrf files [K]
    time: list of times corresponding to var 1st dimension
    ---
    tas: output temperature tas [K]
    atts: attributes of the output variable to be used in the output netcdf 
    """
    if len(time)!=t2.shape[0]:
        sys.exit('ERROR in compute_tas: The lenght of time variable does not correspond to var first dimension')     
    
    #Generating a dictionary with the output attributes of the variable
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="air_temperature",ln="Surface air temperature",un="K",ts=tseconds)

    tas=t2[:]    
    return tas,atts
    
    

def compute_pracc(rainc,rainnc,time):
    """Method to compute precipitation
    The original variable was accumulated throughout the simulation. Accumulation must be removed (But it is not converted to flux)
    rainc: convective rainfall accumulated. Including the timestep previous to the first one in this period [kg m-2]
    rainnc: non-convective rainfall accumulated.Including the timestep previous to the first one in this period [kg m-2]
    time: list of times corresponding to rainc 1st dimension
    ---
    pracc: accumulated rainfall (pracc) since last record [kg m-2]
    atts: attributes of the output variable to be used in the output netcdf
    """
    #rainc and rainnc includes the timestep previous to the first one to remove the accumulation.
    if (len(time)!=rainc.shape[0]+1) or (len(time)!=rainnc.shape[0]+1):
        sys.exit('ERROR in compute_pracc: The lenght of time variable does not correspond to rainc or rainnc first dimension')
    
    #Generating a dictionary with the output attributes of the variable
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="precipitation_amount",ln="Accumulated precipitation",un="Kg m-2",ts=tseconds)


    #Calculating difference between each timestep to remove the accumulation
    pracc=np.zeros((rainc.shape[0]-1,)+rainc.shape[1:],dtype=np.float64)
    pracc[:,:,:]=np.diff(rainc+rainnc,axis=0)
    
    return pracc,atts
    
def compute_huss(q2,time):
    """Method to compute specific humidity
    q2: mixing ratio [kg kg-2]
    time: list of times corresponding to rainc 1st dimension
    ---
    huss: specific humidity 
    atts: attributes of the output variable to be used in the output netcdf
    """
    if len(time)!=q2.shape[0]:
        sys.exit('ERROR in compute_huss: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="specific_humidity",ln="Surface specific humidity",un=" ",ts=tseconds)
    
    huss=q2/(1+q2)
    
    return huss,atts

def compute_wss(u10,v10,time):
    """Method to compute wind speed
    u10: zonal wind [m s-1]
    v10: meridional wind [m s-1]

    time: list of times corresponding to rainc 1st dimension
    ---
    wss: wind speed 
    atts: attributes of the output variable to be used in the output netcdf
    """
    if (len(time)!=u10.shape[0]) or (len(time)!=v10.shape[0]):
        sys.exit('ERROR in compute_wss: The lenght of time variable does not correspond to u10 or v10 first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="air_velocity",ln="Surface wind speed",un="m s-1",ts=tseconds)
    
    #Winds need to be unstagged    
    u10_unstagged=0.5*(u10[:,:,:-1]+u10[:,:,1:])
    v10_unstagged=0.5*(v10[:,:-1,:]+v10[:,1:,:])
    
    wss=(u10_unstagged**2+v10_unstagged**2)**0.5
    
    return wss,atts

def compute_uas(u10,time):
    """Method to compute eastward wind
    u10: zonal wind [m s-1]
    time: list of times corresponding to rainc 1st dimension
    ---
    uas: eastward wind 
    atts: attributes of the output variable to be used in the output netcdf
    """
    if len(time)!=u10.shape[0]:
        sys.exit('ERROR in compute_uas: The lenght of time variable does not correspond to var first dimension')
        
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="eastward_wind",ln="Eastward surface wind",un="m s-1",ts=tseconds)    
    
    #Zonal wind unstagged
    uas=0.5*(u10[:,:,:-1]+u10[:,:,1:])
    
    return uas,atts


def compute_uas(v10,time):
    """Method to compute northward wind
    v10: meridional wind [m s-1]
    time: list of times corresponding to rainc 1st dimension
    ---
    vas: northward wind 
    atts: attributes of the output variable to be used in the output netcdf
    """
    if len(time)!=v10.shape[0]:
        sys.exit('ERROR in compute_uas: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="northward_wind",ln="Northward surface wind",un="m s-1",ts=tseconds)    

    #Zonal wind unstagged
    vas=0.5*(v10[:,:-1,:]+v10[:,1:,:])

    return uas,atts

def compute_evspsbl(sfcevp,time):
    """Method to compute surface evaporation flux
       The original variable was accumulated throughout the simulation. Accumulation must be removed.
       sfcevp: accumulated surfave evaporation. Including the timestep previous to the first one in this period [kg m-2]
       ---
       evspsbl: Surface evaporation flux
       atts: attributes of the output variable to be used in the output netcdf
    """
    #sfcevp includes the timestep previous to the first one to remove the accumulation. 
    if len(time)!sfcevp.shape[0]+1:
        sys.exit('ERROR in compute_evspsbl: The lenght of time variable does not correspond to var first dimension')
    
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="water_evaporation_flux",ln="Surface evaporation",un="kg m-2 s-1",ts=tseconds)
    
    #Calculating difference between each timestep to remove the accumulation
    #Divided by the number of seconds in each timestep to calculate the flux
    evspsbl=np.zeros((sfcevp.shape[0]-1,)+sfcevp.shape[1:],dtype=np.float64)
    evspsbl[:,:,:]=np.diff(sfcevp,axis=0)/tseconds
    
    return evspsbl,atts

def compute_mrso():
    