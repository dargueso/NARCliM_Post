#!/usr/bin/env python

"""compute_vars.py
   Methods to compute the output variables from input variables
"""
import netCDF4 as nc
import numpy as np
import attributes as at
import datetime as dt

def compute_tas(t2,time):
    """Method to compute tas
    t2: T2 from wrf files [K]
    time: list of times corresponding to var 1st dimension
    ---
    tas: output temperature tas [K]
    attributes: attributes of the output variable to be used in the output netcdf 
    """
    if len(time)!=t2.shape[0]:
        sys.exit('ERROR in compute_tas: The lenght of time variable does not correspond to var first dimension')     
    
    #Generating a dictionary with the output attributes of the variable
    tseconds=round(((times[-1]-times[0]).total_seconds()/len(times)))
    atts=pm.get_varatt(sn="air_temperature",ln="Surface air temperature",un="K",ts=tseconds)

    tas=t2[:]    
    return tas,atts
    
    

def compute_pracc(rainc,rainnc,time):
    """Method to compute pracc
    rainc: convective rainfall accumulated [kg m-2]
    rainnc: non-convective rainfall accumulated [kg m-2]
    time: list of times corresponding to rainc 1st dimension
    ---
    pracc: accumulated rainfall (pracc) since last record [kg m-2]
    attributes: attributes of the output variable to be used in the output netcdf
    """
    if (len(time)!=rainc.shape[0]) or (len(time)!=rainnc.shape[0]):
        sys.exit('ERROR in compute_pracc: The lenght of time variable does not correspond to rainc or rainnc first dimension')
    
    #Generating a dictionary with the output attributes of the variable
    tseconds=round(((times[-1]-times[0]).total_seconds()/len(times)))
    atts=pm.get_varatt(sn="precipitation_amount",ln="Accumulated precipitation",un="Kg m-2",ts=tseconds)


    #First hour of each period is set to the same value as second hour
    #### TEMPORARY SOLUTION ######
    # Otherwise the previous file must be loaded too
    pracc=np.zeros(rainc.shape,dtype=np.float64)
    pracc[1:,:,:]=np.diff(rainc+rainnc,axis=0)
    pracc[0,:,:]=pracc[1,:,:]
    
    return pracc,atts
    
def compute_huss(q2,time):
    """Method to compute specific humidity
    q2: mixing ratio [kg kg-2]
    time: list of times corresponding to rainc 1st dimension
    ---
    huss: specific humidity 
    attributes: attributes of the output variable to be used in the output netcdf
    """
    if len(time)!q2.shape[0]:
        sys.exit('ERROR in compute_huss: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((times[-1]-times[0]).total_seconds()/len(times)))
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
    attributes: attributes of the output variable to be used in the output netcdf
    """
    if (len(time)!=u10.shape[0]) or (len(time)!=v10.shape[0]):
        sys.exit('ERROR in compute_wss: The lenght of time variable does not correspond to u10 or v10 first dimension')
    
    tseconds=round(((times[-1]-times[0]).total_seconds()/len(times)))
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
    attributes: attributes of the output variable to be used in the output netcdf
    """
    if len(time)!=u10.shape[0]:
        sys.exit('ERROR in compute_uas: The lenght of time variable does not correspond to var first dimension')
        
    tseconds=round(((times[-1]-times[0]).total_seconds()/len(times)))
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
    attributes: attributes of the output variable to be used in the output netcdf
    """
    if len(time)!=v10.shape[0]:
        sys.exit('ERROR in compute_uas: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((times[-1]-times[0]).total_seconds()/len(times)))
    atts=pm.get_varatt(sn="northward_wind",ln="Northward surface wind",un="m s-1",ts=tseconds)    

    #Zonal wind unstagged
    vas=0.5*(v10[:,:-1,:]+v10[:,1:,:])

    return uas,atts  
    