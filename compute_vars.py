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
    time: list of times corresponding to t2 1st dimension
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
    time: list of times corresponding to q2 1st dimension
    ---
    huss: specific humidity []
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
    time: list of times corresponding to u10 and v10 1st dimensions
    ---
    wss: wind speed [m s-1]
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
    time: list of times corresponding to u10 1st dimension
    ---
    uas: eastward wind [m s-1]
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
    time: list of times corresponding to v10 dimension
    ---
    vas: northward wind [m s-1]
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
       sfcevp: accumulated surface evaporation. Including the timestep previous to the first one in this period [kg m-2]
       time: list of times corresponding to sfcevp 1st dimension
       ---
       evspsbl: Surface evaporation flux [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    #sfcevp includes the timestep previous to the first one to remove the accumulation. 
    if len(time)!=sfcevp.shape[0]+1:
        sys.exit('ERROR in compute_evspsbl: The lenght of time variable does not correspond to var first dimension')
    
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="water_evaporation_flux",ln="Surface evaporation",un="kg m-2 s-1",ts=tseconds)
    
    #Calculating difference between each timestep to remove the accumulation
    #Divided by the number of seconds in each timestep to calculate the flux
    evspsbl=np.zeros((sfcevp.shape[0]-1,)+sfcevp.shape[1:],dtype=np.float64)
    evspsbl[:,:,:]=np.diff(sfcevp,axis=0)/tseconds
    
    return evspsbl,atts

def compute_mrso(smstot,dzs,time):
    """Method to compute the total soil moisture content
       Integrates through all soil layers
       smstot: total soil moisture in each layer [m3 m-3]
       dzs: thickness of each soil layer [m]
       time: list of times corresponding to smstot 1st dimension
       ---
       mrso: total soil moisture content [kg m-2]
       atts: attributes of the output variable to be used in the output netcdf
    """
    if len(time)!=smstot.shape[0]:
        sys.exit('ERROR in compute_mrso: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="soil_moisture_content",ln="Total soil moisture content",un="kg m-2",ts=tseconds)
    
    mrso=smstot*1000*np.sum(dzs)
    
    return mrso,atts

def compute_sst(sst_in,time):
    """Method to compute the sea surface temperature
       sst_in: sea surface temperature [K]
       time: list of times corresponding to sst 1st dimension
       ---
       sst_out: sea surface temperature [K]
       atts: attributes of the output variable to be used in the output netcdf
    """
    if len(time)!=sst_in.shape[0]:
        sys.exit('ERROR in compute_sst: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="sea_surface_temperature",ln="sea surface temperature",un="K",ts=tseconds) 
    
    sst_out=sst_in
    
    return sst_out,atts
    
def compute_potevp(potevp_in,time):
    """Method to compute surface potential evaporation flux
       The original variable was accumulated throughout the simulation. Accumulation must be removed.
       potevp_in: accumulated surface potential evaporation. Including the timestep previous to the first one in this period [kg m-2]
       time: list of times corresponding to potevp 1st dimension
       ---
       potevp_out: Surface evaporation flux [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    #potevp includes the timestep previous to the first one to remove the accumulation. 
    if len(time)!=potevp_in.shape[0]+1:
        sys.exit('ERROR in compute_potevp: The lenght of time variable does not correspond to var first dimension')    
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="water_potential_evaporation_flux",ln="Potential evaporation",un="W m-2",ts=tseconds)
    
    #Calculating difference between each timestep to remove the accumulation
    #Divided by the number of seconds in each timestep to calculate the flux
    potevp_out=np.zeros((potevp_in.shape[0]-1,)+potevp_in.shape[1:],dtype=np.float64)
    potevp_out[:,:,:]=np.diff(potevp_in,axis=0)/tseconds
    
    return potevp_out,atts  
    
      