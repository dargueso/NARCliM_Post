#!/usr/bin/env python

"""compute_vars.py
   Methods to compute the output variables from input variables
"""
import netCDF4 as nc
import numpy as np
import datetime as dt
import sys
import postprocess_modules as pm

def compute_tas(varvals,time,gvars):
    """Method to compute 2-m temperature
    t2: T2 from wrf files [K]
    time: list of times corresponding to t2 1st dimension
    ---
    tas: output temperature tas [K]
    atts: attributes of the output variable to be used in the output netcdf 
    """
    t2=varvals['T2'][:]
    t2=np.ma.masked_equal(t2,pm.const.missingval)
    if len(time)!=t2.shape[0]:
        sys.exit('ERROR in compute_tas: The lenght of time variable does not correspond to t2 first dimension')     
    
    #Generating a dictionary with the output attributes of the variable
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="air_temperature",ln="Surface air temperature",un="K",ts="time: point values %s seconds" %(tseconds),hg="2 m")
    
    tas=t2[:]    
    return tas,atts
    
    
def compute_ps(varvals,time,gvars):
    """Method to compute surface pressure
    psfc: psfc from wrf files [Pa]
    time: list of times corresponding to t2 1st dimension
    ---
    ps: output surface pressure [Pa]
    atts: attributes of the output variable to be used in the output netcdf 
    """
    psfc=varvals['PSFC'][:]
    psfc=np.ma.masked_equal(psfc,pm.const.missingval)
    if len(time)!=psfc.shape[0]:
        sys.exit('ERROR in compute_ps: The lenght of time variable does not correspond to psfc first dimension')     

    #Generating a dictionary with the output attributes of the variable
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="surface_air_pressure",ln="surface pressure",un="Pa",ts="time: point values %s seconds" %(tseconds))
    
    ps=psfc    
    return ps,atts


def compute_pracc(varvals,time,gvars):
    """Method to compute precipitation
    The original variable was accumulated throughout the simulation. Accumulation must be removed (But it is not converted to flux)
    rainc: convective rainfall accumulated. Including the timestep previous to the first one in this period [kg m-2]
    rainnc: non-convective rainfall accumulated.Including the timestep previous to the first one in this period [kg m-2]
    time: list of times corresponding to rainc 1st dimension
    ---
    pracc: accumulated rainfall (pracc) since last record [kg m-2]
    atts: attributes of the output variable to be used in the output netcdf
    """
    rainc=varvals['RAINC'][:]
    rainnc=varvals['RAINNC'][:]
    rainc=np.ma.masked_equal(rainc,pm.const.missingval)
    rainnc=np.ma.masked_equal(rainnc,pm.const.missingval)
    #rainc and rainnc includes the timestep previous to the first one to remove the accumulation.
    if (len(time)!=rainc.shape[0]-1) or (len(time)!=rainnc.shape[0]-1):
        sys.exit('ERROR in compute_pracc: The lenght of time variable does not correspond to rainc or rainnc first dimension')
    
    #Generating a dictionary with the output attributes of the variable
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="precipitation_amount",ln="Accumulated precipitation",un="Kg m-2",ts="time: point values %s seconds" %(tseconds))


    #Calculating difference between each timestep to remove the accumulation
    pracc=np.zeros((rainc.shape[0]-1,)+rainc.shape[1:],dtype=np.float64)
    pracc[:,:,:]=np.diff(rainc+rainnc,axis=0)
    pracc[pracc>pm.const.missingval]=pm.const.missingval
    return pracc,atts
    
def compute_huss(varvals,time,gvars):
    """Method to compute specific humidity
    q2: mixing ratio [kg kg-2]
    time: list of times corresponding to q2 1st dimension
    ---
    huss: specific humidity []
    atts: attributes of the output variable to be used in the output netcdf
    """
    q2=varvals['Q2'][:]
    q2=np.ma.masked_equal(q2,pm.const.missingval)
    if len(time)!=q2.shape[0]:
        sys.exit('ERROR in compute_huss: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="specific_humidity",ln="Surface specific humidity",un=" ",ts="time: point values %s seconds" %(tseconds),hg="2 m")
    
    huss=q2/(1+q2)
    
    return huss,atts

def compute_wss(varvals,time,gvars):
    """Method to compute wind speed
    u10: zonal wind [m s-1]
    v10: meridional wind [m s-1]
    time: list of times corresponding to u10 and v10 1st dimensions
    ---
    wss: wind speed [m s-1]
    atts: attributes of the output variable to be used in the output netcdf
    """
    u10=varvals['U10'][:]
    v10=varvals['V10'][:]
    u10=np.ma.masked_equal(u10,pm.const.missingval)
    v10=np.ma.masked_equal(v10,pm.const.missingval)
    if (len(time)!=u10.shape[0]) or (len(time)!=v10.shape[0]):
        sys.exit('ERROR in compute_wss: The lenght of time variable does not correspond to u10 or v10 first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="air_velocity",ln="Surface wind speed",un="m s-1",ts="time: point values %s seconds" %(tseconds),hg="10 m")
    
    #Winds need to be unstagged    
    #Wu10_unstagged=0.5*(u10[:,:,:-1]+u10[:,:,1:])
    #Wv10_unstagged=0.5*(v10[:,:-1,:]+v10[:,1:,:])
    
    wss=(u10**2+v10**2)**0.5
    
    return wss,atts

def compute_nonrotuas(varvals,time,gvars):
    """Method to compute eastward wind (NOT ROTATED)
    u10: zonal wind [m s-1]
    time: list of times corresponding to u10 1st dimension
    ---
    uas: eastward wind [m s-1]
    atts: attributes of the output variable to be used in the output netcdf
    """
    u10=varvals['U10'][:]
    u10=np.ma.masked_equal(u10,pm.const.missingval)
    if len(time)!=u10.shape[0]:
        sys.exit('ERROR in compute_uas: The lenght of time variable does not correspond to uas first dimension')
        
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="eastward_wind",ln="Eastward near-surface wind (not rotated)",un="m s-1",ts="time: point values %s seconds" %(tseconds),hg="10 m")    
    
    # Zonal wind unstagged
    # uas=0.5*(u10[:,:,:-1]+u10[:,:,1:])
    uas = u10

    return uas,atts


def compute_nonrotvas(varvals,time,gvars):
    """Method to compute northward wind (NOT ROTATED)
    v10: meridional wind [m s-1]
    time: list of times corresponding to v10 dimension
    ---
    vas: northward wind [m s-1]
    atts: attributes of the output variable to be used in the output netcdf
    """
    v10=varvals['V10'][:]
    v10=np.ma.masked_equal(v10,pm.const.missingval)
    if len(time)!=v10.shape[0]:
        sys.exit('ERROR in compute_vas: The lenght of time variable does not correspond to vas first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="northward_wind",ln="Northward near-surface wind (not rotated)",un="m s-1",ts="time: point values %s seconds" %(tseconds),hg="10 m")    

    # Zonal wind unstagged
    # vas=0.5*(v10[:,:-1,:]+v10[:,1:,:])
    vas = v10

    return vas,atts

def compute_uas(varvals,time,gvars):
    """Method to compute eastward wind (ROTATED - EARTH COORDINATES)
    u10: zonal wind [m s-1]
    v10: meridional wind [m s-1]
    time: list of times corresponding to u10 1st dimension
    ---
    uas: earth-coordinates eastward wind [m s-1]
    atts: attributes of the output variable to be used in the output netcdf
    """
    u10=varvals['U10'][:]
    v10=varvals['V10'][:]
    u10=np.ma.masked_equal(u10,pm.const.missingval)
    v10=np.ma.masked_equal(v10,pm.const.missingval)
    if (len(time)!=u10.shape[0]) or (len(time)!=v10.shape[0]):
        sys.exit('ERROR in compute_vas: The lenght of time variable does not correspond to U10 or V10 first dimension')
        
    fileref=nc.Dataset(gvars.fileref_att,'r')
    sina=fileref.variables['SINALPHA'][:]
    cosa=fileref.variables['COSALPHA'][:]
    sina_all=np.tile(sina, (len(time),1,1)) 
    cosa_all=np.tile(cosa, (len(time),1,1))
    
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="eastward_wind",ln="Eastward near-surface wind",un="m s-1",ts="time: point values %s seconds" %(tseconds),hg="10 m")    

    uas = u10[:]*cosa_all[:]-v10[:]*sina_all[:]

    return uas,atts


def compute_vas(varvals,time,gvars):
    """Method to compute northward wind (ROTATED - EARTH COORDINATES)
    u10: zonal wind [m s-1]
    v10: meridional wind [m s-1]
    time: list of times corresponding to v10 dimension
    ---
    vas: earth-coordinates northward wind [m s-1]
    atts: attributes of the output variable to be used in the output netcdf
    """
    u10=varvals['U10'][:]
    v10=varvals['V10'][:]
    u10=np.ma.masked_equal(u10,pm.const.missingval)
    v10=np.ma.masked_equal(v10,pm.const.missingval)
    if (len(time)!=u10.shape[0]) or (len(time)!=v10.shape[0]):
        sys.exit('ERROR in compute_vas: The lenght of time variable does not correspond to U10 or V10 first dimension')
    
    fileref=nc.Dataset(gvars.fileref_att,'r')
    sina=fileref.variables['SINALPHA'][:]
    cosa=fileref.variables['COSALPHA'][:]
    sina_all=np.tile(sina, (len(time),1,1)) 
    cosa_all=np.tile(cosa, (len(time),1,1)) 
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="northward_wind",ln="Northward near-surface wind",un="m s-1",ts="time: point values %s seconds" %(tseconds),hg="10 m")    

    vas = v10[:]*cosa_all[:]+u10[:]*sina_all[:]

    return vas,atts


def compute_evspsbl(varvals,time,gvars):
    """Method to compute surface evaporation flux
       The original variable was accumulated throughout the simulation. Accumulation must be removed.
       sfcevp: accumulated surface evaporation. Including the timestep previous to the first one in this period [kg m-2]
       time: list of times corresponding to sfcevp 1st dimension
       ---
       evspsbl: Surface evaporation flux [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    sfcevp=varvals['SFCEVP'][:]
    sfcevp=np.ma.masked_equal(sfcevp,pm.const.missingval)

   #sfcevp includes the timestep previous to the first one to remove the accumulation. 
    if len(time)!=sfcevp.shape[0]-1:
        sys.exit('ERROR in compute_evspsbl: The lenght of time variable does not correspond to var first dimension')
    
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="water_evaporation_flux",ln="Surface evaporation",un="kg m-2 s-1",ts="time: point values %s seconds" %(tseconds))
    
    #Calculating difference between each timestep to remove the accumulation
    #Divided by the number of seconds in each timestep to calculate the flux
    evspsbl=np.zeros((sfcevp.shape[0]-1,)+sfcevp.shape[1:],dtype=np.float64)
    evspsbl[:,:,:]=np.diff(sfcevp,axis=0)/tseconds
    
    return evspsbl,atts

def compute_mrso(varvals,time,gvars):
    """Method to compute the total soil moisture content
       Integrates through all soil layers
       smstot: total soil moisture in each layer [m3 m-3]
       dzs: thickness of each soil layer [m]
       time: list of times corresponding to smstot 1st dimension
       ---
       mrso: total soil moisture content [kg m-2]
       atts: attributes of the output variable to be used in the output netcdf
    """
    smstot=varvals['SMSTOT'][:]
    smstot=np.ma.masked_equal(smstot,pm.const.missingval)
    filemask=nc.Dataset(gvars.fileref_att,'r')
    mask=filemask.variables['LANDMASK'][:]

    filemask.close()

    if len(time)!=smstot.shape[0]:
        sys.exit('ERROR in compute_mrso: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="soil_moisture_content",ln="Total soil moisture content",un="kg m-2",ts="time: point values %s seconds" %(tseconds))
    
    mrso=smstot
    mrso[:,mask==0]=pm.const.missingval
    return mrso,atts

def compute_sst(varvals,time,gvars):
    """Method to compute the sea surface temperature
       sst_in: sea surface temperature [K]
       time: list of times corresponding to sst 1st dimension
       ---
       sst_out: sea surface temperature [K]
       atts: attributes of the output variable to be used in the output netcdf
    """
    sst_in=varvals['SST'][:]
    sst_in=np.ma.masked_equal(sst_in,pm.const.missingval)
    filemask=nc.Dataset(gvars.fileref_att,'r')
    mask=filemask.variables['LANDMASK'][:]

    filemask.close()

    if len(time)!=sst_in.shape[0]:
        sys.exit('ERROR in compute_sst: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="sea_surface_temperature",ln="sea surface temperature",un="K",ts="time: point values %s seconds" %(tseconds)) 
    
    sst_out=sst_in
    sst_out[:,mask==1]=pm.const.missingval

    return sst_out,atts
    
def compute_potevp(varvals,time,gvars):
    """Method to compute surface potential evaporation flux
       The original variable was accumulated throughout the simulation. Accumulation must be removed.
       potevp_in: accumulated surface potential evaporation. Including the timestep previous to the first one in this period [kg m-2]
       time: list of times corresponding to potevp 1st dimension
       ---
       potevp_out: Surface evaporation flux [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    potevp_in=varvals['POTEVP'][:]
    potevp_in=np.ma.masked_equal(potevp_in,pm.const.missingval)
    #potevp includes the timestep previous to the first one to remove the accumulation. 
    if len(time)!=potevp_in.shape[0]-1:
        sys.exit('ERROR in compute_potevp: The lenght of time variable does not correspond to var first dimension')    
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="water_potential_evaporation_flux",ln="Potential evaporation",un="W m-2",ts="time: point values %s seconds" %(tseconds))
    
    #Calculating difference between each timestep to remove the accumulation
    #Divided by the number of seconds in each timestep to calculate the flux
    potevp_out=np.zeros((potevp_in.shape[0]-1,)+potevp_in.shape[1:],dtype=np.float64)
    potevp_out[:,:,:]=np.diff(potevp_in,axis=0)/tseconds
    
    return potevp_out,atts

def compute_rsds(varvals,time,gvars):
    """Method to compute downward shortwave surface radiation
       swdown: downward short wave flux at ground surface [W m-2]
       time: list of times corresponding to swdown 1st dimension
       ---
       rsds: downward shortwave surface radiation [W m-2]
       atts: attributes of the output variable to be used in the output netcdf
    """ 
    swdown=varvals['SWDOWN'][:]
    swdown=np.ma.masked_equal(swdown,pm.const.missingval)
    if len(time)!=swdown.shape[0]:
        sys.exit('ERROR in compute_rsds: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="surface_downwelling_shortwave_flux_in_air",ln="Downward SW surface radiation",un="W m-2",ts="time: point values %s seconds" %(tseconds))
    
    rsds=swdown
    
    return rsds,atts

def compute_rlds(varvals,time,gvars):
    """Method to compute downward longwave surface radiation
       glw: downward long wave flux at ground surface [W m-2]
       time: list of times corresponding to glw 1st dimension
       ---
       rlds: downward shortwave surface radiation [W m-2]
       atts: attributes of the output variable to be used in the output netcdf
    """
    glw=varvals['GLW'][:]
    glw=np.ma.masked_equal(glw,pm.const.missingval)
    if len(time)!=glw.shape[0]:
        sys.exit('ERROR in compute_rlds: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="surface_downwelling_longwave_flux_in_air",ln="Downward LW surface radiation",un="W m-2",ts="time: point values %s seconds" %(tseconds))

    rlds=glw

    return rlds,atts
    
def compute_hfls(varvals,time,gvars):
    """Method to compute surface latent heat flux
       lh: latent heat flux at the surface [W m-2]
       time: list of times corresponding to lh 1st dimension
       ---
       hfls: upward latent heat flux at the surface [W m-2]
       atts: attributes of the output variable to be used in the output netcdf
    """
    lh=varvals['LH'][:]
    lh=np.ma.masked_equal(lh,pm.const.missingval)
    if len(time)!=lh.shape[0]:
        sys.exit('ERROR in compute_hfls: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="surface_upward_latent_heat_flux",ln="Latent heat flux at surface",un="W m-2",ts="time: point values %s seconds" %(tseconds))

    hfls=lh

    return hfls,atts
    
def compute_hfss(varvals,time,gvars):
    """Method to compute surface sensible heat flux
       hfx: upward heat flux at the surface [W m-2]
       time: list of times corresponding to emiss 1st dimension
       ---
       hfss: upward sensible heat flux at the surface [W m-2]
       atts: attributes of the output variable to be used in the output netcdf
    """
    hfx=varvals['HFX'][:]
    hfx=np.ma.masked_equal(hfx,pm.const.missingval)
    if len(time)!=hfx.shape[0]:
        sys.exit('ERROR in compute_hfss: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="surface_upward_sensible_heat_flux",ln="Heat flux at the surface",un="W m-2",ts="time: point values %s seconds" %(tseconds))

    hfss=hfx

    return hfss,atts

def compute_emiss(varvals,time,gvars):
    """Method to compute surface emissivity
       emiss_in: surface emissivity
       time: list of times corresponding to emiss 1st dimension
       ---
       emiss_out: surface emissivity
       atts: attributes of the output variable to be used in the output netcdf
    """
    emiss_in=varvals['EMISS'][:]
    emiss_in=np.ma.masked_equal(emiss_in,pm.const.missingval)
    if len(time)!=emiss_in.shape[0]:
        sys.exit('ERROR in compute_emiss: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="surface_emissivity",ln="Surface emissivity",un="",ts="time: point values %s seconds" %(tseconds))

    emiss_out=emiss_in

    return emiss_out,atts

def compute_albedo(varvals,time,gvars):
    """Method to compute surface albedo
    albedo_in: surface albedo
    time: list of times corresponding to emiss 1st dimension
    ---
    albedo_out: surface albedo
    atts: attributes of the output variable to be used in the output netcdf
    """
    albedo_in=varvals['ALBEDO'][:]  
    albedo_in=np.ma.masked_equal(albedo_in,pm.const.missingval)
    if len(time)!=albedo_in.shape[0]:
        sys.exit('ERROR in compute_albedo: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="albedo",ln="Surface albedo",un="",ts="time: point values %s seconds" %(tseconds))

    albedo_out=albedo_in

    return albedo_out,atts

def compute_rlus(varvals,time,gvars):
    """Method to compute upward longwave surface radiation

    tsk: surface skin temperature [K]
    emiss: surface emissivity
    time: list of times corresponding to swdown 1st dimension
    ---
    rlus: upward longwave surface radiation [W m-2]
    atts: attributes of the output variable to be used in the output netcdf
    """
    tsk=varvals['TSK'][:]
    emiss=varvals['EMISS'][:]  
    tsk=np.ma.masked_equal(tsk,pm.const.missingval)
    emiss=np.ma.masked_equal(emiss,pm.const.missingval)
    if (len(time)!=tsk.shape[0]) or(len(time)!=emiss.shape[0]) :
        sys.exit('ERROR in compute_rlus: The lenght of time variable does not correspond to emiss or tsk first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="surface_upwelling_longwave_flux_in_air",ln="Upwelling surface LW radiation",un="W m-2",ts="time: point values %s seconds" %(tseconds))

    #calculating with net lw radiation using Stefan-Boltzmann
    rlus=emiss*(pm.const.stefanboltz)*tsk**4

    return rlus,atts
    
def compute_tasmeantstep(varvals,time,gvars):
    """Method to compute the daily mean 2-m temperature using all timesteps of the model
       t2mean: mean 2-m temperature over all timesteps [K]
       time: list of times corresponding to swdown 1st dimension
       ---
       tasmeantstep:mean 2-m temperature over all timesteps [K]
       atts: attributes of the output variable to be used in the output netcdf
    """
    t2mean=varvals['T2MEAN'][:]
    t2mean=np.ma.masked_equal(t2mean,pm.const.missingval)
    if len(time)!=t2mean.shape[0]:
        sys.exit('ERROR in compute_tasmeantstep: The lenght of time variable does not correspond to var first dimension')
    
    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="air_temperature",ln="Mean surface air temperature",un="K",ts="time: mean 24 hour value from point values 60.0 second",hg="2 m")
    
    tasmeantstep=t2mean
    
    return tasmeantstep,atts
    
def compute_tasmintstep(varvals,time,gvars):
    """Method to compute the daily min 2-m temperature using all timesteps of the model
       t2min: min 2-m temperature over all timesteps [K]
       time: list of times corresponding to swdown 1st dimension
       ---
       tasmintstep:min 2-m temperature over all timesteps [K]
       atts: attributes of the output variable to be used in the output netcdf
    """
    
    t2min=varvals['T2MIN'][:]
    print len(time)
    print t2min.shape
    t2min=np.ma.masked_equal(t2min,pm.const.missingval)
    if len(time)!=t2min.shape[0]:
        sys.exit('ERROR in compute_tasmintstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="air_temperature",ln="Minimum surface air temperature",un="K",ts="time: min 24 hour value from point values 60.0 second",hg="2 m")

    tasmintstep=t2min

    return tasmintstep,atts
    
def compute_tasmaxtstep(varvals,time,gvars):
    """Method to compute the daily max 2-m temperature using all timesteps of the model
       t2max: max 2-m temperature over all timesteps [K]
       time: list of times corresponding to swdown 1st dimension
       ---
       tasmaxtstep:max 2-m temperature over all timesteps [K]
       atts: attributes of the output variable to be used in the output netcdf
    """
    t2max=varvals['T2MAX'][:]
    t2max=np.ma.masked_equal(t2max,pm.const.missingval)
    if len(time)!=t2max.shape[0]:
        sys.exit('ERROR in compute_tasmaxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="air_temperature",ln="Maximum surface air temperature",un="K",ts="time: max 24 hour value from point values 60.0 second",hg="2 m")

    tasmaxtstep=t2max

    return tasmaxtstep,atts

def compute_wssmaxtstep(varvals,time,gvars):
    """Method to compute the daily max wind speed using all timesteps of the model
       spduv10max: daily max wind speed over all timesteps [m s-1]
       time: list of times corresponding to uv10max5 1st dimension
       ---
       wssmaxtstep:maximum daily max wind speed over all timesteps [m s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    spduv10max=varvals['SPDUV10MAX'][:]
    spduv10max=np.ma.masked_equal(spduv10max,pm.const.missingval)
    if len(time)!=spduv10max.shape[0]:
        sys.exit('ERROR in compute_wssmaxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="air_velocity",ln="Max. surface wind speed",un="m s-1",ts="time: maximum 24 hour value from point values 60.0 second")

    wssmaxtstep=spduv10max

    return wssmaxtstep,atts

def compute_pr5maxtstep(varvals,time,gvars):
    """Method to compute the max 5-minute precipitationusing all timesteps of the model
       prmax5: maximum 5-minute precipitation using all timesteps [kg m-2 s-1]
       time: list of times corresponding to swdown 1st dimension
       ---
       pr5maxtstep:maximum 5-minute precipitation using all timesteps [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    prmax5=varvals['PRMAX5'][:]
    prmax5=np.ma.masked_equal(prmax5,pm.const.missingval)
    if len(time)!=prmax5.shape[0]:
        sys.exit('ERROR in compute_pr5maxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="5max_precipitation_flux",ln="Max. 5-minute time-window moving averaged precipitation rate",un="kg m-2 s-1",ts="time: maximum 5-minute time-window moving averaged values from point values 60.0 second")

    pr5maxtstep=prmax5

    return pr5maxtstep,atts

def compute_pr10maxtstep(varvals,time,gvars):
    """Method to compute the max 10-minute precipitationusing all timesteps of the model
       prmax10: maximum 10-minute precipitation using all timesteps [kg m-2 s-1]
       time: list of times corresponding to swdown 1st dimension
       ---
       pr10maxtstep:maximum 10-minute precipitation using all timesteps [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    prmax10=varvals['PRMAX10'][:]
    prmax10=np.ma.masked_equal(prmax10,pm.const.missingval)
    if len(time)!=prmax10.shape[0]:
        sys.exit('ERROR in compute_pr10maxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="10max_precipitation_flux",ln="Max. 10-minute time-window moving averaged precipitation rate",un="kg m-2 s-1",ts="time: maximum 10-minute time-window moving averaged values from point values 60.0 second")

    pr10maxtstep=prmax10

    return pr10maxtstep,atts

def compute_pr20maxtstep(varvals,time,gvars):
    """Method to compute the max 20-minute precipitationusing all timesteps of the model
       prmax20: maximum 20-minute precipitation using all timesteps [kg m-2 s-1]
       time: list of times corresponding to swdown 1st dimension
       ---
       pr20maxtstep:maximum 20-minute precipitation using all timesteps [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    prmax20=varvals['PRMAX20'][:]
    prmax20=np.ma.masked_equal(prmax20,pm.const.missingval)
    if len(time)!=prmax20.shape[0]:
        sys.exit('ERROR in compute_pr20maxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="20max_precipitation_flux",ln="Max. 20-minute time-window moving averaged precipitation rate",un="kg m-2 s-1",ts="time: maximum 20-minute time-window moving averaged values from point values 60.0 second")

    pr20maxtstep=prmax20

    return pr20maxtstep,atts
    
def compute_pr30maxtstep(varvals,time,gvars):
    """Method to compute the max 30-minute precipitationusing all timesteps of the model
       prmax30: maximum 30-minute precipitation using all timesteps [kg m-2 s-1]
       time: list of times corresponding to swdown 1st dimension
       ---
       pr30maxtstep:maximum 30-minute precipitation using all timesteps [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    prmax30=varvals['PRMAX30'][:]
    prmax30=np.ma.masked_equal(prmax30,pm.const.missingval)
    if len(time)!=prmax30.shape[0]:
        sys.exit('ERROR in compute_pr30maxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="30max_precipitation_flux",ln="Max. 30-minute time-window moving averaged precipitation rate",un="kg m-2 s-1",ts="time: maximum 30-minute time-window moving averaged values from point values 60.0 second")

    pr30maxtstep=prmax30

    return pr30maxtstep,atts
    

def compute_pr1Hmaxtstep(varvals,time,gvars):
    """Method to compute the max 1-hour precipitationusing all timesteps of the model
       prmax1H: maximum 1-hour precipitation using all timesteps [kg m-2 s-1]
       time: list of times corresponding to swdown 1st dimension
       ---
       pr1Hmaxtstep:maximum 1-hour precipitation using all timesteps [kg m-2 s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    prmax1H=varvals['PRMAX1H'][:]
    prmax1H=np.ma.masked_equal(prmax1H,pm.const.missingval)
    if len(time)!=prmax1H.shape[0]:
        sys.exit('ERROR in compute_pr1Hmaxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="1Hmax_precipitation_flux",ln="Max. 1-hour time-window moving averaged precipitation rate",un="kg m-2 s-1",ts="time: maximum 1-hour time-window moving averaged values from point values 60.0 second")

    pr1Hmaxtstep=prmax1H

    return pr1Hmaxtstep,atts
    
def compute_wss5maxtstep(varvals,time,gvars):
    """Method to compute the max 5-minute wind speed using all timesteps of the model
       uv10max5: maximum 5-mimute wind speed using all timesteps [m s-1]
       time: list of times corresponding to uv10max5 1st dimension
       ---
       wss5maxtstep:maximum 5-mimute wind speed using all timesteps [m s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    uv10max5=varvals['UV10MAX5'][:]
    uv10max5=np.ma.masked_equal(uv10max5,pm.const.missingval)
    if len(time)!=uv10max5.shape[0]:
        sys.exit('ERROR in compute_wss5maxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="5max_air_velocity",ln="Max. 5-minute time-window moving averaged surface wind speed",un="m s-1",ts="time: maximum 5-minute time-window moving averaged values from point values 60.0 second")

    wss5maxtstep=uv10max5

    return wss5maxtstep,atts

def compute_wss10maxtstep(varvals,time,gvars):
    """Method to compute the max 10-minute wind speed using all timesteps of the model
       uv10max10: maximum 10-mimute wind speed using all timesteps [m s-1]
       time: list of times corresponding to uv10max10 1st dimension
       ---
       wss10maxtstep:maximum 10-mimute wind speed using all timesteps [m s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    uv10max10=varvals['UV10MAX10'][:]
    uv10max10=np.ma.masked_equal(uv10max10,pm.const.missingval)
    if len(time)!=uv10max10.shape[0]:
        sys.exit('ERROR in compute_wss10maxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="10max_air_velocity",ln="Max. 10-minute time-window moving averaged surface wind speed",un="m s-1",ts="time: maximum 10-minute time-window moving averaged values from point values 60.0 second")

    wss10maxtstep=uv10max10

    return wss10maxtstep,atts

def compute_wss20maxtstep(varvals,time,gvars):
    """Method to compute the max 20-minute wind speed using all timesteps of the model
       uv10max20: maximum 20-mimute wind speed using all timesteps [m s-1]
       time: list of times corresponding to uv10max20 1st dimension
       ---
       wss20maxtstep:maximum 20-mimute wind speed using all timesteps [m s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    uv10max20=varvals['UV10MAX20'][:]
    uv10max20=np.ma.masked_equal(uv10max20,pm.const.missingval)
    if len(time)!=uv10max20.shape[0]:
        sys.exit('ERROR in compute_wss20maxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="20max_air_velocity",ln="Max. 20 minute time-window moving averaged surface wind speed",un="m s-1",ts="time: maximum 20-minute time-window moving averaged values from point values 60.0 second")

    wss20maxtstep=uv10max20

    return wss20maxtstep,atts

def compute_wss30maxtstep(varvals,time,gvars):
    """Method to compute the max 30-minute wind speed using all timesteps of the model
       uv10max30: maximum 30-mimute wind speed using all timesteps [m s-1]
       time: list of times corresponding to uv10max30 1st dimension
       ---
       wss30maxtstep:maximum 30-mimute wind speed using all timesteps [m s-1]
       atts: attributes of the output variable to be used in the output netcdf
    """
    uv10max30=varvals['UV10MAX30'][:]
    uv10max30=np.ma.masked_equal(uv10max30,pm.const.missingval)
    if len(time)!=uv10max30.shape[0]:
        sys.exit('ERROR in compute_wss30maxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="30max_air_velocity",ln="Max. 30 minute time-window moving averaged surface wind speed",un="m s-1",ts="time: maximum 30-minute time-window moving averaged values from point values 60.0 second")

    wss30maxtstep=uv10max30

    return wss30maxtstep,atts

def compute_wss1Hmaxtstep(varvals,time,gvars):
    """Method to compute the max 1-hour wind speed using all timesteps of the model
    uv10max1H: maximum 1-hour wind speed using all timesteps [m s-1]
    time: list of times corresponding to uv10max1H 1st dimension
    ---
    wss1Hmaxtstep:maximum 1-hour wind speed using all timesteps [m s-1]
    atts: attributes of the output variable to be used in the output netcdf
    """
    uv10max1H=varvals['UV10MAX1H'][:]
    uv10max1H=np.ma.masked_equal(uv10max1H,pm.const.missingval)
    if len(time)!=uv10max1H.shape[0]:
        sys.exit('ERROR in compute_wss1Hmaxtstep: The lenght of time variable does not correspond to var first dimension')

    tseconds=round(((time[-1]-time[0]).total_seconds()/len(time)))
    atts=pm.get_varatt(sn="1Hmax_air_velocity",ln="Max. 1-hour time-window moving averaged surface wind speed",un="m s-1",ts="time: maximum 1-hour time-window moving averaged values from point values 60.0 second")

    wss1Hmaxtstep=uv10max1H

    return wss1Hmaxtstep,atts
