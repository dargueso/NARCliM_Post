import re
import sys
import os
import netCDF4 as nc
import numpy as np
import datetime as dt
import glob as glob
import compute_stats as coms
import compute_vars as comv
from dateutil.relativedelta import relativedelta
from collections import OrderedDict
class const:
  """Class that contains most used atmospheric constant values
  """
  Rd = 287.04
  Rv = 461.5
  RdRv = Rd / Rv
  cp = 7.*Rd/2.
  epsilon_gamma = 0.62197
  es_base_bolton = 0.6112
  es_Abolton = 17.67
  es_Bbolton = 243.5
  es_base_tetens = 6.1078
  es_Atetens_vapor = 7.5
  es_Btetens_vapor = 237.3
  es_Atetens_ice = 9.5
  es_Btetens_ice = 265.5
  g = 9.81
  stefanboltz = 5.670373e-8
  p1000mb = 100000.
  rcp = Rd/cp
  tkelvin = 273.15
  missingval = 1.e+20
  
  
# *************************************************************************************
class gvar:
  def __init__(self,inputinf):
    if inputinf['GCM']=='CCCMA3.1':
      self.GCM_calendar='no_leap'
    else:
      self.GCM_calendar='standard'
    self.ref_date=dt.datetime(1949,12,1,00,00,00)
    self.pathin=inputinf['pathin']
    self.pathout=inputinf['pathout']
    self.GCM=inputinf['GCM']
    self.RCM=inputinf['RCM']
    self.syear=int(inputinf['syear'])
    self.eyear=int(inputinf['eyear'])
    self.domain=inputinf['domain']
    self.overwrite=inputinf['overwrite']
    self.outfile_patt=inputinf['outfile_patt']
    self.fileref_att='%s/wrfout_%s_%s-01-01_00:00:00' %(self.pathin,self.domain,self.syear)


# *************************************************************************************
def create_outdir(gvars):
  fullpathout='%spostprocess/%s-%s/%s/%s/%s/' %(gvars.pathout,gvars.syear,gvars.eyear,gvars.GCM,gvars.RCM,gvars.domain)
  
  # CREATE OUTPUT DIR IF IT DOESN'T EXIST
  if not os.path.exists(fullpathout):
    os.makedirs(fullpathout)
  
  # CREATE A TEMPORAL DIR WITHIN THE OUTPUT DIR IF IT DOESN'T EXIST
  if not os.path.exists("%stemp/" %(fullpathout)):
    os.makedirs("%stemp/" %(fullpathout))
  
  return fullpathout
  

# *************************************************************************************
def create_outtime(dates,gvars):
  datehours=date2hours(dates,gvars.ref_date)
  time=np.zeros(len(datehours),dtype=np.float64)
  time[:-1]=datehours[:-1]+np.diff(datehours)/2
  time[-1]=datehours[-1]+(datehours[-1]-datehours[-2])/2
  return time


# *************************************************************************************
def create_timebnds(time):
  time_bnds=np.zeros((len(time),2),dtype=np.float64)

  time_bnds[:-1,1]=time[:-1]+np.diff(time)/2
  time_bnds[:-1,0]=time[:-1]-np.diff(time)/2

  time_bnds[-1,0]=time[-1]-(time[-1]-time[-2])/2
  time_bnds[-1,1]=time[-1]+(time[-1]-time[-2])/2
  
  return time_bnds


# *************************************************************************************
def add_timestep_acc(wrfvar,varvals,year,gvars,filet):
  accvar={}
  for wrfv in wrfvar:
    if year<gvars.eyear:
      next_file='%s%s_%s_%s-01-01_00:00:00' % (gvars.pathin,filet,gvars.domain,year+1)
      ncfile=nc.Dataset(next_file,'r')
      next_tstep=np.squeeze(ncfile.variables[wrfv][:])
      accvar[wrfv]=np.concatenate((varvals[wrfv],next_tstep[0:1,:,:]),axis=0)
    else:
      fillvar=np.ones((1,)+varvals[wrfv].shape[1:],dtype=np.float64)*const.missingval
      accvar[wrfv]=np.concatenate((varvals[wrfv][:],fillvar),axis=0)

  return accvar


# *************************************************************************************
def mv_timestep(wrfvar,varvals,year,gvars,filet):
  accvar={}
  for wrfv in wrfvar:
    if year<gvars.eyear:
      next_file='%s%s_%s_%s-01-01_00:00:00' % (gvars.pathin,filet,gvars.domain,year+1)
      ncfile=nc.Dataset(next_file,'r')
      print 'READ ONE MORE TIME STEP: ', next_file
      next_tstep=np.squeeze(ncfile.variables[wrfv][:])
      accvar[wrfv]=np.concatenate((varvals[wrfv],next_tstep[0:1,:,:]),axis=0)
      accvar[wrfv]=accvar[wrfv][1:,:,:] 
    else:
      fillvar=np.ones((1,)+varvals[wrfv].shape[1:],dtype=np.float64)*const.missingval
      accvar[wrfv]=np.concatenate((varvals[wrfv][:],fillvar),axis=0)
      accvar[wrfv]=accvar[wrfv][1:,:,:]
  return accvar


# *************************************************************************************
def get_wrfvars(wrfvar,fin):
  variabs={}
  for wrfv in wrfvar:
    variabs[wrfv]=fin.variables[wrfv][:].astype('float64')
  return  variabs
  

# *************************************************************************************
def read_input(filename):
  """
  Read input file with input arguments:
  pathin: path to input wrf files
  pathout: path to output postprocessed files
  GCM: Name of the GCM (e.g. MIROC3.2)
  RCM: Name of the RCM (e.g. R1)
  syear: First year to postprocess (e.g. 1990)
  eyear: Last year to postprocess (e.g. 2009)
  domain: domain to postprocess (e.g. 'd02')
  
  """
  filein=open(filename,'r')

  
  options,sentinel,varnames=filein.read().partition('#### Requested output variables (DO NOT CHANGE THIS LINE) ####')
  fileopts=open('fileopts.tmp','w')
  fileopts.write(options)
  filevarnames=open('filevarnames.tmp','w')
  filevarnames.write(varnames)
  
  
  fileopts=open('fileopts.tmp','r')
  lines=fileopts.readlines()

  inputinf={}
  entryname=[]
  entryvalue=[]
  for line in lines:
    line=re.sub('\s+',' ',line)
    line=line.replace(" ", "")
    li=line.strip()
    #Ignore empty lines
    if li:
      #Ignore commented lines
      if not li.startswith("#"): 
        values=li.split('=')
        entryname.append(values[0])
        entryvalue.append(values[1])
  fileopts.close()
  
  for ii in xrange(len(entryname)):
    inputinf[entryname[ii]]=entryvalue[ii]
    
  filevarnames=open('filevarnames.tmp','r')
  lines=filevarnames.readlines()
  varnames=[]
  for line in lines:
    line=re.sub('\s+',' ',line)
    line=line.replace(" ", "")
    li=line.strip()
    #Ignore empty lines
    if li:
      #Ignore commented lines
      if not li.startswith("#"):
        values=li.split()
        varnames.append(values[0])
  filevarnames.close()     
  
  os.remove('filevarnames.tmp')
  os.remove('fileopts.tmp') 
  
  print 'Variables that will be obtained from postprocessing:',varnames
  return inputinf,varnames

# *************************************************************************************

def read_varinfo(filename):
  """
  Read file (filename) with information regarding the output variables.
  varname: Name of the variable
  filetype: Type of file where the original variables are stored (wrfhrly, wrfout, wrfxtrm, wrfdly)
  freqreq: Frequency requested (now it set to DAY and MON by default to all variables)
  statsreq: Statistics requested for each variable (acc, mean, max, min)
  ---
  varinfo: dictionary with filetype,varnames and stats
  """
  filein=open(filename,"r")
  lines=filein.readlines()
  varinfo={}
   
  varname=[]
  filetype=[]
  freqreq=[]
  statsreq=[]
  for line in lines:
    line=re.sub('\s+',' ',line)
    li=line.strip()
    #Ignore empty lines
    if li:
      #Ignore commented lines 
      if not li.startswith("#"):
        values=li.split(' ')
        varname.append(values[0])
        filetype.append(values[1])
        freqreq.append(values[2])
        statsreq.append(values[3])

  filein.close()


  varinfo=dictionary2entries(filetype,varname,statsreq)

  return varinfo

 
# *************************************************************************************
def getwrfname(varname):
  """ Dictionary containing the equivalence between CF variable names and WRF 
  variable names.
  """

  dic={'tas':['T2'],'pracc':['RAINC-RAINNC'],'ps':['PSFC'], 'uas':['U10-V10'], 'vas':['U10-V10'],\
       'huss':['Q2'], 'wss':['U10-V10'], 'sst':['SST'],'rsds':['SWDOWN'], 'rlds':['GLW'],\
       'emiss':['EMISS'], 'albedo':['ALBEDO'], 'hfls':['LH'], 'hfss':['HFX'],'evspsbl':['SFCEVP'],\
       'mrso':['SMSTOT'], 'potevp':['POTEVP'], 'rlus':['TSK-EMISS'],'tasmeantstep':['T2MEAN'],\
       'tasmintstep':['T2MIN'],'tasmaxtstep':['T2MAX'],'wssmaxtstep':['SPDUV10MAX'], 'pr5maxtstep':['PRMAX5'],\
       'pr10maxtstep':['PRMAX10'],'pr20maxtstep':['PRMAX20'],'pr30maxtstep':['PRMAX30'], 'pr1Hmaxtstep':['PRMAX1H'],\
       'wss5maxtstep':['UV10MAX5'], 'wss10maxtstep':['UV10MAX10'], 'wss20maxtstep':['UV10MAX20'],\
       'wss30maxtstep':['UV10MAX30'], 'wss1Hmaxtstep':['UV10MAX1H']}
  
  return dic[varname]


# *************************************************************************************
def read_schemes(filename):
  import netCDF4 as nc

  """Method that reads a sample WRF file and gets the Physics option used in the simulation
     Searches in WRF_schemes.inf file that contains the names and references of all schemes
     Outputs a dictionary with the schemes types and the schemes names (and references used in the simulation)
     To be used in get_globatt
     ---
     sch_info: dictionary with the schemes options information
  """
  filein=nc.Dataset(filename,'r')
  schemes={}
  schemes['ra_lw_physics']=filein.RA_LW_PHYSICS
  schemes['sf_sfclay_physics']=filein.SF_SFCLAY_PHYSICS
  schemes['cu_physics']=filein.CU_PHYSICS
  schemes['bl_pbl_physics']=filein.BL_PBL_PHYSICS
  schemes['ra_sw_physics']=filein.RA_SW_PHYSICS
  schemes['sf_surface_physics']=filein.SF_SURFACE_PHYSICS
  
  fileschemes=open("./WRF_schemes.inf")
  lines=fileschemes.readlines()
  sch_infoall={}
  sch_type=[]
  sch_number=[]
  sch_name=[]
  sch_ref=[]
  for line in lines:
    line=re.sub('\s+',' ',line)
    li=line.strip()
    #Ignore empty lines
    if li:
      #Ignore commented lines 
      if not li.startswith("#"):
        values=li.split(' ')
        sch_type.append(values[0])
        sch_number.append(values[1])
        sch_name.append('%s; %s' %(values[2],values[3]))
        
  fileschemes.close()
  
  sch_infoall=dictionary2entries(sch_type,sch_number,sch_name)
  
  
  sch_info={}
  for ii in xrange(len(schemes)):
    sch=schemes.keys()[ii]
    num=str(schemes[schemes.keys()[ii]])
    sch_info[sch]=sch_infoall[sch][num].replace("!"," ")


  return sch_info


# *************************************************************************************
def get_varatt(sn,ln,un,ts,hg=None):
  att=OrderedDict()
  att['standard_name'] = sn
  att['long_name'] = ln
  att['units']=un
  if hg!=None:
    att['height']=hg
  
  att['coordinates'] = "lon lat"
  att['cell_method'] = "%s" %(ts)
  att['grid_mapping'] = "Rotated_pole"
  att['_FillValue']=const.missingval

  return att


# *************************************************************************************
def get_wrfdate(time):
  """Method to extract year, month, day and hour from a particular timestep of the
  time variable of the wrf output.

  time: string from the time variable of a wrf output.
  """

  time=np.squeeze(time)
  year=int(time[0])*1000+int(time[1])*100+int(time[2])*10+int(time[3]) 
  month=int(time[5])*10+int(time[6])
  day=int(time[8])*10+int(time[9])
  hour=int(time[11])*10+int(time[12])

  return year, month, day, hour


# *************************************************************************************
def add_leapdays(varvals,wrfvar,date,time_step):
  """Method that add missing values to the 29th February in leap year for
  those simulations that are run without leap years.
  
  varval: 3-D matrix with time, latitude and longitude in the first, second and third dimensions. 
  The time dimension has one year of data without 29th february in a leap-year case.
  
  index: is the index of the 29th feburary in the time dimension.
  """
  
  months_all=np.asarray([date[i].month for i in xrange(len(date))]) 
  days_all=np.asarray([date[i].day for i in xrange(len(date))])
  index=np.where(np.logical_and(months_all==2,days_all==29))
  tt=24/time_step
  
  for wrfv in wrfvar:
    temp=np.ones((varvals[wrfv].shape[0]+tt,varvals[wrfv].shape[1],varvals[wrfv].shape[2]))*const.missingval
    
    temp[0:index[0][0]+1,:,:]=varvals[wrfv][0:index[0][0]+1,:,:]
    temp[index[0][0]+tt+1:,:,:]=varvals[wrfv][index[0][0]+1:,:,:]
    
    varvals[wrfv]=temp
    
  return varvals
  

# *************************************************************************************
def get_globatt(GCM,RCM,sch_info,perturb=None):
  import datetime as dt
  """Method that generates the global attributes
  Defines the name of the schemes and the references to include in the global attributes

  GCM: Name of the GCM (e.g.:MIROC3.2)
  RCM: Name of the RCM (e.g.: R1)
  Period: Name of the perturbation of the GCM
  sch_info: dictionary containing the kind of schemes as keys, and the name and references as values
  e.g.: global_attributes=ga.globalatt("MIROC3.2","R1",sch_info,"d1")
  """
  glatt=OrderedDict()
  glatt['institution']        = "University of New South Wales (Australia)" 
  glatt['contact']          = "jason.evans@unsw.edu.au"
  glatt['institute_id']       = "CCRC"
  glatt['institute']          = "Climate Change Research Centre"
  glatt['references']         = "http://www.ccrc.unsw.edu.au"
  glatt['model_id']         = "WRF"
  glatt['RCM_version_id']       = "v3.3"
  glatt['RCM_source']         = "WRF 3.3 modified at the Climate Change Research Centre"  
  glatt['CORDEX_domain']        = "AUS44"
  glatt['project_id']         = "NARCliM"
  glatt['title']            = "Projection run"
  glatt['experiment_id']        = "projection"
  glatt['experiment']         = "Projection run with %s" %(GCM)
  glatt['driving_experiment']     = "%s, projection, %s%s%s" %(GCM,RCM,GCM,perturb)
  glatt['driving_model_id']     = "%s%s" %(GCM,perturb)
  glatt['driving_model_ensemble_member'] = "%s%s%s" %(RCM,GCM,perturb)
  glatt['wrf_options']    = "sst_update & tmn_update"
  glatt['driving_experiment_name']  = "projection_%s_%s" %(GCM,RCM) 
  glatt['product']          = "output"
  glatt['creation_date'] = dt.datetime.utcnow().strftime("%Y/%m/%d %H:%M:%S UTC")
  glatt['wrf_schemes_ra_lw_physics']        = "%s" %(sch_info['ra_lw_physics']) 
  glatt['wrf_schemes_sf_sfclay_physics']      = "%s" %(sch_info['sf_sfclay_physics']) 
  glatt['wrf_schemes_cu_physics']       = "%s" %(sch_info['cu_physics'])  
  glatt['wrf_schemes_bl_pbl_physics']     = "%s" %(sch_info['bl_pbl_physics']) 
  glatt['wrf_schemes_ra_sw_physics']        = "%s" %(sch_info['ra_sw_physics']) 
  glatt['wrf_schemes_sf_surface_physics']   = "%s" %(sch_info['sf_surface_physics'])
    
  return glatt


# *************************************************************************************
def dictionary2entries(vals1, vals2, vals3):
  """ Function to create a dictionary with 3 entries (thanks to Jeff Exbrayat, CoECSCC-CCRC)
  >>>country = ['India', 'Germany', 'Guatemala', 'Burma', 'South Africa', 'South Africa']
  >>>person = ['Mohandas Karamchand Ghandi', 'Albert Einstein', 'Rigoberta Menchu', 'Aung San Suu Kyi', 'Nelson Mandela', 'Desmond Tutu']
  >>>fight = ['Independence', 'Nuclear bomb', 'Indigenous Rights', 'Freedom', 'anti-Apartheid', 'Reconciliation']

  >>>dict3 = dictionary3entries(country, person, fight)

  >>>print dict3['South Africa']['Nelson Mandela']
  anti-Apartheid
  >>>print dict3['South Africa']
  {'Nelson Mandela': 'anti-Apartheid', 'Desmond Tutu': 'Reconciliation'}
  """

  dicjeff={}

  for ii in range(len(vals1)):
    if vals1[ii] not in dicjeff.keys():
      dicjeff[vals1[ii]]={}
      if vals2[ii] not in dicjeff[vals1[ii]].keys():
        dicjeff[vals1[ii]][vals2[ii]]={}
        dicjeff[vals1[ii]][vals2[ii]]=vals3[ii]
      else:
        dicjeff[vals1[ii]][vals2[ii]]=vals3[ii]
    else:
      if vals2[ii] not in dicjeff[vals1[ii]].keys():
        dicjeff[vals1[ii]][vals2[ii]]={}
        dicjeff[vals1[ii]][vals2[ii]]=vals3[ii]
      else:
        dicjeff[vals1[ii]][vals2[ii]]=vals3[ii]
  return dicjeff


# *************************************************************************************
def create_netcdf(info,gvars, varval, time, time_bnds):
    

  """ Create a netcdf file for the post-processed variables of NARCliM simulations
       
  By default, the module do not overwrite the file so if there is one with the same name then 
  the script does not do anything.

  Input: global  attributres from pre-defined classes:   
  Output: a netcdf file
  Author: Alejanro Di Luca, Daniel Argueso
  Created: 14/09/2013
  Last Modification: 07/08/2013

  """
  print '\n', ' CALLING CREATE_NETCDF MODULE ','\n'
  import numpy as np
  import netCDF4 as nc
  import sys
  

  file_out=info[0]
  varname=info[1]
  varatt=info[2]
  time_bounds=info[3]

  # **********************************************************************
  # Read attributes from the geo_file of the corresponding domain
  fin1=nc.Dataset(gvars.fileref_att,mode='r')
  lon=np.squeeze(fin1.variables['XLONG'][:]) # Getting longitude
  lat=np.squeeze(fin1.variables['XLAT'][:]) # Getting latitude
  dx=getattr(fin1, 'DX')
  dy=getattr(fin1, 'DY')
  cen_lat=getattr(fin1, 'CEN_LAT')
  cen_lon=getattr(fin1, 'CEN_LON')
  pole_lat=getattr(fin1, 'POLE_LAT')
  pole_lon=getattr(fin1, 'POLE_LON')
  stand_lon=getattr(fin1, 'STAND_LON')
  fin1.close()
  sch_info=read_schemes(gvars.fileref_att) 
  #**********************************************************************
  # CREATING NETCDF FILE
  # Create output file
  fout=nc.Dataset(file_out,mode='w', format='NETCDF4_CLASSIC')

  # ------------------------
  # Create dimensions
  print '   CREATING AND WRITING DIMENSIONS: '
  print '              TIME, X, Y(, TIME_BNDS)'
  fout.createDimension('x',varval.shape[2])
  fout.createDimension('y',varval.shape[1])
  fout.createDimension('time',None)
  fout.createDimension('bnds', 2)

  # ------------------------
  # Create and assign values to variables
  print "\n"
  print '   CREATING AND WRITING VARIABLES:'

  # VARIABLE: longitude
  print '    ---   LONGITUDE VARIABLE CREATED ' 
  varout=fout.createVariable('lon','f',['y', 'x'])
  varout[:]=lon[:]
  setattr(varout, 'standard_name','longitude')
  setattr(varout, 'long_name','Longitude')
  setattr(varout, 'units','degrees_east')
  setattr(varout, '_CoordinateAxisType','Lon')
  
  # VARIABLE: latitude
  print '    ---   LATITUDE VARIABLE CREATED ' 
  varout=fout.createVariable('lat','f',['y', 'x'])
  varout[:]=lat[:]
  setattr(varout, 'standard_name','latitude')
  setattr(varout, 'long_name','Latitude')
  setattr(varout, 'units','degrees_north')
  setattr(varout, '_CoordinateAxisType','Lat')
     
  # VARIABLE: time 
  print '    ---   TIME VARIABLE CREATED ' 
  varout=fout.createVariable('time','f',['time'])
  varout[:]=time[:]
  setattr(varout, 'standard_name','time')
  setattr(varout, 'long_name','time')
  setattr(varout, 'bounds','time_bnds')
  setattr(varout, 'units','hours since %s' %(gvars.ref_date.strftime("%Y-%m-%d %H:%M:%S")))
  setattr(varout, 'calendar','standard')

  # VARIABLE: time_bnds 
  if time_bounds==True:
    print '  ---   TIME_BNDS VARIABLE CREATED ' 
    varout=fout.createVariable('time_bnds','f8',['time', 'bnds'])
    varout[:]=time_bnds[:]
    setattr(varout, 'units','hours since %s' %(gvars.ref_date.strftime("%Y-%m-%d %H:%M:%S")))
    setattr(varout, 'calendar','standard')
    
  # VARIABLE: variable
  print '    ---   ',varname, ' VARIABLE CREATED ' 
  varout=fout.createVariable(varname,'f',['time', 'y', 'x'], fill_value=varatt['_FillValue'])
  varout[:]=varval[:]
  for att in varatt.keys():
    if att!='_FillValue':
      if varatt[att]!=None:
        setattr(varout, att, varatt[att])
      
  # VARIABLE: Rotated_Pole 
  print '    ---   Rotated_pole VARIABLE CREATED ' 
  varout=fout.createVariable('Rotated_pole','c',[])
  setattr(varout, 'grid_mapping_name', 'rotated_latitude_longitude')
  setattr(varout, 'dx_m', dx)
  setattr(varout, 'dy_m', dy)
  setattr(varout, 'latitude_of_projection_origin', cen_lat)
  setattr(varout, 'longitude_of_central_meridian',cen_lon)
  setattr(varout, 'true_longitude_of_projection',stand_lon)
  setattr(varout, 'grid_north_pole_latitude',  pole_lat)
  setattr(varout, 'grid_north_pole_longitude', pole_lon)
    
  # WRITE GLOBAL ATTRIBUTES
  print '\n', ' CREATING AND WRITING GLOBAL ATTRIBUTES:'
  gblatt = get_globatt(gvars.GCM,gvars.RCM,sch_info)
  for att in gblatt.keys():
    setattr(fout, att, gblatt[att])
  fout.close()
  print '  ===> FILE: ', file_out
    
  print '     ------------  SUCCESFULLY CREATED!!!  ------------ '


#**************************************************************************************

def checkpoint(ctime):
  import time

  """ Computes the spent time from the last checkpoint
      
  Input: a given time
  Output: present time
  Print: the difference between the given time and present time
  Author: Alejandro Di Luca
  Created: 07/08/2013
  Last Modification: 14/08/2013
  
  """
  if ctime==0:
    ctime=time.time()
    dtime=0
  else:
    dtime=time.time()-ctime
    ctime=time.time()
    print '======> DONE in ',float('%.2g' %(dtime)),' seconds',"\n"
  return ctime

#**************************************************************************************

def checkfile(file_out,overwrite):
  """Checks if the output file exist and whether it should be written or not
  """

  # ***********************************************************
  # BEFORE READING AND PROCESSING THE VARIABLE OF INTEREST CHECK 
  # IF THE FILE ALREADY EXISTS
  # If it does then go to the next one...
  fileexist=os.path.exists(file_out)
  filewrite=False
  if overwrite=='False':
    overwrite=False

  print '  --> OUTPUT FILE:'
  print '         ', file_out
  if fileexist==True:
    if overwrite==False:
      print '          +++ FILE ALREADY EXISTS +++'
      filewrite=False
    else:
      print '           +++ FILE EXISTS AND WILL BE OVERWRITTEN +++'
      filewrite=True
  else:
    print '         +++ FILE DOES NOT EXISTS YET +++'
    filewrite=True
  return filewrite

#**************************************************************************************
def get_dates(year,month,day,hour,mins,time_step,n_timesteps):
  import datetime as dt
  """ Gives a dates vector starting on year/month/day/time with a total 
  of n_timesteps each time_steps in hours.
  """
  dates=[dt.datetime(year,month,day,hour,mins)+ \
     dt.timedelta(hours=x) for x in xrange(0,n_timesteps*time_step,time_step)]

  return dates


# ***********************************************************
def date2hours(datelist,ref_date):
  hourssince=[(datelist[i]-ref_date).days*24. + (datelist[i]-ref_date).seconds/3600. for i in xrange(len(datelist))]
  return hourssince
  

# ***********************************************************
def get_yearsfile(fileall,varname):
  eyfile=np.asarray([int(fileall[i].split('_%s.nc' %(varname))[0][-4:]) for i in xrange(len(fileall))])
  syfile=np.asarray([int(fileall[i].split('_%s.nc' %(varname))[0][-9:-5]) for i in xrange(len(fileall))])
  return syfile,eyfile


# ***********************************************************
def create_dailyfiles(gvars,varname,stat_all):
  syp=gvars.syear
  fullpathout=create_outdir(gvars)
  fileall=sorted(glob.glob('%s/%s0?H_*_%s.nc' %(fullpathout,gvars.outfile_patt,varname)))
  fileref=nc.Dataset(fileall[0],'r')
  syfile,eyfile=get_yearsfile(fileall,varname)
  while syp<gvars.eyear:
    ctime_var=checkpoint(0)
    eyp=((int(syp)/5)+1)*5
    loadfile=False
    for stat in stat_all:
      if varname=='pracc':
        varstat=varname
      else:
        varstat=varname+stat
      
      file_out='%s/%sDAY_%s-%s_%s.nc' % (fullpathout,gvars.outfile_patt,syp,eyp-1,varstat) # Specify output file
      filewrite=checkfile(file_out,gvars.overwrite)
      if filewrite==True:
        loadfile=True
    if loadfile==True:

      print 'PROCESSING PERIOD %s-%s for variable %s' %(syp,eyp,varname)
      sel_files=[fileall[i] for i in xrange(len(syfile)) if ((syfile[i]>=syp) & (eyfile[i]<eyp))]
      files=nc.MFDataset(sel_files)
      time=nc.num2date(files.variables['time'][:],units=files.variables['time'].units)
      var=files.variables[varname][:]
      ctime=checkpoint(ctime_var)

      for stat in stat_all:
        ctime_var=checkpoint(0)
        if varname=='pracc':
          varstat=varname
        else:
          varstat=varname+stat
        file_out='%s/%sDAY_%s-%s_%s.nc' % (fullpathout,gvars.outfile_patt,syp,eyp-1,varstat) # Specify output file
        filewrite=checkfile(file_out,gvars.overwrite)
        if filewrite==True:
          dvar,dtime=coms.compute_daily(var,time,stat)
          dtime_nc=date2hours(dtime,gvars.ref_date)
          time_bnds=create_timebnds(dtime_nc)
          varatt={}
          for att in fileref.variables[varname].ncattrs():
            varatt[att]=getattr(fileref.variables[varname],att)

          # INFO NEEDED TO WRITE THE OUTPUT NETCDF
          netcdf_info=[file_out, varstat, varatt, True]

          # CREATE NETCDF FILE
          print 'yes'
          create_netcdf(netcdf_info,gvars, dvar, dtime_nc, time_bnds)
          ctime=checkpoint(ctime_var)
          print '=====================================================', '\n', '\n', '\n'
    syp=eyp


# ***********************************************************
def create_monthlyfiles(gvars,varname,stat_all):
  syp=gvars.syear
  fullpathout=create_outdir(gvars)
  for stat in stat_all:
    if varname=='pracc' or varname[-4:]=='step':
      varstat=varname
    else:
      varstat=varname+stat
    fileall=sorted(glob.glob('%s/%sDAY_*_%s.nc' %(fullpathout,gvars.outfile_patt,varstat)))
    print '%s/%sDAY_*_%s.nc' %(fullpathout,gvars.outfile_patt,varstat)
    fileref=nc.Dataset(fileall[0],'r')
    syfile,eyfile=get_yearsfile(fileall,varstat)
    syp=gvars.syear
    while syp<gvars.eyear:
      print 'start period year:', syp
      ctime_var=checkpoint(0)
      eyp=((int(syp)/10)+1)*10
      print 'PROCESSING PERIOD %s-%s for variable %s' %(syp,eyp,varstat)
      sel_files=[fileall[i] for i in xrange(len(syfile)) if ((syfile[i]>=syp) & (eyfile[i]<eyp))]
      files=nc.MFDataset(sel_files)
      time=nc.num2date(files.variables['time'][:],units=files.variables['time'].units)
      var=files.variables[varstat][:]
      ctime=checkpoint(ctime_var)
      file_out=fullpathout+'/%sMON_%s-%s_%s.nc' % (gvars.outfile_patt,syp,eyp-1,varstat) # Specify output file
      filewrite=checkfile(file_out,gvars.overwrite)
      if filewrite==True:
        mvar,mtime=coms.compute_monthly(var,time,stat)
        mtime_nc=date2hours(mtime,gvars.ref_date)
        mtime_bnds_inf=date2hours([dt.datetime(mtime[i].year,mtime[i].month,1,0,0,0) for i in xrange(len(mtime))],gvars.ref_date)
        mtime_bnds_sup=date2hours([dt.datetime(mtime[i].year,mtime[i].month,1,0,0,0)+relativedelta(months=1) for i in xrange(len(mtime))],gvars.ref_date)
        mtime_bnds=np.reshape(np.concatenate([mtime_bnds_inf,mtime_bnds_sup],axis=0), (len(mtime),2),order='F')
        varatt={}

        for att in fileref.variables[varstat].ncattrs():
          varatt[att]=getattr(fileref.variables[varstat],att)
        
        # INFO NEEDED TO WRITE THE OUTPUT NETCDF
        netcdf_info=[file_out, varstat, varatt, True]

        # CREATE NETCDF FILE
        create_netcdf(netcdf_info,gvars,mvar, mtime_nc, mtime_bnds)
        ctime=checkpoint(ctime_var)
        print '=====================================================', '\n', '\n', '\n'
      syp=eyp



# ***********************************************************
def intersect(a, b):
   return list(set(a) & set(b))


# ***********************************************************
def file_list(gvars, per, per_f, filet,n_files):
  tot_files=0
  print '\n', ' -> PROCESSING PERIOD: ', str(per)+' - '+str(per_f)
  #  files_in=[]
  #  for pp,pper in enumerate(np.arange(per,per_f+1)):
  #  loadfiles = gvars.pathin+'%s_%s_%s*' % (filet,gvars.domain,pper) # Specify path
  #  files_in.append(sorted(glob.glob(loadfiles)))
  #  tot_files=len(files_in[pp])+tot_files
  #  print '  -->  Number of files to read:', tot_files 

            
  loadfiles = sorted(glob.glob(gvars.pathin+'%s_%s_%s*' % (filet,gvars.domain,per))) # Specify path
  files_in=loadfiles
  for pp,pper in enumerate(np.arange(per+1,per_f+1)):
    loadfiles = sorted(glob.glob(gvars.pathin+'%s_%s_%s*' % (filet,gvars.domain,pper))) # Specify path
    files_in=np.concatenate((files_in,loadfiles))
  
  print '  -->  Number of files to read:', len(files_in)  ,n_files
  
  # CHECKING: Check if the number of files is right
  if len(files_in)!=n_files:
    print '\n', 'ERROR: the number of ',filet, ' files in period ', per,' is INCORRECT'
    print ' ---- SOME FILES ARE MISSING ---'
    print 'SCRIPT stops running ','\n' 
    sys.exit(0)
  
  return list(files_in)
  
# ***********************************************************
def read_list(files_list,var):
  from joblib import Parallel, delayed
  ctime_read=checkpoint(0)

  print '  -->  READING FILES '
  wrfvar=(getwrfname(var)[0]).split('-')

           
  if len(files_list)<=15:
    method='MFDataset'
  else:
    method='Dataset'

  
  # ---------------------
  if method=='MFDataset':
    print files_list
    fin=nc.MFDataset(files_list) # Read all files
    print '   -->   EXTRACTING VARIABLE Time'
    time = fin.variables['Times'][:] # Get time variable
  
    print '   -->   EXTRACTING VARIABLE ',var
    varvals=get_wrfvars(wrfvar,fin)
    fin.close()
   # ---------------------


  # ---------------------
  if method=='Dataset':
    varvals={}

    njobs=10
    nlen = len(files_list)/njobs #time step block length
    a=len(files_list)-njobs*nlen 
    nt_v=np.zeros(njobs)
    nt_v[:]=nlen
    nt_v[njobs-1]=nlen+a #block length for each job
    nt_v=nt_v.cumsum()
  
    files_in = [(files_list[0:int(nt_v[0])],wrfvar)]
    for tt in np.arange(1,njobs):
      files_in.append((files_list[int(nt_v[tt-1]):int(nt_v[tt])], wrfvar))
  
  
    var_v = Parallel(n_jobs=njobs)(delayed(read_block)(*files_in[i]) for i in xrange(len(files_in)))

    for i in np.arange(0,njobs):
      if i==0:
        time=var_v[i][0]
        for ii,wrfv in enumerate(wrfvar):
          varvals[wrfv]=var_v[i][1][wrfv]
      else:
        time=np.concatenate((time,var_v[i][0]))
        for ii,wrfv in enumerate(wrfvar):
          varvals[wrfv]=np.concatenate((varvals[wrfv],var_v[i][1][wrfv]))
  # ---------------------

  ctime=checkpoint(ctime_read)
  return np.asarray(time), varvals

#**************************************************************************************
def read_block(files_in,wrfvar):
  """ Extract wrfvar and time variables from the list of files files_in.
      The output is an array with times and a dictionary caonting arrays with
      the different variables.
  """  
  fin=nc.MFDataset(files_in) # Read all files
  temptime = fin.variables['Times'][:] # Get time variable
  tempvar=get_wrfvars(wrfvar,fin)
  fin.close()

  return temptime, tempvar
  
#**************************************************************************************
def read_block_Dataset(files_in,wrfvar):
  """ Extract wrfvar and time variables from the list of files files_in.
      The output is an array with times and a dictionary caonting arrays with
      the different variables.
      It does the same as read_block but using Dataset instead of MFDatset. It might
      be slower than read_block but it has the advantage to work even if some files 
      of the list do not have the same number of variables.
  """  
  tempvar={}
  for ff,ifile in enumerate(files_in):
    fin=nc.Dataset(ifile) # Read all files
    temptimet = fin.variables['Times'][:] # Get time variable
    tempvart=get_wrfvars(wrfvar,fin)
    fin.close()  

    if ff==0:
      temptime=temptimet
      tempvar=tempvart

    else:
      temptime=np.concatenate((temptime,temptimet))
      for ii,wrfv in enumerate(wrfvar):
        tempvar[wrfv]=np.concatenate((tempvar[wrfv],tempvart[wrfv]))


  return temptime, tempvar

  
#**************************************************************************************
def get_filefreq(filet):
  """ Method to get information about the type of file
      filet: type of file (wrfhrly, wrfout,wrfxtrm, wrfdly)
      ---
      nfiles: number of files per period
      time_step: hours between two time steps
      file_freq: string defining the frequency
      tbounds: whether the output file should include time bounds.
  """
  if filet=='wrfhrly':
    n_files=12      
    time_step=1 #hours between two time steps
    file_freq='01H'
    tbounds=False
    perstep=1

  elif filet=='wrfout':
    n_files=0
    time_step=3 #hours between two time steps
    file_freq='03H'
    tbounds=False
    perstep=5

  elif filet=='wrfxtrm':
    n_files=60      
    time_step=24 #hours between two time steps
    file_freq='DAY'
    tbounds=True
    perstep=5

  elif filet=='wrfdly':
    n_files=60      
    time_step=24 #hours between two time steps
    file_freq='DAY'
    tbounds=True
    perstep=5
  
  else:  
    sys.exit('The file tipe %s is erroneus' %(filet))
    
    
  return n_files,time_step,file_freq,tbounds,perstep
