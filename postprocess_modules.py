import re
import sys
import os
import netCDF4 as nc

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
  p1000mb = 100000.
  rcp = Rd/cp
  tkelvin = 273.15
  missingval = 1.e+20


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
    
    print 'Variables that will be obtained from postprocessing:'
    print varnames   
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

    dic={'tas':['T2'],'pracc':['RAINC-RAINNC'],'ps':['PSFC'], 'uas':['U10'], 'vas':['V10'],\
           'huss':['Q2'], 'wss':['U10-V10'], 'sst':['SST'],'rsds':['SWDOWN'], 'rlds':['GLW'],\
           'emiss':['EMISS'], 'albedo':['ALBEDO'], 'hfls':['LH'], 'hfss':['HFX'],'evspsbl':['SFCEVP'],\
           'mrso':['SMSTOT-DZS'], 'potevp':['POTEVP'], 'rlus':['TSK-EMISS'],'tasmeantstep':['T2MEAN'],\
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
    
    fileschemes=open("./info_files/WRF_schemes.inf")
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
    att={}
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
  import numpy as np
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
def add_leap(varval,index):
  import numpy as np
  """Method that add missing values to the 29th February in leap year for
  those simulations that are run without leap years.
  
  varval: 3-D matrix with time, latitude and longitude in the first, second and third dimensions. 
  The time dimension has one year of data without 29th february in a leap-year case.
  
  index: is the index of the 29th feburary in the time dimension.
  """
  
  temp=np.zeros((varval.shape[0]+1,varval.shape[1],varval.shape[2]))
  temp[0:index[0][0]-1,:,:]=varval[0:index[0][0]-1,:,:]
  temp[index[0][0]:index[0][23],:,:]=const.missingval
  temp[index[0][23]+1:,:,:]=varval[0:index[0][0]:,:,:]
  del varval
  varval=temp
  del temp 

  return varval



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
    
    glatt={}
    glatt['institution']				= "University of New South Wales (Australia)" 
    glatt['contact']					= "jason.evans@unsw.edu.au"
    glatt['institute_id']				= "CCRC"
    glatt['institute']					= "Climate Change Research Centre"
    glatt['references']					= "http://www.ccrc.unsw.edu.au"
    glatt['model_id']					= "WRF"
    glatt['RCM_version_id'] 			= "v3.3"
    glatt['RCM_source']					= "WRF 3.3 modified at the Climate Change Research Centre"	
    glatt['CORDEX_domain']              = "AUS44"
    glatt['project_id']					= "NARCliM"
    glatt['title']						= "Projection run"
    glatt['experiment_id']				= "projection"
    glatt['experiment']					= "Projection run with %s" %(GCM)
    glatt['driving_experiment'] 		= "%s, projection, %s%s%s" %(GCM,RCM,GCM,perturb)
    glatt['driving_model_id']			= "%s%s" %(GCM,perturb)
    glatt['driving_model_ensemble_member'] = "%s%s%s" %(RCM,GCM,perturb)
    glatt['wrf_options'] 		= "sst_update & tmn_update"
    glatt['driving_experiment_name'] 	= "projection_%s_%s" %(GCM,RCM) 
    glatt['product']					= "output"
    glatt['creation_date'] = dt.datetime.utcnow().strftime("%Y/%m/%d %H:%M:%S UTC")
    glatt['wrf_schemes_ra_lw_physics']            = "%s" %(sch_info['ra_lw_physics']) 
    glatt['wrf_schemes_sf_sfclay_physics']        = "%s" %(sch_info['sf_sfclay_physics']) 
    glatt['wrf_schemes_cu_physics']           = "%s" %(sch_info['cu_physics'])  
    glatt['wrf_schemes_bl_pbl_physics']       = "%s" %(sch_info['bl_pbl_physics']) 
    glatt['wrf_schemes_ra_sw_physics']            = "%s" %(sch_info['ra_sw_physics']) 
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
def create_netcdf(info, varval, time, time_bnds):
        

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
        calendar=info[3]
        domain=info[4]
        wrf_file_eg=info[5]
        GCM=info[6]
        RCM=info[7]
        time_bounds=info[8]

        # **********************************************************************
        # Read attributes from the geo_file of the corresponding domain
        file10='/srv/ccrc/data18/z3393242/studies/domains/NARCliM/geo_em.'+domain+'.nc'
        fin1=nc.Dataset(file10,mode='r')
        lon=np.squeeze(fin1.variables['XLONG_M'][:,:]) # Getting longitude
        lat=np.squeeze(fin1.variables['XLAT_M'][:,:]) # Getting latitude
        dx=getattr(fin1, 'DX')
        dy=getattr(fin1, 'DY')
        cen_lat=getattr(fin1, 'CEN_LAT')
        cen_lon=getattr(fin1, 'CEN_LON')
        pole_lat=getattr(fin1, 'POLE_LAT')
        pole_lon=getattr(fin1, 'POLE_LON')
        stand_lon=getattr(fin1, 'STAND_LON')
        fin1.close()

      
	#**********************************************************************
	# CREATING NETCDF FILE
	# Create output file
	fout=nc.Dataset(file_out,mode='w', format='NETCDF4_CLASSIC')

	# ------------------------
        # Create dimensions
        print '   CREATING AND WRITING DIMENSIONS: '
        print '                        TIME, X, Y(, TIME_BNDS)'
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
        setattr(varout, 'units','hours since 1949-12-01 00:00:00')
        setattr(varout, 'calendar',calendar)

        # VARIABLE: variable
        print '    ---   ',varname, ' VARIABLE CREATED ' 
        varout=fout.createVariable(varname,'f',['time', 'y', 'x'], fill_value=varatt['_FillValue'])
        varout[:]=varval[:]
        for att in varatt.keys():
          if att!='_FillValue':
            if varatt[att]!=None:
              setattr(varout, att, varatt[att])
            
        # VARIABLE: time_bnds 
        if time_bounds==True:
          print '    ---   TIME_BNDS VARIABLE CREATED ' 
          varout=fout.createVariable('time_bnds','f',['time', 'bnds'])
          varout[:]=time_bnds[:]
          setattr(varout, 'units','hours since 1949-12-01 00:00:00')
          setattr(varout, 'calendar',calendar)
        
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
        print '\n', '   CREATING AND WRITING GLOBAL ATTRIBUTES:'
        sch_info=read_schemes(wrf_file_eg)
        gblatt = get_globatt(GCM,RCM,sch_info)
        for att in gblatt.keys():
          setattr(fout, att, gblatt[att])
        
   
	fout.close()
        print '  ===> FILE: ', file_out
				
        return '         ------------  SUCCESFULLY CREATED!!!    ------------ '


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

