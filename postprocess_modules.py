import re
import sys
import os


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
<<<<<<< .merge_file_v6RVqS
    
    options,sentinel,varnames=filein.read().partition('#### Requested output variables (DO NOT CHANGE THIS LINE) ####')
    fileopts=open('fileopts.tmp','w')
    fileopts.write(options)
    filevarnames=open('filevarnames.tmp','w')
    filevarnames.write(varnames)
    
    
    fileopts=open('fileopts.tmp','r')
    lines=fileopts.readlines()
=======
    lines=filein.readlines()
>>>>>>> .merge_file_lBdPWd
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
    for ii in xrange(len(entryname)):
        inputinf[entryname[ii]]=entryvalue[ii]
<<<<<<< .merge_file_v6RVqS
        
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
    os.remove('filevarnames.tmp')
    os.remove('fileopts.tmp') 
    
    print 'Variables that will be obtained from postprocessing:'
    print varnames   
    return inputinf,varnames
=======
    return inputinf
>>>>>>> .merge_file_lBdPWd
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
    	        print li
                values=li.split(' ')
                varname.append(values[0])
                filetype.append(values[1])
                freqreq.append(values[2])
                statsreq.append(values[3])

    filein.close()


    varinfo=dictionary2entries(filetype,varname,statsreq)
    return varinfo

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

def create_netcdf(info, varval, lat, lon, time, overwrite=None):
        

        """ Create a netcdf file for the post-processed variables of NARCliM simulations
                   
	By default, the module do not overwrite the file so if there is one with the same name then 
	the script does not do anything.

        Input: global  attributres from pre-defined classes:   
        Output: a netcdf file
        Author: Alejandro Di Luca, Daniel Argueso
        Created: 14/09/2013
        Last Modification: 07/08/2013
        
        """
        import numpy as np
        import netCDF4 as nc

        filename=info[0]
        varname=info[1]
        varatt=info[2]

        # **********************************************************************
        # Read attributes from the geo_file of the corresponding domain
        file10='/srv/ccrc/data18/z3393242/studies/domains/NARCliM/geo_em.d02.nc'
        fin1=nc.Dataset(file10,mode='r')
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
        print '   Creating DIMENSIONS '
        fout.createDimension('bnds', 2)
        fout.createDimension('time',None)
        fout.createDimension('y',varval.shape[2])
        fout.createDimension('x',varval.shape[1])

        # ------------------------
        # Create and assign values to variables
        print "\n"
        print 'Creating variables: ', fin.variables.keys() 
	
        # VARIABLE 1: Rotated_Pole 
        varout=fout.createVariable('Rotated_pole','c',[])
        setattr(varout, 'grid_mapping_name', 'rotated_latitude_longitude')
        setattr(varout, 'dx_m', dx)
        setattr(varout, 'dy_m', dy)
        setattr(varout, 'latitude_of_projection_origin', cen_lat)
        setattr(varout, 'longitude_of_central_meridian',cen_lon)
        setattr(varout, 'true_longitude_of_projection',stand_lon)
        setattr(varout, 'grid_north_pole_latitude',  pole_lat)
        setattr(varout, 'grid_north_pole_longitude', pole_lon)

        # VARIABLE 2: time_bnds 
        varout=fout.createVariable('time_bnds','d',['time', 'bnds'])
        varout[:]=time_bnds[0:time_dim]

        # VARIABLE 3: time 
        varout=fout.createVariable('time','d',['time'])
        varout[:]=time[:]
        setattr(varout, 'standard_name','time')
        setattr(varout, 'long_name','time')
        setattr(varout, 'bounds','time_bnds')
        setattr(varout, 'units','hours since 1949-12-01 00:00:00')
        setattr(varout, 'calendar','standard')


        # VARIABLE 4: lon
        varout=fout.createVariable('lon','f',['y', 'x'])
        varout[:]=lon[:]
        setattr(varout, 'standard_name','longitude')
        setattr(varout, 'long_name','Longitude')
        setattr(varout, 'units','degrees_east')
        setattr(varout, '_CoordinateAxisType','Lon')
        
        # VARIABLE 5: lat
        varout=fout.createVariable('lat','f',['y', 'x'])
        varout[:]=lat[:]
        setattr(varout, 'standard_name','latitude')
        setattr(varout, 'long_name','Latitude')
        setattr(varout, 'units','degrees_north')
        setattr(varout, '_CoordinateAxisType','Lat')
           
        # VARIABLE 6: variable
        varout=fout.createVariable(varname,'f',['time', 'y', 'x'])
        varout[:]=varval
        # Copy attributes from varatt:
        for att in varatt:
            setattr(varout, att, getattr(fin.variables[var],att))


	# Add global attributes
	for att in fin.ncattrs():
		if hasattr(fout, att): 
			pass
		else:
			setattr(fout, att, getattr(fin,att))

		setattr(fout,"temporal_scale_hours", ts)
		setattr(fout,"spatial_scale_km", ss)
		setattr(fout,"orig_temporal_scale_hours", d['TempRes'])				
		setattr(fout,"orig_spatial_scale_km", d['SpatRes'])
		setattr(fout,"dataset", sname)
		setattr(fout,"input_dir", d['path'])	
		setattr(fout,"year", year)
		setattr(fout,"interpolation_method", 'griddata-cubic')
		setattr(fout,"Projection", Projection)
		setattr(fout,"Geofile",'/srv/ccrc/data23/z3444417/studies/Data/RefGridMeshes/GeoFiles/geo_em_D10'+Projection+'.nc')
		setattr(fout,"region_limits", file10)			
		setattr(fout,"author","Alejandro Di Luca @ CCRC, UNSW, Australia")
		setattr(fout,"comment","Variable extracted using extract_psl.py.")
		setattr(fout,"date",datetime.today().strftime('%Y-%m-%d'))
		setattr(fout,"contact","a.diluca@unsw.edu.au")

	fout.close()				
	#****************************************************************************************

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

