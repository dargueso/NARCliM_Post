# *************************************************************************************
def create_netcdf(filename, var, lat, lon, times, rotpole, varatt, overwrite=None):
        import sys

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

	#***************************************************************************************
	# CREATING NETCDF FILE
	# Create output file
	file_out='%s%s-%s-%s_%s_%s%s%s'%(path_out,sname,ss,ts,varname_out[var],str(year),Projection,'.nc')
	fout=ncdf(file_out,mode='w', format='NETCDF4_CLASSIC')

	# Create dimensions
	ldims_out=temp.shape			
	jj=1
	for dim in dims_out:
		print '      Processing ', '%s' %(dim),' dimension: ', ldims_out[jj]
		fout.createDimension('%s' %(dim),ldims_out[jj] )
		jj=jj+1
	fout.createDimension('time',None)

	# Create variables
	outvar=fout.createVariable(varname_out[var],vartypes[var],['time','latitude','longitude'],\
					   fill_value=-32767)
	outvar[:]=temp


	# Copy attributes from invar and create some new ones:
	setattr(outvar, 'units', 'hPa')
	setattr(outvar, 'scale_factor', 1.0)
	setattr(outvar, 'add_offset', 0)
	print '     Adding Attributes: ', fin.variables[d['varname_in']].ncattrs()
	for att in fin.variables[d['varname_in']].ncattrs():

		if hasattr(outvar, att):
			pass
		else:
			setattr(outvar, att, getattr(fin.variables[d['varname_in']],att))

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

