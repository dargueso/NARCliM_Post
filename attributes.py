"""attributes.py
   Class that contains the global and variable attributes.
   To be used in the NARCliM postprocessing
   Created:13 September 2013
"""
import datetime as dt

class globalatt:
    """Class that generates the global attributes
    Defines the name of the schemes and the references to include in the global attributes

    GCM: Name of the GCM (e.g.:MIROC3.2)
    RCM: Name of the RCM (e.g.: R1)
    Period: Name of the perturbation of the GCM
    sch_info: dictionary containing the kind of schemes as keys, and the name and references as values
    e.g.: global_attributes=ga.globalatt("MIROC3.2","R1","p1",sch_info)
    """
    def __init__(self,GCM,RCM,sch_info,perturb=None):
        self.institution				= "University of New South Wales (Australia)" 
        self.contact					= "jason.evans@unsw.edu.au"
        self.institute_id				= "CCRC"
        self.institute					= "Climate Change Research Centre"
        self.references					= "http://www.ccrc.unsw.edu.au"
        self.model_id					= "WRF"
        self.RCM_version_id 			= "v3.3"
        self.RCM_source					= "WRF 3.3 modified at the Climate Change Research Centre"	
        self.CORDEX_domain              = "AUS44"
        self.project_id					= "NARCliM"
        self.title						= "Projection run"
        self.experiment_id				= "projection"
        self.experiment					= "Projection run with %s" %(GCM)
        self.driving_experiment 		= "%s, projection, %s%s%s" %(GCM,RCM,GCM,perturb)
        self.driving_model_id			= "%s%s" %(GCM,period)
        self.driving_model_ensemble_member = "%s%s%s" %(RCM,GCM,perturb)
        self.wrf_options 		= "sst_update & tmn_update"
        self.driving_experiment_name 	= "projection_%s_%s" (GCM,RCM) 
        self.product					= "output"
        self.creation_date = dt.datetime.utcnow().strftime("%Y/%m/%d %H:%M:%S UTC")
        self.wrf_schemes_ra_lw_physics 			= "%s" %s(sch_info['ra_lw_physics']) 
        self.wrf_schemes_sf_sfclay_physics 		= "%s" %s(sch_info['sf_sfclay_physics']) 
        self.wrf_schemes_cu_physics 			= "%s" %s(sch_info['cu_physics'])  
        self.wrf_schemes_bl_pbl_physics 		= "%s" %s(sch_info['bl_pbl_physics']) 
        self.wrf_schemes_ra_sw_physics 			= "%s" %s(sch_info['ra_sw_physics ']) 
        self.wrf_schemes_sf_surface_physics 	= "%s" %s(sch_info['sf_surface_physics']) 

		