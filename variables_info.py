

class VariablesInfo(object):

    RAW_SOURCE = {
        'wrfhrly': {'tas'          :{'D': ('mean','min','max'), 'M': ('mean','min','max', 'minmean', 'maxmean')}, 
                    'wss'          :{'D': ('mean','max'), 'M': ('mean','max')},
                    'pracc'        :{'D': ('acc' ,), 'M': ('acc' ,)}, 
                    'prcacc'       :{'D': ('acc' ,), 'M': ('acc' ,)},
                    'prncacc'      :{'D': ('acc' ,), 'M': ('acc' ,)},
                    'ps'           :{'D': ('mean',), 'M': ('mean',)}, 
                    'huss'         :{'D': ('mean',), 'M': ('mean',)},
                    'hurs'         :{'D': ('mean',), 'M': ('mean',)},  
                    'uas'          :{'D': ('mean',), 'M': ('mean',)},  
                    'vas'          :{'D': ('mean',), 'M': ('mean',)},
                    },
        'wrfout' : {'evspsbl'      :{'D': ('mean',), 'M': ('mean',)}, 
                    'mrso'         :{'D': ('mean',), 'M': ('mean',)}, 
                    'sst'          :{'D': ('mean',), 'M': ('mean',)}, 
                    'potevp'       :{'D': ('mean',), 'M': ('mean',)}, 
                    'rsds'         :{'D': ('mean',), 'M': ('mean',)}, 
                    'rlds'         :{'D': ('mean',), 'M': ('mean',)}, 
                    'hfls'         :{'D': ('mean',), 'M': ('mean',)},
                    'hfss'         :{'D': ('mean',), 'M': ('mean',)},  
                    'albedo'       :{'D': ('mean',), 'M': ('mean',)}, 
                    'emiss'        :{'D': ('mean',), 'M': ('mean',)}, 
                    'rlus'         :{'D': ('mean',), 'M': ('mean',)},
                    'clt'          :{'D': ('mean',), 'M': ('mean',)},
                    'snm'          ::{'D': ('mean',), 'M': ('mean',)},
                    'snc'          ::{'D': ('mean',), 'M': ('mean',)},
                    'snw'          ::{'D': ('mean',), 'M': ('mean',)},
                    'snd'          ::{'D': ('mean',), 'M': ('mean',)},
                    },
        'wrfdly' : {'pr5maxtstep'  :{'D': ('max',), 'M': ('max',)}, 
                    'pr10maxtstep' :{'D': ('max',), 'M': ('max',)}, 
                    'pr20maxtstep' :{'D': ('max',), 'M': ('max',)}, 
                    'pr30maxtstep' :{'D': ('max',), 'M': ('max',)}, 
                    'pr1Hmaxtstep' :{'D': ('max',), 'M': ('max',)}, 
                    'wss5maxtstep' :{'D': ('max',), 'M': ('max',)}, 
                    'wss10maxtstep':{'D': ('max',), 'M': ('max',)}, 
                    'wss20maxtstep':{'D': ('max',), 'M': ('max',)}, 
                    'wss30maxtstep':{'D': ('max',), 'M': ('max',)}, 
                    'wss1Hmaxtstep':{'D': ('max',), 'M': ('max',)}
                    },
        'wrfxtrm': {'tasmeantstep' :{'D': ('mean',), 'M': ('mean',)}, 
                    'tasmintstep'  :{'D': ('min' ,), 'M': ('min', 'minmean')}, 
                    'tasmaxtstep'  :{'D': ('max' ,), 'M': ('max', 'maxmean')}
                    },
        }
        
        
        
        
    def get_wrf_file_types(self):
        return self.RAW_SOURCE.keys()
        
    def get_variables(self, wrf_file_type):
        """ get the list of configured variables
        """
        return self.RAW_SOURCE[wrf_file_type].keys()
       
    def get_monthly_variables(self, wrf_file_type):
        return self.get_variables(wrf_file_type)
        
    def get_daily_variables(self, wrf_file_type):
        return self.get_variables(wrf_file_type)
       
    def is_supported(self, vname):

        for raw_type, vnames in self.RAW_SOURCE.items():
            if vname in vnames.keys():
                return True

        return False

    def get_daily_variable_stats(self, wrf_file_type, vname):
        """get the configured daily statistics to compute for a given variable 'vname'
        """
        return self.RAW_SOURCE[wrf_file_type][vname]['D']
        
    def get_monthly_variable_stats(self, wrf_file_type, vname):
        """get the configured monthly statistics to compute for a given variable 'vname'
        """
        return self.RAW_SOURCE[wrf_file_type][vname]['M']       


    def get_source_variables_for_daily_stat(self, vname, stat):
        """ get the name of the variable from which to compute the daily stat
            
            returns: the name of the source sub daily stat from which to process the variable
                     and the name of the daily target stat
        """
        if vname in ['pracc','prcacc','prncacc']:
            return vname,vname 
            

        return vname+stat, vname+stat 


    def get_source_variables_for_monthly_stat(self, vname, stat):
        """ get the name of the variable from which to compute the monthly stat
            
            returns: the name of the source daily stat from which to process the variable
                     and the name of the monthly target stat
         """
        if vname[-4:]=='step':
          if stat in ['maxmean','minmean']:
            return vname, vname + 'mean'
          else:
            return vname, vname + stat
        
        if stat == 'maxmean':
            return vname + 'max', vname + 'maxmean'
        
        if stat == 'minmean':
            return vname + 'min', vname + 'minmean'

        if vname in ['pracc','prcacc','prncacc']:
            return vname,vname
            
        
            

        return vname+stat, vname+stat



    ##############################################################################
    #
    #
    #   PRIVATE IMPL.
    #
    #
    ##############################################################################


    


