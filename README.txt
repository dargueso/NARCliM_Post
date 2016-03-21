INSTRUCTIONS FOR THE NEW POSTPROCESSING
Authors: Alejandro Di Luca and Daniel Argueso
Created: 27 September 2013
Institution: CCRC, UNSW. (Australia)

NOTE: Updated WRF_schemes.inf is good for v3.6

In order to use the new post-process, please read the following instructions.

1. Make sure your environment is adequate: It has been tested in python 2.7.5, so better use that version. It also needs joblib (which is installed in the new storm servers) for parallelisation of the code.


2. Get the latest version from the git repository: https://dargueso@bitbucket.org/dargueso/narclim_scripts.git
    It includes various NARCliM-related scripts and a directory named 'NARCliM_postprocess' which host all necessary files for the new postproces:
    postprocess_NARCliM.py: The main program
    postprocess_modules.py: File containing the postprocess methods called in the main program.
    compute_vars.py: File containing methods to compute each of the required variables from WRF variables
    compute_stats.py: File containing methods to calculate daily and monthly statistics
    WRF_schemes.inf: File with information regarding the available schemes in WRF (used to write the global attributes)
    variables.inf: File with information of the variables that are processed (name, frequency, WRF file from where it is retrieved…)
    
    And the most relevant to the user: NARCliM_post.input - The input namelist with all options to run the postprocessing.

3. Create your own NARCliM_post.input for your particular simulation. This is a sample, you can name it as you like (e.g. NARCliM_post_MIROC3.2-R3-1990-2010.input). The parameters you have to edit within the namelist are:
         pathin: where the original files are located
         pathout: where the postprocessed files will be written. Please note that this is the root where the necessary directories will be created depending on the GCM, RCM, period, domain. So the final output directory will be: pathout/period/GCM/RCM/domain/postprocess/
         GCM: The GCM driving the simulation (MIROC3.2, ECHAM5, CCCMA3.1, MK3.0)
         RCM: The RCM of the simulation (R1, R2, R3)
         syear: first year to be post processed
         eyear: last year to be post processed
         outfile_patt: the pattern of the output files. For NARCliM: CCRC_NARCliM_
         overwrite: Whether the existing files will be overwritten or not.

          A line that separates the options from the variables and should not be modified: #### Requested output variables (DO NOT CHANGE THIS LINE) ####

         The variables that will be post processed. Do not remove variables, it is enough to comment them using #.
 
4. Run the script from the folder where the scripts live: python postprocess_NARCliM.py -i [namelist_input]
    For example:  python postprocess_NARCliM.py -i NARCliM_post_MIROC3.2-R3-1990-2010.input

5. The postprocess is ready to generate all NARCliM variables. It provides yearly files for the 01H, 5-year files for the 03H and daily statistics, and 10-year files for monthly statistics. It also generates a log in the same output folder as the postprocessed files named as:
postprocess_[GCM]_[RCM]_[SYEAR]-[EYEAR]_[DOMAIN]_[DATE_OF_CREATION].log 
     The date of creation is the time when the postprocessing started, so it doesn't overwrite previous log files.
