# ======================================================================
# PURPOSE
# =========
# For all runs, create tasmean_bc files from decadal tasminmean_bc and
# tasmaxmean_bc files, for the NARCliM project

# STATUS
# ========

# INPUT
# =======
# inputdir   => Full input directory where bias-corrected files are
#               located. String. Must end with a slash, and contain
#               bias-corrected files in
#               NARCliM format as [GCM]/[RCM]/[period]/[d01/d02]. Example:
#              NNRP/R2/1950-2010/d01/CCRC_NARCliM_MON_1950-1959_tasminmean_bc.nc
# outputdur  => Full output directory where to place the
#               file. String. Can be the
#               same as input directory. If different, has to be an
#               empty directory structure of the same form as the
#               input directory structure. 
# author     => Author for the globals of output files. String. 
# email      => Author email for the globals of output files. String. 

# CALLS
# =======
# make_tasmean_bc

# OUTPUT
# ========
# tasmean_bc files in the output directory structure. For example:  
# NNRP/R2/1950-2010/d01/CCRC_NARCliM_MON_1950-1959_tasminmean_bc.nc

# REQUIRES
# =========

# HISTORY
# ========
# Jan 30 2015 => Written by Roman Olson, CCRC, UNSW.
# ======================================================================

def make_tasmean_bc_all(inputdir, outputdir, author, email):

   # IMPORT MODULES
   import os
   import glob
   from make_tasmean_bc import make_tasmean_bc

   # WALK THROUGH DIRECTORIES
   # Note that dirpath does not end with a slash
   for dirpath, dirnames, filenames in os.walk(inputdir):

      print dirpath

      # Only look into bottom directories ending in d01 or d02
      if (dirpath.endswith(("d01", "d02"))):

        reldirpath = os.path.relpath(dirpath, inputdir)
        print reldirpath
        minfiles = glob.glob(dirpath + "/CCRC_NARCliM_MON_*_tasminmean_bc.nc")
        maxfiles = glob.glob(dirpath + "/CCRC_NARCliM_MON_*_tasmaxmean_bc.nc")
        
        # Iterate across files for each directory. These are full names
        for minfile, maxfile in zip(minfiles, maxfiles):

           baseminfile = os.path.basename(minfile)
           basemaxfile = os.path.basename(maxfile)

           # Check if files are for the same period
           if (baseminfile[:26] == basemaxfile[:26]):

              basemeanfile = baseminfile[:26] + "_tasmean_bc.nc"
              meanfile = os.path.join(outputdir, reldirpath, basemeanfile)

              make_tasmean_bc(minfile, maxfile, meanfile, 
                 author, email)

           else:
              raise ValueError, "Min and max files have different times"


      
