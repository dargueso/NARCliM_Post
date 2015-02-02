# To test make_tasmean_bc_all
# Feb 02 2014 => Roman Olson, CCRC, UNSW.
# Program runs during set-up, and the results are compared to
# individual program calls in the tests. The individual file check is
# passed no matter what. However, the output goes to the screen and
# needs to be checked visually by the user. The good result is
# differences in creation date, etc., which is what to be expected. 
import pdb
import unittest
import numpy as np
from make_tasmean_bc_all import make_tasmean_bc_all
from make_tasmean_bc import make_tasmean_bc
from ccrc_utils import filediff
import subprocess
import os.path

class test_make_tasmean_bc_all(unittest.TestCase):
   
   # Set-up. This is done prior to each test.
   def setUp(self):
      self.biasdir = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected/"
      self.outdir = "/srv/ccrc/data30/z3481416/narclim_scripts_git/\
narclim_scripts/NARCliM_postprocess/tests/bcor/"
      make_tasmean_bc_all(self.biasdir, self.outdir, "Roman Olson", 
   "privatemail@email.privatemail")

      
   # Tear-down. This is done after each test. 
   def tearDown(self):
      pass
      print "Removing test files!"
      subprocess.call(['bash', 'rmtestfiles'])
      subprocess.call(['rm', 'make_tasmean_bc_all_test1_1.nc'])
      subprocess.call(['rm', 'make_tasmean_bc_all_test1_2.nc'])
      subprocess.call(['rm', 'make_tasmean_bc_all_test1_3.nc'])

   # Test1.
   def test1(self):
      pass
      # CHECK TOTAL NUMBER OF FILES AND CHECK THE FILENAMES
      GCMs = ('MIROC3.2', 'ECHAM5', 'CCCMA3.1', 'CSIRO-MK3.0')
      RCMs = ('R1', 'R2', 'R3')
      domains = ('d01', 'd02')
      periods = ('1990-2010', '2020-2040', '2060-2080')
      modelfilesok = []
      rnlsfilesok = []

      # GCM files
      print "================="
      print "Testing GCM files"
      print "================="
      for GCM in GCMs:
         for RCM in RCMs:
            for period in periods:
               for domain in domains:
                  curdir = os.path.join("bcor/", GCM, RCM, period, domain, "")
                  pstart = int(period[:4])

                  file1 = "CCRC_NARCliM_MON_" + str(pstart) + "-" +\
 str(pstart+9) + "_tasmean_bc.nc"
                  file2 = "CCRC_NARCliM_MON_" + str(pstart+10) + "-" +\
 str(pstart+19) + "_tasmean_bc.nc"
                  fullfile1 = os.path.join(curdir, file1)
                  fullfile2 = os.path.join(curdir, file2)
                  print fullfile1
                  print fullfile2
                  modelfilesok.append(os.path.isfile(fullfile1))
                  modelfilesok.append(os.path.isfile(fullfile2))
      gcmfilesallok = np.all(modelfilesok)
      gcmnumok = len(modelfilesok) == 144 

      # Reanalysis files
      rp_start = [1950, 1960, 1970, 1980, 1990, 2000]
      print "========================"
      print "Testing reanalysis files"
      print "========================"
      for RCM in RCMs:
         for domain in domains:
            curdir = os.path.join("bcor/NNRP", RCM, "1950-2010", domain, "")

            for rp in rp_start:
               myfile = "CCRC_NARCliM_MON_"  + str(rp) + "-" +\
 str(rp+9) + "_tasmean_bc.nc"
               myfullfile = os.path.join(curdir, myfile)
               rnlsfilesok.append(os.path.isfile(myfullfile))
               print myfullfile
      rnlsfilesallok = np.all(rnlsfilesok)
      rnlsnumok = len(rnlsfilesok) == 36 

      filesok = gcmfilesallok and gcmnumok and rnlsfilesallok and rnlsnumok
      self.assertTrue(filesok)#!+


      #Check 1: CCCMA3.1 R1 1990-2000 d01
      print "========================="
      print "Differences for CCCMA3.1 R1 1990-2000 d01"
      print "========================="
      minfile = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected/CCCMA3.1/\
R1/1990-2010/d01/CCRC_NARCliM_MON_1990-1999_tasminmean_bc.nc"
      maxfile = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected/CCCMA3.1/\
R1/1990-2010/d01/CCRC_NARCliM_MON_1990-1999_tasmaxmean_bc.nc"
      outfile = "bcor/CCCMA3.1/R1/1990-2010/d01/CCRC_NARCliM_MON_1990-1999_\
tasmean_bc.nc"
      outfile_save = "make_tasmean_bc_all_test1_1.nc"
      make_tasmean_bc(minfile, maxfile, outfile_save, "Roman Olson", 
"privatemail@email.privatemail")
      subprocess.call(["ncdiff", outfile, outfile_save])#!+

      #Check 2: CSIRO-MK3.0 R2 2070-2080 d02
      print "========================="
      print "Differences for CSIRO-MK3.0 R2 2070-2080 d02"
      print "========================="
      minfile = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected/CSIRO-MK3.0/\
R2/2060-2080/d02/CCRC_NARCliM_MON_2070-2079_tasminmean_bc.nc"
      maxfile = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected/CSIRO-MK3.0/\
R2/2060-2080/d02/CCRC_NARCliM_MON_2070-2079_tasmaxmean_bc.nc"
      outfile = "bcor/CSIRO-MK3.0/R2/2060-2080/d02/CCRC_NARCliM_MON_2070-2079_\
tasmean_bc.nc"
      outfile_save = "make_tasmean_bc_all_test1_2.nc"
      make_tasmean_bc(minfile, maxfile, outfile_save, "Roman Olson", 
"privatemail@email.privatemail")
      subprocess.call(["ncdiff", outfile, outfile_save])#!+


      #Check 3: NNRP R3 1970-1980 d02
      print "========================="
      print "Differences for NNRP R3 1970-1980 d02"
      print "========================="
      minfile = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected/NNRP/\
R3/1950-2010/d02/CCRC_NARCliM_MON_1970-1979_tasminmean_bc.nc"
      maxfile = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected/NNRP/\
R3/1950-2010/d02/CCRC_NARCliM_MON_1970-1979_tasmaxmean_bc.nc"
      outfile = "bcor/NNRP/R3/1950-2010/d02/CCRC_NARCliM_MON_1970-1979_\
tasmean_bc.nc"
      outfile_save = "make_tasmean_bc_all_test1_3.nc"
      make_tasmean_bc(minfile, maxfile, outfile_save, "Roman Olson", 
"privatemail@email.privatemail")
      subprocess.call(["ncdiff", outfile, outfile_save])#!+

   # # Test run
   # def test_run(self):
   #    make_tasmean_bc_all(self.biasdir, self.outdir, "Roman Olson", 
   #       "roman.olson@unsw.edu.au")
      
suite = unittest.TestLoader().loadTestsFromTestCase(test_make_tasmean_bc_all)
unittest.TextTestRunner(verbosity=2).run(suite)
