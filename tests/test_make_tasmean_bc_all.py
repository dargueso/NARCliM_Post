# To test make_tasmean_bc_all

import pdb
import unittest
import numpy as np
from make_tasmean_bc_all import make_tasmean_bc_all
from ccrc_utils import filediff

class test_make_tasmean_bc_all(unittest.TestCase):
   
   # Set-up. This is done prior to each test.
   def setUp(self):
      self.biasdir = "/srv/ccrc/data30/z3393020/NARCliM/Bias_corrected/"
      self.outdir = "/home/z3481416/hdrive/work/narclim_scripts_git/\
narclim_scripts/NARCliM_postprocess/tests/bcor/"
      
   # Tear-down. This is done after each test
   def TearDown(self):
      pass

   # Test run
   def test_run(self):
      make_tasmean_bc_all(self.biasdir, self.outdir, "Roman Olson", 
         "roman.olson@unsw.edu.au")
      







suite = unittest.TestLoader().loadTestsFromTestCase(test_make_tasmean_bc_all)
unittest.TextTestRunner(verbosity=2).run(suite)
