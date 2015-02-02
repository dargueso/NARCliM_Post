# To test make_tasmean_bc
# Feb 02, 2015. Roman Olson, CCRC, UNSW.

import pdb
import unittest
import numpy as np
import numpy.ma as ma
import netCDF4 as ncdf
from make_tasmean_bc import make_tasmean_bc
from ccrc_utils import filediff
import subprocess

class test_make_tasmean_bc(unittest.TestCase):
   
   # Set-up. This is done prior to each test.
   def setUp(self):
      self.testfiledir = "/home/z3481416/hdrive/work/test_files/bias_corrected/"
      self.CCCMA_d01_tmin = self.testfiledir + "CCCMA3.1_d01_tasminmean.nc"
      self.CCCMA_d01_tmax = self.testfiledir + "CCCMA3.1_d01_tasmaxmean.nc"
      self.CCCMA_2030tmax = self.testfiledir + "CCCMA3.1_d01_2030-2039\
_tasmaxmean.nc"
      self.CCCMA_d02_tmax = self.testfiledir + "CCCMA3.1_d02_tasmaxmean.nc"
      self.NNRP_d02_tmin = self.testfiledir + "NNRP_d02_tasminmean.nc"
      self.NNRP_d02_tmax = self.testfiledir + "NNRP_d02_tasmaxmean.nc"
      
   # Tear-down. This is done after each test
   def tearDown(self):
      print "Running TearDown"
      subprocess.call(["rm", "bcor/out.nc"])
      subprocess.call(["rm", "bcor/outfile.nc"])

   # Failure cases
   def test1_except(self):
      # Different domains
      self.assertRaisesRegexp(ValueError, "Latitudes don't match", 
         make_tasmean_bc, self.CCCMA_d01_tmin, self.CCCMA_d02_tmax, "out.nc",
         "John Smith", "email@email.email")
      self.assertRaisesRegexp(ValueError, "Times don't match", 
         make_tasmean_bc, self.CCCMA_d01_tmin, self.CCCMA_2030tmax, "out.nc",
         "John Smith", "email@email.email")#!+
      

   # TEST WITH CCCMA3.1 R1 2020-2029 d01
   def test2(self):
      make_tasmean_bc(self.CCCMA_d01_tmin, self.CCCMA_d01_tmax, "bcor/out.nc", 
         "John Smith", "email@email.email")

      # Test coordinates
      ncnew    = ncdf.Dataset("bcor/out.nc")
      ncminold = ncdf.Dataset(self.CCCMA_d01_tmin)

      # Latitudes
      nclat    = ncnew.variables["lat"]
      nclatold = ncminold.variables["lat"]

      latdataok  = np.all(nclat[:] == nclatold[:])
      latshapeok = nclat[:].shape == (144, 215)
      latdictnew = nclat.__dict__.copy()
      latdictold = nclatold.__dict__
      latfillok  = latdictnew["_FillValue"] == 1e+20
      latdictnew.pop("_FillValue")
      latdictok = latdictold == latdictnew

      latok = latdataok and latshapeok and latfillok and latdictok
      self.assertTrue(latok) #!+

      # Longitudes
      nclon    = ncnew.variables["lon"]
      nclonold = ncminold.variables["lon"]

      londataok  = np.all(nclon[:] == nclonold[:])
      lonshapeok = nclon[:].shape == (144, 215)
      londictnew = nclon.__dict__.copy()
      londictold = nclonold.__dict__
      lonfillok  = londictnew["_FillValue"] == 1e+20
      londictnew.pop("_FillValue")
      londictok = londictold == londictnew

      lonok = londataok and lonshapeok and lonfillok and londictok
      self.assertTrue(lonok)
      
      # Time
      nctime    = ncnew.variables["time"]
      nctimeold = ncminold.variables["time"]

      timedataok  = np.all(nctime[:] == nctimeold[:])
      timeshapeok = nctime[:].shape == (120,)
      timedictnew = nctime.__dict__.copy()
      timedictold = nctimeold.__dict__
      timefillok  = timedictnew["_FillValue"] == 1e+20
      timedictnew.pop("_FillValue")
      timedictold.pop("bounds")
      timedictok = timedictold == timedictnew

      timeok = timedataok and timeshapeok and timefillok and timedictok
      self.assertTrue(timeok) 

      # Test data
      ncdatavar  = ncnew.variables["tasmean_bc"]
      ncdata     = ncdatavar[:]
      ncoldtmin  = ncminold.variables["tasminmean_bc"]
      ncdatamask = ma.getmaskarray(ncdata)

      self.assertTrue(ncdatamask[0,0,0] == True)
      self.assertTrue(ncdatamask[0, 39, 48] == True)
      self.assertTrue(ncdata.data[0,0,0] == 1E+20)
      self.assertTrue(ncdata.data[0, 39, 48] == 1E+20)
      self.assertAlmostEqual(ncdata[0, 40, 48], 295.214, places=3)
      self.assertAlmostEqual(ncdata[4, 44, 52], 288.413, places=3)
      self.assertAlmostEqual(ncdata[119, 101, 95], 302.404, places=3)
      datafillok  = ncdata.fill_value == 1E+20

      datashapeok = ncdata.shape == (120, 144, 215)
      datatypeok  = type(ncdata) == np.ma.core.MaskedArray
      datadtypeok = ncdata.dtype == np.float64

      datadictnew = ncdatavar.__dict__
      datadictold = ncoldtmin.__dict__
      dataattrshistok = datadictnew["history"] == "Created by averaging\
 tasminmean_bc (monthly\
 mean minimum temperature) and tasmaxmean_bc (monthly mean maximum temperature)\
."
      datacomok = datadictnew["comment"] == "tasmean_bc is not an original\
 NARCliM variable. It is calculated as an average of tasminmean_bc and\
 tasmaxmean_bc variables"
      datadictnew.pop("history")
      datadictnew.pop("comment")
      # Can't compare fill values since they have changed type
      self.assertEqual(datadictnew["_FillValue"], 1e+20)
      datadictnew.pop("_FillValue")
      datadictold.pop("_FillValue")
      dataattrsok = dict(datadictnew) == dict(datadictold)

      dataok = datafillok and datashapeok and datatypeok and datadtypeok \
 and dataattrshistok and datacomok and dataattrsok
      self.assertTrue(dataok) 

      
      # Test global attributes
      globnew = ncnew.__dict__.copy()
      globold = ncminold.__dict__.copy()

      histok = str(globnew["history"]).endswith("UTC created mean temperature\
 from tasminmean_bc and tasmaxmean_bc.")
      authorok = globnew["author"] == "John Smith"
      emailok  = globnew["contact"] == "email@email.email"
      comok = globnew["comments"] == "tasmean_bc is not an original NARCliM\
 variable. It is calculated as an average of tasminmean_bc and tasmaxmean_bc\
 variables"
      globnew.pop("comments")
      globold.pop("comments")
      globnew.pop("history")
      globold.pop("history")
      globnew.pop("author")
      globold.pop("author")
      globnew.pop("contact")
      globold.pop("contact")
      globnew.pop("creation_date")
      globold.pop("creation_date")

      # New file has no "Conventions" or "date"
      globold.pop("Conventions")
      globold.pop("date")

      globattrsok = dict(globnew) == dict(globold)
      globok = histok and authorok and emailok and comok and globattrsok

      self.assertTrue(globok) #!+
     

   # TEST WITH NNRP R3 1950-1959 d02
   def test3(self):
      make_tasmean_bc(self.NNRP_d02_tmin, self.NNRP_d02_tmax,"bcor/outfile.nc", 
         "Roman Olson", "privatemail@email.privatemail")

      # Test coordinates
      ncnew    = ncdf.Dataset("bcor/outfile.nc")
      ncminold = ncdf.Dataset(self.NNRP_d02_tmin)

      # Latitudes
      nclat    = ncnew.variables["lat"]
      nclatold = ncminold.variables["lat"]

      latdataok  = np.all(nclat[:] == nclatold[:])
      latshapeok = nclat[:].shape == (200, 325)
      latdictnew = nclat.__dict__.copy()
      latdictold = nclatold.__dict__
      latfillok  = latdictnew["_FillValue"] == 1e+20
      latdictnew.pop("_FillValue")
      latdictok = latdictold == latdictnew

      latok = latdataok and latshapeok and latfillok and latdictok
      self.assertTrue(latok) 

      # Longitudes
      nclon    = ncnew.variables["lon"]
      nclonold = ncminold.variables["lon"]

      londataok  = np.all(nclon[:] == nclonold[:])
      lonshapeok = nclon[:].shape == (200, 325)
      londictnew = nclon.__dict__.copy()
      londictold = nclonold.__dict__
      lonfillok  = londictnew["_FillValue"] == 1e+20
      londictnew.pop("_FillValue")
      londictok = londictold == londictnew

      lonok = londataok and lonshapeok and lonfillok and londictok
      self.assertTrue(lonok) #!+


      
      # Time
      nctime    = ncnew.variables["time"]
      nctimeold = ncminold.variables["time"]

      timedataok  = np.all(nctime[:] == nctimeold[:])
      timeshapeok = nctime[:].shape == (120,)
      timedictnew = nctime.__dict__.copy()
      timedictold = nctimeold.__dict__
      timefillok  = timedictnew["_FillValue"] == 1e+20
      timedictnew.pop("_FillValue")
      timedictold.pop("bounds")
      timedictok = timedictold == timedictnew

      timeok = timedataok and timeshapeok and timefillok and timedictok
      self.assertTrue(timeok) #!+

      # Test data
      ncdatavar  = ncnew.variables["tasmean_bc"]
      ncdata     = ncdatavar[:]
      ncoldtmin  = ncminold.variables["tasminmean_bc"]
      ncdatamask = ma.getmaskarray(ncdata)

      self.assertTrue(ncdatamask[0,0,0] == True)
      self.assertTrue(ncdatamask[0, 30, 144] == True)
      self.assertTrue(ncdata.data[0,0,0] == 1E+20)
      self.assertTrue(ncdata.data[0, 30, 144] == 1E+20)
      self.assertAlmostEqual(ncdata[0, 31, 144], 290.694, places=3)
      self.assertAlmostEqual(ncdata[4, 100, 99], 287.235, places=3)
      self.assertAlmostEqual(ncdata[119, 199, 0], 302.023, places=3)

      # Test for means
      ncmean = np.mean(ncdata, dtype=np.float64)
      self.assertAlmostEqual(ncmean, 292.101, places=3) 

      datafillok  = ncdata.fill_value == 1E+20

      datashapeok = ncdata.shape == (120, 200, 325)
      datatypeok  = type(ncdata) == np.ma.core.MaskedArray
      datadtypeok = ncdata.dtype == np.float64

      datadictnew = ncdatavar.__dict__
      datadictold = ncoldtmin.__dict__
      dataattrshistok = datadictnew["history"] == "Created by averaging\
 tasminmean_bc (monthly\
 mean minimum temperature) and tasmaxmean_bc (monthly mean maximum temperature)\
."
      datacomok = datadictnew["comment"] == "tasmean_bc is not an original\
 NARCliM variable. It is calculated as an average of tasminmean_bc and\
 tasmaxmean_bc variables"
      datadictnew.pop("comment")
      datadictnew.pop("history")
      # Can't compare fill values since they have changed type
      self.assertEqual(datadictnew["_FillValue"], 1e+20)
      datadictnew.pop("_FillValue")
      datadictold.pop("_FillValue")
      dataattrsok = dict(datadictnew) == dict(datadictold)

      dataok = datafillok and datashapeok and datatypeok and datadtypeok \
 and dataattrshistok and datacomok and dataattrsok
      self.assertTrue(dataok) #!+

      
      # Test global attributes
      globnew = ncnew.__dict__.copy()
      globold = ncminold.__dict__.copy()

      histok = str(globnew["history"]).endswith("UTC created mean temperature\
 from tasminmean_bc and tasmaxmean_bc.")
      authorok = globnew["author"] == "Roman Olson"
      emailok  = globnew["contact"] == "privatemail@email.privatemail"
      comok = globnew["comments"] == "tasmean_bc is not an original NARCliM\
 variable. It is calculated as an average of tasminmean_bc and tasmaxmean_bc\
 variables"
 
      globnew.pop("comments")
      globold.pop("comments")
      globnew.pop("history")
      globold.pop("history")
      globnew.pop("author")
      globold.pop("author")
      globnew.pop("contact")
      globold.pop("contact")
      globnew.pop("creation_date")
      globold.pop("creation_date")

      # New file has no "Conventions" or "date"
      globold.pop("Conventions")
      globold.pop("date")

      globattrsok = dict(globnew) == dict(globold)
      globok = histok and authorok and emailok and comok and globattrsok

      self.assertTrue(globok) #!+





suite = unittest.TestLoader().loadTestsFromTestCase(test_make_tasmean_bc)
unittest.TextTestRunner(verbosity=2).run(suite)
