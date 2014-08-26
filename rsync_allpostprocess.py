#!/usr/bin/env python

""" rsync_allpostprocess.py

Author: Daniel Argues @ CCRC, UNSW. Sydney (Australia)
email: d.argueso@ unsw.edu.au
Created: Fri May 16 10:45:12 EST 2014

"""

import subprocess as subprocess
import ccrc_utils as cu
import os


GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5','CSIRO-MK3.0']
RCM_names=['R1','R2','R3']
Period_names=['1990-2010','2020-2040','2060-2080']
Domain_names=['d01','d02']
Period_covers=['1990-2009','2020-2039','2060-2079']


for gind,gname in enumerate(GCM_names):
  for rind,rname in enumerate(RCM_names):
    for pind,pname in enumerate(Period_names):
      for dind,dname in enumerate(Domain_names):
        fullpath_out=cu.get_postproc_location(gname,rname,pname)[0]
        if not os.path.exists(fullpath_out):
          os.makedirs(fullpath_out)
        fullpath_in="/srv/ccrc/data13/z3393020/NARCliM_newpost/postprocess/%s/%s/%s/%s" %(Period_covers[pind],gname,rname,dname)
        print "rsync -avz --stats %s/* %s%s/" %(fullpath_in,fullpath_out,dname)
        subprocess.call("rsync -avz --stats %s/* %s%s/" %(fullpath_in,fullpath_out,dname),shell=True)