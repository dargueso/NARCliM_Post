#!/usr/bin/env python
"""runpostprocess_allsim.py script
   Script to run postprocess for all NARCliM simulations
   It reads an input file (NARClIM_post_template.input), where the input arguments are provided
   Authors: Daniel Argueso (d.argueso@unsw.edu.au), Alejandro Di Luca (a.diluca@unsw.edu.au)
   Institution: CoECSS, UNSW. CCRC, UNSW.
   Created: 09 May 2014
"""
import subprocess
import numpy as np
GCM_names=['MIROC3.2','CCCMA3.1','ECHAM5','CSIRO-MK30']
RCM_names=['R1','R2','R3']
Period_names=['1990-2010','2020-2040','2060-2080']
Domain_names=['d01','d02']
indeck="NARClIM_post.input.deck"



for gind,gcm in enumerate(GCM_names):
  for rind,rcm in enumerate(RCM_names):
    for pind,period in enumerate(Period_names):
      for dind,domain in enumerate(Domain_names):
        pathin = "/home/z3393020//WRFouts/NARCliM/%s/%s/%s/out/" %(gcm,rcm,period)
        pathout="/srv/ccrc/data13/z3393020/NARCliM_newpost/"
        fin = open (indeck,"r")
        fout = open ("NARCliM_repost_%s_%s_%s.input"%(gcm,rcm,period),"w")
        syear=np.int(period[0:4])
        eyear=np.int(period[5:])-1
        print "Processing GCM: %s; RCM: %s; Period: %s; Domain: %s" %(gcm,rcm,period,domain)
        
        namelist_dic={'%pathin%'  : pathin,
                      '%pathout%' : pathout,
                      '%GCM%'     : gcm,
                      '%RCM%'     : rcm,
                      '%syear'    : str(syear),
                      '%eyear%'   : str(eyear),
                      '%domain%'  : domain,
                      }
        
        for line in fin.readlines():
          for linerep in namelist_dic.keys():
            line=line.replace(linerep,namelist_dic[linerep])

          fout.write(lines)




        fin.close()
        fout.close()
        
        subprocess.call("python ./postprocess_NARCliM.py -i NARCliM_repost_%s_%s_%s.input" %(gcm,rcm,period), shell=True)