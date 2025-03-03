# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 09:44:58 2024

@author: dowel
"""
from analysis_funs.regression import fci_regmodel

import numpy as np
import pandas as pd
import src.utilities.funcs as func
from analysis_funs.optogenetics import opto 
import os
import matplotlib.pyplot as plt
from analysis_funs.CX_imaging import CX
from analysis_funs.CX_analysis_tan import CX_tan
plt.rcParams['pdf.fonttype'] = 42 
savedir = "Y:\\Applications\\HHW\\Figures\\Figure3"
#%% FB4X panel
# Get mean trajectory and mean activity

datadirs = ["Y:\Data\FCI\Hedwig\\SS70711_FB4X\\240307\\f1\\Trial3",
            "Y:\Data\FCI\Hedwig\\SS70711_FB4X\\240313\\f1\\Trial3",
            "Y:\Data\FCI\Hedwig\\SS70711_FB4X\\240531\\f1\\Trial3"]
xoffset= 0
plt.figure()
for i,datadir in enumerate(datadirs):
    cxt = CX_tan(datadir)
    trj,ca = cxt.mean_traj_nF()
    cxt.mean_traj_heat(xoffset=xoffset,set_cmx=False)
    xoffset = xoffset+30
    
    if i==0:
        T = np.expand_dims(trj,2)
        C = np.expand_dims(ca,1)
    else:
        T = np.append(T,np.expand_dims(trj,2),axis=2)
        C = np.append(C,np.expand_dims(ca,1),axis=1)
plt.xlim([-10, 110])
savename = os.path.join(savedir,'FB4X_Ca_traj.pdf')
plt.savefig(savename)
C_x = C
#%% FB4R panel
datadirs = ['Y:\\Data\FCI\\AndyData\\47H09_GDrive\\47H09\\20220128_47H09jGCaMP7f_Fly1-002',
         'Y:\\Data\FCI\\AndyData\\47H09_GDrive\\47H09\\20220130_47H09sytjGCaMP7f_Fly1-001',
         'Y:\\Data\FCI\\AndyData\\47H09_GDrive\\47H09\\20220128_47H09jGCaMP7f_Fly2-001',
         'Y:\\Data\FCI\\AndyData\\47H09_GDrive\\47H09\\20220127_47H09jGCaMP7f_Fly1-001']
plt.figure()
xoffset= 0

for i,datadir in enumerate(datadirs):
    cxt = CX_tan(datadir,tnstring='fb4r')
    trj,ca = cxt.mean_traj_nF()
    if i==0:
        T = np.expand_dims(trj,2)
        C = np.expand_dims(ca,1)
    else:
        T = np.append(T,np.expand_dims(trj,2),axis=2)
        C = np.append(C,np.expand_dims(ca,1),axis=1)
    cxt.mean_traj_heat(xoffset=xoffset,set_cmx=False)
    xoffset = xoffset+30
    
plt.xlim([-10, 110])
savename = os.path.join(savedir,'FB4R_Ca_traj.pdf')
plt.savefig(savename)
C_R = C
#%% Mean fluor
plt.close('all')
cnx = np.max(C_x,axis=0)
cm_x = np.mean(C_x/cnx,axis=1)
cnr = np.max(C_R,axis=0)
cm_r = np.mean(C_R/cnr,axis=1)
colours = np.array([[0,127,10],[0,43,127]])/256
plt.fill([49,99,99,49],[-1,-1,1,1],color=[0.7,0.7,0.7])
plt.plot([0,99],[0,0],color='k',linestyle='--')
plt.plot([49,49],[-1,1],color='k',linestyle='--')
plt.plot([0,0],[-1,1],color='k',linestyle='--')
plt.plot([99,99],[-1,1],color='k',linestyle='--')
plt.plot(cm_x,color=colours[0,:])
plt.plot(cm_r,color=colours[1,:])
plt.xticks([49,99],labels=['Plume entry','Plume exit'])
plt.xlim([0,101])
plt.ylabel('Normalised fluor')
savename = os.path.join(savedir,'Mean_Fluor_ent_ext.pdf')
plt.savefig(savename)