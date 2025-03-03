# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 17:08:04 2024

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
    trj,ca = cxt.fc.mean_traj_nF()
    #cxt.mean_traj_heat(xoffset=xoffset,set_cmx=False)
    xoffset = xoffset+30
    
    if i==0:
        T = np.expand_dims(trj,2)
        C = np.expand_dims(ca,1)
    else:
        T = np.append(T,np.expand_dims(trj,2),axis=2)
        C = np.append(C,np.expand_dims(ca,1),axis=1)
plt.xlim([-10, 110])
savename = os.path.join(savedir,'FB4X_Ca_traj.pdf')
#plt.savefig(savename)
C_x = C
#%% FB4P_b panel
plt.close('all')
savedir = "Y:\\Applications\\HHW\\Figures\\Figure3"
datadirs = [
    "Y:\\Data\\FCI\\Hedwig\\FB4P_b_SS60296\\240912\\f2\\Trial3",
    "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240809\\f2\\Trial2",
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240806\\f1\\Trial5",#New mask
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240805\\f1\\Trial3",#New mask
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240909\\f1\\Trial5",#New mask
       "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240906\\f2\\Trial2",#New mask
      ]




for i,datadir in enumerate(datadirs):
    plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    trj,ca = cxt.fc.mean_traj_nF()
    if i==0:
        T = np.expand_dims(trj,2)
        C = np.expand_dims(ca,1)
    else:
        T = np.append(T,np.expand_dims(trj,2),axis=2)
        C = np.append(C,np.expand_dims(ca,1),axis=1)
    #cxt.mean_traj_heat(xoffset=xoffset,set_cmx=False)
    cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=xoffset)
    xoffset = xoffset+30
    savename = os.path.join(savedir,'FB4P_b_Ca_traj' +str(i)+'.pdf')
    plt.savefig(savename)
#plt.xlim([-10, 110])
#savename = os.path.join(savedir,'FB4P_b_Ca_traj.pdf')

C_PB = C
#%% FB4R
plt.close('all')
datadirs = ["Y:\Data\FCI\\Hedwig\\SS61646_FB4R\\240828\\f3\\Trial1",
            "Y:\Data\FCI\\Hedwig\\SS61646_FB4R\\240910\\f1\\Trial1",
            "Y:\Data\FCI\\Hedwig\\SS61645_FB4R\\240911\\f1\\Trial3" # Only 2 jumps
            ]


for i,datadir in enumerate(datadirs):
    plt.figure()
    xoffset= 0
    
    cxt = CX_tan(datadir)
    trj,ca = cxt.fc.mean_traj_nF()
    if i==0:
        T = np.expand_dims(trj,2)
        C = np.expand_dims(ca,1)
    else:
        T = np.append(T,np.expand_dims(trj,2),axis=2)
        C = np.append(C,np.expand_dims(ca,1),axis=1)
    #cxt.mean_traj_heat(xoffset=xoffset,set_cmx=False)
    cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=xoffset)
    savename = os.path.join(savedir,'FB4R_traj' +str(i) + '.pdf')
    plt.savefig(savename)
    #xoffset = xoffset+30
C_R = C
#%% FB5I
datadirs = [
   "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240628\\f1\\Trial2",#Nice
            "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240917\\f1\\Trial2",#Best for this fly
           "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240917\\f3\\Trial3",
           ]

for i,datadir in enumerate(datadirs):
    plt.figure()
    xoffset= 0
    
    cxt = CX_tan(datadir)
    trj,ca = cxt.fc.mean_traj_nF()
    if i==0:
        T = np.expand_dims(trj,2)
        C = np.expand_dims(ca,1)
    else:
        T = np.append(T,np.expand_dims(trj,2),axis=2)
        C = np.append(C,np.expand_dims(ca,1),axis=1)
    #cxt.mean_traj_heat(xoffset=xoffset,set_cmx=False)
    cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=xoffset)
    savename = os.path.join(savedir,'FB5I_traj' +str(i) + '.pdf')
    plt.savefig(savename)
    #xoffset = xoffset+30
C_I = C
#%% 

plt.close('all')
cnx = np.max(C_x,axis=0)
cm_x = np.mean(C_x/cnx,axis=1)
cnpb = np.max(C_PB,axis=0)
cm_pb = np.mean(C_PB/cnpb,axis=1)

cnr = np.max(C_R,axis=0)
cm_r = np.mean(C_R/cnr,axis=1)

cni = np.max(C_I,axis=0)
cm_i = np.mean(C_I/cni,axis=1)

colours = np.array([[0,127,10],[0,43,127],[127,0,43],[200,50,200]])/256
plt.fill([49,99,99,49],[-1,-1,1,1],color=[0.7,0.7,0.7])
plt.plot([0,99],[0,0],color='k',linestyle='--')
plt.plot([49,49],[-1,1],color='k',linestyle='--')
plt.plot([0,0],[-1,1],color='k',linestyle='--')
plt.plot([99,99],[-1,1],color='k',linestyle='--')
plt.plot(cm_x,color=colours[0,:])
plt.plot(cm_pb,color=colours[1,:])
plt.plot(cm_r,color=colours[2,:])
plt.plot(cm_i,color=colours[3,:])
plt.xticks([49,99],labels=['Plume entry','Plume exit'])
plt.xlim([0,101])
plt.ylabel('Normalised fluor')
savename = os.path.join(savedir,'Mean_Fluor_ent_ext.pdf')
plt.savefig(savename)

#%%
datadirs = [
    #"Y:\\Data\\FCI\\Hedwig\\FB4P_b_SS67631\\240912\\f2\\Trial3",
    #"Y:\Data\FCI\\Hedwig\\FB4P_b_SS67631_sytGC7f\\240806\\f1\\Trial5",#New mask
    #"Y:\Data\FCI\\Hedwig\\FB4P_b_SS67631_sytGC7f\\240805\\f1\\Trial3",#New mask
    #"Y:\Data\FCI\\Hedwig\\FB4P_b_SS67631_sytGC7f\\240909\\f1\\Trial5",#New mask
    #"Y:\Data\FCI\\Hedwig\\FB4P_b_SS67631_sytGC7f\\240906\\f2\\Trial2",#New mask
    "Y:\Data\FCI\\Hedwig\\FB4P_b_SS67631_sytGC7f\\240809\\f2\\Trial2",
      ]
for i,datadir in enumerate(datadirs):
    
    d = datadir.split("\\")
    name = d[-3] + '_' + d[-2] + '_' + d[-1]
    cx = CX(name,['fsbTN'],datadir)
    # save preprocessing, consolidates behavioural data
    cx.save_preprocessing()
    # Process ROIs and saves csv
    cx.process_rois()
    # Post processing, saves data as h5
    cx.crop = False
    cx.save_postprocessing()










