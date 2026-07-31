# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 16:36:16 2026

@author: dowel
"""

from analysis_funs.regression import fci_regmodel

import numpy as np
import pandas as pd
import src.utilities.funcs as fc
from analysis_funs.optogenetics import opto 
import os
import matplotlib.pyplot as plt 
from analysis_funs.utilities import imaging as im
from skimage import io, data, registration, filters, measure
from scipy import signal as sg
from analysis_funs.CX_imaging import CX
from analysis_funs.CX_analysis_col import CX_a
from src.utilities import funcs as fn
from scipy import stats
from Utilities.utils_general import utils_general as ug
from Utilities.utils_plotting import uplt
plt.rcParams['font.sans-serif'] = 'Arial'
from EdgeTrackingOriginal.ETpap_plots.ET_paper import ET_paper

plt.rcParams['pdf.fonttype'] = 42 
import matplotlib as mpl
mpl.rcParams['path.simplify'] = True
plt.rcParams['path.simplify'] = True
#%% Load fly data

groups = {
    "FC2": ["Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f1\\Trial3",
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f2\\Trial3",
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\240502\\f1\\Trial2",
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\240514\\f1\\Trial2",
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\241104\\f1\\Trial5",
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\241106\\f1\\Trial2",
    r"Y:\Data\FCI\Hedwig\FC2_maimon2\250128\f1\Trial4"
    ],
    "hDeltaC": [r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250411\f1\Trial1",# Phase recording is not the best - strong pointer
    #r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f1\Trial3",# Not many jumps, weak pointer
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f2\Trial2", # Strong pointer
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250714\f1\Trial2",# Strong pointer, just like FC2
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f1\Trial1", # Points away ******
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f2\Trial3", # Strong pointer
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250721\f1\Trial2",# Strong pointer
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250907\f1\Trial2",# Good pointer
       
    ],
    "hDeltaJ": [r"Y:\Data\FCI\Hedwig\hDeltaJ\240529\f1\Trial3",#Good pointer
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251011\f1\Trial1",
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251022\f1\Trial2",
    r'Y:\Data\FCI\Hedwig\hDeltaJ\251121\f2\Trial1',
    r'Y:\Data\FCI\Hedwig\hDeltaJ\260115\f1\Trial2',    
    r'Y:\Data\FCI\Hedwig\hDeltaJ\260127\f1\Trial1',
    r'Y:\Data\FCI\Hedwig\hDeltaJ\260203\f1\Trial2'
    ],
    
    "hDeltaH":[ 
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f1\Trial2",# Varied recovery with ET and plume distance
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f2\Trial3",# Long returns and plume traversals, fast recovery dynamics
        r'Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251209\f1\Trial1',# Ok behaviour
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f1\Trial2",# Decent tracking, some traversing also
               r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f2\Trial1"], #Nice dataset. Some ET bouts followed by plume traversal]
    
    
    "hDeltaC_FC2": [r'Y:\Data\FCI\Hedwig\hDeltaC_68A10_FC2_GCaMP_RCaMP\260131\f1\Trial2',
        r'Y:\Data\FCI\Hedwig\hDeltaC_68A10_FC2_GCaMP_RCaMP\260201\f1\Trial3',
        r'Y:\Data\FCI\Hedwig\hDeltaC_68A10_FC2_GCaMP_RCaMP\260211\f1\Trial2'],
    
    
    "hDeltaC_ACV_pulses": [
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250320\f2\Trial5",#ACV pulses
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250328\f1\Trial5",#'ACV
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f1\Trial6", # ACV pulses
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250714\f1\Trial6", # ACV pulses
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250716\f1\Trial6", # ACV pulses
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f1\Trial5",#ACV pulses
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250721\f1\Trial6", # ACV pulses
        
    ],
    "hDeltaC_Oct_pulses": [
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250320\f2\Trial4",#Octanol [?] pulses - neuron is inhibited
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250328\f1\Trial4",#'Octanol'
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f1\Trial7", # octanol pulses
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250714\f1\Trial7", # Octanol pulses
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250716\f1\Trial5", #  oct pulses, imaging artefact in pre air period
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f1\Trial7",#OCT pulses
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250721\f1\Trial5", # Oct pulses
       
    ],
    "hDeltaJ_ACV_pulses": [ 
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251011\f1\Trial4", # ACV pulses
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f1\Trial6", # ACV pulses
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f2\Trial5", # ACV pulses
    r'Y:\Data\FCI\Hedwig\hDeltaJ\251121\f2\Trial3', # ACV pulses
     r'Y:\Data\FCI\Hedwig\hDeltaJ\260115\f1\Trial3'
        
    ],
    "hDeltaJ_Oct_pulses": [
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251011\f1\Trial5", # Oct pulses
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f1\Trial5", # Oct pulses
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f2\Trial4", # Oct pulses
    r'Y:\Data\FCI\Hedwig\hDeltaJ\251121\f2\Trial4',# Oct pulses
    r'Y:\Data\FCI\Hedwig\hDeltaJ\260115\f1\Trial4'   
    ],
    
    "hDeltaH_ACV_pulses": [r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f1\Trial4",
                           r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f2\Trial5",
                           r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f1\Trial5",
                           r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f2\Trial4"],
    "hDeltaH_Oct_pulses":[
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f1\Trial5",
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f2\Trial4",
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f1\Trial4",
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f2\Trial5"]
}
all_flies = {}

for group_name, datadirs in groups.items():
    all_flies[group_name] = {}
    
    for i, datadir in enumerate(datadirs):
        print(f"[{group_name}] {datadir}")
        if group_name=='hDeltaC_FC2':
            etp = ET_paper(datadir,regions=['fsb1_ch1','fsb2_ch2'])
        else:
            etp = ET_paper(datadir)
        all_flies[group_name][i] = etp
        
#%% Plot example arrows for hdeltaJ
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC']
mastercolours = np.array([[0,0,0],[81,61,204],[49,99,125],[81,156,205]])/255
c_combos = np.array([[0,1],[0,2]])
savedir = r'Y:\Posters\IWNC\panels'
savedir = r'D:\Posters\IWNC'
for di,d in enumerate(dchoices):
    data = all_flies[d]
    for i in data:
        etp = data[i]
        plt.figure()
        colours = mastercolours[c_combos[di,:]]
        etp.plt_tmp(scalar=8,phase_num=6,colours=colours)
        plt.xlim([-100,100])
        plt.ylim([-100,200])
        plt.savefig(os.path.join(savedir,d+'_'+str(i)+ '_side_'+ str(etp.cxa.side) +'.pdf'))
        plt.savefig(os.path.join(savedir,d+'_'+str(i)+ '_side_'+ str(etp.cxa.side) +'.png'))
    
#%% Plot example mean arrows for hDeltaC
#%% Pseudo time plots

achoice=1
dchoices = ['hDeltaJ','hDeltaC','FC2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
numpoints = 12  
plt.close('all')
grandphase = np.zeros((200,2,len(dchoices)))
grandtraj = np.zeros((200,2))
grandamp = np.zeros((200,2,2,len(dchoices)))
for i,d1 in enumerate(dchoices):
    data=all_flies[d1]
        
    phases = np.zeros((200,2,len(data)))
    trajs = np.zeros((200,2,len(data)))
    amps = np.zeros((200,2,2,len(data)))
    for d in data:
        print(d1)
        etp = data[d]
        phase,traj,amp = etp.trajectory_mean(time_threshold2=0,regions=['eb','fsb_upper'])
        print(traj.shape)
        phases[:,:,int(d)] = stats.circmean(phase,axis=2,low=-np.pi,high=np.pi)
        trajs[:,:,int(d)] = np.mean(traj,axis=2)
        amps[:,:,:,int(d)] = np.mean(amp,axis=3)
    pmean = stats.circmean(phases,axis=2,low=-np.pi,high=np.pi)
    tmean = np.mean(trajs,axis=2)
    grandphase[:,:,i] = pmean
    grandtraj = grandtraj+tmean*(1/len(dchoices))
    grandamp[:,:,:,i] =  np.mean(amps,axis=3)    

pdx = np.linspace(0,199,numpoints,dtype= 'int')
offset = -1
for p in pdx:
    offset +=1
    for i in range(len(dchoices)):
        for ip in range(2):
            tphase = grandphase[p,ip,i]
            x = np.array([np.sin(tphase)*grandamp[p,achoice,ip,i],0])+offset
            y = np.array([np.cos(tphase)*grandamp[p,achoice,ip,i],0])
            if ip>0:
                
                plt.plot(x,y,color=colours[i,:])
            else:
                plt.plot(x,y,color=[0.5,0.5,0.5])
            
                
                
g = plt.gca()
g.set_aspect('equal')
plt.fill([-.5,5.5,5.5,-.5],[-.5,-.5,.5,.5],color=[0.8,0.8,0.8],zorder=-1)

#%%
savedir = r'Y:\Presentations\2026\03_CSHL\Panels'
plt.close('all')
x = np.linspace(0,1,100)
x  = np.append(x,np.linspace(x[-1]+x[1],5+x[1],100))
plt.fill([0,1,1,0],[0,0,.8,.8],color=[0.8,.8,.8])
for i,d1 in enumerate(dchoices):
    data=all_flies[d1]
        
    phases = np.zeros((200,2,len(data)))
    trajs = np.zeros((200,2,len(data)))
    amps = np.zeros((200,2,2,len(data)))
    for d in data:
        print(d1)
        etp = data[d]
        phase,traj,amp = etp.trajectory_mean(time_threshold2=0,regions=['eb','fsb_upper'])
        print(traj.shape)
        phases[:,:,int(d)] = stats.circmean(phase,axis=2,low=-np.pi,high=np.pi)
        trajs[:,:,int(d)] = np.mean(traj,axis=2)
        amps[:,:,:,int(d)] = np.mean(amp,axis=3)
        
    plt.plot(x,amps[:,achoice,1,:],color=colours[i,:],alpha=0.4)

for i in range(3):
    plt.plot(x,grandamp[:,achoice,1,i],color=colours[i,:],linewidth=2)


#timepoints = [.75,1,1.5,1.75,2,2.5,3,3.5,4]
timepoints = [.75,1.75,4.5]

for r in timepoints:
# In plume timepoint
    plt.plot([r,r],[0,.8],color='k')

plt.xlim([0,5])
plt.yticks([0,.2,.4,.6,.8])
plt.ylim([0,.65])
plt.ylabel('mean dF/F')
plt.xlabel('pseudotime')
plt.savefig(os.path.join(savedir,'Recovery_pseudotime.pdf'))

plt.figure()
for r in timepoints:
    dx = ug.find_nearest(r,x)
    tamps = grandamp[dx,achoice,1,:]
    tphase = grandphase[dx,1,:]
    for i in range(3):
        tx = np.array([0,np.sin(tphase[i])])*(tamps[i]-0.15)
        ty = np.array([0,np.cos(tphase[i])])*(tamps[i]-0.15)
        plt.plot(tx+r,ty,color=colours[i,:])
plt.xlim([0,5])
plt.gca().set_aspect('equal')
plt.savefig(os.path.join(savedir,'Phase_3_timepoints.pdf'))



plt.figure()
ymn = np.min(grandtraj[:,1])
plt.fill([-5,0,0,-5],[ymn,ymn,0,0],color=[0.8,0.8,0.8])
ymn = np.max(grandtraj[:,1])
plt.fill([-8,-3,-3,-8],[ymn,ymn,0,0],color=[0.8,0.8,0.8])
plt.plot(grandtraj[:,0],grandtraj[:,1],color='k')
plt.gca().set_aspect('equal')
scale = 20
for r in timepoints:
    dx = ug.find_nearest(r,x)
    tamps = grandamp[dx,achoice,1,:]
    tphase = grandphase[dx,1,:]
    for i in range(3):
        tx = np.array([0,np.sin(tphase[i])])*(tamps[i]-.15)*scale
        ty = np.array([0,np.cos(tphase[i])])*(tamps[i]-.15)*scale
        plt.plot(tx+grandtraj[dx,0],ty+grandtraj[dx,1],color=colours[i,:])
        
plt.savefig(os.path.join(savedir,'Phase_3_on_traj.pdf'))
#%%
savedir = r'Y:\Posters\IWNC\panels'
plt.close('all')
plt.figure(figsize=(10,3.5))
x = np.linspace(0,1,100)
x  = np.append(x,np.linspace(x[-1]+x[1],5+x[1],100))
plt.fill([0,1,1,0],[-np.pi,-np.pi,np.pi,np.pi],color=[0.8,.8,.8])
for i,d1 in enumerate(dchoices):
    data=all_flies[d1]
        
    phases = np.zeros((200,2,len(data)))
    trajs = np.zeros((200,2,len(data)))
    amps = np.zeros((200,2,2,len(data)))
    for di,d in enumerate(data):
        print(d1)
        etp = data[d]
        phase,traj,amp = etp.trajectory_mean(time_threshold2=0,regions=['eb','fsb_upper'])
        print(traj.shape)
        phases[:,:,int(d)] = stats.circmean(phase,axis=2,low=-np.pi,high=np.pi)
        trajs[:,:,int(d)] = np.mean(traj,axis=2)
        amps[:,:,:,int(d)] = np.mean(amp,axis=3)
        
        #plt.scatter(x,phases[:,1,di],color=colours[i,:],alpha=0.3,s=10,linewidth=0)

for i in range(3):
    #plt.plot(x,grandphase[:,1,i],color=colours[i,:])
    tamp = grandamp[:,1,1,i]
    tamp = (tamp-.2)*2
    tamp[tamp<0] = 0
    plt.scatter(x,grandphase[:,1,i],color=colours[i,:],linewidth=0,alpha = tamp,s=tamp*100)
    
    
plt.plot([1,5],[0,0],color='k',linestyle='--')
plt.plot([1,5],[-np.pi/2,-np.pi/2],color='r',linestyle='--')

plt.xlim([0,5])
plt.ylim([-np.pi,np.pi])
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])

timepoints = [.75,1.75,4.5]

for r in timepoints:
# In plume timepoint
    plt.plot([r,r],[-np.pi,np.pi],color='k')

plt.xlabel('pseudotime')
plt.savefig(os.path.join(savedir,'PhaseProgression.pdf'))






