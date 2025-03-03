# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 14:24:05 2024

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
savedir = "Y:\\Data\\FCI\\FCI_summaries\\TangentialSummaries"
from Utils.utils_general import utils_general as ugn
#%% FB4P_b
plt.close('all')
datadirs = [
    "Y:\\Data\\FCI\\Hedwig\\FB4P_b_SS60296\\240912\\f2\\Trial3",#27 jumps
    "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240809\\f2\\Trial2",#24 jumps
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240806\\f1\\Trial5",#2 jumps
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240805\\f1\\Trial3",# 2 jumps
    # "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240909\\f1\\Trial5",# 1 jump
       "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240906\\f2\\Trial2",# 2 jumps
      ]

bins = 100
timecourse = 15
traj_fb4pb = np.empty((bins*2,2,len(datadirs)))
ca_fb4pb = np.empty((bins*2,len(datadirs)))
cart_fb4pb = np.empty((timecourse*10,len(datadirs)))
# for i,datadir in enumerate(datadirs):
#     cxt = CX_tan(datadir)
#     j = cxt.get_jumps()
#     print(len(j))

for i,datadir in enumerate(datadirs):
    
    print(datadir)
    #plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    #trj,ca = cxt.fc.mean_traj_nF()
    
    ca,trj = cxt.mean_traj_jump(bins=bins)
    
    
    traj_fb4pb[:,:,i] = np.mean(trj,axis=2)
    ca_fb4pb[:,i] = np.mean(ca,axis=1)
    ct = cxt.mean_ca_realtime_jump(timecourse,norm=False)
    
    x = np.nanmean(ct,axis=1)
    cart_fb4pb[:,i] = x
    #cart_fb4pb[:,i] = np.nanmean(ct,axis=1)
    
   # plt.figure()
    #cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=0,set_cmx =True,cmx=0.25)
    #savename = os.path.join(savedir,'FB4P_b_Ca_traj' +str(i)+'.pdf')
   # plt.savefig(savename)

savedict = {'ca_fb4pb': ca_fb4pb,'traj_fb4pb':traj_fb4pb,'ca_fb4pb_realtime':cart_fb4pb}

# FB4R
plt.close('all')
datadirs = ["Y:\Data\FCI\\Hedwig\\SS61646_FB4R\\240828\\f3\\Trial1",
            "Y:\Data\FCI\\Hedwig\\SS61646_FB4R\\240910\\f1\\Trial1",
            "Y:\Data\FCI\\Hedwig\\SS61645_FB4R\\240911\\f1\\Trial3" # Only 2 jumps
            ]


traj_fb4r = np.empty((bins*2,2,len(datadirs)))
ca_fb4r = np.empty((bins*2,len(datadirs)))
cart_fb4r = np.empty((timecourse*10,len(datadirs)))
for i,datadir in enumerate(datadirs):
    #plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    
    ca,trj = cxt.mean_traj_jump(bins=bins)
    traj_fb4r[:,:,i] = np.mean(trj,axis=2)
    ca_fb4r[:,i] = np.mean(ca,axis=1)
    
    ct = cxt.mean_ca_realtime_jump(timecourse)
    
    #ct[np.abs(ct)==0.0] = np.nan
    x = np.nanmean(ct,axis=1)
    cart_fb4r[:,i] = x
    
    # plt.figure()
    # cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=0,set_cmx =True,cmx=0.25)
    # savename = os.path.join(savedir,'FB4R_Ca_traj' +str(i)+'.pdf')
    # plt.savefig(savename)

savedict.update({'ca_fb4r': ca_fb4r,'traj_fb4r':traj_fb4r,'ca_fb4r_realtime':cart_fb4r})


# FB5I
datadirs = [
   "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240628\\f1\\Trial2",#Nice
            "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240917\\f1\\Trial2",#Best for this fly
           "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240917\\f3\\Trial3",
           ]

traj_fb5i = np.empty((bins*2,2,len(datadirs)))
ca_fb5i = np.empty((bins*2,len(datadirs)))
cart_fb5I = np.empty((timecourse*10,len(datadirs)))

for i,datadir in enumerate(datadirs):
   # plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    
    ca,trj = cxt.mean_traj_jump(bins=bins)
    #cxt.fc.jump_heat(np.mean(trj,axis=2),np.mean(ca,axis=1),xoffset=0)
    
    traj_fb5i[:,:,i] = np.mean(trj,axis=2)
    ca_fb5i[:,i] = np.mean(ca,axis=1)
    
    ct = cxt.mean_ca_realtime_jump(timecourse)
    #ct[np.abs(ct)==0.0] = np.nan
    x = np.nanmean(ct,axis=1)
    cart_fb5I[:,i] = x
    
    # plt.figure()
    # cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=0,set_cmx =True,cmx=0.25)
    # savename = os.path.join(savedir,'FB5I_Ca_traj' +str(i)+'.pdf')
    # plt.savefig(savename)

savedict.update({'ca_fb5i': ca_fb5i,'traj_fb5i':traj_fb5i,'ca_fb5i_realtime':cart_fb5I})


# FB4X
datadirs = [
   "Y:\Data\\FCI\\Hedwig\\SS70711_FB4X\\241030\\f3\\Trial3",
   "Y:\\Data\\FCI\\Hedwig\\SS70711_FB4X\\241031\\f1\\Trial3"
           ]

traj_fb4x = np.empty((bins*2,2,len(datadirs)))
ca_fb4x = np.empty((bins*2,len(datadirs)))
cart_fb4x = np.empty((timecourse*10,len(datadirs)))

for i,datadir in enumerate(datadirs):
    #plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    
    ca,trj = cxt.mean_traj_jump(bins=bins)
    #cxt.fc.jump_heat(np.mean(trj,axis=2),np.mean(ca,axis=1),xoffset=0)
    
    traj_fb4x[:,:,i] = np.mean(trj,axis=2)
    ca_fb4x[:,i] = np.mean(ca,axis=1)
    
    ct = cxt.mean_ca_realtime_jump(timecourse)
    #ct[np.abs(ct)==0.0] = np.nan
    x = np.nanmean(ct,axis=1)
    cart_fb4x[:,i] = x
    
    # plt.figure()
    # cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=0,set_cmx =True,cmx=0.25)
    # savename = os.path.join(savedir,'FB4X_Ca_traj' +str(i)+'.pdf')
    # plt.savefig(savename)

savedict.update({'ca_fb4x': ca_fb4x,'traj_fb4x':traj_fb4x,'ca_fb4x_realtime': cart_fb4x})

# Save data as pickle

u = ugn
u.save_pick(savedict,os.path.join(savedir,'TN_mean_jump.pickle'))
data = u.load_pick(os.path.join(savedir,'TN_mean_jump.pickle'))

#%% plot max
u = ugn
plt.close('all')
colours2 = np.array([[106,207,246],[237,30,36],[168,170,173],[6,149,207]])/255

data = u.load_pick(os.path.join(savedir,'TN_mean_jump.pickle'))
regions =['fb5i','fb4pb','fb4r']
#regions = ['fb4pb']
rcount = 0
for i,r in enumerate(regions):
    traj = data['traj_' +r]
    if i==0:
        tmean = np.sum(traj,axis=2)
    else:
        tmean = tmean+np.sum(traj,axis=2)
    rcount = rcount+traj.shape[2]
    
tmean = tmean/rcount
ymn = np.min(tmean[:,1])
ymx= np.max(tmean[:,1])
x = np.array([-5,5,5,-5])-5
y = [ymn,ymn,0,0]
y2 = [0,0,ymx,ymx]
plt.fill(x,y,color=[0.8,0.8,0.8])
plt.fill(x-3,y2,color=[0.8,0.8,0.8])

plt.plot(tmean[:bins,0],tmean[:bins,1],color=colours2[1,:])
plt.plot(tmean[bins-1:,0],tmean[bins-1:,1],color=colours2[0,:])
colours = np.array([[11,15,157],
[0,175,171],
[0,131,201],
[0,68,27]])/255

for i,r in enumerate(regions):
    ca = data['ca_' + r]
    im = np.argmax(ca,axis=0)
    im_mean = np.round(np.mean(im)).astype('int')
    #im_mean = np.argmax(np.mean(ca,axis=1))
    #plt.scatter(tmean[im,0],tmean[im,1],s=10,color=colours[i,:],zorder=10,alpha=1)
    plt.scatter(tmean[im_mean,0],tmean[im_mean,1],s=100,color=colours[i,:],zorder=10,alpha=1)
        
ax = plt.gca()
    
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show()
plt.savefig(os.path.join(savedir,'SummaryPeak.pdf'))
#%% plot ca
plt.close('all')
cells = ['fb4pb','fb4r','fb5i']
colours = np.array([
    [155,116,255],
    [63,184,68],
    [0,104,56],[0,0,0]])/255
newtime = np.arange(0,200,2)
oldtime = np.arange(0,200)
for i,c in enumerate(cells):
    tdat = data['ca_'+c]
    tdat = tdat/np.max(tdat,axis=0)
    plt.fill([0,100,100,0],[-1,-1,1,1],color=[0.8,0.8,0.8])
    plt.plot([100,100],[-1,1],color='k',linestyle='--')
   #plt.plot(tdat,color=colours[i,:],alpha=0.2)
    plt.plot([200,200],[-1,1],color='k',linestyle='--')
    y = np.mean(tdat,axis=1)
    ynew = np.interp(newtime,oldtime,y)
    plt.plot(y,color=colours[i,:])
    
plt.ylim([0,1])
plt.xlim([0,200])
plt.xlabel('Normalised time')
plt.ylabel('Norm dF/F')
plt.savefig(os.path.join(savedir,'SummaryTNpseudotime.pdf'))
#%% FC2 comparison
from EdgeTrackingOriginal.ETpap_plots.ET_paper import ET_paper
from scipy.stats import circmean, circstd

FC2_dirs = [
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f1\\Trial3",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f2\\Trial3",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240502\\f1\\Trial2",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240514\\f1\\Trial2",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\241104\\f1\\Trial5",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\241106\\f1\\Trial2"]
timecourse = 15
t = np.arange(0,timecourse,0.1)
for i, d in enumerate(FC2_dirs):
    etp = ET_paper(d)
    p = etp.phase_time(timecourse)
    p[np.abs(p)==0.0] = np.nan
    pmean = circmean(p,high=np.pi, low=-np.pi,axis=2,nan_policy='omit')
    if i==0:
        pltmean = np.zeros((len(pmean),3,len(FC2_dirs)))
    pltmean[:,:,i] = pmean
zorder = [1,3,0]
pm = circmean(pltmean,high=np.pi, low=-np.pi,axis=2)

#%%
x = np.arange(0,15,0.1)
fig, ax1 = plt.subplots()

colours = np.array([[81,156,204],[84,39,143],[0,0,0],[6,149,207]])/255
colours = np.array([[0,0,0],[81,156,204]])/255
pm = circmean(pltmean,high=180,low=-180,axis=2)
zorder = [1,3,0]
for i in range(2):
   # plt.plot(pltmean[:,i,:],t,color = colours[i,:],alpha=0.3,zorder=zorder[i])
    ax1.plot(x,pm[:,i],color=colours[i,:],zorder=zorder[i])

ax1.plot([0,15],[0,0],color='k',linestyle='--',zorder=0)

#ax1.plot(pm[:,1],x)
ax1.set_ylim([np.pi,-np.pi])
ax2 = ax1.twinx() 
#ax2.plot(-cart_fb4pb,x,color=[0.2,0.6,0.4],alpha=0.3)

ax2.plot(x,np.nanmean(cart_fb4pb,axis=1),color=[0.2,0.6,0.4])
ax2.plot(x,np.nanmean(cart_fb5I,axis=1),color=[0.2,0.6,0.6])
ax2.set_ylim([0,0.6])

#%% All data trajectories
datadirs = [
    "Y:\\Data\\FCI\\Hedwig\\FB4P_b_SS60296\\240912\\f2\\Trial3",#27 jumps
    "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240809\\f2\\Trial2",#24 jumps
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240806\\f1\\Trial5",#2 jumps
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240805\\f1\\Trial3",# 2 jumps
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240909\\f1\\Trial5",# 1 jump
       "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240906\\f2\\Trial2",# 2 jumps
      ]
offset =0
plt.close('all')
savedir = r'Y:\Presentations\2025\02_BSV\Tangentials'
for d in datadirs:
    cxt = CX_tan(d)
    # cxt.fc.example_trajectory_jump(cmin=-0.4,cmax =0.4)
    # plt.title('FB4P_b '+ cxt.name)
    # plt.figure()
    # cxt.fc.mean_traj_nF_jump(cxt.fc.ca,plotjumps=True,cmx=0.25)
    # plt.title('FB4R '+ cxt.name)
    
    # plt.figure()
    cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=offset,set_cmx=True,cmx=0.25)
    offset = offset+20
    plt.savefig(os.path.join(savedir,'FB4P_b_jump_mean.pdf'))
#%%
datadirs = ["Y:\Data\FCI\\Hedwig\\SS61646_FB4R\\240828\\f3\\Trial1",
            "Y:\Data\FCI\\Hedwig\\SS61646_FB4R\\240910\\f1\\Trial1",
            "Y:\Data\FCI\\Hedwig\\SS61645_FB4R\\240911\\f1\\Trial3" # Only 2 jumps
            ]
offset = 0
savedir = r'Y:\Presentations\2025\02_BSV\Tangentials'
plt.close('all')
for d in datadirs:
    cxt = CX_tan(d)
    # cxt.fc.example_trajectory_jump(cmin=-0.4,cmax =0.4)
    # plt.title('FB4R '+ cxt.name)
    # plt.figure()
    # cxt.fc.mean_traj_nF_jump(cxt.fc.ca,plotjumps=True,cmx=0.25)
    # plt.title('FB4R '+ cxt.name)
    
    plt.figure()
    cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=offset,set_cmx=True,cmx=0.25)
    offset = offset+20
    plt.savefig(os.path.join(savedir,'FB4R_jump_mean.pdf'))
#%%    
datadirs = [
   "Y:\Data\\FCI\\Hedwig\\SS70711_FB4X\\241030\\f3\\Trial3",
   "Y:\\Data\\FCI\\Hedwig\\SS70711_FB4X\\241031\\f1\\Trial3"
           ]
for d in datadirs:
    cxt = CX_tan(d)
    cxt.fc.example_trajectory_jump(cmin=-0.4,cmax =0.4)
    plt.title('FB4X '+ cxt.name)
#%%
plt.close('all')
datadirs = [
   "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240628\\f1\\Trial2",#Nice
            "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240917\\f1\\Trial2",#Best for this fly
           "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240917\\f3\\Trial3",
           ]

savedir = r'Y:\Presentations\2025\02_BSV\Tangentials'
for d in datadirs:
    cxt = CX_tan(d)
    cxt.fc.example_trajectory_jump(cmin=-0.4,cmax =0.4)
    plt.title('FB5I '+ cxt.name)
    plt.figure()
    cxt.fc.mean_traj_nF_jump(cxt.fc.ca,plotjumps=True,cmx=0.25)
    plt.title('FB5I '+ cxt.name)
    
    # cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=offset,set_cmx=True,cmx=0.25)
    # offset = offset+20
    # plt.savefig(os.path.join(savedir,'FB5I_jump_mean.pdf'))
    
    
    
    
    