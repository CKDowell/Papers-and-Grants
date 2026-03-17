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
from Utilities.utils_general import utils_general as ugn
data = ugn.load_pick(os.path.join(savedir,'TN_mean_jump.pickle'))

#%% FB4P_b
plt.close('all')
datadirs = [
    "Y:\\Data\\FCI\\Hedwig\\FB4P_b_SS60296\\240912\\f2\\Trial3",#27 jumps
    "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240809\\f2\\Trial2",#24 jumps
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240806\\f1\\Trial5",#2 jumps
     "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240805\\f1\\Trial3",# 2 jumps
    # "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240909\\f1\\Trial5",# 1 jump
       "Y:\Data\FCI\\Hedwig\\FB4P_b_SS60296_sytGC7f\\240906\\f2\\Trial2",# 2 jumps
       r"Y:\Data\FCI\Hedwig\FB4P_b_SS60296\250311\f1\Trial1"
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
alldat = {}
for i,datadir in enumerate(datadirs):
    
    print(datadir)
    #plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    #trj,ca = cxt.fc.mean_traj_nF()
    
    ca,trj = cxt.mean_traj_jump(bins=bins)
    alldat.update({'ca_' +str(i):ca,
                   'traj_'+str(i):trj})
    
    
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

savedict = {'ca_fb4pb': ca_fb4pb,'traj_fb4pb':traj_fb4pb,'ca_fb4pb_realtime':cart_fb4pb,'fb4pb_alldat':alldat}

# FB4R
plt.close('all')
datadirs = ["Y:\Data\FCI\\Hedwig\\SS61646_FB4R\\240828\\f3\\Trial1",
            "Y:\Data\FCI\\Hedwig\\SS61646_FB4R\\240910\\f1\\Trial1",
            "Y:\Data\FCI\\Hedwig\\SS61645_FB4R\\240911\\f1\\Trial3" # Only 2 jumps
            ]


traj_fb4r = np.empty((bins*2,2,len(datadirs)))
ca_fb4r = np.empty((bins*2,len(datadirs)))
cart_fb4r = np.empty((timecourse*10,len(datadirs)))
alldat = {}
for i,datadir in enumerate(datadirs):
    #plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    
    ca,trj = cxt.mean_traj_jump(bins=bins)
    traj_fb4r[:,:,i] = np.mean(trj,axis=2)
    ca_fb4r[:,i] = np.mean(ca,axis=1)
    alldat.update({'ca_' +str(i):ca,
                   'traj_'+str(i):trj})
    
    ct = cxt.mean_ca_realtime_jump(timecourse)
    
    #ct[np.abs(ct)==0.0] = np.nan
    x = np.nanmean(ct,axis=1)
    cart_fb4r[:,i] = x
    
    # plt.figure()
    # cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=0,set_cmx =True,cmx=0.25)
    # savename = os.path.join(savedir,'FB4R_Ca_traj' +str(i)+'.pdf')
    # plt.savefig(savename)

savedict.update({'ca_fb4r': ca_fb4r,'traj_fb4r':traj_fb4r,'ca_fb4r_realtime':cart_fb4r,'fb4r_alldat':alldat})


# FB5I
datadirs = [
   "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240628\\f1\\Trial2",#Nice
            "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240917\\f1\\Trial2",#Best for this fly
           "Y:\\Data\\FCI\\Hedwig\\FB5I_SS100553\\240917\\f3\\Trial3",
           ]

traj_fb5i = np.empty((bins*2,2,len(datadirs)))
ca_fb5i = np.empty((bins*2,len(datadirs)))
cart_fb5I = np.empty((timecourse*10,len(datadirs)))
alldat ={}
for i,datadir in enumerate(datadirs):
   # plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    
    ca,trj = cxt.mean_traj_jump(bins=bins)
    #cxt.fc.jump_heat(np.mean(trj,axis=2),np.mean(ca,axis=1),xoffset=0)
    alldat.update({'ca_' +str(i):ca,
                   'traj_'+str(i):trj})
    
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

savedict.update({'ca_fb5i': ca_fb5i,'traj_fb5i':traj_fb5i,'ca_fb5i_realtime':cart_fb5I,'fb5i_alldat':alldat})


# FB4X
datadirs = [
   "Y:\Data\\FCI\\Hedwig\\SS70711_FB4X\\241030\\f3\\Trial3",
   "Y:\\Data\\FCI\\Hedwig\\SS70711_FB4X\\241031\\f1\\Trial3"
           ]

traj_fb4x = np.empty((bins*2,2,len(datadirs)))
ca_fb4x = np.empty((bins*2,len(datadirs)))
cart_fb4x = np.empty((timecourse*10,len(datadirs)))
alldat = {}
for i,datadir in enumerate(datadirs):
    #plt.figure()
    xoffset= 0
    cxt = CX_tan(datadir)
    
    ca,trj = cxt.mean_traj_jump(bins=bins)
    #cxt.fc.jump_heat(np.mean(trj,axis=2),np.mean(ca,axis=1),xoffset=0)
    
    traj_fb4x[:,:,i] = np.mean(trj,axis=2)
    ca_fb4x[:,i] = np.mean(ca,axis=1)
    
    alldat.update({'ca_' +str(i):ca,
                   'traj_'+str(i):trj})
    
    ct = cxt.mean_ca_realtime_jump(timecourse)
    #ct[np.abs(ct)==0.0] = np.nan
    x = np.nanmean(ct,axis=1)
    cart_fb4x[:,i] = x
    
    # plt.figure()
    # cxt.fc.mean_traj_heat_jump(cxt.fc.ca,xoffset=0,set_cmx =True,cmx=0.25)
    # savename = os.path.join(savedir,'FB4X_Ca_traj' +str(i)+'.pdf')
    # plt.savefig(savename)

savedict.update({'ca_fb4x': ca_fb4x,'traj_fb4x':traj_fb4x,'ca_fb4x_realtime': cart_fb4x,'fb4x_alldat':alldat})

# Save data as pickle

u = ugn
u.save_pick(savedict,os.path.join(savedir,'TN_mean_jump.pickle'))
data = u.load_pick(os.path.join(savedir,'TN_mean_jump.pickle'))

#%% plot max
bins = 100
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
newtime = np.arange(0,bins*2,2)
oldtime = np.arange(0,bins*2)
for i,c in enumerate(cells):
    tdat = data['ca_'+c]
    tdat = tdat/np.max(tdat,axis=0)
    plt.fill([0,bins,bins,0],[-1,-1,1,1],color=[0.8,0.8,0.8])
    plt.plot([bins,bins],[-1,1],color='k',linestyle='--')
   #plt.plot(tdat,color=colours[i,:],alpha=0.2)
    plt.plot([bins*2,bins*2],[-1,1],color='k',linestyle='--')
    y = np.mean(tdat,axis=1)
    ynew = np.interp(newtime,oldtime,y)
    plt.plot(y,color=colours[i,:])
    
#plt.ylim([0,1])
plt.xlim([0,bins*2])
plt.xlabel('Normalised time')
plt.ylabel('Norm dF/F')
plt.savefig(os.path.join(savedir,'SummaryTNpseudotime.pdf'))
#%% plot all ca
plt.close('all')

from Utilities.utils_plotting import uplt
tneuron = 'fb5i'
alldat = data[tneuron+'_alldat']
camean = data['ca_'+tneuron]
trajmean = data['traj_'+tneuron]
numflies = camean.shape[1]
cmin= -0.25
cmax = 0.25
offset = 50
for i in range(numflies):
    aca = alldat['ca_'+str(i)]
    atraj = alldat['traj_'+str(i)]
    calen= aca.shape[1]
    fig,ax=plt.subplots()
    
    tmean = trajmean[:,:,i]
    pedge = tmean[-1,0]
    plt.fill(np.array([-10,0,0,-10])+pedge,np.array([0,0,50,50]),color=[0.7,0.7,0.7])
    plt.fill(np.array([-10,0,0,-10]),np.array([0,0,-50,-50]),color=[0.7,0.7,0.7])
    
    for ic in range(calen):
        ca =aca[:,0]
        x = atraj[:,0,ic]
        y = atraj[:,1,ic]
        x = np.append([np.nan,np.nan],x)
        y = np.append([np.nan,np.nan],y)
        ca = np.append([cmin,cmax],ca)
        # y[0] = cmin
        # y[1] = cmax
        ca[ca<cmin] =cmin
        ca[ca>cmax] = cmax
        #uplt.coloured_line_simple(atraj[:,0,ic],atraj[:,1,ic],ca,'coolwarm',-0.25,0.25,alpha=0.5)
        uplt.coloured_line(x,y,aca[:,ic],ax,cmap='coolwarm',alpha=0.3)
    
        
    plt.fill(np.array([-10,0,0,-10])+pedge+offset,np.array([0,0,50,50]),color=[0.7,0.7,0.7])
    plt.fill(np.array([-10,0,0,-10])+offset,np.array([0,0,-50,-50]),color=[0.7,0.7,0.7])
    
    tca = camean[:,i]
    tca[tca<cmin] =cmin
    tca[tca>cmax] = cmax
    x = tmean[:,0]
    y = tmean[:,1]
    x = np.append([np.nan,np.nan],x)
    y = np.append([np.nan,np.nan],y)
    tca = np.append([cmin,cmax],tca)
    uplt.coloured_line(x+offset,y,tca,ax,cmap='coolwarm',linewidth=2)
    ax.set_aspect('equal')  
    ax.set_xlim([-25,80])
    ax.set_ylim([-50,60])
    plt.savefig(os.path.join(savedir,'ReturnMean_'+tneuron+'_'+str(i) +'.pdf'))
    plt.savefig(os.path.join(savedir,'ReturnMean_'+tneuron+'_'+str(i) +'.png'))
#%%

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
    
    
    
    
    