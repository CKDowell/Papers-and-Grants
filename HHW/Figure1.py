# -*- coding: utf-8 -*-
"""
Created on Fri May 17 13:46:34 2024

@author: dowel

This script will produce the panels for Figure 1 of the Helen Hay Whitney grant
application.

Panels: 
    Example ET with FSB activity. 
    Histogram of orientation

"""

#%%
from analysis_funs.regression import fci_regmodel

import numpy as np
import pandas as pd
import src.utilities.funcs as fc
from analysis_funs.optogenetics import opto 
import os
import matplotlib.pyplot as plt 
from src.utilities import imaging as im
from skimage import io, data, registration, filters, measure
from scipy import signal as sg
from analysis_funs.CX_imaging import CX
from analysis_funs.CX_analysis_col import CX_a
from src.utilities import funcs as fn
plt.rcParams['pdf.fonttype'] = 42 
#%%
savedir = "Y:\\Applications\\HHW\\Figures\\Figure1"
datadirs = [
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f1\\Trial3",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f2\\Trial3",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240502\\f1\\Trial2",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240514\\f1\\Trial2"]

#%% Panel: Example data - choose carefully
#%% Load example data
td = 3
datadir = datadirs[td]
d = datadir.split("\\")
name = d[-3] + '_' + d[-2] + '_' + d[-1]
cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
#%% Plot it
plt.close('all')

cxa.point2point_heat(2550,3000,toffset=0,arrowpoint=np.array([20,80,122,193,227,335]),regions=['eb','fsb_upper','fsb_lower'])
plt.savefig(os.path.join(savedir,'ReturnLongTrace2.pdf'))
plt.savefig(os.path.join(savedir,'ReturnLongTrace2.png'))
#%%
cxa.point2point_heat(2550,3000,toffset=0,arrowpoint=np.array([20,122,227]),regions=['eb','fsb_upper','fsb_lower'])
plt.savefig(os.path.join(savedir,'LeaveLongTrace2.pdf'))
plt.savefig(os.path.join(savedir,'LeaveLongTrace2.png'))
#%%
cxa.point2point_heat(2550,3000,toffset=0,arrowpoint=np.array([80,193,335]),regions=['eb','fsb_upper','fsb_lower'])
plt.savefig(os.path.join(savedir,'ReturnLongTrace2.pdf'))
plt.savefig(os.path.join(savedir,'ReturnLongTrace2.png'))
#%%
plt.close('all')
cxa.point2point_heat(3950,4450,toffset=0,arrowpoint=np.array([20,75,145,193,400,470]),regions=['eb','fsb_upper','fsb_lower'])
plt.savefig(os.path.join(savedir,'ReturnLongExample.pdf'))
plt.savefig(os.path.join(savedir,'ReturnLongExample.png'))
#%%

cxa.point2point_heat(8500,9000,toffset=0,arrowpoint=np.array([20,75,145,193,400,470]),regions=['eb','fsb_upper','fsb_lower'])
#%%
cxa.point2point_heat(0,500,toffset=0,arrowpoint=np.array([100,250,390]))
plt.savefig(os.path.join(savedir,'Amenotaxis.png'))
#%% Histograms of returns
#%% Load data

for i,datadir  in enumerate(datadirs):
    print(i)
    d = datadir.split("\\")
    name = d[-3] + '_' + d[-2] + '_' + d[-1]
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    meno_array,et_array,et_array_ex,pltbins = cxa.plume_meno_comp(diff_phase=True,diff_val='eb')
    meno_array = np.expand_dims(meno_array,2)
    et_array = np.expand_dims(et_array,2)
    et_array_ex = np.expand_dims(et_array_ex,2)
    if i==0:
        plt_meno = meno_array
        plt_et = et_array
        plt_et_ex = et_array_ex
    else:
        plt_meno = np.append(plt_meno,meno_array,axis=2)
        plt_et = np.append(plt_et,et_array,axis=2)
        plt_et_ex = np.append(plt_et_ex,et_array,axis=2)

#%% Plot it
plt.close('all')
colours = np.array([[0,0,0],[0.3,0.3,1],[1,0.3,1]])
i = 1
fig, ax1 = plt.subplots()

tdat = plt_et[:,i,:]
pltmean = np.mean(tdat,axis=1)
pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[i,:],alpha=0.5)
ax1.plot(pltbins,pltmean,color=colours[i,:],alpha=1)
ax1.set_xlabel('FC2-EPG phase (deg)',fontsize=15)
ax1.set_xticks([-180,-90,0,90,180],labels = [-180,-90,0,90,180],fontsize=12)
ax1.set_ylabel('Probability',fontsize=15)
ax1.set_title('Plume jump returns',fontsize=15)
ax1.set_yticks(np.arange(0,0.16,0.05),labels =[0,0.05,0.1,0.15], fontsize=12)

ax1.set_ylim([0,0.175])
ax1.plot([0,0],[0,0.2],color='k',linestyle='--')
ax1.set_xlim([-180,180])
plt.savefig(os.path.join(savedir,'PlumeJumpReturnHisto_FC2upper-EPG.pdf'))

i = 1
fig, ax1 = plt.subplots()

tdat = plt_meno[:,i,:]
pltmean = np.mean(tdat,axis=1)
pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[i,:],alpha=0.5)
ax1.plot(pltbins,pltmean,color=colours[i,:],alpha=1)
ax1.set_xlabel('FC2-EPG phase (deg)',fontsize=15)
ax1.set_xticks([-180,-90,0,90,180],labels = [-180,-90,0,90,180],fontsize=12)

ax1.set_ylabel('Probability - EB',fontsize=15)
ax1.set_title('Amenotaxis',fontsize=15)
ax1.set_ylim([0,0.175])
ax1.plot([0,0],[0,0.2],color='k',linestyle='--')
ax1.set_yticks(np.arange(0,0.16,0.05),labels =[0,0.05,0.1,0.15], fontsize=12)

ax1.set_xlim([-180,180])
plt.savefig(os.path.join(savedir,'AmenotaxisHisto_FC2upper-EPG.pdf'))
#plm_mn = np.sum(pltmean[pltbins<0])
#ax1.plot([-180,0],[plm_mn,plm_mn],color='k',linestyle='--')


i = 1
fig, ax1 = plt.subplots()

tdat = plt_et_ex[:,i,:]
pltmean = np.mean(tdat,axis=1)
pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[0,:],alpha=0.5)
ax1.plot(pltbins,pltmean,color=colours[0,:],alpha=1)

tdat = plt_et[:,i,:]
pltmean = np.mean(tdat,axis=1)
pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[i,:],alpha=0.5)
ax1.plot(pltbins,pltmean,color=colours[i,:],alpha=1)
ax1.set_xlabel('FC2-EPG phase (deg)')
ax1.set_xticks([-180,-90,0,90,180])
ax1.set_ylabel('Probability - EB')
ax1.set_title('Plume jump returns')
ax1.set_ylim([0,0.175])
ax1.plot([0,0],[0,0.2],color='k',linestyle='--')
ax1.set_xlim([-180,180])
#%% Plot it
plt.close('all')
colours = np.array([[0,0,0],[0.3,0.3,1],[1,0.3,1]])

fig, ax1 = plt.subplots()
i = 0
tdat = plt_et[:,i,:]
pltmean = np.mean(tdat,axis=1)
pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
# ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[i,:],alpha=0.5)
# ax1.plot(pltbins,pltmean,color=colours[i,:],alpha=1)
ax1.plot([0,0],[0,max(pltmean+pltse)],linestyle='--',color='k',zorder=3)
ax1.set_xlabel('EPG-FC2 phase (deg)')
ax1.set_xticks([-180,-90,0,90,180])
ax1.set_ylabel('Probability')
ax1.set_title('Plume jump returns')
ax1.set_ylim([0,0.175])

for i in range(2):
    #plt.plot(pltbins,plt_et[:,i,:],color=colours[i,:],alpha=0.5)
    tdat = plt_et[:,i+1,:]
    pltmean = np.mean(tdat,axis=1)
    pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
    ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[i+1,:],alpha=0.5)
    ax1.plot(pltbins,pltmean,color=colours[i+1,:],alpha=1)
    #
fig.tight_layout()  
plt.show()
plt.savefig(os.path.join(savedir,'PlumeJumpReturnHisto.png'))

fig, ax1 = plt.subplots()
i = 0
tdat = plt_et_ex[:,i,:]
pltmean = np.mean(tdat,axis=1)
pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
# ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[i,:],alpha=0.5)
# ax1.plot(pltbins,pltmean,color=colours[i,:],alpha=1)
ax1.plot([0,0],[0,max(pltmean+pltse)],linestyle='--',color='k',zorder=3)
ax1.set_xlabel('EPG-FC2 phase (deg)')
ax1.set_xticks([-180,-90,0,90,180])
ax1.set_ylabel('Probability')
ax1.set_title('Plume exits')
ax1.set_ylim([0,0.22])
for i in range(2):
    #plt.plot(pltbins,plt_et[:,i,:],color=colours[i,:],alpha=0.5)
    tdat = plt_et_ex[:,i+1,:]
    pltmean = np.mean(tdat,axis=1)
    pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
    ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[i+1,:],alpha=0.5)
    ax1.plot(pltbins,pltmean,color=colours[i+1,:],alpha=1)
    #
#ax2.set_ylim([0,0.15])
fig.tight_layout()  
plt.show()
plt.yticks([0,0.05,0.1,0.15,0.2])
plt.savefig(os.path.join(savedir,'PlumeExitHisto.png'))
plt.savefig(os.path.join(savedir,'PlumeExitHisto.pdf'))

fig, ax1 = plt.subplots()
i = 0
tdat = plt_meno[:,i,:]
pltmean = np.mean(tdat,axis=1)
pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
ax1.plot([0,0],[0,max(pltmean+pltse)],linestyle='--',color='k',zorder=3)
ax1.set_xlabel('EPG-FC2 phase (deg)')
ax1.set_xticks([-180,-90,0,90,180])
ax1.set_ylabel('Probability')
ax1.set_title('Amenotaxis')
ax1.set_ylim([0,0.175])

for i in range(2):
    #plt.plot(pltbins,plt_et[:,i,:],color=colours[i,:],alpha=0.5)
    tdat = plt_meno[:,i+1,:]
    pltmean = np.mean(tdat,axis=1)
    pltse  = np.std(tdat,axis=1)/np.sqrt(len(datadirs))
    ax1.fill_between(pltbins,pltmean-pltse,pltmean+pltse,color=colours[i+1,:],alpha=0.5)
    ax1.plot(pltbins,pltmean,color=colours[i+1,:],alpha=1)
    #
#ax2.set_ylim([0,0.15])
fig.tight_layout()  
plt.show()
plt.savefig(os.path.join(savedir,'AmenotaxisHisto.png'))
        
        