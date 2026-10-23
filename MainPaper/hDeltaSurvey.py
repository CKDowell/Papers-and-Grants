# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 13:46:48 2026

@author: dowel

Script does a comparison of hDelta and other FC2 input activity for
the first figure of the paper



"""
#%% Import packages
from analysis_funs.regression import fci_regmodel

import numpy as np
import pandas as pd
import src.utilities.funcs as fc
from analysis_funs.optogenetics import opto 
import os
import matplotlib.pyplot as plt 
from analysis_funs.utilities import imaging as im
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

#%% Load data - pulses
groups = {
    
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
    "hDeltaH_ACV_pulses": [
    r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f1\Trial4",
    r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f2\Trial5",
    r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f1\Trial5",
    r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f2\Trial4"
    ],
    "hDeltaH_Oct_pulses":[
    r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f1\Trial5",
    r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f2\Trial4",
    r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f1\Trial4",
    r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f2\Trial5"
    ],
    "hDeltaF_ACV_pulses": [
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260805\f1\Trial4',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260811\f2\Trial5',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260812\f1\Trial3',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260813\f1\Trial5',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260813\f2\Trial3',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260818\f1\Trial4'
    ],
    "hDeltaF_Oct_pulses": [
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260805\f1\Trial5',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260811\f2\Trial4',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260812\f1\Trial4',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260813\f1\Trial4',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260813\f2\Trial4',
    r'Y:\Data\FCI\Hedwig\hDeltaF_SS53683_sytGCaMP7f\260818\f1\Trial3'
    ],
    
    "hDeltaD_ACV_pulses": [ 
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260819\f1\Trial3', # ACV pulses 
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260819\f2\Trial4',# ACV pulses 
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260820\f1\Trial3', #ACV pulses
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260820\f2\Trial4', #ACV pulses
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260825\f1\Trial3', #ACV pulses
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260826\f1\Trial4', # ACV pulses
        
        ],
    "hDeltaD_Oct_pulses": [
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260819\f1\Trial4', # Octanol pulses
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260819\f2\Trial3',# Octanol pulses
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260820\f1\Trial4', # Oct pulses
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260820\f2\Trial3',# Oct pulses
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260825\f1\Trial4', #Oct pulses
        r'Y:\Data\FCI\Hedwig\hDeltaD_SS65202_sytGCaMP7f\260826\f1\Trial3', #Oct pulses
        ],
    
    "PFNpm_ACV_pulses": [r'Y:\Data\FCI\Hedwig\PFNpm_SS00191\260826\f1\Trial4', # ACV pulses
        ],
    
    "PFNpm_Oct_pulses": [ r'Y:\Data\FCI\Hedwig\PFNpm_SS00191\260826\f1\Trial5', # Oct pulses
        ]
}
pulse_data = {}
regions = {'hDeltaC': np.array(['eb','fsb_upper','fsb_lower']),
           'hDeltaJ': np.array(['eb','fsb_upper','fsb_lower']),
           'hDeltaH': np.array(['eb','fsb_upper','fsb_lower']),
           'hDeltaF': np.array(['fsb']),
           'hDeltaD': np.array(['fsb']),
           'PFNpm': np.array(['pfn_l_16','pfn_r_16'])
           }

savedir = r'Y:\Papers, Review, Theses\EarlyDraft26\Figure1_panels'
for group_name, datadirs in groups.items():
    pulse_data[group_name] = {}
    gname = group_name.split('_')[0]
    for i, datadir in enumerate(datadirs):
        print(f"[{group_name}] {datadir}")
        etp = ET_paper(datadir,regions=regions[gname])
        if gname=='PFNpm':
            etp.cxa.pdat['wedges_pb_all'] = np.append(etp.cxa.pdat['wedges_'+regions[gname][0]],etp.cxa.pdat['wedges_'+regions[gname][1]],axis=1)
        
        
        
        pulse_data[group_name][i] = etp




#%% ACV pulses - mean


plt.close('all')
minlen = 50
tbef = 70
x = np.arange(-tbef,600)/10
dchoices = [
    'hDeltaJ','hDeltaC','hDeltaH','hDeltaF','hDeltaD','PFNpm'
    ]
regions_plot = {'hDeltaC': np.array(['fsb_upper']),
           'hDeltaJ': np.array(['fsb_upper']),
           'hDeltaH': np.array(['fsb_upper']),
           'hDeltaF': np.array(['fsb']),
           'hDeltaD': np.array(['fsb']),
           'PFNpm': np.array(['pb_all']),
           }


colours = uplt.hDelta_colours()

for di,d in enumerate(dchoices):
    dname = d + '_ACV_pulses'
    data = pulse_data[dname]
    print(len(data))
    meandat = np.zeros((len(x),len(data)))
    for f in data:
        cxa = data[f].cxa
        
        amp = np.mean(cxa.pdat['wedges_'+regions_plot[d][0]],axis=1)
        jumps = cxa.get_entries_exits_like_jumps()
        jumplen = jumps[:,2]-jumps[:,1]
        jumps = jumps[jumplen>minlen]
        jumps  = jumps[:15]
        all_js = np.zeros((len(x),len(jumps)))
        for i,j in enumerate(jumps):
            dx = np.arange(j[1]-tbef,j[2])
            dxnorm = np.arange(j[0],j[1])
            dxnorm = np.arange(j[0]-10,j[0])
            nrm = np.mean(amp[dxnorm])
           # plt.plot(amp[dx],color='k',alpha=0.3)
            tamp = amp[dx]-nrm
            amplen  = np.min([600+tbef,len(tamp)])
            all_js[:amplen,i] = tamp[:amplen]
            
        all_js[all_js==0] = np.nan    
        meandat[:,int(f)] = np.nanmean(all_js,axis=1)

    
    plt.plot(x[:200+tbef],meandat[:200+tbef,:],color=colours[d],alpha=0.3)
    plt.plot(x[:200+tbef],np.nanmean(meandat[:200+tbef,:],1),color=colours[d],alpha=1,linewidth=3,label=d)
plt.legend()
plt.ylabel('mean dF/F0',fontsize=15)
plt.xlabel('time from odour offset (s)',fontsize=15)

plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(-.6,1,.2),fontsize=15)
plt.ylim([-.6,.5])
plt.fill([-2,0,0,-2],[-.5,-.5,.2,.2],linewidth=0,color=[.8,.8,.8])
plt.xlim([-5,15])
plt.pause(.1)
plt.savefig(os.path.join(savedir,'ACV_pulse_mean.pdf'))

#%% Tau determination
from scipy.optimize import curve_fit

def sat_exp(t,A,tau,y0):
    return A*(1-np.exp(-t/tau))+y0

tau_start = 70
tau_end = 170

dchoices = ['hDeltaC']
for di,d in enumerate(dchoices):
    dname = d + '_ACV_pulses'
    data = pulse_data[dname]
    print(len(data))
    meandat = np.zeros((len(x),len(data)))
    for f in data:
        cxa = data[f].cxa

        amp = np.mean(cxa.pdat['wedges_'+regions_plot[d][0]],axis=1)
        jumps = cxa.get_entries_exits_like_jumps()
        jumplen = jumps[:,2]-jumps[:,1]
        jumps = jumps[jumplen>minlen]
        jumps  = jumps[:15]
        all_js = np.zeros((len(x),len(jumps)))
        for i,j in enumerate(jumps):
            dx = np.arange(j[1]-tbef,j[2])
            dxnorm = np.arange(j[0],j[1])
            dxnorm = np.arange(j[0]-10,j[0])
            nrm = np.mean(amp[dxnorm])
           # plt.plot(amp[dx],color='k',alpha=0.3)
            tamp = amp[dx]-nrm
            amplen  = np.min([600+tbef,len(tamp)])
            all_js[:amplen,i] = tamp[:amplen]

        all_js[all_js==0] = np.nan
        meandat[:,int(f)] = np.nanmean(all_js,axis=1)

    t_fit = x[tau_start:tau_end]-x[tau_start]
    taus = np.full(meandat.shape[1],np.nan)
    fits = np.full((tau_end-tau_start,meandat.shape[1]),np.nan)
    for c in range(meandat.shape[1]):
        y_fit = meandat[tau_start:tau_end,c]
        if np.any(np.isnan(y_fit)):
            continue
        p0 = [y_fit[-1]-y_fit[0],1.0,y_fit[0]]
        try:
            popt,_ = curve_fit(sat_exp,t_fit,y_fit,p0=p0,maxfev=5000)
            taus[c] = popt[1]
            fits[:,c] = sat_exp(t_fit,*popt)
        except RuntimeError:
            pass

    print(d,'tau (s): mean=',np.nanmean(taus),'per fly=',taus)

    plt.figure()
    plt.plot(x[:200+tbef],meandat[:200+tbef,:],color=colours[d],alpha=0.3)
    plt.plot(x[:200+tbef],np.nanmean(meandat[:200+tbef,:],1),color=colours[d],alpha=1,linewidth=3,label=d)

    # Per-fly fit quality plots
    ncols = int(np.ceil(np.sqrt(meandat.shape[1])))
    nrows = int(np.ceil(meandat.shape[1]/ncols))
    fig,axs = plt.subplots(nrows,ncols,figsize=(3*ncols,2.5*nrows),squeeze=False)
    axs = axs.flatten()
    for c in range(meandat.shape[1]):
        ax = axs[c]
        ax.plot(x[:200+tbef],meandat[:200+tbef,c],color='k',alpha=0.5,label='data')
        ax.plot(t_fit+x[tau_start],meandat[tau_start:tau_end,c],color='k',linewidth=2,label='fit window')
        if not np.isnan(taus[c]):
            ax.plot(t_fit+x[tau_start],fits[:,c],color='r',linewidth=2,label='fit')
            ax.set_title('fly %d, tau=%.2f s'%(c,taus[c]),fontsize=9)
        else:
            ax.set_title('fly %d, fit failed'%c,fontsize=9)
        ax.axvline(0,color='grey',linestyle='--',linewidth=0.5)
    for c in range(meandat.shape[1],len(axs)):
        axs[c].axis('off')
    axs[0].legend(fontsize=7)
    fig.suptitle(d+' tau fit quality')
    fig.tight_layout()



#%% ACV pulses - first pulse
plt.close('all')
minlen = 50
tbef = 70
x = np.arange(-tbef,600)/10
dchoices = [
    'hDeltaJ','hDeltaC','hDeltaH','hDeltaF','hDeltaD','PFNpm'
    ]
regions_plot = {'hDeltaC': np.array(['fsb_upper']),
           'hDeltaJ': np.array(['fsb_upper']),
           'hDeltaH': np.array(['fsb_upper']),
           'hDeltaF': np.array(['fsb']),
           'hDeltaD': np.array(['fsb']),
           'PFNpm': np.array(['pb_all']),
           }


colours = uplt.hDelta_colours()

for di,d in enumerate(dchoices):
    dname = d + '_ACV_pulses'
    data = pulse_data[dname]
    print(len(data))
    meandat = np.zeros((len(x),len(data)))
    for i,f in enumerate(data):
        cxa = data[f].cxa
        
        amp = np.mean(cxa.pdat['wedges_'+regions_plot[d][0]],axis=1)
        jumps = cxa.get_entries_exits_like_jumps()
        jumplen = jumps[:,2]-jumps[:,1]
        jumps = jumps[jumplen>minlen]
        jumps  = jumps[:15]
        all_js = np.zeros((len(x),len(jumps)))
        
        j = jumps[0,:]
        dx = np.arange(j[1]-tbef,j[2])
        dxnorm = np.arange(j[0],j[1])
        dxnorm = np.arange(j[0]-10,j[0])
        nrm = np.mean(amp[dxnorm])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]-nrm
        amplen  = np.min([600+tbef,len(tamp)])
        meandat[:amplen,i] = tamp[:amplen]
        
        
    
    plt.plot(x[:200+tbef],meandat[:200+tbef,:],color=colours[d],alpha=0.3)
    plt.plot(x[:200+tbef],np.nanmean(meandat[:200+tbef,:],1),color=colours[d],alpha=1,linewidth=3,label=d)
plt.legend()
plt.ylabel('mean dF/F0',fontsize=15)
plt.xlabel('time from odour offset (s)',fontsize=15)

plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(-.6,1,.2),fontsize=15)
#plt.ylim([-.5,.2])
plt.ylim([-.6,.5])

plt.fill([-2,0,0,-2],[-.5,-.5,.2,.2],linewidth=0,color=[.8,.8,.8])
plt.xlim([-5,15])
plt.pause(.1)
plt.savefig(os.path.join(savedir,'ACV_pulse_first.pdf'))
#%% Octanol pulses

plt.close('all')
minlen = 50
tbef = 70
x = np.arange(-tbef,600)/10
dchoices = [
    'hDeltaJ','hDeltaC','hDeltaH','hDeltaF','hDeltaD','PFNpm'
    ]
regions_plot = {'hDeltaC': np.array(['fsb_upper']),
           'hDeltaJ': np.array(['fsb_upper']),
           'hDeltaH': np.array(['fsb_upper']),
           'hDeltaF': np.array(['fsb']),
           'hDeltaD': np.array(['fsb']),
           'PFNpm': np.array(['pb_all']),
           }

colours = uplt.hDelta_colours()

for di,d in enumerate(dchoices):
    dname = d + '_Oct_pulses'
    data = pulse_data[dname]
    print(len(data))
    meandat = np.zeros((len(x),len(data)))
    for f in data:
        cxa = data[f].cxa
        
        amp = np.mean(cxa.pdat['wedges_'+regions_plot[d][0]],axis=1)
        jumps = cxa.get_entries_exits_like_jumps(odour='universal')
        jumplen = jumps[:,2]-jumps[:,1]
        jumps = jumps[jumplen>minlen]
        jumps  = jumps[:15]
        print(jumps[:,1]-jumps[:,0])
        all_js = np.zeros((len(x),len(jumps)))
        for i,j in enumerate(jumps):
            dx = np.arange(j[1]-tbef,j[2])
            dxnorm = np.arange(j[0],j[1])
            dxnorm = np.arange(j[0]-10,j[0])
            nrm = np.mean(amp[dxnorm])
           # plt.plot(amp[dx],color='k',alpha=0.3)
            tamp = amp[dx]-nrm
            amplen  = np.min([600+tbef,len(tamp)])
            all_js[:amplen,i] = tamp[:amplen]
            
        all_js[all_js==0] = np.nan    
        meandat[:,int(f)] = np.nanmean(all_js,axis=1)

    
    plt.plot(x[:200+tbef],meandat[:200+tbef,:],color=colours[d],alpha=0.3)
    plt.plot(x[:200+tbef],np.nanmean(meandat[:200+tbef,:],1),color=colours[d],alpha=1,linewidth=3,label=d)
plt.legend()
plt.ylabel('mean dF/F0',fontsize=15)
plt.xlabel('time from odour offset (s)',fontsize=15)

plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(-.6,1,.2),fontsize=15)
plt.ylim([-.6,.5])
plt.fill([-2,0,0,-2],[-.5,-.5,.2,.2],linewidth=0,color=[.8,.8,.8])
plt.xlim([-5,15])
plt.pause(.1)
plt.savefig(os.path.join(savedir,'Oct_pulse_mean.pdf'))
#%% Oct pulses - first pulse
plt.close('all')
minlen = 50
tbef = 70
x = np.arange(-tbef,600)/10
dchoices = [
    'hDeltaJ','hDeltaC','hDeltaH','hDeltaF','hDeltaD','PFNpm'
    ]
regions_plot = {'hDeltaC': np.array(['fsb_upper']),
           'hDeltaJ': np.array(['fsb_upper']),
           'hDeltaH': np.array(['fsb_upper']),
           'hDeltaF': np.array(['fsb']),
           'hDeltaD': np.array(['fsb']),
           'PFNpm': np.array(['pb_all']),
           }


colours = uplt.hDelta_colours()

for di,d in enumerate(dchoices):
    dname = d + '_Oct_pulses'
    data = pulse_data[dname]
    print(len(data))
    meandat = np.zeros((len(x),len(data)))
    for i,f in enumerate(data):
        cxa = data[f].cxa
        
        amp = np.mean(cxa.pdat['wedges_'+regions_plot[d][0]],axis=1)
        jumps = cxa.get_entries_exits_like_jumps(odour='universal')
        jumplen = jumps[:,2]-jumps[:,1]
        jumps = jumps[jumplen>minlen]
        jumps  = jumps[:15]
        all_js = np.zeros((len(x),len(jumps)))
        
        j = jumps[0,:]
        dx = np.arange(j[1]-tbef,j[2])
        dxnorm = np.arange(j[0],j[1])
        dxnorm = np.arange(j[0]-10,j[0])
        nrm = np.mean(amp[dxnorm])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]-nrm
        amplen  = np.min([600+tbef,len(tamp)])
        meandat[:amplen,i] = tamp[:amplen]
        
        
    
    plt.plot(x[:200+tbef],meandat[:200+tbef,:],color=colours[d],alpha=0.3)
    plt.plot(x[:200+tbef],np.nanmean(meandat[:200+tbef,:],1),color=colours[d],alpha=1,linewidth=3,label=d)
plt.legend()
plt.ylabel('mean dF/F0',fontsize=15)
plt.xlabel('time from odour offset (s)',fontsize=15)

plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(-.6,1,.2),fontsize=15)
plt.ylim([-.6,.5])
plt.fill([-2,0,0,-2],[-.5,-.5,.2,.2],linewidth=0,color=[.8,.8,1])
plt.xlim([-5,15])
plt.pause(.1)
plt.savefig(os.path.join(savedir,'Oct_pulse_first.pdf'))
#%% Stacked heatmaps of wedge activity from first ACV odour onset
plt.close('all')
tplot = 150 # seconds of recording to show after the first ACV onset - tune if the plot gets too crowded
fs = 10 # imaging sample rate (Hz)
nplot = int(tplot*fs)

dchoices = [
    'hDeltaJ','hDeltaC','hDeltaH','hDeltaF','hDeltaD','PFNpm'
    ]
regions_plot = {'hDeltaC': np.array(['fsb_upper']),
           'hDeltaJ': np.array(['fsb_upper']),
           'hDeltaH': np.array(['fsb_upper']),
           'hDeltaF': np.array(['fsb']),
           'hDeltaD': np.array(['fsb']),
           'PFNpm': np.array(['pb_all']),
           }

for d in dchoices:
    dname = d+'_ACV_pulses'
    data = pulse_data[dname]
    nflies = len(data)
    fig,axs = plt.subplots(nflies,1,figsize=(8,1.5*nflies),sharex=True)
    axs = np.atleast_1d(axs)
    for i,f in enumerate(data):
        cxa = data[f].cxa
        wedges = cxa.pdat['wedges_'+regions_plot[d][0]]
        jumps = cxa.get_entries_exits_like_jumps()
        onset0 = jumps[0,0]
        edx = np.arange(onset0,min(onset0+nplot,len(wedges)))
        plotdat = wedges[edx,:].T
        t = (edx-onset0)/fs

        ax = axs[i]
        ax.imshow(plotdat,aspect='auto',cmap='Blues',
                  extent=[t[0],t[-1],plotdat.shape[0],0])

        onsets = jumps[:,0]
        onsets = onsets[(onsets>=onset0) & (onsets<onset0+nplot)]
        for o in onsets:
            ax.axvline((o-onset0)/fs,color='r',linewidth=2)
        ax.set_ylabel('Fly '+str(f))

    axs[-1].set_xlabel('Time from first ACV onset (s)')
    fig.suptitle(d+' ACV wedge activity')
    plt.tight_layout()
    plt.pause(.1)
    plt.savefig(os.path.join(savedir,d+'_ACV_wedge_heatmaps.pdf'))
#%% Stacked heatmaps - one example fly per neuron type
plt.close('all')
tplot = 150 # seconds of recording to show after the first ACV onset - tune if the plot gets too crowded
fs = 10 # imaging sample rate (Hz)
nplot = int(tplot*fs)

dchoices = [
    'hDeltaJ','hDeltaC','hDeltaH','hDeltaF','hDeltaD',#'PFNpm'
    ]
regions_plot = {'hDeltaC': np.array(['fsb_upper']),
           'hDeltaJ': np.array(['fsb_upper']),
           'hDeltaH': np.array(['fsb_upper']),
           'hDeltaF': np.array(['fsb']),
           'hDeltaD': np.array(['fsb']),
          # 'PFNpm': np.array(['pb_all']),
           }

example_fly = {'hDeltaC': 3,
               'hDeltaJ': 1,
               'hDeltaH': 3,
               'hDeltaF': 0,
               'hDeltaD': 0,
              # 'PFNpm': 0
               } # tune which fly index is shown as the example for each neuron type

example_cmap = {'hDeltaC': 'Blues',
                'hDeltaJ': 'Blues',
                'hDeltaH': 'Blues',
                'hDeltaF': 'Blues',
                'hDeltaD': 'Blues',
               # 'PFNpm': 'Greys'
                } # tune the colormap used for each neuron type

fig,axs = plt.subplots(len(dchoices),1,figsize=(8,1.5*len(dchoices)),sharex=True)
axs = np.atleast_1d(axs)
for i,d in enumerate(dchoices):
    dname = d+'_ACV_pulses'
    data = pulse_data[dname]
    f = example_fly[d]
    cxa = data[f].cxa
    wedges = cxa.pdat['wedges_'+regions_plot[d][0]]
    jumps = cxa.get_entries_exits_like_jumps()
    onset0 = jumps[0,0]-50
    edx = np.arange(onset0,min(onset0+nplot,len(wedges)))
    plotdat = wedges[edx,:].T
    t = (edx-onset0)/fs

    ax = axs[i]
    ax.imshow(plotdat,aspect='auto',cmap=example_cmap[d],
              extent=[t[0],t[-1],plotdat.shape[0],0],vmin=0,vmax=1)

    onsets = jumps[:,0]
    onsets = onsets[(onsets>=onset0) & (onsets<onset0+nplot)]
    offsets = jumps[:,1]
    offsets = offsets[(offsets>=onset0) & (offsets<onset0+nplot)]
    for o in onsets:
        ax.axvline((o-onset0)/fs,color='r',linewidth=2,linestyle='--')
    for o in offsets:
        ax.axvline((o-onset0)/fs,color='r',linewidth=2,linestyle='--')
    ax.set_ylabel(d)

axs[-1].set_xlabel('Time from first ACV onset (s)')
fig.suptitle('Example ACV wedge activity')
plt.tight_layout()
plt.pause(.1)
plt.savefig(os.path.join(savedir,'Example_ACV_wedge_heatmaps.pdf'))































