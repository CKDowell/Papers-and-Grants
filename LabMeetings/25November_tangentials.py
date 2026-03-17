# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 12:24:10 2025

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
from analysis_funs.CX_analysis_tan import CX_tan

plt.rcParams['pdf.fonttype'] = 42 
import matplotlib as mpl
mpl.rcParams['path.simplify'] = True
plt.rcParams['path.simplify'] = True
#%% Load data
groups = {
    "FB5I": [
        r"Y:\Data\FCI\Hedwig\FB5I_SS100553\240628\f1\Trial2",
        r"Y:\Data\FCI\Hedwig\FB5I_SS100553\240917\f1\Trial2",
        r"Y:\Data\FCI\Hedwig\FB5I_SS100553\240917\f3\Trial3",
        r"Y:\Data\FCI\Hedwig\FB5I_SS100553_sytGC7f\251104\f1\Trial2",
    ],
    "FB6H": [
        r"Y:\Data\FCI\Hedwig\FB6H_SS95649\250318\f1\Trial2",
        r"Y:\Data\FCI\Hedwig\FB6H_SS95649\250321\f1\Trial2",
        r"Y:\Data\FCI\Hedwig\FB6H_SS95649\250328\f2\Trial3",
        r"Y:\Data\FCI\Hedwig\FB6H_SS95649\250507\f2\Trial2",
        r"Y:\Data\FCI\Hedwig\FB6H_SS95649\250729\f1\Trial4",
    ],
    "FB4R": [
        r"Y:\Data\FCI\Hedwig\SS61646_FB4R\240828\f3\Trial1",
        r"Y:\Data\FCI\Hedwig\SS61646_FB4R\240910\f1\Trial1",
        r"Y:\Data\FCI\Hedwig\SS61645_FB4R\240911\f1\Trial3",
    ],
    "FB5AB": [
        r"Y:\Data\FCI\Hedwig\FB5AB_SS53640\241205\f2\Trial3",
        
        r'Y:\Data\FCI\AndyData\21D07\20220103_21D07_sytjGCaMP7f_Fly1_001\processed',
        r'Y:\Data\FCI\AndyData\21D07\20220107_21D07_sytjGCaMP7f_Fly1_001\processed',
        r'Y:\Data\FCI\AndyData\21D07\20220110_21D07_sytjGCaMP7f_Fly1_001\processed',
        r'Y:\Data\FCI\AndyData\21D07\20220111_21D07_sytjGCaMP7f_Fly1_001\processed',
    ],
    "FB6A": [
        r"Y:\Data\FCI\Hedwig\FB6A_SS95731\250821\f2\Trial2",
        r"Y:\Data\FCI\Hedwig\FB6A_SS95731\250813\f1\Trial3",
    ],
    "FB4X": [
        r"Y:\Data\FCI\Hedwig\SS70711_FB4X\241030\f3\Trial3",
        r"Y:\Data\FCI\Hedwig\SS70711_FB4X\241031\f1\Trial3",
    ],
    "FB4P_b": [
        r"Y:\Data\FCI\Hedwig\FB4P_b_SS60296\240912\f2\Trial3",
        r"Y:\Data\FCI\Hedwig\FB4P_b_SS60296_sytGC7f\240809\f2\Trial2",
        r"Y:\Data\FCI\Hedwig\FB4P_b_SS60296_sytGC7f\240806\f1\Trial5",
        r"Y:\Data\FCI\Hedwig\FB4P_b_SS60296_sytGC7f\240805\f1\Trial3",
        r"Y:\Data\FCI\Hedwig\FB4P_b_SS60296_sytGC7f\240906\f2\Trial2",
        r"Y:\Data\FCI\Hedwig\FB4P_b_SS60296\250311\f1\Trial1",
    ],
}
all_flies = {}

for group_name, datadirs in groups.items():
    all_flies[group_name] = {}
    
    for i, datadir in enumerate(datadirs):
        print(f"[{group_name}] {datadir}")
        try:
            cxt = CX_tan(datadir)
        except:
            cxt = CX_tan(datadir,Andy=True)
        all_flies[group_name][i] = cxt

#%%
plt.close('all')
super_group = {'group_0':['FB5AB','FB4R'],
               'group_1':['FB5I','FB4P_b'],
               'group_2':['FB6H','FB6A'],
               'group_3':['FB4X']
               }


savedir = r'Y:\Presentations\2025\11_LabMeeting\tangentials'
minlen=50
x = np.arange(0,600)/10
cmap = plt.get_cmap('Accent')
colours = np.array( [[109, 212, 117],[217    , 137, 169]])/255
for si,s in enumerate(super_group):
    group_names = super_group['group_'+str(si)]
    plt.figure()
    
    for gi,g in enumerate(group_names):
        
        colour = colours[gi,:]
        this_group = all_flies[g]
        glen= len(this_group)
        print(g,glen)
        meandat = np.zeros((600,glen))
        for i in this_group:
            cxt = this_group[i]
            jumps = cxt.get_entries_exits_like_jumps()
            amp = cxt.ca
            jumplen = jumps[:,2]-jumps[:,1]
            jumps = jumps[jumplen>minlen]
            all_js = np.zeros((600,len(jumps)))
            for i2,j in enumerate(jumps):
                dx = np.arange(j[1],j[2])
               # plt.plot(amp[dx],color='k',alpha=0.3)
                tamp = amp[dx].flatten()
                #tamp = tamp-np.mean(tamp[:5])
                amplen  = np.min([600,len(tamp)])
                all_js[:amplen,i] = tamp[:amplen]
                
            all_js[all_js==0] = np.nan  
            meandat[:,i] = np.nanmean(all_js,axis=1)
            
            
        pmean = np.nanmean(meandat,axis=1)
        
        
        plt.plot(x,meandat,color=colour,alpha=0.3)
        plt.plot(x,pmean,color=colour,alpha=1,linewidth=3,label=g)
            
    plt.xlabel("time from plume exit (s)",fontsize=15)
    plt.ylabel("dF/F0",fontsize=15)
    plt.legend()
    plt.tight_layout()
    plt.show()
    plt.xlim([0,20])
    plt.ylim([-0.1,1])
    plt.xticks(np.arange(0,25,5),fontsize=15) 
    plt.yticks(np.arange(0,1.2,.2),fontsize=15)
    savename= 'AverageExitFluorescence'+"_".join(group_names) +'.png'
    plt.savefig(os.path.join(savedir,savename))