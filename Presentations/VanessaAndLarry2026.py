# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:45:44 2026

@author: dowel

Script plots figures for Vanessa and Larry's flatiron talk in June 2026

Plots:
    1) Comparison of in vs out of plume activity
    PFNp/m, hDeltaJ, hDeltaC, hDeltaH
    
    2) Ramping of hDeltaC - relation to FC2 during stops
    
    3) PFL3 model, and its prediction of turns


"""
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
    
    "hDeltaC_FC2_EPG": [
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260331\f1\Trial3', # 6 jumps
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260426\f1\Trial1', # 8 jumps
        #r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260425\f2\Trial3', # 3 jumps
        #r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260401\f1\Trial2', # Some nice tracking 3 jumps
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260429\f2\Trial1',# 17 jumps
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260616\f1\Trial2', # 6 jumps
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260623\f1\Trial2', # 6 jumps
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260707\f1\Trial3' # 7 jumps
        ],
    
    "PFN_pm":[ r'Y:\Data\FCI\Hedwig\PFNpm_SS00191\251210\f1\Trial1'],
    
    "hDeltaA":[r"Y:\Data\FCI\Hedwig\hDeltaA_SS64464\241121\f2\Trial1"],
    
    "hDeltaM":[r"Y:\Data\FCI\Hedwig\hDeltaM_SS59766\250829\f1\Trial2" ],
    
    "hDeltaK": [r"Y:\Data\FCI\Hedwig\hDeltaK_SS63089\250807\f2\Trial4"],
    
}
all_flies = {}

for group_name, datadirs in groups.items():
    all_flies[group_name] = {}
    
    for i, datadir in enumerate(datadirs):
        print(f"[{group_name}] {datadir}")
        if group_name=='hDeltaC_FC2':
            etp = ET_paper(datadir,regions=['fsb1_ch1','fsb2_ch2'])
        elif group_name=='hDeltaC_FC2_EPG':
            etp = ET_paper(datadir,regions=['eb_ch1','fsb1_ch1','fsb2_ch2'])
        elif group_name=='PFN_pm':
            etp = ET_paper(datadir,regions= ['pfn_l_16','pfn_r_16','pfn_l','pfn_r'] )
        else:
            etp = ET_paper(datadir)
        all_flies[group_name][i] = etp
#%% Bump amplitude comparison
savedir = r'Y:\Data\FCI\FCI_summaries\bump_analysis'
plt.close('all')

dchoice = ['hDeltaJ','hDeltaC','FC2']

pltall = np.zeros((16,2,2,len(dchoice),7))
for ids,d in enumerate(dchoice):
    tdata = all_flies[d]
    for t in tdata:
        etp = tdata[t]
        pltall[:,:,:,ids,t] = etp.cxa.phase_null_in_out_jump(fsb_names=['eb','fsb_upper'])

colours = uplt.columnar_colours()

pmean = np.mean(pltall,axis=4)

offset = 0
x= np.arange(0,16) 
for i1 in range(2):
    for i in range(len(dchoice)):
        plt.plot(x+offset,pmean[:,i1,1,i],color = colours[i,:],linewidth=2)
        plt.plot(x+offset,pltall[:,i1,1,i,:],color = colours[i,:],alpha=.25)
        offset = offset+17
        if i==2:
            plt.plot(x+offset,pltall[:,i1,0,i,:],color = 'k',alpha=.25)
            plt.plot(x+offset,pmean[:,i1,0,i],color = 'k',linewidth=2)
            offset = offset+17
    offset = offset+16
plt.xticks([40,120],labels= ['plume exit','plume return'])
plt.ylabel('dF/F0')
plt.savefig(os.path.join(savedir,'hDelta_FC2_bump.pdf'))

#%% Do single offset on phase  to avoid data twists
eg_fly= 4
etp = all_flies['hDeltaC_FC2_EPG'][eg_fly]
heading = etp.cxa.ft2['ft_heading'].to_numpy()
eb = etp.cxa.pdat['phase_eb_ch1']
fc2 = etp.cxa.pdat['phase_fsb1_ch1']
hdc = etp.cxa.pdat['phase_fsb2_ch2']

ebw = etp.cxa.pdat['wedges_eb_ch1']
fc2w = etp.cxa.pdat['wedges_fsb1_ch1']
hdcw = etp.cxa.pdat['wedges_fsb2_ch2']

offset = stats.circmean(ug.circ_subtract(eb,heading),low=-np.pi,high=np.pi)
etp.cxa.pdat['offset_eb_ch1_phase'] = pd.Series(ug.circ_subtract(eb,offset))
etp.cxa.pdat['offset_fsb1_ch1_phase'] = pd.Series(ug.circ_subtract(fc2,offset))
etp.cxa.pdat['offset_fsb2_ch2_phase'] = pd.Series(ug.circ_subtract(hdc,offset))
offset = np.zeros_like(eb)+offset
etp.cxa.pdat['wedges_offset_eb_ch1'] = etp.cxa.cx.rotate_wedges(ebw,phase_offset =offset)
etp.cxa.pdat['wedges_offset_fsb1_ch1'] = etp.cxa.cx.rotate_wedges(fc2w,phase_offset =offset)
etp.cxa.pdat['wedges_offset_fsb2_ch2'] = etp.cxa.cx.rotate_wedges(hdcw,phase_offset =offset)
#%% FC2 hDeltaC co-imaging arrows and heatmaps
savedir = r'Y:\Presentations\2026\05_ForVanessaAndLarry'
colours = uplt.columnar_colours()
colours[0,:] = 0
colours = np.append(colours,np.array([6,149,207])[np.newaxis,:]/255,axis=0)
colours[[1,2],:] = colours[[2,1],:]

plt.close('all')
#Full trajectory with boxes
# tstart = 0
# tend = 100000
# etp.example_trajectory_arrows(tstart,tend,colours,samp_period = 5,boxes=np.array([[60,110],[60*4.2,60*5.5]]),regions=['eb_ch1','fsb1_ch1','fsb2_ch2'])
# #plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Amenotaxis trajectory
plt.close('all')
tstart = 32
tend = 110
etp.example_trajectory_arrows(tstart,tend,colours,samp_period = 5,regions=['eb_ch1','fsb1_ch1','fsb2_ch2'])
plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Amenotaxis heats

etp.example_heatmaps(tstart,tend,regions = ['eb_ch1','fsb1_ch1','fsb2_ch2'],cmaps=['Greys','Blues','Greens'],colours= colours)
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'hDeltaC_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Phase_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(6)
plt.savefig(os.path.join(savedir,'Plume_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))


# Jump trajectory
plt.close('all')
tstart = 60*15.06
#tend = 60*7.5
tend = 60*16.06
etp.example_trajectory_arrows(tstart,tend,colours,samp_period = 5,regions=['eb_ch1','fsb1_ch1','fsb2_ch2'])
#plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Jump heats
etp.example_heatmaps(tstart,tend,regions = ['eb_ch1','fsb1_ch1','fsb2_ch2'],cmaps=['Greys','Blues','Greens'],colours= colours)
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'hDeltaC_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Phase_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(6)
plt.savefig(os.path.join(savedir,'Plume_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

#%% example trajectory at the end of the experiment
savedir = r'Y:\Presentations\2026\05_ForVanessaAndLarry'
colours = uplt.columnar_colours()
colours[0,:] = 0
colours = np.append(colours,np.array([6,149,207])[np.newaxis,:]/255,axis=0)
colours[[1,2],:] = colours[[2,1],:]

plt.close('all')
#Full trajectory with boxes
# tstart = 0
# tend = 100000
# etp.example_trajectory_arrows(tstart,tend,colours,samp_period = 5,boxes=np.array([[60,110],[60*4.2,60*5.5]]),regions=['eb_ch1','fsb1_ch1','fsb2_ch2'])
# #plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Amenotaxis trajectory
plt.close('all')
tstart = 60*22.78
tend = 60*23.2
etp.example_trajectory_arrows(tstart,tend,colours,scale=4,samp_period = 1,regions=['eb_ch1','fsb1_ch1','fsb2_ch2'])
plt.savefig(os.path.join(savedir,'EndOfPlumeSearch.pdf'))
#%% Figure 1
"""
Raster plot of entries and exits ranked by length and ordered by entry and exit

# This can be a colour plot

"""


dchoice = ['FC2','hDeltaJ','hDeltaC','hDeltaH','PFN_pm']

example = [2,      1,        1,        4,          0]
minlen = 4
padding = 100
plt.close('all')
alignment = ['entry','exit']
for i,d in enumerate(dchoice):
    etp = all_flies[d][example[i]]
    e_e = etp.cxa.get_entries_exits_like_jumps()
    e_len = e_e[:,1]-e_e[:,0]
    e_e = e_e[e_len>=minlen]
    e_len = e_len[e_len>=minlen]
    isort = np.argsort(e_len)
    
    e_e = e_e[isort,:]
    e_len = e_len[isort]
    plot_array = np.zeros((len(e_e),600,2))*np.nan
    # 3rd dim is for stim off/stim on ordering
    if d == 'PFN_pm':
        mnflur = np.mean(np.append(etp.cxa.pdat['wedges_pfn_l_16'],etp.cxa.pdat['wedges_pfn_r_16'],axis=1),axis=1)
    else:
        mnflur = np.mean(etp.cxa.pdat['wedges_fsb_upper'],axis=1)
    
    
    for ie,e in enumerate(e_e):
        dx = np.arange(e[0]-padding,e[1]+padding)
        tflur = mnflur[dx]
        plot_array[ie,:len(dx),0] = tflur #align to onset
        plot_array[ie,-len(dx):,1] = tflur#align to offset
                          
    plt.figure()
    for i2 in range(2):
        plt.subplot(1,2,i2+1)
        plt.imshow(plot_array[:,:,i2],interpolation='None',aspect='auto',cmap='Blues',vmin=0,vmax=1)
        plt.title(d+' '+ alignment[i2]+ ' aligned')
    plt.subplot(1,2,1)
    plt.plot([padding,padding],[-.5,plot_array.shape[0]-.5],color='r',linewidth=.5)
    plt.subplot(1,2,2)
    plt.plot([plot_array.shape[1]-padding,plot_array.shape[1]-padding],[-.5,plot_array.shape[0]-.5],color='r',linewidth=.5)

    for ie_e,e in enumerate(e_len):
        plt.subplot(1,2,1)
        plt.plot([padding+e,padding+e],[ie_e-.5,ie_e+.5],color='r',linewidth=.5)
        plt.subplot(1,2,2)
        plt.plot([plot_array.shape[1]-padding-e,plot_array.shape[1]-padding-e],[ie_e-.5,ie_e+.5],color='r',linewidth=.5)


#%% Figure data on mean fluorescence for 5 minutes. 1 before first odour presentation
"""
Plot is line plot, can do colour raster if better

"""
dchoice = ['FC2','PFN_pm','hDeltaA','hDeltaK','hDeltaM', 'hDeltaJ','hDeltaC','hDeltaH']

example = [4,      0,          0,       0,        0,        0,       1,        0  ]
minlen = 5

plt.close('all')
alignment = ['entry','exit']
yoffset = 0
maxlen = 200
for i,d in enumerate(dchoice):
    etp = all_flies[d][example[i]]
    e_e = etp.cxa.get_entries_exits_like_jumps()
    
   
    # 3rd dim is for stim off/stim on ordering
    if d == 'PFN_pm':
        mnflur = np.mean(np.append(etp.cxa.pdat['wedges_pfn_l_16'],etp.cxa.pdat['wedges_pfn_r_16'],axis=1),axis=1)
    else:
        mnflur = np.mean(etp.cxa.pdat['wedges_fsb_upper'],axis=1)

    dx = np.arange(e_e[0,0]-100,e_e[0,0]+300*1)
    ins = etp.cxa.ft2['instrip'].to_numpy()
    tins = ins[dx]
    tt = (dx-e_e[0,0])*.1 
    tmnflur = mnflur[dx]
    
    plt.plot(tt,tmnflur+yoffset,color='k')
    ttn = tt.copy()
    ttn[tins==0] = np.nan
    plt.plot(ttn,np.zeros(len(ttn))+yoffset-.1,color='r',linewidth=3)
    yoffset = yoffset-1.5
yt = np.linspace(0,yoffset+1.5,len(dchoice))
plt.yticks(yt,labels=dchoice)
plt.xlabel('time (s)')


# Same but stacked



minlen = 5

alignment = ['entry','exit']
yoffset = 0
colours = uplt.columnar_colours()
extra_colours = np.array([[93, 128, 75],[224, 43, 40]])/255
colours = np.append(colours,extra_colours,axis=0)
padding= 50

yoffset= 0
xoffset = max(tt)+3
tt = np.arange(0,600)/10
for i,d in enumerate(dchoice):
    etp = all_flies[d][example[i]]
    e_e = etp.cxa.get_entries_exits_like_jumps()
    e_len = e_e[:,1]-e_e[:,0]
    e_e = e_e[e_len>=minlen]
    e_len = e_len[e_len>=minlen]
    
    e_e = e_e[e_len<=maxlen]
    e_len = e_len[e_len<=maxlen]
   
    # 3rd dim is for stim off/stim on ordering
    if d == 'PFN_pm':
        mnflur = np.mean(np.append(etp.cxa.pdat['wedges_pfn_l_16'],etp.cxa.pdat['wedges_pfn_r_16'],axis=1),axis=1)
    else:
        mnflur = np.mean(etp.cxa.pdat['wedges_fsb_upper'],axis=1)
        
    plot_array = np.zeros((len(e_e),600,2))*np.nan
    # 3rd dim is for stim off/stim on ordering
    if d == 'PFN_pm':
        mnflur = np.mean(np.append(etp.cxa.pdat['wedges_pfn_l_16'],etp.cxa.pdat['wedges_pfn_r_16'],axis=1),axis=1)
    else:
        mnflur = np.mean(etp.cxa.pdat['wedges_fsb_upper'],axis=1)
    
    normaliser = np.median(mnflur[e_e[0,0]-200:e_e[0,0]])
    for ie,e in enumerate(e_e):
        dx = np.arange(e[0]-padding,e[1])
        tflur = mnflur[dx]
        plot_array[ie,:len(dx),0] = tflur #align to onset
        dx = np.arange(e[0],e[1]+padding)
        tflur = mnflur[dx]
        plot_array[ie,-len(dx):,1] = tflur#align to offset
    
    plot_array = plot_array-normaliser
    
    pmean = np.nanmean(plot_array,axis=0)
    print(i)
    plt.figure(1)
    yoffset2 = yoffset+normaliser
    plt.plot(tt[:padding]+xoffset,   pmean[30:padding+30,0]+yoffset2,color='k')
    plt.plot([xoffset,tt[padding]+xoffset],[yoffset2,yoffset2],color='k',linestyle='--')
    # plt.xlim([3,tt[padding]+3])
    plt.fill([tt[padding-30]+xoffset,tt[padding]+xoffset,tt[padding]+xoffset,tt[padding-30]+xoffset],np.array([-.4,-.4,.4,.4])+yoffset2,color=[.8,.8,.8])
    
    # plt.figure(2)
    plt.plot(tt[:padding+10]+xoffset+9,pmean[-padding-10:,1]+yoffset2,color='k')
    plt.plot([xoffset+9,tt[padding+10]+xoffset+9],[yoffset2,yoffset2],color='k',linestyle='--')
    plt.fill([xoffset+9,xoffset+9+1,xoffset+9+1,xoffset+9,],np.array([-.4,-.4,.4,.4,])+yoffset2,color=[.8,.8,.8])
    # plt.xlim([tt[-padding]-1,tt[-1]])
    yoffset= yoffset-1.5


#%%

#%% Figure data with stacked images on left, mean fluorescence on right
"""
Left panels: wedge images stacked vertically
Right panels: mean fluorescence entry/exit aligned (as before)
"""
dchoice = ['FC2','PFN_pm','hDeltaA','hDeltaK','hDeltaM', 'hDeltaJ','hDeltaC','hDeltaH']

example = [3,       0,          0,         0,        0,      0,       1,        4,        ]
minlen = 5
plt.close('all')
alignment = ['entry','exit']
yoffset = 0
maxlen = 200
image_height = 1.0   # vertical space each image occupies
row_spacing = 1.5    # vertical spacing between rows (matches original 1.5 offset)

fig, ax = plt.subplots()

for i,d in enumerate(dchoice):
    etp = all_flies[d][example[i]]
    e_e = etp.cxa.get_entries_exits_like_jumps()
    
    # Get 2D wedge data instead of mean
    if d == 'PFN_pm':
        # wedges = np.append(etp.cxa.pdat['wedges_pfn_l_16'],
        #                    etp.cxa.pdat['wedges_pfn_r_16'], axis=1)
        wedges = np.append(etp.cxa.pv2.filter(regex='pfn_l').to_numpy()[:,:8],etp.cxa.pv2.filter(regex='pfn_r').to_numpy()[:,:8],axis=1)
        
    else:
        wedges = etp.cxa.pdat['wedges_fsb_upper']
    
    dx = np.arange(e_e[0,0]-100, e_e[0,0]+300*1)
    ins = etp.cxa.ft2['instrip'].to_numpy()
    tins = ins[dx]
    tt = (dx-e_e[0,0])*.1
    twedges = wedges[dx, :]  # shape: (time, wedges)
    
    # Plot image: transpose so wedges are on y-axis, time on x-axis
    ax.imshow(twedges.T,
              aspect='auto',
              extent=[tt[0], tt[-1],
                      yoffset - image_height/2,
                      yoffset + image_height/2],
              origin='lower',
              cmap='Blues',vmin=0,vmax=1)
    
    # Keep the red instrip indicator line
    ttn = tt.copy()+.5#adding to account for MFC delay as with Vanessa's suggestion
    ttn[tins==0] = np.nan
    ax.plot(ttn, np.zeros(len(ttn))+yoffset-image_height/2-0.1,
            color='r', linewidth=3)
    
    yoffset = yoffset - row_spacing

yt = np.linspace(0, yoffset+row_spacing, len(dchoice))
ax.set_yticks(yt)
ax.set_yticklabels(dchoice)
ax.set_xlabel('time (s)')

# Right-hand panels: mean fluorescence entry/exit aligned (as in original)
minlen = 5
alignment = ['entry','exit']
yoffset = 0
colours = uplt.columnar_colours()
extra_colours = np.array([[93, 128, 75],[224, 43, 40]])/255
colours = np.append(colours, extra_colours, axis=0)
padding = 50
yoffset = 0
xoffset = max(tt)+3
tt = np.arange(0,600)/10

for i,d in enumerate(dchoice):
    etp = all_flies[d][example[i]]
    e_e = etp.cxa.get_entries_exits_like_jumps()
    e_len = e_e[:,1]-e_e[:,0]
    e_e = e_e[e_len>=minlen]
    e_len = e_len[e_len>=minlen]
    
    e_e = e_e[e_len<=maxlen]
    e_len = e_len[e_len<=maxlen]
   
    # 3rd dim is for stim off/stim on ordering
    if d == 'PFN_pm':
        mnflur = np.mean(np.append(etp.cxa.pdat['wedges_pfn_l_16'],
                                   etp.cxa.pdat['wedges_pfn_r_16'], axis=1), axis=1)
    else:
        mnflur = np.mean(etp.cxa.pdat['wedges_fsb_upper'], axis=1)
        
    plot_array = np.zeros((len(e_e),600,2))*np.nan
    
    normaliser = np.median(mnflur[e_e[0,0]-200:e_e[0,0]])
    for ie,e in enumerate(e_e):
        dx = np.arange(e[0]-padding, e[1])
        tflur = mnflur[dx]
        plot_array[ie,:len(dx),0] = tflur  # align to onset
        dx = np.arange(e[0], e[1]+padding)
        tflur = mnflur[dx]
        plot_array[ie,-len(dx):,1] = tflur  # align to offset
    
    plot_array = plot_array - normaliser
    
    pmean = np.nanmean(plot_array, axis=0)
    print(i)
    
    yoffset2 = yoffset #+ normaliser
    ax.plot(tt[:padding]+xoffset, pmean[30:padding+30,0]+yoffset2, color='k')
    ax.plot([xoffset, tt[padding]+xoffset], [yoffset2, yoffset2],
            color='k', linestyle='--')
    ax.fill([tt[padding-30]+xoffset, tt[padding]+xoffset,
             tt[padding]+xoffset, tt[padding-30]+xoffset],
            np.array([-.4,-.4,.4,.4])+yoffset2, color=[.8,.8,.8])
    
    ax.plot(tt[:padding+10]+xoffset+9, pmean[-padding-10:,1]+yoffset2, color='k')
    ax.plot([xoffset+9, tt[padding+10]+xoffset+9], [yoffset2, yoffset2],
            color='k', linestyle='--')
    ax.fill([xoffset+9, xoffset+9+1, xoffset+9+1, xoffset+9],
            np.array([-.4,-.4,.4,.4])+yoffset2, color=[.8,.8,.8])
    
    yoffset = yoffset - 1.5
plt.plot([49,49],[-5.5,-6],color='k')
plt.text(49,-5.75,'0.5 dF/F0')
plt.savefig(os.path.join(savedir,'Heat_Deltas_FC2_Example.pdf'))

#%% Average plot of odour onset/offset with pre odour baseline


dchoice = ['hDeltaJ','hDeltaC','FC2','hDeltaH','PFN_pm']

example = [0,      1,        2,        4,          0]
minlen = 5

plt.close('all')
alignment = ['entry','exit']
yoffset = 0
colours = uplt.columnar_colours()
extra_colours = np.array([[93, 128, 75],[224, 43, 40]])/255
colours = np.append(colours,extra_colours,axis=0)
padding= 50
tt = np.arange(0,600)/10
for i,d in enumerate(dchoice):
    etp = all_flies[d][example[i]]
    e_e = etp.cxa.get_entries_exits_like_jumps()
    e_len = e_e[:,1]-e_e[:,0]
    e_e = e_e[e_len>=minlen]
    e_len = e_len[e_len>=minlen]
    
    e_e = e_e[e_len<=maxlen]
    e_len = e_len[e_len<=maxlen]
   
    # 3rd dim is for stim off/stim on ordering
    if d == 'PFN_pm':
        mnflur = np.mean(np.append(etp.cxa.pdat['wedges_pfn_l_16'],etp.cxa.pdat['wedges_pfn_r_16'],axis=1),axis=1)
    else:
        mnflur = np.mean(etp.cxa.pdat['wedges_fsb_upper'],axis=1)
        
    plot_array = np.zeros((len(e_e),600,2))*np.nan
    # 3rd dim is for stim off/stim on ordering
    if d == 'PFN_pm':
        mnflur = np.mean(np.append(etp.cxa.pdat['wedges_pfn_l_16'],etp.cxa.pdat['wedges_pfn_r_16'],axis=1),axis=1)
    else:
        mnflur = np.mean(etp.cxa.pdat['wedges_fsb_upper'],axis=1)
    
    normaliser = np.median(mnflur[e_e[0,0]-200:e_e[0,0]])
    for ie,e in enumerate(e_e):
        dx = np.arange(e[0]-padding,e[1])
        tflur = mnflur[dx]
        plot_array[ie,:len(dx),0] = tflur #align to onset
        dx = np.arange(e[0],e[1]+padding)
        tflur = mnflur[dx]
        plot_array[ie,-len(dx):,1] = tflur#align to offset
    
    plot_array = plot_array-normaliser
    
    pmean = np.nanmean(plot_array,axis=0)
    print(i)
    plt.figure(1)
    plt.plot(tt,pmean[:,0],color=colours[i,:])
    plt.figure(2)
    plt.plot(tt,pmean[:,1],color=colours[i,:])
    
plt.figure(1)
plt.title()
plt.plot([0,tt[-1]],[0,0],color='k',linestyle='--')
plt.xlim([3,tt[padding]+3])
plt.fill([tt[padding],tt[padding]+3,tt[padding]+3,tt[padding]],[-.4,-.4,.4,.4],color=[.8,.8,.8])

plt.figure(2)
plt.plot([0,tt[-1]],[0,0],color='k',linestyle='--')
plt.fill([tt[-padding],tt[-padding],tt[-padding]-5,tt[-padding]-5],[-.4,.4,.4,-.4],color=[.8,.8,.8])
plt.xlim([tt[-padding]-1,tt[-1]])
#%% hDeltaC and FC2 correlation on trajectory

dchoice = 'hDeltaC_FC2_EPG'
tdata = all_flies[dchoice]
colours = uplt.defined_cmap('coolwarm', 5)
regions = ['eb_ch1','fsb2_ch2']
titles = ['EPG-FC2','hDeltaC-FC2']
for ir, reg in enumerate(regions):
    plt.figure()
    plt.title(titles[ir])
    for fly in tdata:
        etp = tdata[fly]
        jumps = etp.cxa.get_jumps()
        print(fly,len(jumps))
        
        whdc = etp.cxa.pdat['wedges_'+reg]
        wfc2 = etp.cxa.pdat['wedges_fsb1_ch1']
        r = ug.rowwise_pearson(whdc,wfc2)
        jdx = np.zeros(len(r))
        for j in jumps:
            jdx[j[1]:j[2]] = 1
        jdx = jdx>0
        #etp.cxa.plot_traj_arrow_heat(['fsb1_ch1'],r,cmin=-1,cmax=1,a_sep=5)
        offsets = np.arange(1,100)
        for o in offsets:
            a = whdc[o:,:]
            b = wfc2[:-o,:]
            r = ug.rowwise_pearson(a,b)
            #r = ug.rowwise_cov(a,b)
            plt.scatter(-o*.1,np.mean(r),color=colours[fly,:])
            #plt.scatter(-o*.1,np.mean(r[jdx[:-o]]),color=colours[fly,:])
        
        offsets = np.arange(1,100)
        for o in offsets:
            a = whdc[:-o,:]
            b = wfc2[o:,:]
            r = ug.rowwise_pearson(a,b)
            #r = ug.rowwise_cov(a,b)
            plt.scatter(o*.1,np.mean(r),color=colours[fly,:])
            #plt.scatter(o*.1,np.mean(r[jdx[o:]]),color=colours[fly,:])
        #plt.text(o*.1,np.mean(r),'fly: ' +str(fly) + ' ' + 'jumps: '+ str(len(jumps)) )
    plt.xlabel('offset in time -ve FC2 leads +ve FC2 lags')
    plt.ylabel('Mean pearson rho')
#%% plotting R variation
plt.close('all')
import scipy.signal as sg
lag = 25
fly = 4
tdata = all_flies[dchoice]
etp = tdata[fly]
whdc = etp.cxa.pdat['wedges_fsb2_ch2']
wfc2 = etp.cxa.pdat['wedges_fsb1_ch1']
cdiff = np.mean((wfc2[lag:]-whdc[:-lag])**2,axis=1)
shdc = np.std(whdc[:-lag],axis=1)
sfc2 = np.std(whdc[lag:],axis=1)
#ccc = (2*r*shdc*sfc2)/(shdc**2+sfc2**2+(np.mean(whdc[lag:]))**2)
r = ug.rowwise_pearson(whdc[:-lag],wfc2[lag:])
ro = ug.rowwise_pearson(whdc,wfc2)
cv = ug.rowwise_cov(whdc[:-lag],wfc2[lag:])
rk = ug.kalman_filter_1d(r, process_var=1e-3, measurement_var=2e-1)
rk2 = ug.kalman_filter_1d(np.flipud(r), process_var=1e-4, measurement_var=1e-1)
rko = (rk+np.flipud(rk2))/2

ins = etp.cxa.ft2['instrip'][lag:].to_numpy()
plt.plot(r,color='b')

#plt.plot(ro,color='k')

rsmth = sg.savgol_filter(r,20,3)
rpeaks = rsmth.copy()
rpeaks[rpeaks<0] = 0
pks,pkmeta = sg.find_peaks(rpeaks,height=.5,width=10,)
plt.plot(rsmth,color=[0.2,0.8,1])
plt.plot(rko,color='k')
plt.plot(cv*10,color='g')
#plt.plot(ccc,color='m')
plt.scatter(pks,rsmth[pks],color='g',zorder=5)

plt.plot(ins,color='r')

cxa = etp.cxa
x = cxa.ft2['ft_posx'].to_numpy()
y = cxa.ft2['ft_posy'].to_numpy()
x = x[lag:]
y = y[lag:]
plt.figure()
plt.plot(x,y,color='k')

plt.gca().set_aspect('equal')
plt.scatter(x[ins>0],y[ins>0],color=[.8,.8,.8])
plt.scatter(x[pks],y[pks],color='r')
fc2p = cxa.pdat['offset_fsb1_ch1_phase'].to_numpy()[lag:]
hdcp = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()[:-lag]

pfc2 = fc2p[pks]
pfc22 = fc2p[pks-lag]
phdc = hdcp[pks]
xpk = x[pks]
ypk = y[pks]

xpk2 = x[pks-lag]
ypk2 = y[pks-lag]
scale = 20
colours = uplt.columnar_colours()
for i,p in enumerate(pfc2):
    x = [xpk[i],xpk[i]+np.sin(pfc2[i])*scale]
    y = [ypk[i],ypk[i]+np.cos(pfc2[i])*scale]
    plt.plot(x,y,color=colours[2,:])
    
    x = [xpk[i],xpk[i]+np.sin(phdc[i])*scale]
    y = [ypk[i],ypk[i]+np.cos(phdc[i])*scale]
    plt.plot(x,y,color=colours[1,:])

    x = [xpk2[i],xpk2[i]+np.sin(pfc22[i])*scale]
    y = [ypk2[i],ypk2[i]+np.cos(pfc22[i])*scale]
    plt.plot(x,y,color=colours[2,:])
    

etp.cxa.plot_traj_arrow_heat(['fsb1_ch1'],cv*10,cmin=-.2,cmax=.2,a_sep=5)
plt.figure()
plt.plot(rko,color='k')
plt.plot(ins-1,color='r')
#%% Plot covariance variation
lag = 25
fly = 4
tdata = all_flies[dchoice]
etp = tdata[fly]
whdc = etp.cxa.pdat['wedges_fsb2_ch2']
wfc2 = etp.cxa.pdat['wedges_fsb1_ch1']
cdiff = np.mean((wfc2[lag:]-whdc[:-lag])**2,axis=1)
shdc = np.std(whdc[:-lag],axis=1)
sfc2 = np.std(whdc[lag:],axis=1)
ins = etp.cxa.ft2['instrip'][lag:].to_numpy()

#ccc = (2*r*shdc*sfc2)/(shdc**2+sfc2**2+(np.mean(whdc[lag:]))**2)
cv = ug.rowwise_cov(whdc[:-lag],wfc2[lag:])
cvo = ug.rowwise_cov(whdc,wfc2)

plt.plot(cv*20,color='b')
plt.plot(cvo*20,color='k')
plt.plot(ins,color='r')

cv = np.append(np.zeros(lag),cv)
etp.cxa.plot_traj_arrow_heat(['fsb1_ch1','fsb2_ch2'],cv*10,cmin=-.2,cmax=.2,a_sep=5)


#%% hDeltaC ramps
plt.close('all')

d = 'hDeltaC'
tdata = all_flies[d]
padding = 100
minsize = 10 
maxsize = 100
for t in tdata:
    etp = tdata[t]
    u = ug()
    
    
    cxa = etp.cxa
    dx,dy,dd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'].to_numpy())

    allstps = dd<.5
    bstart,bsize = ug.find_blocks(allstps)
    bstart = bstart[:-1]
    bsize = bsize[:-1]
    
    bstart = bstart[np.logical_and(bsize>=minsize,bsize<=maxsize)]
    bsize = bsize[np.logical_and(bsize>=minsize,bsize<=maxsize)]
    
    
    
    brank = np.argsort(bsize)
    bstart = bstart[brank]
    bsize = bsize[brank]
    if d=='hDeltaC_FC2_EPG':
        mnflur = np.mean(cxa.pdat['wedges_fsb2_ch2'],axis=1)
        pva = ug.get_pvas(cxa.pdat['wedges_fsb2_ch2'])
    else:
        mnflur = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        pva = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    
    plotdata = np.zeros((len(bstart),np.max(bsize)+padding+10,2))*np.nan
    
    for i,b in enumerate(bstart):
        
        dx = np.arange(b+5,b+bsize[i]+padding)
        if dx[-1]>=len(pva):
            continue
            #plotdata = plotdata[:-1,:,:]
            
        tflur = mnflur[dx]
        plotdata[i,-len(dx):,0] = tflur
        
        tflur = pva[dx]
        plotdata[i,-len(dx):,1] = tflur
        
    pmn = np.median(plotdata[:,-padding:,:],axis=1)
    pt = np.tile(pmn[:,np.newaxis,:],(1,plotdata.shape[1],1))
    plotdata = plotdata-np.tile(pmn[:,np.newaxis,:],(1,plotdata.shape[1],1))
    plt.figure()
    plt.subplot(2,2,1)
    plt.imshow(plotdata[:,:,0],aspect='auto',interpolation='none',vmin=-.2,vmax=.2,cmap='coolwarm')
    plt.plot([plotdata.shape[1]-padding,plotdata.shape[1]-padding],[-.5,plotdata.shape[0]-.5],color='k')
    
    plt.subplot(2,2,2)
    plt.imshow(plotdata[:,:,1],aspect='auto',interpolation='none',vmin=-.02,vmax=.02,cmap='coolwarm')
    plt.plot([plotdata.shape[1]-padding,plotdata.shape[1]-padding],[-.5,plotdata.shape[0]-.5],color='k')
    
    plt.subplot(2,2,3)
    plotmean = np.nanmedian(plotdata,axis=0)
    plt.plot(plotmean[:,0])
    plt.plot([plotdata.shape[1]-padding,plotdata.shape[1]-padding],[-.2,.2],color='k')
    
    plt.subplot(2,2,4)
    plt.plot(plotmean[:,1])
    plt.plot([plotdata.shape[1]-padding,plotdata.shape[1]-padding],[-.02,.02],color='k')
    
    
#%% Stop start PCA

plt.close('all')
from sklearn.decomposition import PCA
d = 'hDeltaC_FC2_EPG'
tdata = all_flies[d]
padding = 1
minsize = 50 
maxsize = 100
scrall = np.zeros((minsize-5+1,len(tdata)))
total_rho = np.array([])
for t in tdata:
    etp = tdata[t]
    u = ug()
    
    
    cxa = etp.cxa
    dx,dy,dd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'].to_numpy())

    allstps = dd<.5
    bstart,bsize = ug.find_blocks(allstps)
    bstart = bstart[:-1]
    bsize = bsize[:-1]
    
    bstart = bstart[np.logical_and(bsize>=minsize,bsize<=maxsize)]
    bsize = bsize[np.logical_and(bsize>=minsize,bsize<=maxsize)]
    
    
    
    brank = np.argsort(bsize)
    bstart = bstart[brank]
    bsize = bsize[brank]
    if d=='hDeltaC_FC2_EPG':
        mnflur = np.mean(cxa.pdat['wedges_fsb2_ch2'],axis=1)
        pva = ug.get_pvas(cxa.pdat['wedges_fsb2_ch2'])
    else:
        mnflur = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        pva = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    
    plotdata = np.zeros((len(bstart),np.max(bsize)+padding+10,2))*np.nan
    
    for i,b in enumerate(bstart):
        
        dx = np.arange(b+5,b+bsize[i]+padding)
        if dx[-1]>=len(pva):
            continue
            #plotdata = plotdata[:-1,:,:]
            
        tflur = mnflur[dx]
        plotdata[i,-len(dx):,0] = tflur
        
        tflur = pva[dx]
        plotdata[i,-len(dx):,1] = tflur
        
    pmn = np.median(plotdata[:,-padding:,:],axis=1)
    pt = np.tile(pmn[:,np.newaxis,:],(1,plotdata.shape[1],1))
    plotdata = plotdata-np.tile(pmn[:,np.newaxis,:],(1,plotdata.shape[1],1))
    plt.figure(t)
    plt.subplot(2,2,1)
    plt.imshow(plotdata[:,:,0],aspect='auto',interpolation='none',vmin=-.2,vmax=.2,cmap='coolwarm')
    plt.plot([plotdata.shape[1]-padding,plotdata.shape[1]-padding],[-.5,plotdata.shape[0]-.5],color='k')
    
    plt.subplot(2,2,2)
    plt.imshow(plotdata[:,:,1],aspect='auto',interpolation='none',vmin=-.02,vmax=.02,cmap='coolwarm')
    plt.plot([plotdata.shape[1]-padding,plotdata.shape[1]-padding],[-.5,plotdata.shape[0]-.5],color='k')
    
    plt.subplot(2,2,3)
    plotmean = np.nanmedian(plotdata,axis=0)
    plt.plot(plotmean[:,0])
    plt.plot([plotdata.shape[1]-padding,plotdata.shape[1]-padding],[-.2,.2],color='k')
    
    plt.subplot(2,2,4)
    plt.plot(plotmean[:,1])
    plt.plot([plotdata.shape[1]-padding,plotdata.shape[1]-padding],[-.02,.02],color='k')

    
    pca_matrix = plotdata[:,-padding-minsize+5:,]
    
    ndx = ~np.isnan(pca_matrix[:,10,0])
    pca_matrix = pca_matrix[ndx,:,:]
    pcm = pca_matrix[:,:,0].T
    pcm = (pcm-np.mean(pcm,axis=0))/np.std(pcm,axis=0)
    pca = PCA()
    pca.fit(pcm)
    
    expl_var = np.cumsum(pca.explained_variance_ratio_)
    scr = np.matmul(pcm,pca.components_[0,:].T)
    scrall[:,t] = scr
    
    rho = np.matmul((scr[np.newaxis,:]-np.mean(scr)),(pcm-np.mean(pcm,axis=0)[np.newaxis,:]))/((np.std(scr)*np.std(pcm,axis=0))*(minsize-5+1-1))
    plt.figure(101)
    plt.plot(scr,color='k')
    total_rho = np.append(total_rho,pca.components_[0,:]*np.sqrt(pca.explained_variance_[0]))
    
    #plt.scatter(rho[0],pca.components_[0,:]*np.sqrt(pca.explained_variance_[0]))
    # scr2 = np.matmul(pcm,pca.components_[1,:].T)
    # plt.plot(scr2)

    # scr = np.matmul(pcm,pca.components_[2,:].T)
    # plt.plot(scr)
    


plt.plot(np.mean(scrall,axis=1),color='r',linewidth=3)


#%% Plot offset between hDeltaC and FC2 phases during jump returns - dual
data = all_flies['hDeltaC_FC2_EPG']
colours = uplt.columnar_colours()
bins= 50
t = np.arange(0,bins*2)
t = np.append(np.linspace(0,1,bins),np.linspace(1,5,bins))
grandmean = np.zeros((len(t),len(data),3))
plt.close('all')
plt.fill_between([0,1],[-np.pi,-np.pi],[np.pi,np.pi],color=[.8,.8,.8])
for d in data:
    etp = data[d]
    phases,traj,amp= etp.trajectory_mean(bins = bins,regions=['eb_ch1','fsb1_ch1','fsb2_ch2'],time_threshold=60,time_threshold2=0,min_inplume=5)
    diff = np.angle(np.exp(1j*phases[:,1,:] - 1j*phases[:,2,:]))
    dwind = np.angle(np.exp(1j*phases+1j*np.pi))
    
    pmean = stats.circmean(diff,axis=1,high=np.pi,low=-np.pi)
    dmean = stats.circmean(dwind,axis=2,high=np.pi,low=-np.pi)
    grandmean[:,d,:] = stats.circmean(phases,axis=2,high=np.pi,low=-np.pi)
    #plt.scatter(t,dmean[:,1],color=colours[2,:],s=5,alpha=.5)
    #plt.scatter(t,dmean[:,2],color=colours[1,:],s=5,alpha=.5)
    plt.scatter(t,grandmean[:,d,0],color='k',s=15,alpha=.5,linewidth=0)
    plt.scatter(t,grandmean[:,d,1],color=colours[2,:],s=15,alpha=.5,linewidth=0)
    plt.scatter(t,grandmean[:,d,2],color=colours[1,:],s=15,alpha=.5,linewidth=0)
    #plt.plot(pmean,color='k')
pltmean = stats.circmean(grandmean,axis=1,low=-np.pi,high=np.pi)
plt.scatter(t,pltmean[:,0],color='k',s=50)
plt.scatter(t,pltmean[:,1],color=colours[2,:],s=50)
plt.scatter(t,pltmean[:,2],color=colours[1,:],s=50)
plt.plot([1,5],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
plt.plot([0,5],[0,0],color='k',linestyle='--')
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,90,0,90,180])
plt.ylabel('phase (deg)')
plt.xlabel('pseudo time')
plt.ylim([-np.pi,np.pi])
#plt.xlim([t[0],t[-1]])
savedir = r'Y:\Data\FCI\FCI_summaries\FC2_hDeltaC_comparison'
plt.savefig(os.path.join(savedir,'pseudotime_phase_hDC_FC2_Dual.pdf'))
plt.savefig(os.path.join(savedir,'pseudotime_phase_hDC_FC2_Dual.png'))


#%% Plot offset between hDeltaC and FC2 phases during jump returns -single
data = all_flies['hDeltaC_FC2_EPG']
colours = uplt.columnar_colours()
bins= 50
t = np.arange(0,bins*2)
t = np.append(np.linspace(0,1,bins),np.linspace(1,5,bins))
grandmean = np.zeros((len(t),7,2,2))
plt.close('all')
plt.fill_between([0,1],[-np.pi,-np.pi],[np.pi,np.pi],color=[.8,.8,.8])

dchoice = ['hDeltaC','FC2']
for i,di in enumerate(dchoice):
    data = all_flies[di]
    for d in data:
        etp = data[d]
        phases,traj,amp= etp.trajectory_mean(bins = bins,regions=['eb','fsb_upper'],time_threshold=60,time_threshold2=0,min_inplume=5)
        # diff = np.angle(np.exp(1j*phases[:,1,:] - 1j*phases[:,2,:]))
        # dwind = np.angle(np.exp(1j*phases+1j*np.pi))
        
        # pmean = stats.circmean(diff,axis=1,high=np.pi,low=-np.pi)
        # dmean = stats.circmean(dwind,axis=2,high=np.pi,low=-np.pi)
        grandmean[:,d,:,i] = stats.circmean(phases,axis=2,high=np.pi,low=-np.pi)
        #plt.scatter(t,dmean[:,1],color=colours[2,:],s=5,alpha=.5)
        #plt.scatter(t,dmean[:,2],color=colours[1,:],s=5,alpha=.5)
        plt.scatter(t,grandmean[:,d,0,i],color='k',s=15,alpha=.5,linewidth=0)
        plt.scatter(t,grandmean[:,d,1,i],color=colours[i+1,:],s=15,alpha=.5,linewidth=0)
       
        #plt.plot(pmean,color='k')
pltmean = stats.circmean(grandmean,axis=1,low=-np.pi,high=np.pi)

epg = grandmean[:,:,0,:]
epg = epg.reshape([len(t),7*2])
epgm = stats.circmean(epg,axis=1,low=-np.pi,high=np.pi)





plt.scatter(t,epgm,color='k',s=50)
plt.scatter(t,pltmean[:,1,1],color=colours[2,:],s=50)
plt.scatter(t,pltmean[:,1,0],color=colours[1,:],s=50)
plt.plot([1,5],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
plt.plot([0,5],[0,0],color='k',linestyle='--')
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,90,0,90,180])
plt.ylabel('phase (deg)')
plt.xlabel('pseudo time')
plt.ylim([-np.pi,np.pi])

savedir = r'Y:\Data\FCI\FCI_summaries\FC2_hDeltaC_comparison'
plt.savefig(os.path.join(savedir,'pseudotime_phase_hDC_FC2_Individual.pdf'))
plt.savefig(os.path.join(savedir,'pseudotime_phase_hDC_FC2_Individual.png'))


















































