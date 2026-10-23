# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 11:10:12 2026

@author: dowel
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
               r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f2\Trial1"], #Nice dataset. Some ET bouts followed by plume traversal]
    
    
    # "hDeltaC_FC2": [r'Y:\Data\FCI\Hedwig\hDeltaC_68A10_FC2_GCaMP_RCaMP\260131\f1\Trial2',
    #     r'Y:\Data\FCI\Hedwig\hDeltaC_68A10_FC2_GCaMP_RCaMP\260201\f1\Trial3',
    #     r'Y:\Data\FCI\Hedwig\hDeltaC_68A10_FC2_GCaMP_RCaMP\260211\f1\Trial2'],
    
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
#%% Average return duration
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2',"hDeltaC_FC2_EPG"]
all_jumps = np.array([])
n = 0
for d in dchoices:
    data = all_flies[d]
    for f in data:
        cxa = data[f].cxa
        jumps = cxa.get_jumps()
        jlen = (jumps[:,2]-jumps[:,1])/10
        all_jumps = np.append(all_jumps,jlen)
        n+=1
        
plt.hist(all_jumps,np.arange(0,60,1),color='k',density=True)
medjump = np.median(all_jumps)
iqr = [np.percentile(all_jumps,75),np.percentile(all_jumps,25)]
modjump = stats.mode(np.round(all_jumps))[0]
meanjump = np.mean(all_jumps)

plt.plot([medjump,medjump],[0,.1],color='r')
plt.plot([modjump,modjump],[0,.1],color='b')
plt.plot([meanjump,meanjump],[0,.1],color='g')
plt.text(40,.1,'Median: '+str(np.round(medjump,decimals=1))+ ' s',color='r',fontsize=15)
plt.text(40,.09,'Mode: '+str(np.round(modjump,decimals=1))+' s',color='b',fontsize=15)
plt.text(40,.08,'Mean: '+str(np.round(meanjump,decimals=1))+ ' s',color='g',fontsize=15)

plt.xlabel('Jump return time (s)',fontsize=15)
plt.ylabel('Probability',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim([0.1,60])
savedir = 'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'
plt.savefig(os.path.join(savedir,'jump_return_time.png'))

#%% J epg ofset over time
d = 'hDeltaJ'
data = all_flies[d]
colours = uplt.columnar_colours()
plt.close('all')
for f in data:
    cxa = data[f].cxa
    hdb = cxa.pdat['phase_eb']
    hdj = cxa.pdat['phase_fsb_upper']
    ins = cxa.ft2['instrip'].to_numpy()
    t = np.arange(0,len(ins))/10
    # plt.plot(t,ins,color='r')
    # plt.scatter(t,hdb,color='k',s=3)
    # plt.scatter(t,hdj,color=colours[0,:],s=3)
    u =  ug()
    dx,dy,dd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),t)
    plt.figure()
    da = 1+ug.circ_subtract(hdj,hdb)/np.pi
    #plt.plot(t,da,color='k')
    plt.scatter(t,da,color='k',s=5)
    plt.plot(t,ins*2,color='r')
    plt.plot(t[1:],-2+dd/10,color=[.5,.5,.5])
    jumps = cxa.get_entries_exits_like_jumps()
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    pva = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    plt.plot(t,amp*2-1,color=colours[0,:])
    for ij,j in enumerate(jumps):
        plt.text(t[j[0]],np.pi+.3,ij)
    
    plt.plot([0,t[-1]],[1,1],color='g')
    cxa.get_jumps()
    plt.title(cxa.side)
    plt.plot([0,t[-1]],[1+cxa.side/2,1+cxa.side/2],color='g',linestyle='--')
    plt.plot(t,2*pva+2.1,color = colours[0,:])
    
    cxa.plot_traj_arrow_new(regions=['eb','fsb_upper'],a_sep=5)
    
    

#%% post odoour histogram
pbins = [0,4,8,16,32,60]

savedir = r'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison\histograms2D'
plt.close('all')
dchoice = 'hDeltaC_FC2_EPG'
dsecond = 'hDeltaC'
data = all_flies[dchoice]
# Scatter pre
epgall = np.array([])
hdcall = np.array([])
colours = uplt.hDelta_colours()

for f in data:
    cxa = data[f].cxa
    cxa.get_jumps()
    j = cxa.get_entries_exits_like_jumps()
    if dchoice =='hDeltaC' or dchoice=='hDeltaJ' or d=='FC2':
        thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        epg = cxa.pdat['offset_eb_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        dsecond = ''
        tc = colours[dchoice]
    elif dchoice =='hDeltaC_FC2_EPG':
        if dsecond=='FC2':
            thdc = cxa.pdat['offset_fsb1_ch1_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        elif dsecond == 'hDeltaC':
            thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        tc = colours[dsecond]
    hdcall = np.append(hdcall,thdc)
    
    epgall = np.append(epgall,epg)
    
epgall = 180*epgall/np.pi
hdcall = 180*hdcall/np.pi
plotbins = np.linspace(-180,180,17)
h,xb,yb = np.histogram2d(epgall,hdcall,plotbins)
medians, edges, bin_numbers = stats.binned_statistic(
     epgall,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
)

pedge = plotbins[1:]-np.mean(np.diff(plotbins))/2
hs = np.sum(h,axis=1)[:,np.newaxis]
hsT = np.tile(hs,(1,len(hs)))
hnorm = h/hsT
extent = [xb[0], xb[-1], yb[0], yb[-1]]
plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=.2)
plt.scatter(pedge,medians,s=40,color=tc)

plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='--')
plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='--')
plt.yticks([-180,-90,0,90,180])
plt.xticks([-180,-90,0,90,180])
plt.xlabel('EPG phase (deg)')
plt.ylabel(dchoice+ ' phase (deg)')
#plt.colorbar()
plt.title('Pre odour')

plt.savefig(os.path.join(savedir,dchoice+dsecond+'_PreAirAll_hist_.png'))
plt.savefig(os.path.join(savedir,dchoice+dsecond+'_PreAirAll_hist_.pdf'))
# post bins all not just jumps


for ip,p in enumerate(pbins):
    epgall = np.array([])
    hdcall = np.array([])
    if ip<len(pbins)-1:
        b1 = pbins[ip]*10
        b2 = pbins[ip+1]*10
    else:
        b1 = pbins[ip]*10
        b2 = np.inf
    
    for f in data:
        cxa = data[f].cxa
        cxa.get_jumps()
        jumps = cxa.get_entries_exits_like_jumps()
        if  dchoice =='hDeltaC' or dchoice=='hDeltaJ' or d=='FC2':
            thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
            epg = cxa.pdat['offset_eb_phase'].to_numpy()*cxa.side*-1
        elif dchoice =='hDeltaC_FC2_EPG':
            if dsecond=='FC2':
                thdc = cxa.pdat['offset_fsb1_ch1_phase'].to_numpy()*cxa.side*-1
            elif dsecond == 'hDeltaC':
                thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()*cxa.side*-1
            epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()*cxa.side*-1
            
            
            
        dx = np.array([],dtype='int')
        
        for ij,j in enumerate(jumps):
            d1 = min(j[1]+b1,j[2])
            d2 = min(j[1]+b2,j[2])
            if d1==d2:
                continue
            dx = np.append(dx,np.arange(d1,d2))
            
        
        
        hdcall = np.append(hdcall,thdc[dx])
        epgall = np.append(epgall,epg[dx])

    epgall = 180*epgall/np.pi
    hdcall = 180*hdcall/np.pi
    #plotbins = np.linspace(-180,180,20)
    plt.figure()
    h,xb,yb = np.histogram2d(epgall,hdcall,plotbins)
    hs = np.sum(h,axis=1)[:,np.newaxis]
    hsT = np.tile(hs,(1,len(hs)))
    hnorm = h/hsT
    extent = [xb[0], xb[-1], yb[0], yb[-1]]
    plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=.2)
    
    h,xb,yb = np.histogram2d(epgall,hdcall,plotbins)
    medians, edges, bin_numbers = stats.binned_statistic(
         epgall,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
    )

    pedge = plotbins[1:]-np.mean(np.diff(plotbins))/2
    plt.scatter(pedge,medians,s=40,color=tc)


    
    hs = np.sum(h,axis=1)[:,np.newaxis]
    
    plt.yticks([-180,-90,0,90,180])
    plt.xticks([-180,-90,0,90,180])
    plt.xlabel('EPG phase (deg)')
    plt.ylabel(dchoice+ ' phase (deg)')
    plt.plot([-180,180],[-90,-90],color='r',linestyle='--')
    plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='-')
    plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='-')
    plt.plot([-180,180],[-180,180],color=[.2,.2,1],linestyle='--')
    #plt.colorbar()
    plt.title(dchoice+' Post odour '+str(b1/10) +' - '+str(b2/10))
    # post bins
    plt.savefig(os.path.join(savedir,dchoice+dsecond+'_PostAirAll_hist_'+ str(b1/10) +' - '+str(b2/10) + '.png'))
    plt.savefig(os.path.join(savedir,dchoice+dsecond+'_PostAirAll_hist_'+ str(b1/10) +' - '+str(b2/10)   + '.pdf'))

# for just jump returns

for ip,p in enumerate(pbins):
    epgall = np.array([])
    hdcall = np.array([])
    if ip<len(pbins)-1:
        b1 = pbins[ip]*10
        b2 = pbins[ip+1]*10
    else:
        b1 = pbins[ip]*10
        b2 = np.inf
    
    for f in data:
        cxa = data[f].cxa
        cxa.get_jumps()
        jumps = cxa.get_jumps()
        if  dchoice =='hDeltaC' or dchoice=='hDeltaJ' or d== 'FC2':
            thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
            epg = cxa.pdat['offset_eb_phase'].to_numpy()*cxa.side*-1
        elif dchoice =='hDeltaC_FC2_EPG':
            if dsecond=='FC2':
                thdc = cxa.pdat['offset_fsb1_ch1_phase'].to_numpy()*cxa.side*-1
            elif dsecond == 'hDeltaC':
                thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()*cxa.side*-1
            epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()*cxa.side*-1
        dx = np.array([],dtype='int')
        
        for ij,j in enumerate(jumps):
            d1 = min(j[1]+b1,j[2])
            d2 = min(j[1]+b2,j[2])
            if d1==d2:
                continue
            dx = np.append(dx,np.arange(d1,d2))
            
        
        
        hdcall = np.append(hdcall,thdc[dx])
        epgall = np.append(epgall,epg[dx])

    epgall = 180*epgall/np.pi
    hdcall = 180*hdcall/np.pi
   # plotbins = np.linspace(-180,180,20)
    plt.figure()
    
    h,xb,yb = np.histogram2d(epgall,hdcall,plotbins)
    hs = np.sum(h,axis=1)[:,np.newaxis]
    hsT = np.tile(hs,(1,len(hs)))
    hnorm = h/hsT
    extent = [xb[0], xb[-1], yb[0], yb[-1]]
    plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=.2)
    medians, edges, bin_numbers = stats.binned_statistic(
         epgall,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
    )

    pedge = plotbins[1:]-np.mean(np.diff(plotbins))/2
    plt.scatter(pedge,medians,s=40,color=tc)

    
    plt.yticks([-180,-90,0,90,180])
    plt.xticks([-180,-90,0,90,180])
    plt.xlabel('EPG phase (deg)')
    plt.ylabel(dchoice+ ' phase (deg)')
    plt.plot([-180,180],[-90,-90],color='r',linestyle='--')
    plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='-')
    plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='-')
    plt.plot([-180,180],[-180,180],color=[.2,.2,1],linestyle='--')
   # plt.colorbar()
    plt.title(dchoice+' Post odour '+str(b1/10) +' - '+str(b2/10))
    # post bins
    plt.savefig(os.path.join(savedir,dchoice+dsecond+'_PostAirJump_hist_'+ str(b1/10) +' - '+str(b2/10) + '.png'))
    plt.savefig(os.path.join(savedir,dchoice+dsecond+'_PostAirJump_hist_'+ str(b1/10) +' - '+str(b2/10)   + '.pdf'))



#%% Do comparison of binned means with vector sums, fit weights
plotbins = np.linspace(-np.pi,np.pi,17)
plt.close('all')
pbins = [0,4,8,12,16,32,60]
dchoices = ['FC2','hDeltaC','hDeltaJ']
datamatrix = np.zeros((len(plotbins)-1,len(pbins),len(dchoices)))
amps = np.zeros((len(pbins),len(dchoices)))
for ip,p in enumerate(pbins):
    for di,d in enumerate(dchoices):
        data = all_flies[d]
        epgall = np.array([])
        hdcall = np.array([])
        Aall = np.array([])
        if ip<len(pbins)-1:
            b1 = pbins[ip]*10
            b2 = pbins[ip+1]*10
        else:
            b1 = pbins[ip]*10
            b2 = np.inf
        
        for f in data:
            cxa = data[f].cxa
            cxa.get_jumps()
            jumps = cxa.get_entries_exits_like_jumps()
            if  d =='hDeltaC' or d=='hDeltaJ' or 'FC2':
                thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
                epg = cxa.pdat['offset_eb_phase'].to_numpy()*cxa.side*-1
                A = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
            elif d =='hDeltaC_FC2_EPG':
                thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()*cxa.side*-1
                epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()*cxa.side*-1
                A = np.mean(cxa.pdat['offset_fsb2_ch2_phase'],axis=1)
            dx = np.array([],dtype='int')
            
            for ij,j in enumerate(jumps):
                d1 = min(j[1]+b1,j[2])
                d2 = min(j[1]+b2,j[2])
                if d1==d2:
                    continue
                dx = np.append(dx,np.arange(d1,d2))
                
            
            
            hdcall = np.append(hdcall,thdc[dx])
            epgall = np.append(epgall,epg[dx])
            Aall  = np.append(Aall,A[dx])
        #epgall = 180*epgall/np.pi
       # hdcall = 180*hdcall/np.pi
       
        medians, edges, bin_numbers = stats.binned_statistic(
             epgall,hdcall, statistic=lambda a: stats.circmean(a,high=np.pi,low=-np.pi), bins=plotbins
        )
        
        datamatrix[:,ip,di] = medians
        amps[ip,di] = np.mean(Aall)

pedge = plotbins[1:]-np.mean(np.diff(plotbins))/2
px = np.sin(pedge);py=np.cos(pedge)
colours = uplt.columnar_colours()
for ip,p in enumerate(pbins):
    plt.figure()
    tdata = datamatrix[:,ip,:].squeeze()
    
    #hDeltaJ vector sum
    
    tx = np.sin(tdata[:,2]);ty = np.cos(tdata[:,2])
    hsum = np.arctan2(tx+px,ty+py)
    
    
    plt.subplot(2,2,1)
    plt.scatter(pedge,tdata[:,0],color=colours[2,:])
    plt.scatter(pedge,tdata[:,2],color=colours[0,:])
    plt.scatter(pedge,hsum,color='k')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    #hDeltaC vector sum
    tx = np.sin(tdata[:,1]);ty = np.cos(tdata[:,1])
    hsum = np.arctan2(tx+px,ty+py)
    plt.subplot(2,2,2)
    plt.scatter(pedge,tdata[:,0],color=colours[2,:])
    plt.scatter(pedge,tdata[:,1],color=colours[1,:])
    plt.scatter(pedge,hsum,color='k')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    #hDeltaC and hDeltaJ vector sum
    tx = np.sin(tdata[:,1])+np.sin(tdata[:,2]);ty = np.cos(tdata[:,1])+np.cos(tdata[:,2])
    hsum = np.arctan2(tx+px,ty+py)
    plt.subplot(2,2,3)
    plt.scatter(pedge,tdata[:,0],color=colours[2,:])
    plt.scatter(pedge,tdata[:,1],color=colours[1,:])
    plt.scatter(pedge,tdata[:,2],color=colours[0,:])
    plt.scatter(pedge,hsum,color='k')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    
#%% Fit vector sum
colours = uplt.columnar_colours()
savedir = r'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'
from scipy.optimize import minimize 
from scipy.special import i0e, softmax


def wrap(a):
    return (a+np.pi)% (2 *np.pi)-np.pi

def predicted_direction(phi,w):
    mx = np.cos(phi)@ w
    my = np.sin(phi)@ w
    
    return np.arctan2(my,mx)

def vm_logpdf(theta,mu,kappa):
    return kappa* (np.cos(theta-mu)-1.0) - np.log(2*np.pi*i0e(kappa))

def neg_log_lik(params,theta,phi):
    p = phi.shape[1]
    w = softmax(np.append(params[: p-1],0.0))
    kappa = np.exp(params[-1])
    mu = predicted_direction(phi,w)
    return -np.sum(vm_logpdf(theta,mu,kappa))


rng = np.random.default_rng(42)
def fit(theta,phi,n_starts=10):
    p = phi.shape[1]
    best = None
    for _ in range(n_starts):
        x0 = np.append(rng.normal(0,1,p-1),np.log(2.0))
        res = minimize(neg_log_lik,x0,args = (theta,phi),method="BFGS")
        if best is None or res.fun<best.fun:
            best = res
        
    w = softmax(np.append(best.x[: p-1],0.0))
    return w,np.exp(best.x[-1]),best.fun

def aic(nll_val, k_params):
    return 2 * k_params + 2 * nll_val

def bootstrap_ci(theta,phi,n_boot=300,alpha=0.05):
    n = len(theta)
    ws, ks = [],[]
    for _ in range(n_boot):
        idx = rng.integers(0,n,n)
        w, k, _ = fit(theta[idx],phi[idx],n_starts=3)
        ws.append(w)
        ks.append(k)
        
    ws, ks = np.array(ws),np.array(ks)
    q = [100 * alpha /2 , 100*(1-alpha/2)]
    return np.percentile(ws,q,axis=0),np.percentile(ks,q)



plotbins = np.linspace(-np.pi,np.pi,17)
pbx = plotbins[1:]-np.mean(np.diff(plotbins))/2
plt.close('all')
aic_values = np.zeros((datamatrix.shape[1],4))
all_w = np.zeros((datamatrix.shape[1],3))
pdex = np.abs(pbx)<3*np.pi/4
groups = ["0-4","4-8","8-12","12-16","16-32","32-60","60 <"]
plt.figure(2,figsize=(15,15))
for i in range(datamatrix.shape[1]):
    y = datamatrix[:,i,0]
    x = np.column_stack([datamatrix[:, i, 1:], pbx])
    assert y.ndim == 1 and x.ndim == 2 and len(y) == len(x)
    
    wh,kh,n11 = fit(y[pdex],x[pdex,:],20)
    all_w[i,:] =wh
    plt.figure(1)
    
    plt.plot(wh)
    plt.plot(amps[i,1:])
    plt.ylim([0,1])
    plt.figure(2)
    plt.subplot(1,len(groups),i+1)
    plt.scatter(pbx,y,color=colours[2,:])
    plt.title(groups[i])
    yp = predicted_direction(x,wh)
    plt.scatter(pbx,yp,color='r')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    plt.plot([0,0],[-np.pi,np.pi],color='k')
    plt.plot([-np.pi,np.pi],[0,0],color='k')
    plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
    plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
    plt.gca().set_aspect('equal')
    #(w_lo, w_hi), (k_lo, k_hi) = bootstrap_ci(y,x)
    
    
    # plt.figure()
    # plt.subplot(1,2,1)
    # plt.plot(w_lo)
    # plt.plot(w_hi)
    # plt.subplot(1,2,2)
    # plt.scatter(pbx,y,color='b')
    # yp = predicted_direction(x,w_hi)
    # plt.scatter(pbx,yp,color='r')
    # plt.scatter(pbx,yp,color='r')
    # plt.xlim([-np.pi,np.pi])
    # plt.ylim([-np.pi,np.pi])
    # plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    
    aic_values[i,0] = aic(n11,3)
    
    for j in range(3):
        keep = [i for i in range(3) if i != j]
        xp = x[:,keep]
        _, _, nll_red = fit(y[pdex], xp[pdex, :],20)
        aic_values[i,j+1] = aic(nll_red, 2)

plt.figure(2)
plt.gca().set_aspect('equal')
plt.savefig(os.path.join(savedir,'VectorSumExamples.pdf'))
plt.savefig(os.path.join(savedir,'VectorSumExamples.png'))

best_mods = np.argmin(aic_values,axis=1)
plt.figure(figsize=(10,5))
plt.plot(pbins,best_mods,color='k')
plt.scatter(pbins,best_mods,color='r',zorder=5)
plt.yticks(np.arange(0,4),labels=['hdJ hdC heading','hdJ heading','hdC heading','hdC hDj'])
plt.subplots_adjust(left=.3)
plt.ylabel('AIC best model')
plt.xticks(pbins,labels=groups,rotation=45)
plt.savefig(os.path.join(savedir,'AIC_best_model.png'))
plt.savefig(os.path.join(savedir,'AIC_best_model.pdf'))

#%%        
fig,ax = plt.subplots(layout="constrained")
width =.8
offset = 0
xbar = np.arange(3)
colours = uplt.columnar_colours()
colours[2,:] = 0.0
colours[[0,1],:] = colours[[1,0],:]
pbins = [0,4,8,12,16,32,60]
groups = ["0-4","4-8","8-12","12-16","16-32","32-60","60 <"]
for i in range(all_w.shape[0]):
    offset =width*i*5 
    rects = ax.bar(xbar+offset,all_w[i,:],width,color=colours,edgecolor=None)
    
    #ax.bar_label(rects, padding=2)
xtick = np.arange(1,29,4)
ax.set_xticks(xtick, groups)
plt.xlabel('timebins (s)')
plt.ylabel('weights')
plt.savefig(os.path.join(savedir,'VectorSumWeights.pdf'))
plt.savefig(os.path.join(savedir,'VectorSumWeights.png'))

#%% Similar phase analysis but for dual colour imaging

plotbins = np.linspace(-np.pi,np.pi,17)
plt.close('all')
pbins = [0,4,8,12,16,32,60]
dchoices = ['hDeltaC_FC2_EPG','hDeltaC_FC2_EPG','hDeltaJ']
datamatrix = np.zeros((len(plotbins)-1,len(pbins),len(dchoices)))
amps = np.zeros((len(pbins),len(dchoices)))
for ip,p in enumerate(pbins):
    for di,d in enumerate(dchoices):
        data = all_flies[d]
        epgall = np.array([])
        hdcall = np.array([])
        Aall = np.array([])
        if ip<len(pbins)-1:
            b1 = pbins[ip]*10
            b2 = pbins[ip+1]*10
        else:
            b1 = pbins[ip]*10
            b2 = np.inf
        
        for f in data:
            cxa = data[f].cxa
            cxa.get_jumps()
            jumps = cxa.get_entries_exits_like_jumps()
            
            if d =='hDeltaC_FC2_EPG':
                if di==0:
                    thdc = cxa.pdat['offset_fsb1_ch1_phase'].to_numpy().squeeze()*cxa.side*-1
                    A = np.mean(cxa.pdat['wedges_fsb2_ch2'],axis=1)
                else:
                    thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy().squeeze()*cxa.side*-1
                    A = np.mean(cxa.pdat['wedges_fsb1_ch1'],axis=1)
                epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()*cxa.side*-1
                
                
            else: 
                thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
                epg = cxa.pdat['offset_eb_phase'].to_numpy()*cxa.side*-1
                A = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
            dx = np.array([],dtype='int')
            
            for ij,j in enumerate(jumps):
                d1 = min(j[1]+b1,j[2])
                d2 = min(j[1]+b2,j[2])
                if d1==d2:
                    continue
                dx = np.append(dx,np.arange(d1,d2))
                
            
            
            hdcall = np.append(hdcall,thdc[dx])
            epgall = np.append(epgall,epg[dx])
            Aall  = np.append(Aall,A[dx])
        #epgall = 180*epgall/np.pi
       # hdcall = 180*hdcall/np.pi
       
        medians, edges, bin_numbers = stats.binned_statistic(
             epgall,hdcall, statistic=lambda a: stats.circmean(a,high=np.pi,low=-np.pi), bins=plotbins
        )
        
        datamatrix[:,ip,di] = medians
        amps[ip,di] = np.mean(Aall)

pedge = plotbins[1:]-np.mean(np.diff(plotbins))/2
px = np.sin(pedge);py=np.cos(pedge)
colours = uplt.columnar_colours()
for ip,p in enumerate(pbins):
    plt.figure()
    tdata = datamatrix[:,ip,:].squeeze()
    
    #hDeltaJ vector sum
    
    tx = np.sin(tdata[:,2]);ty = np.cos(tdata[:,2])
    hsum = np.arctan2(tx+px,ty+py)
    
    
    plt.subplot(2,2,1)
    plt.scatter(pedge,tdata[:,0],color=colours[2,:])
    plt.scatter(pedge,tdata[:,2],color=colours[0,:])
    plt.scatter(pedge,hsum,color='k')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    #hDeltaC vector sum
    tx = np.sin(tdata[:,1]);ty = np.cos(tdata[:,1])
    hsum = np.arctan2(tx+px,ty+py)
    plt.subplot(2,2,2)
    plt.scatter(pedge,tdata[:,0],color=colours[2,:])
    plt.scatter(pedge,tdata[:,1],color=colours[1,:])
    plt.scatter(pedge,hsum,color='k')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    #hDeltaC and hDeltaJ vector sum
    tx = np.sin(tdata[:,1])+np.sin(tdata[:,2]);ty = np.cos(tdata[:,1])+np.cos(tdata[:,2])
    hsum = np.arctan2(tx+px,ty+py)
    plt.subplot(2,2,3)
    plt.scatter(pedge,tdata[:,0],color=colours[2,:])
    plt.scatter(pedge,tdata[:,1],color=colours[1,:])
    plt.scatter(pedge,tdata[:,2],color=colours[0,:])
    plt.scatter(pedge,hsum,color='k')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    


#%%
plotbins = np.linspace(-np.pi,np.pi,17)
pbx = plotbins[1:]-np.mean(np.diff(plotbins))/2
plt.close('all')
aic_values = np.zeros((datamatrix.shape[1],4))
all_w = np.zeros((datamatrix.shape[1],3))
pdex = np.abs(pbx)<3*np.pi/4
groups = ["0-4","4-8","8-12","12-16","16-32","32-60","60 <"]
plt.figure(2,figsize=(15,15))
colours = uplt.columnar_colours()
for i in range(datamatrix.shape[1]):
    y = datamatrix[:,i,0]
    x = np.column_stack([datamatrix[:, i, 1:], pbx])
    assert y.ndim == 1 and x.ndim == 2 and len(y) == len(x)
    
    wh,kh,n11 = fit(y[pdex],x[pdex,:],20)
    all_w[i,:] =wh
    plt.figure(1)
    
    plt.plot(wh)
    plt.plot(amps[i,1:])
    plt.ylim([0,1])
    plt.figure(2)
    plt.subplot(1,len(groups),i+1)
    plt.scatter(pbx,y,color=colours[2,:])
    plt.title(groups[i])
    yp = predicted_direction(x,wh)
    plt.scatter(pbx,yp,color='r')
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    plt.plot([0,0],[-np.pi,np.pi],color='k')
    plt.plot([-np.pi,np.pi],[0,0],color='k')
    plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
    plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
    plt.gca().set_aspect('equal')
    #(w_lo, w_hi), (k_lo, k_hi) = bootstrap_ci(y,x)
    
    
    # plt.figure()
    # plt.subplot(1,2,1)
    # plt.plot(w_lo)
    # plt.plot(w_hi)
    # plt.subplot(1,2,2)
    # plt.scatter(pbx,y,color='b')
    # yp = predicted_direction(x,w_hi)
    # plt.scatter(pbx,yp,color='r')
    # plt.scatter(pbx,yp,color='r')
    # plt.xlim([-np.pi,np.pi])
    # plt.ylim([-np.pi,np.pi])
    # plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='k',linestyle='--')
    
    aic_values[i,0] = aic(n11,3)
    
    for j in range(3):
        keep = [i for i in range(3) if i != j]
        xp = x[:,keep]
        _, _, nll_red = fit(y[pdex], xp[pdex, :],20)
        aic_values[i,j+1] = aic(nll_red, 2)

plt.figure(2)
plt.gca().set_aspect('equal')
plt.savefig(os.path.join(savedir,'VectorSumExamples_hDC_dual_colour.pdf'))
plt.savefig(os.path.join(savedir,'VectorSumExamples_hDC_dual_colour.png'))

best_mods = np.argmin(aic_values,axis=1)
plt.figure(figsize=(10,5))
plt.plot(pbins,best_mods,color='k')
plt.scatter(pbins,best_mods,color='r',zorder=5)
plt.yticks(np.arange(0,4),labels=['hdJ hdC heading','hdJ heading','hdC heading','hdC hDj'])
plt.subplots_adjust(left=.3)
plt.ylabel('AIC best model')
plt.xticks(pbins,labels=groups,rotation=45)
plt.savefig(os.path.join(savedir,'AIC_best_model_hDC_dual_colour.png'))
plt.savefig(os.path.join(savedir,'AIC_best_model_hDC_dual_colour.pdf'))

fig,ax = plt.subplots(layout="constrained")
width =.8
offset = 0
xbar = np.arange(3)
colours = uplt.columnar_colours()
colours[2,:] = 0.0
colours[[0,1],:] = colours[[1,0],:]
pbins = [0,4,8,12,16,32,60]
groups = ["0-4","4-8","8-12","12-16","16-32","32-60","60 <"]
for i in range(all_w.shape[0]):
    offset =width*i*5 
    rects = ax.bar(xbar+offset,all_w[i,:],width,color=colours,edgecolor=None)
    
    #ax.bar_label(rects, padding=2)
xtick = np.arange(1,29,4)
ax.set_xticks(xtick, groups)
plt.xlabel('timebins (s)')
plt.ylabel('weights')
plt.savefig(os.path.join(savedir,'VectorSumWeights_hDC_dual_colour.pdf'))
plt.savefig(os.path.join(savedir,'VectorSumWeights_hDC_dual_colour.png'))