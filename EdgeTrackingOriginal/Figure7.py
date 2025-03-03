# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 11:53:52 2024

@author: dowel
"""

#%%from analysis_funs.regression import fci_regmodel

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
from EdgeTrackingOriginal.ETpap_plots.ET_paper import ET_paper
from scipy.stats import circmean, circstd
plt.rcParams['pdf.fonttype'] = 42 
colours2 = np.array([[106,207,246],[237,30,36],[168,170,173],[6,149,207]])/255
colours = np.array([[81,156,204],[84,39,143],[237,30,36],[6,149,207]])/255
#%%
savedir = "Y:\\Papers, Review, Theses\\RutaLabPapers\\ET_2024\\MyContribution\\Panels"
datadirs = [
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f1\\Trial3",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f2\\Trial3",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240502\\f1\\Trial2",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240514\\f1\\Trial2",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\241104\\f1\\Trial5",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\241106\\f1\\Trial2",
#r"Y:\Data\FCI\Hedwig\FC2_maimon2\250128\f1\Trial4" left out of original six to be consistent
]
#%% Load data
eg_fly =3
datadir = datadirs[eg_fly]
d = datadir.split("\\")
name = d[-3] + '_' + d[-2] + '_' + d[-1]
etp = ET_paper(datadir)

#%% Panel B heatmaps, phase, behaviour
plt.close('all')
#Full trajectory with boxes
tstart = 0
tend = 100000
etp.example_trajectory(tstart,tend,boxes=np.array([[60,110],[60*4.2,60*5.5]]))
plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Amenotaxis trajectory
plt.close('all')
tstart = 32
tend = 110
etp.example_trajectory(tstart,tend)
plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Amenotaxis heats
etp.example_heatmaps(tstart,tend)
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'Phase_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Plume_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Jump trajectory
plt.close('all')
tstart = 60*4.2
tend = 60*5.5
etp.example_trajectory(tstart,tend)
plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Jump heats
etp.example_heatmaps(tstart,tend)
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'Phase_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Plume_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
#%% Colour bars
x = np.zeros((10,10))
plt.close('all')
plt.figure()
plt.imshow(x,cmap='Blues')
plt.colorbar()
plt.savefig(os.path.join(savedir,'BlueBar.pdf'))

plt.figure()
plt.imshow(x,cmap='Purples')
plt.colorbar()
plt.savefig(os.path.join(savedir,'PurpleBar.pdf'))

plt.figure()
plt.imshow(x,cmap='Greys')
plt.colorbar()
plt.savefig(os.path.join(savedir,'GreyBar.pdf'))

#%% Panel C Segment amplitude - does not look amazing
plt.close('all')
colours = np.array([[81,156,204],[84,39,143],[237,30,36],[6,149,207]])/255
fig, ax1 = plt.subplots()
ax1_2 = ax1.twinx()
fig2, ax2 = plt.subplots()
ax2_2 = ax2.twinx()
fig2,ax3 = plt.subplots()
ax3_2 = ax3.twinx()

plt.figure()

total_wedges = np.zeros((16,2,3,len(datadirs)))
for i,d in enumerate(datadirs):
    print(d)
    etp = ET_paper(d)
    plt.figure(1)
    m_wedge = etp.mean_wedge_rotate('jumps',1)
    total_wedges[:,:,0,i] = m_wedge
    ax1.plot(m_wedge[:,0],color=colours[0,:],alpha=0.3)
    ax1_2.plot(m_wedge[:,1],color=colours[1,:],alpha=0.3)
    plt.figure(2)
    m_wedge = etp.mean_wedge_rotate('pre air',1)
    ax2.plot(m_wedge[:,0],color=colours[0,:],alpha=0.3)
    ax2_2.plot(m_wedge[:,1],color=colours[1,:],alpha=0.3)
    total_wedges[:,:,1,i] = m_wedge
    
    plt.figure(3)
    m_wedge = etp.mean_wedge_rotate('pre jumps',1)
    ax3.plot(m_wedge[:,0],color=colours[0,:],alpha=0.3)
    ax3_2.plot(m_wedge[:,1],color=colours[1,:],alpha=0.3)
    total_wedges[:,:,2,i] = m_wedge
    
    plt.figure(4)
    plt.plot(m_wedge[:,1],color=colours2[1,:],alpha=0.3)
    plt.plot(total_wedges[:,1,0,i],color=colours2[0,:],alpha=0.3)
    plt.figure(5)
    plt.plot(m_wedge[:,0],color=colours2[1,:],alpha=0.3)
    plt.plot(total_wedges[:,0,0,i],color=colours2[0,:],alpha=0.3)
plt_wedges = np.mean(total_wedges,axis=3)



plt.figure(1)
ax1.plot(plt_wedges[:,0,0],color=colours[0,:],marker='s')
ax1_2.plot(plt_wedges[:,1,0],color=colours[1,:],marker='s')
ax1.set_ylabel('dF/F0 EPG',color=colours[0,:])
ax1_2.set_ylabel('dF/F0 FC2',color=colours[1,:])
ax1.set_ylim([0,0.7])
ax1_2.set_ylim([0,0.7])
ax1.plot([3,3],[0,1],color='r',linestyle='--')
ax1.set_xticks([0,7,15])
ax1.set_xlabel('EB/FSB wedge number')

plt.figure(2)
ax2.plot(plt_wedges[:,0,1],color=colours[0,:],marker='s')
ax2_2.plot(plt_wedges[:,1,1],color=colours[1,:],marker='s')
ax2.set_ylabel('dF/F0 EPG',color=colours[0,:])
ax2_2.set_ylabel('dF/F0 FC2',color=colours[1,:])
ax2.set_ylim([0,0.7])
ax2_2.set_ylim([0,0.7])
ax2.set_xticks([0,7,15])
ax1.set_xlabel('EB/FSB wedge number')

plt.figure(3)
ax3.plot(plt_wedges[:,0,2],color=colours[0,:],marker='s')
ax3_2.plot(plt_wedges[:,1,2],color=colours[1,:],marker='s')
ax3.set_ylabel('dF/F0 EPG',color=colours[0,:])
ax3_2.set_ylabel('dF/F0 FC2',color=colours[1,:])
ax3.set_ylim([0,0.7])
ax3_2.set_ylim([0,0.7])
ax3.set_xticks([0,7,15])
ax3.set_xlabel('EB/FSB wedge number')

plt.figure(4)
plt.plot(plt_wedges[:,1,2],color=colours2[1,:],marker='s')
plt.plot(plt_wedges[:,1,0],color=colours2[0,:],marker='s')

plt.figure(5)
plt.plot(plt_wedges[:,0,2],color=colours2[1,:],marker='s')
plt.plot(plt_wedges[:,0,0],color=colours2[0,:],marker='s')
#%% Panel C Histograms
# 1. Phase epg vs FC2 scatter/2D histogram
plt.close('all')
bins = np.linspace(-np.pi,np.pi,10)
pltbins =180* (bins[1:]-np.mean(np.diff(bins))/2)/np.pi
plt.figure()
plt.figure()
plt.figure()
plt.figure()
colours2 = np.array([[106,207,246],[237,30,36],[168,170,173],[6,149,207]])/255
colours = np.array([[0,0,0,],[81,156,204],[237,30,36],[6,149,207],[84,39,143]])/255
colours_leavereturn = np.array([[247,115,17],[48,120,198]])/255
phall = np.zeros((len(pltbins),3,3,len(datadirs)))
for i,d in enumerate(datadirs):
    etp = ET_paper(d)
    ph = etp.phase_histogram('jumps',1,bins)
    phall[:,:,0,i] = ph
    plt.figure(1)
    plt.plot(pltbins,ph[:,1],color=colours[0,:],alpha=0.3)
    plt.plot(pltbins,ph[:,2],color=colours[1,:],alpha=0.3)
    
    plt.figure(4)
    plt.plot(pltbins,ph[:,0],color=colours_leavereturn[1,:],alpha=0.3)
    
    #pham = etp.phase_histogram('pre air',1,bins)
    pham = etp.phase_histogram('pre air anemotaxis',1,bins)
    phall[:,:,1,i] = pham
    plt.figure(2)
    plt.plot(pltbins,pham[:,1],color=colours[0,:],alpha=0.3)
    plt.plot(pltbins,pham[:,2],color=colours[1,:],alpha=0.3)
    
    plt.figure(5)
    plt.plot(pltbins,pham[:,0],color=colours[1,:],alpha=0.3)
    
    
    pham = etp.phase_histogram('pre jumps',0.5,bins)
    phall[:,:,2,i] = pham
    plt.figure(3)
    plt.plot(pltbins,pham[:,1],color=colours[0,:],alpha=0.3)
    plt.plot(pltbins,pham[:,2],color=colours[1,:],alpha=0.3)
    
    plt.figure(4)
    plt.plot(pltbins,pham[:,0],color=colours_leavereturn[0,:],alpha=0.3)
pltmean = np.mean(phall,axis=3)

for r in range(3):
    plt.figure(r+1)
    for i in range(2):
        plt.plot(pltbins,pltmean[:,i+1,r],color=colours[i,:],marker='s')
    plt.xticks([-180,-90,0,90,180])
    
    if r==0:
        plt.plot([-90,-90],[0,0.4],color=colours2[1,:],linestyle='--')
        plt.plot([0,0],[0,0.4],color='k',linestyle='--')
    if r==2:
        plt.plot([90,90],[0,0.4],color=colours2[1,:],linestyle='--')
        plt.plot([0,0],[0,0.4],color='k',linestyle='--')
    plt.xticks([-180,-90,0,90,180])    
    plt.xlabel('phase (degrees)')
    plt.ylabel('probability')
    s = plt.gcf()
    s.set_size_inches([5,4])
plt.figure(4)
plt.plot(pltbins,pltmean[:,0,0],color=colours_leavereturn[1,:],marker='s')
plt.plot(pltbins,pltmean[:,0,2],color=colours_leavereturn[0,:],marker='s')
plt.plot([-90,-90],[0,0.4],color=colours2[1,:],linestyle='--')
plt.plot([90,90],[0,0.4],color=colours2[1,:],linestyle='--')
plt.plot([0,0],[0,0.4],color='k',linestyle='--')

plt.xticks([-180,-90,0,90,180])
plt.xlabel('EPG-FC2 phase (degrees)')
plt.ylabel('probability')
s = plt.gcf()
s.set_size_inches([5,4])


plt.figure(5)
plt.plot(pltbins,pltmean[:,0,1],color=colours_leavereturn[1,:],marker='s')
#plt.plot(pltbins,pltmean[:,0,1],color=colours_leavereturn[0,:],marker='s')
plt.plot([-90,-90],[0,0.4],color=colours2[1,:],linestyle='--')
plt.plot([90,90],[0,0.4],color=colours2[1,:],linestyle='--')
plt.plot([0,0],[0,0.4],color='k',linestyle='--')

plt.xticks([-180,-90,0,90,180])
plt.xlabel('EPG-FC2 phase (degrees)')
plt.ylabel('probability')
s = plt.gcf()
s.set_size_inches([5,4])

savenames = ['PhaseHist_Jump','PhaseHist_PreAir','PhaseHist_PreJump','PhaseHist_Diff','PhaseHist_Diff_preair']
for i in range(5):
    plt.figure(i+1)
    plt.savefig(os.path.join(savedir,savenames[i]+'.pdf'))

#%%
plt.close('all')
plt.figure()
plt.figure()
tphase_jump = np.empty((1,3))
tphase_air = np.empty((1,3))
for i, d in enumerate(datadirs):
    etp = ET_paper(d)
    plt.figure(1)
    tphase = etp.output_phases('jumps',1,i)
    tphase_jump = np.append(tphase_jump,tphase,axis=0)
    plt.scatter(tphase[:,0],tphase[:,1],s=2,color='k',alpha=0.2)
    # 1. Phase and probability (as for grant)
    plt.figure(2)
    tphase = etp.output_phases('pre air',1,i)
    tphase_air = np.append(tphase_air,tphase,axis=0)
    plt.scatter(tphase[:,0],tphase[:,1],s=2,color='k',alpha=0.2)
    
plt.figure()
h = plt.hist2d(tphase_jump[:,0],tphase_jump[:,1],bins=30,density=True,cmap='bwr',vmin=-0.3,vmax=0.3)
plt.figure()
h = plt.hist2d(tphase_air[:,0],tphase_air[:,1],bins=30,density=True,cmap='Blues',vmin=0,vmax=0.3)
#%% Panel D return timecourse
# 1. Pseudotime

timecourse = 15
# 2. Realtime
plt.close('all')
for i, d in enumerate(datadirs):
    
    
    etp = ET_paper(d)
    p = etp.phase_time(timecourse)
    
    p[np.abs(p)==0.0] = np.nan
    pmean = circmean(p,high=np.pi, low=-np.pi,axis=2,nan_policy='omit') #%need to do circ mean
    t = np.arange(0,timecourse,0.1)
    #plt.plot(t,p[:,0,:],color=colours[0,:],alpha=0.1,linestyle=':')
    #plt.plot(t,p[:,1,:],color=colours[1,:],alpha=0.1,linestyle=':')
    tt = np.tile(t,(np.shape(p)[2],1))
    #plt.scatter(tt,p[:,0,:],alpha=0.4,s=5,color=colours[0,:])
    #plt.scatter(tt,p[:,1,:],alpha=0.4,s=5,color=colours[1,:])
    #plt.plot(t,pmean[:,0],color=colours[0,:],alpha=1)
    #plt.plot(t,pmean[:,1],color=colours[1,:],alpha=1)
    if i==0:
        pltmean = np.zeros((len(pmean),4,len(datadirs)))
    pltmean[:,:,i] = 180*pmean/np.pi


colours = np.array([[81,156,204],[84,39,143],[0,0,0],[6,149,207]])/255
colours = np.array([[0,0,0],[81,156,204]])/255
plt.figure()
pm = circmean(pltmean,high=180,low=-180,axis=2,nan_policy='omit')
zorder = [1,3,0]
for i in range(2):
    plt.plot(pltmean[:,i,:],t,color = colours[i,:],alpha=0.3,zorder=zorder[i])
    plt.plot(pm[:,i],t,color=colours[i,:],zorder=zorder[i])
 
    
    
plt.plot([-180,180],[timecourse,timecourse],color='k',linestyle='--')
plt.plot([0,0],[0,timecourse],color='k',linestyle='--')
plt.plot([-90,-90],[0,timecourse],color=colours2[1,:],linestyle='--')
plt.xlim([-180,180])
plt.xlabel('angle (deg)')
plt.ylabel('time (s)')

plt.ylim([0,timecourse+0.2])
plt.xticks([-180,0,180])
plt.savefig(os.path.join(savedir,'Phase_return_jump_realtime.pdf'))

plt.figure()
plt.plot(pltmean[:,3,:],t,color = colours[1,:],alpha=0.3,zorder=zorder[i])
plt.plot(pm[:,3],t,color=colours[1,:],zorder=zorder[i])
plt.plot([-180,180],[timecourse,timecourse],color='k',linestyle='--')
plt.plot([0,0],[0,timecourse],color='k',linestyle='--')
plt.plot([-90,-90],[0,timecourse],color=colours2[1,:],linestyle='--')
plt.xlim([-180,180])
plt.xlabel('FC2-EPG phase (deg)')
plt.ylabel('time (s)')
plt.ylim([0,timecourse+0.2])
plt.xticks([-180,0,180])
plt.savefig(os.path.join(savedir,'Phase_return_diff_jump_realtime.pdf'))
#%% Panel D bump amplitude
timecourse = 15
t = np.arange(0,timecourse,0.1)
for i, d in enumerate(datadirs):
    etp = ET_paper(d)
    p = etp.amp_time(timecourse)
    p[p<0.00000000001] = np.nan
    p[np.abs(p)==0.0] = np.nan
    pmean = np.nanmean(p,axis=3)
    if i==0:
        pltmean = np.zeros((len(pmean),2,3,len(datadirs)))
    pltmean[:,:,:,i] = pmean
zorder = [1,3,0]
pm = np.mean(pltmean,axis=3)
for i in range(3):
    plt.figure()
    for ir in range(2):
        plt.plot(pltmean[:,ir,i,:],t,color= colours[ir,:],alpha=0.3,zorder = zorder[ir])
        plt.plot(pm[:,ir,i],t,color=colours[ir,:],zorder=zorder[i])
#%% Trajectory mean
#phases,trajs= etp.trajectory_mean()
plt.close('all')
for i, d in enumerate(datadirs):
    plt.figure()
    
    etp = ET_paper(d)
    etp.plt_tmp(scalar=4,phase_num=6)
    plt.ylim([-35,35])
    plt.title(str(etp.cxa.side))
    if etp.cxa.side==1:
        plt.gca().invert_xaxis()
    plt.ylim([-50,50])
    plt.savefig(os.path.join(savedir,'Arrows_'+str(i)+ 'side_' +str(etp.cxa.side) + '.pdf'))
    
    
#etp.cxa.mean_jump_arrows()
#%% Trajectory all
plt.close('all')
for i,d in enumerate(datadirs):
    plt.figure()
    etp = ET_paper(d)
    etp.plt_tmp_all(scalar =7.5,phase_num=5)
    plt.savefig(os.path.join(savedir,'AllArrows'+str(i) + 'side_' +str(etp.cxa.side) + '.pdf'))

#%% Phase transition plot
plt.close('all')
plt.figure()
plt.plot([0,0],[-180,180],color='k',linestyle='--')
plt.plot([-180,180],[0,0],color='k',linestyle='--')
bins = 10
for i, d in enumerate(datadirs):
    
    
    etp = ET_paper(d)
    phases,trajs = etp.trajectory_mean(bins=bins)
    pmean = 180*circmean(phases,low=-np.pi,high=np.pi,axis=2)/np.pi
    plt.plot(pmean[:bins,1],pmean[:bins,0],color=colours2[1,:],alpha=0.3)
    plt.plot(pmean[bins-1:,1],pmean[bins-1:,0],color=colours2[3,:],alpha=0.3)
    plt.scatter(pmean[0,1],pmean[0,0],color=colours2[1,:],alpha=0.5)
    plt.scatter(pmean[-1,1],pmean[-1,0],color=colours2[3,:],alpha=0.5,marker='x')
    if i==0:
        pall = np.zeros((pmean.shape[0],3,len(datadirs)))
    pall[:,:,i] = pmean
pm = circmean(pall,low=-180,high=180,axis=2)
plt.plot(pm[:bins,1],pm[:bins,0],color=colours2[1,:])
plt.scatter(pm[-1,1],pm[-1,0],color=colours2[3,:],marker='x')
plt.scatter(pm[0,1],pm[0,0],color=colours2[1,:])
plt.plot(pm[bins-1:,1],pm[bins-1:,0],color=colours2[3,:])
plt.xlim([-180,180])
plt.ylim([-180,180])
plt.plot([-180,180],[-180,180],color='k',linestyle='--')
plt.xticks([-180,0,180])
plt.yticks([-180,0,180])
plt.ylabel('EPG phase (deg)')
plt.xlabel('FC2 phase (deg)')
plt.savefig(os.path.join(savedir,'PhaseRelationshipFC2_EPG.pdf'))
# %% distance from plume boundary

for i,d in enumerate(datadirs):
    etp = ET_paper(d)
    td,tdin = etp.distance_plume_boundary()
    tdat = np.ones((len(td),2))*i
    tdat[:,1] = td
    tdatin = np.ones((len(tdin),2))*i
    tdatin[:,1] = tdin
    if i==0:
        dists = tdat
        distsin = tdatin
    else:
        dists = np.append(dists,tdat,axis=0)
        distsin = np.append(distsin,tdatin,axis=0)
#%%
dmean = np.mean(dists[:,1])
dstd = np.std(dists[:,1])
n = len(dists)
print('D out mean: ',dmean ,' +- ', dstd,' N ',n)

dmedian = np.median(dists[:,1])
iqr = np.percentile(dists[:,1],[25,75])
print('D out median: ',dmedian ,' +- ', iqr,' N ',n)



#%%


plt.close('all')
c,h = np.histogram(dists[:,1],bins = np.arange(0,100,1))
plt.hist(-dists[:,1],np.arange(-150,1,1),color=colours2[0,:])
plt.hist(distsin[:,1],np.arange(0,11,1),color=colours2[1,:])
plt.xticks(np.append(np.arange(-120,1,20),10))
plt.xlim([-125,10])
dm = np.median(dists[:,1])
plt.plot([-dm,-dm],[0,100],color='k',linestyle='--')

dmin = np.median(distsin[:,1])
plt.plot([dmin,dmin],[0,100],color='k',linestyle='--')
plt.text(-dm-4,100,str(np.round(dm,decimals=1)))
plt.text(dmin+1,100,str(np.round(dmin,decimals=1)))
plt.xlabel('Crosswind range (mm)')
plt.ylabel('Number of bouts')
plt.savefig(os.path.join(savedir,'CrosswindPlumeDistance.pdf'))
#%% Extra: bump amplitude on returns and inside plume
# Bump amplitude and regression.





