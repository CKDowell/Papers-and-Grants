# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 14:37:35 2026

@author: dowel

Script plots panels for Figure 5 of K99 application

Includes combination of sine functions. Feed


"""
import numpy as np
import matplotlib.pyplot as plt 

#%% Model

angle1 = -np.pi/2
angle2 = np.pi/2
column_number = 16
gain1 =.5
gain2 = .25
angles = np.linspace(-np.pi,np.pi,column_number)

s1 = (np.cos(angles-angle1)+1)*gain1
s2 = (np.cos(angles-angle2)+1)*gain2

plt.plot(angles,s1)
plt.plot(angles,s2)


plt.plot(angles,s1+s2)


#%%
output = s1+s2
output[output<0] = 0


plt.figure()
plt.plot(angles,s1)
plt.plot(angles,s2)
plt.plot(angles,output)


plt.figure()
plt.plot(angles,s1)
plt.plot(angles,s2)
for i in range(15):
    output = s1+s2-i/10
    output[output<0] = 0



    plt.plot(angles,output)
#%% Toy network - feed4forward inputs + global feedback inhibition
# Two feedforward populations (s1,s2) connect one-to-one onto the output
# population according to angular position (identity connectivity - neuron
# at a given angle only drives the output neuron at that same angle).
# A single inhibitory neuron pools activity across all output neurons and
# feeds back uniformly (untuned) onto every output neuron.

gain1 =.5
gain2 = .5
s1 = (np.cos(angles-angle1)+1)*gain1
s2 = (np.cos(angles-angle2)+1)*gain2
ff_drive = s1+s2 # one-to-one angular connectivity means the drive to output neuron i is just s1[i]+s2[i]

tau = 1. # output neuron time constant
tau_i = .2 # inhibitory neuron time constant (faster than the output population)
dt = .01
T = 5
nsteps = int(T/dt)

def relu(x):
    return np.where(x>0,x,0)

def run_network(ff_drive,gain_inh):
    r = np.zeros_like(ff_drive)
    r_inh = 0
    for t in range(nsteps):
        r_inh = (r_inh+dt/tau_i*(-r_inh+r))*gain_inh
        r_inh[r_inh>1] = 1
        #r_inh = r_inh+dt/tau_i*(-r_inh+np.sum(r))
        drive = ff_drive-r_inh
        r = r+dt/tau*(-r+relu(drive))
    return r,r_inh

def decode_angle(r,angles):
    if np.sum(r)==0:
        return np.nan
    return np.angle(np.sum(r*np.exp(1j*angles)))

#%% Output activity profile for a range of inhibitory gains
inh_gains = np.linspace(0,3,20)
plt.figure()
for g in inh_gains:
    r,r_inh = run_network(ff_drive,g)
    plt.plot(angles,r,label=f'gain={g:.2f}')
    #plt.plot(angles,r_inh,label=f'gain={g:.2f}')
plt.plot(angles,ff_drive,'k--',label='FF drive')
plt.xlabel('Angle (rad)')
plt.ylabel('Output activity')
plt.legend(fontsize=6)
plt.title('Output profile vs inhibitory gain')

#%% Effect of inhibitory gain on output amplitude and decoded angle
inh_gains = np.linspace(0,2,40)
amps = np.zeros(len(inh_gains))
decoded = np.zeros(len(inh_gains))
for i,g in enumerate(inh_gains):
    r,r_inh = run_network(ff_drive,g)
    amps[i] = np.max(r)
    decoded[i] = decode_angle(r,angles)

plt.figure()
plt.subplot(2,1,1)
plt.plot(inh_gains,amps)
plt.ylabel('Peak output activity')

plt.subplot(2,1,2)
plt.plot(inh_gains,np.degrees(decoded))
plt.ylabel('Decoded angle (deg)')
plt.xlabel('Inhibitory gain')



#%% Change drive
inh_gain = .75
ffgains = np.linspace(0,1,20)
for g in ffgains:
    
    r,r_inh = run_network(ff_drive*g,inh_gain)
    plt.plot(angles,r,label=f'gain={g:.2f}')
    #plt.plot(angles,r_inh,label=f'gain={g:.2f}')
plt.plot(angles,ff_drive,'k--',label='FF drive')
plt.xlabel('Angle (rad)')
plt.ylabel('Output activity')
plt.legend(fontsize=6)
plt.title('Output profile vs inhibitory gain')


#%% Dual colour imaging phase and pseudotime plots


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



#%%

achoice=1
dchoices = ['hDeltaC_FC2_EPG']
#colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
colours = uplt.columnar_colours()
numpoints = 12  
plt.close('all')
grandphase = np.zeros((200,3,len(dchoices)))
grandtraj = np.zeros((200,2))
grandamp = np.zeros((200,2,3,len(dchoices)))
for i,d1 in enumerate(dchoices):
    data=all_flies[d1]
        
    phases = np.zeros((200,3,len(data)))
    trajs = np.zeros((200,2,len(data)))
    amps = np.zeros((200,2,3,len(data)))
    for d in data:
        print(d1)
        etp = data[d]
        phase,traj,amp = etp.trajectory_mean(time_threshold2=0,regions=['eb_ch1','fsb2_ch2','fsb1_ch1'])
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
        for ip in range(3):
            tphase = grandphase[p,ip,i]
            x = np.array([np.sin(tphase)*grandamp[p,achoice,ip,i],0])+offset
            y = np.array([np.cos(tphase)*grandamp[p,achoice,ip,i],0])
            if ip>0:
                
                plt.plot(x,y,color=colours[ip,:])
            else:
                plt.plot(x,y,color=[0.5,0.5,0.5])
            
                
                
g = plt.gca()
g.set_aspect('equal')
plt.fill([-.5,5.5,5.5,-.5],[-.5,-.5,.5,.5],color=[0.8,0.8,0.8],zorder=-1)


#%%
plt.close('all')
savedir = r'Y:\Applications\K99\Contents\PHS 398 Career Development Award Supplemental Form\Figures\Figure4'
scale=10
colours[0,:] = 0.0
timepoints = [.75,2.0,4.5]
t = np.linspace(0,1,100)
t  = np.append(t,np.linspace(t[-1]+t[1],5+t[1],100))

# Trajectory
plt.fill([-10,0,0,-10],[-50,-50,0,0],color=[.8,.8,.8],linewidth=0)
plt.fill([-13,-3,-3,-13],[50,50,0,0],color=[.8,.8,.8],linewidth=0)
plt.plot(grandtraj[:,0],grandtraj[:,1],color='k')
plt.gca().set_aspect('equal')

for r in timepoints:
    dx = ug.find_nearest(r,t)
    tamps = grandamp[dx,achoice,:,:].squeeze()*0.0+1
    tphase = grandphase[dx,:,:].squeeze()
    for i in range(3):
        tx = np.array([0,np.sin(tphase[i])])*(tamps[i]-.15)*scale
        ty = np.array([0,np.cos(tphase[i])])*(tamps[i]-.15)*scale
        
        plt.plot(tx+grandtraj[dx,0],ty+grandtraj[dx,1],color=colours[i,:])

plt.savefig(os.path.join(savedir,'Dual_col_traj.pdf'))

# Pseudotime
plt.figure()
plt.fill([0,0,1,1],[-np.pi,np.pi,np.pi,-np.pi],color=[.8,.8,.8])
for i in range(3):
    for i2 in range(phases.shape[2]):
        plt.scatter(t,phases[:,i,i2],color=colours[i,:],s=10,alpha=.25)
    plt.scatter(t,grandphase[:,i,0],color=colours[i,:],s=40)
    
for r in timepoints:
    plt.plot([r,r],[-np.pi,np.pi],color='k')

plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])

plt.savefig(os.path.join(savedir,'Dual_col_pseudo.pdf'))


















#%%
dchoices = ['hDeltaC_FC2_EPG']
#savedir = r'Y:\Presentations\2026\03_CSHL\Panels'
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































