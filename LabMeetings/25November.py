# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 10:39:34 2025

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
#%% Load data
datadirs_hdj = [
                r"Y:\Data\FCI\Hedwig\hDeltaJ\240529\f1\Trial3",#Good pointer
                r"Y:\Data\FCI\Hedwig\hDeltaJ\251011\f1\Trial1",
                r"Y:\Data\FCI\Hedwig\hDeltaJ\251022\f1\Trial2",
                r'Y:\Data\FCI\Hedwig\hDeltaJ\251121\f2\Trial1']

all_flies_hdj = {}
etp_hdj = {}
hdj_names = np.array([])

for i,datadir in enumerate(datadirs_hdj):
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    hdj_names = np.append(hdj_names,cxa.name)
    all_flies_hdj.update({str(i):cxa})
    etp = ET_paper(datadir)
    etp_hdj.update({str(i):etp})
    
    
datadirs_hdj_bad = [
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f1\Trial4",
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f2\Trial2",#Actually not that bad a tracker
    r"Y:\Data\FCI\Hedwig\hDeltaJ\251028\f1\Trial2",
                ]

all_flies_hdj_b = {}
etp_hdj_b = {}

for i,datadir in enumerate(datadirs_hdj_bad):
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    all_flies_hdj_b.update({str(i):cxa})
    etp = ET_paper(datadir)
    etp_hdj_b.update({str(i):etp})
    
datadirs_hdc = [
                r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250411\f1\Trial1",# Phase recording is not the best - strong pointer
                #r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f1\Trial3",# Not many jumps, weak pointer
                r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f2\Trial2", # Strong pointer
                r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250714\f1\Trial2",# Strong pointer, just like FC2
                r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f1\Trial1", # Points away ******
                r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f2\Trial3", # Strong pointer
                r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250721\f1\Trial2",# Strong pointer
                r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250907\f1\Trial2",# Good pointer
                ]

all_flies_hdc = {}
etp_hdc = {}
hdc_names = np.array([])
for i,datadir in enumerate(datadirs_hdc):
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    hdc_names = np.append(hdc_names,cxa.name)
    all_flies_hdc.update({str(i):cxa})
    etp = ET_paper(datadir)
    etp_hdc.update({str(i):etp})
datadirs_fc2 = [
    "Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f1\\Trial3",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240418\\f2\\Trial3",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240502\\f1\\Trial2",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\240514\\f1\\Trial2",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\241104\\f1\\Trial5",
"Y:\Data\FCI\Hedwig\FC2_maimon2\\241106\\f1\\Trial2",
r"Y:\Data\FCI\Hedwig\FC2_maimon2\250128\f1\Trial4"]
all_flies_fc2 = {}
etp_fc2 = {}
for i,datadir in enumerate(datadirs_fc2):
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    all_flies_fc2.update({str(i):cxa})
    etp = ET_paper(datadir)
    etp_fc2.update({str(i):etp})    
    
datadirs_hdc_bad = [r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250328\f1\Trial3",
                 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250313\f1\Trial2",
                 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250716\f1\Trial3"]
all_flies_hdc_b = {
                 }
for i,datadir in enumerate(datadirs_hdc_bad):
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    all_flies_hdc_b.update({str(i):cxa})
  
all_flies_hdc_pulses =[
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250320\f2\Trial5",#ACV pulses
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250320\f2\Trial4",#Octanol [?] pulses - neuron is inhibited
 
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250328\f1\Trial5",#'ACV
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250328\f1\Trial4",#'Octanol'
 
 
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f1\Trial6", # ACV pulses
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f1\Trial7", # octanol pulses
 
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250714\f1\Trial6", # ACV pulses
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250714\f1\Trial7", # Octanol pulses
 
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250716\f1\Trial6", # ACV pulses
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250716\f1\Trial5", #  oct pulses, imaging artefact in pre air period
 
 
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f1\Trial5",#ACV pulses
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f1\Trial7",#OCT pulses
 
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250721\f1\Trial6", # ACV pulses
 r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250721\f1\Trial5", # Oct pulses


 ]
  
all_flies_hdc_p = {
                 }
pulse_names = np.array([])
in_imaging = np.zeros(int(len(all_flies_hdc_pulses)/2))
for i,datadir in enumerate(all_flies_hdc_pulses[0:-1:2]):
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    all_flies_hdc_p.update({str(i)+'_ACV':cxa})
    pulse_names = np.append(pulse_names,cxa.name)
    
    name_check = [1 for name in hdc_names if name[:9] ==cxa.name[:9]]
    if len(name_check)>0:
        in_imaging[i] = name_check[0]
    
    datadir = all_flies_hdc_pulses[i*2+1]
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    all_flies_hdc_p.update({str(i)+'_Oct':cxa})
all_flies_hdc_p.update({'in_imaging': in_imaging})

all_flies_hdj_pulses =[
r"Y:\Data\FCI\Hedwig\hDeltaJ\251011\f1\Trial4", # ACV pulses
 r"Y:\Data\FCI\Hedwig\hDeltaJ\251011\f1\Trial5", # Oct pulses

r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f1\Trial6", # ACV pulses
r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f1\Trial5", # Oct pulses


r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f2\Trial5", # ACV pulses
r"Y:\Data\FCI\Hedwig\hDeltaJ\251023\f2\Trial4", # Oct pulses

r'Y:\Data\FCI\Hedwig\hDeltaJ\251121\f2\Trial3', # ACV pulses
r'Y:\Data\FCI\Hedwig\hDeltaJ\251121\f2\Trial4'# Oct pulses
 ]
  
all_flies_hdj_p = {
                 }
pulse_names = np.array([])
in_imaging = np.zeros(int(len(all_flies_hdj_pulses)/2))
for i,datadir in enumerate(all_flies_hdj_pulses[0:-1:2]):
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    all_flies_hdj_p.update({str(i)+'_ACV':cxa})
    pulse_names = np.append(pulse_names,cxa.name)
    
    name_check = [1 for name in hdj_names if name[:9] ==cxa.name[:9]]
    if len(name_check)>0:
        in_imaging[i] = name_check[0]
    
    datadir = all_flies_hdj_pulses[i*2+1]
    print(datadir)
    cxa = CX_a(datadir,regions=['eb','fsb_upper','fsb_lower'],denovo=False)
    all_flies_hdj_p.update({str(i)+'_Oct':cxa})
    
all_flies_hdj_p.update({'in_imaging': in_imaging})
#%% fit tau recovery  compare across cell types
plt.close('all')


dchoices  =['hdc','hdj']
offset = 0
for dchoice in dchoices:
    if dchoice=='hdc':
        data = all_flies_hdc
    elif dchoice=='hdj':
        data= all_flies_hdj
    output= np.zeros(len(data))
    for d in data:
        cxa = data[d]
        
        jumps = cxa.get_entries_exits_like_jumps()
        amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        from scipy.optimize import curve_fit
        def exp_rise(x,A,K,C):
            return A*(1-np.exp(K*x)) + C
        bounds = ([0,-10,-.01],
                  [np.max(amp),-0.01,.01])
        taus = np.zeros(len(jumps))
        for ij,j in enumerate(jumps):
            dx = np.arange(j[1],np.min([j[1]+100,j[2]]))
            if len(dx)<30:
                continue
            
            
            
            y = amp[dx]
            x = (dx-dx[0])/10
            
            # Clip x and y so that start is minimum
            ymn = np.argmin(y[:20])
            
            
            y = y[ymn:]
            x = x[ymn:]
            y = y-np.min(y) 
            x = x-x[0]
            
            
            p0 = [np.max(y)-np.min(y),-0.5,np.min(y)]
            try:
                popt,c = curve_fit(exp_rise,x,y,p0=p0,bounds=bounds)
            except:
                continue
            
            ypred = exp_rise(x,popt[0],popt[1],popt[2])
           
            taus[ij] =-1/ popt[1]
        print(np.median(taus[taus>0]))
        output[int(d)] = np.median(taus[taus>0])
    xplt =offset+np.random.rand(len(output))*.25
    plt.scatter(xplt,output,color='k')
    offset = offset+1
#%% Look at tau deviation for ACV versus octanol
#%% look at tau deviation across each experiment
plt.close('all')
cxa = all_flies_hdc['3']

jumps = cxa.get_entries_exits_like_jumps()
amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
amp = sg.savgol_filter(amp,20,3)
u = ug()
_,_,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
from scipy.optimize import curve_fit
def exp_rise(x,A,K,C):
    return A*(1-np.exp(K*x)) + C

# def exp_rise3(x,A,K,A2,K2):
#     return A*(1-np.exp(K*x)) + A2*(1-np.exp(K2*x))
bounds = ([0,-10,-.01],
          [np.max(amp),-0.01,.01])
bounds2 = ([0,     -10,   -.01,  -2,  -10],
           [2,       -0.01,  .01,   0,   -0.01])

def exp_rise2(x,A,K,C,A2,K2):
    return A*(1-np.exp(K*x)) + A2*(1-np.exp(K2*x)) + C


#bounds3= ([0,-np.inf,-np.inf,0],[np.inf,np.inf,np.inf,.1])
taus = np.zeros(len(jumps))
plt.figure()
offset = 0
for ij,j in enumerate(jumps):
    dx = np.arange(j[1],np.min([j[1]+100,j[2]]))
    if len(dx)<20:
        continue
    if len(dx)<50:
        dx = np.arange(j[1],j[2]+10)
    
    y = amp[dx]
    x = (dx-dx[0])/10
    v = vd[dx]
    
    # Clip x and y so that start is minimum
    ymn = np.argmin(y[1:20])
    
    
    y = y[ymn:]
    x = x[ymn:]
    v = v[ymn:]
    y = y-np.min(y) 
    x = x-x[0]
    p0 = [np.max(y)-np.min(y),-0.5,np.min(y)]
    p02 = [np.max(y)-np.min(y),-0.5,0,-1,-0.5]
    
    try:
        popt,c = curve_fit(exp_rise,x,y,p0=p0,bounds=bounds,maxfev=20000)
    except:
        plt.plot(x+offset,y,color='k')
        offset = offset+11
        
       
        print('fail')
        break
        continue
    
    ypred = exp_rise(x,popt[0],popt[1],popt[2])
    #ypred = exp_rise2(x,popt[0],popt[1],popt[2],popt[3],popt[4])
    taus[ij] =-1/ popt[1]
    if taus[ij]>100:
        break
    
    plt.plot(x+offset,y,color='k')
    plt.plot(x+offset,ypred,color='r')
    plt.plot(x+offset,v/20-.5,color=[0.5,0.5,0.5])
    plt.text(offset,.5,str(ij) +' Tau: '+ str(np.round(taus[ij],2)))
    #plt.title(str(ij) +' Tau: '+ str(taus[ij]))
    #plt.xlim([0,10])
    offset = offset+11
    print(popt)

print(np.median(taus[taus>0]))
#%%
taus[taus>100] = 0
gdx = np.where(taus>0)[0]
x = cxa.ft2['ft_posx'].to_numpy()
y = cxa.ft2['ft_posy'].to_numpy()
x,y = cxa.fictrac_repair(x,y)
ins = cxa.ft2['instrip'].to_numpy()
plt.figure()
plt.plot(x,y,color='k')
plt.scatter(x[ins>0],y[ins>0],color=[0.7,0.7,0.7])
for i,g in enumerate(gdx):
    t_tau = ug.round_to_sig_figs(taus[g],2)
    plt.text(x[jumps[g,1]],y[jumps[g,1]],str(t_tau))
    recp = int(t_tau*10)
    try :
        plt.scatter(x[jumps[g,1]+recp],y[jumps[g,1]+recp],color='r')
    except:
        print('end of line')
gtaus = taus[taus>0]
g= plt.gca()
g.set_aspect('equal')
#%% Transition behavioural raster
cxa = all_flies_hdj['3']
plt.figure()
jumps = cxa.get_jumps()
amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
heading = cxa.ft2['ft_heading'].to_numpy()*cxa.side*-1
_,_,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
stopping = vd<0.5
offset=0
for ij,j in enumerate(jumps):
    dx = np.arange(j[1],j[2])
    theading = heading[dx]
    tamp = amp[dx]
    tstop = stopping[dx]
    tphase = phase[dx]
    x = dx-dx[0]
    x= x/10
    y = np.ones_like(dx)-offset
    uplt.coloured_line_simple(x,y,theading,'coolwarm',cmin=-np.pi/2,cmax=np.pi/2,linewidth=5)
    uplt.coloured_line_simple(x,y-.1,tphase,'coolwarm',cmin=-np.pi,cmax=np.pi,linewidth=5)
    uplt.coloured_line_simple(x,y-.2,tamp,'viridis',cmin=0,cmax=0.75,linewidth=5)
    plt.scatter(x[tstop],y[tstop]-.3,color='k',s=20,zorder=5)
    offset= offset+.5
plt.xlim([0,30])
#%% Mean arrows overlaid
dchoices = ['hdj','hdc','fc2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255


grandphase = np.zeros((200,len(dchoices)))
grandtraj = np.zeros((200,2))
grandamp = np.zeros((200,2,len(dchoices)))
for i,d in enumerate(dchoices):
    if d=='hdj':
        data = etp_hdj
    elif d =='fc2':
        data = etp_fc2
    elif d=='hdc':
        data = etp_hdc
        
    phases = np.zeros((200,len(data)))
    trajs = np.zeros((200,2,len(data)))
    amps = np.zeros((200,2,len(data)))
    for d in data:
        etp = data[d]
        phase,traj,amp = etp.trajectory_mean(time_threshold2=10,regions=['fsb_upper'])
        print(traj.shape)
        phases[:,int(d)] = stats.circmean(phase[:,0,:],axis=1,low=-np.pi,high=np.pi)
        trajs[:,:,int(d)] = np.mean(traj,axis=2)
        amps[:,:,int(d)] = np.mean(amp[:,:,0,:],axis=2)
    pmean = stats.circmean(phases,axis=1,low=-np.pi,high=np.pi)
    tmean = np.mean(trajs,axis=2)
    grandphase[:,i] = pmean
    grandtraj = grandtraj+tmean*(1/len(dchoices))
    grandamp[:,:,i] =  np.mean(amps,axis=2)    
    
    
plt.plot(grandtraj[:,0],grandtraj[:,1],color='k')
pdx = np.linspace(0,199,20,dtype= 'int')

fx = np.array([-10,0,0,-10])
fy = np.array([-50,-50,0,0])
plt.fill(fx,fy,color=[0.7,0.7,.7])
plt.fill(fx-3,fy+50,color=[.7,.7,.7])
achoice=1
scale=15
for i in range(len(dchoices)):
    for p in pdx:
        tphase = grandphase[p,i]
        x = [grandtraj[p,0]+np.sin(tphase)*scale*grandamp[p,achoice,i],grandtraj[p,0]]
        y = [grandtraj[p,1]+np.cos(tphase)*scale*grandamp[p,achoice,i],grandtraj[p,1]]
        plt.plot(x,y,color=colours[i,:])

g = plt.gca()
g.set_aspect('equal')
plt.ylim([-40,40])
savedir= r'Y:\Presentations\2025\11_LabMeeting'
#plt.savefig(os.path.join(savedir,'hDeltaFC2CompTraj.pdf'))




#%% Example of hDeltaJ plume timeseries

savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'

hdeltaj_colour = np.array([81,61,204])/255
fly = 0
lims = [3500,5500]
dx = np.arange(lims[0],lims[1])
cxa = all_flies_hdj[str(fly)]
plt.figure(figsize=(9,3))
phase = cxa.pdat['phase_fsb_upper']
phase_eb = cxa.pdat['phase_eb']
ins = cxa.ft2['instrip'].to_numpy().astype(int)
u = ug()
_,_,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'].to_numpy())

x = np.arange(0,len(phase_eb))/10
plt.scatter(x[dx],phase[dx],color=hdeltaj_colour,s=3,zorder=10)
plt.scatter(x[dx],phase_eb[dx],color='k',s=3)
plt.scatter(x[dx],ug.circ_subtract(phase_eb[dx],np.pi),s=3,color=[0.6,0.6,0.6])

plt.fill_between(x[dx],np.pi*(2*ins[dx]-1),ins[dx]*0-np.pi,color='r',zorder=-1)
vd  =vd/np.std(vd)
vd = np.append(0,vd)
plt.plot(x[dx],vd[dx]/4-4.5,color=[0.3,.3,.3])
plt.xlim([350,550])
plt.xlabel('time (s)')
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
plt.ylim([-4.6,np.pi])
savename = 'Example_trace_'+str(fly)+'.pdf'
plt.savefig(os.path.join(savedir,savename))
#%% Example of hDeltaC plume timeseries
hdeltac_colour = np.array([49,99,125])/255
plt.figure(figsize=(9,3))
fly =4

cxa = all_flies_hdc[str(fly)]
dx = np.arange(8000,11000)

phase = cxa.pdat['phase_fsb_upper']
phase_eb = cxa.pdat['phase_eb']
dx = np.arange(0,len(phase))
#phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
#phase_eb =cxa.pdat['offset_eb_phase'].to_numpy()
heading_sign = -np.sign(cxa.ft2['ft_heading'].to_numpy())
#heading_sign= -np.sign(phase_eb)
phase_pred = ug.circ_subtract(phase_eb,3*(np.pi/4)*heading_sign)

ins = cxa.ft2['instrip'].to_numpy().astype(int)
u = ug()
_,_,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'].to_numpy())

x = np.arange(0,len(phase_eb))/10
x = x-x[dx[0]]
plt.scatter(x[dx],phase[dx],color=hdeltac_colour,s=3,zorder=10)
plt.scatter(x[dx],phase_eb[dx],color='k',s=3)
plt.scatter(x[dx],phase_pred[dx],s=3,color=[0.8,0.5,0.5])
plt.fill_between(x[dx],np.pi*(2*ins[dx]-1),ins[dx]*0-np.pi,color='r',zorder=-1)
vd  =vd/np.std(vd)
vd = np.append(0,vd)
plt.plot(x[dx],vd[dx]/4-4.5,color='k')
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdc'

plt.xlim([x[dx][0],x[dx][-1]])
plt.xlabel('time (s)')
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
plt.ylim([-4.6,np.pi])
savename = 'Example_trace_'+str(fly)+'.pdf'
#plt.savefig(os.path.join(savedir,savename))
#%% Plot hDeltaC mean fluorescence over experiment
hdeltac_colour = np.array([49,99,125])/255
plt.close('all')
for i in all_flies_hdc:
    plt.figure()
    cxa = all_flies_hdc[i]
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    ins = cxa.ft2['instrip'].to_numpy().astype(int)
    u = ug()
    _,_,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
    
    x = np.arange(0,len(amp))
    plt.plot(x,amp,color=hdeltac_colour,zorder=10)
    
    plt.fill_between(x,ins,ins*0,color='r',zorder=-1)
    vd  =vd/np.std(vd)
    vd = np.append(0,vd)
    plt.plot(x,vd/4-1.5,color=[0.3,0.3,0.3])
#%% Sign of phase offset with EPG phase
minlen = 10
plt.close('all')
for i in all_flies_hdc:
    cxa = all_flies_hdc[i]
    phase = cxa.pdat['phase_fsb_upper']
    phase_eb = cxa.pdat['phase_eb']
    heading = cxa.ft2['ft_heading'].to_numpy()
    heading = ug.savgol_circ(heading,20,3)
    pdiff = ug.circ_subtract(phase,phase_eb)
    
    jumps = cxa.get_jumps()
    
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    plt.figure()
    for ij,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
        plt.scatter(pdiff[dx],heading[dx],color='k',s=1,alpha=.2)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
#%%
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
minlen = 5
nbins = 100
hist_sum = None
hist_count = 0

# Predefine bin edges (same for all histograms)
xedges = np.linspace(-np.pi, np.pi, nbins+1)   # new X axis = heading
yedges = np.linspace(-np.pi, np.pi, nbins+1)   # new Y axis = pdiff
fed = 135
fed2 = 45
datlab = 'fc2'
if datlab=='fc2':
    data = all_flies_fc2
elif datlab=='hdc':
    data = all_flies_hdc
elif datlab =='hdj':
    data = all_flies_hdj
    

for i in data:

    cxa = data[i]
    
    
    # jumps
    jumps = cxa.get_jumps()
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2] - jumps[:,1]
    jumps = jumps[jumplen > minlen]
    
    phase = cxa.pdat['phase_fsb_upper']*-cxa.side
    phase_eb = cxa.pdat['phase_eb']*-cxa.side
    heading = ug.savgol_circ(cxa.ft2['ft_heading'].to_numpy(), 20, 3)*-cxa.side

    pdiff = ug.circ_subtract(phase, phase_eb)

    

    pdiff_all = []
    heading_all = []

    for j in jumps:
        dx = np.arange(j[1], j[2])
        pdiff_all.append(pdiff[dx])
        heading_all.append(heading[dx])

    if len(pdiff_all) == 0:
        continue

    pdiff_all = np.concatenate(pdiff_all)
    heading_all = np.concatenate(heading_all)

    # -------------------------------
    #     FLIPPED AXIS HISTOGRAM
    # -------------------------------
    # NOTE: we swap the order: (x = heading, y = pdiff)
    H, _, _ = np.histogram2d(
        heading_all,   # x-axis
        pdiff_all,     # y-axis
        bins=[xedges, yedges]
    )

    # Plot (no transpose needed now — histogram2d aligns directly)
    plt.figure(figsize=(5,5))
    plt.imshow(
        H.T,
        origin='lower',
        extent=[-np.pi, np.pi, -np.pi, np.pi],
        aspect='auto',
        cmap='viridis'
    )
    plt.title(f"Fly {i} — 2D histogram (axes flipped)")
    plt.xlabel("heading")
    plt.ylabel("pdiff")
    plt.colorbar(label="Counts")

    # Accumulate
    if hist_sum is None:
        hist_sum = H.astype(float)
    else:
        hist_sum += H
    hist_count += 1


# --------------------------------------
#            MEAN HISTOGRAM
# --------------------------------------

if hist_count > 0:
    Hmean = hist_sum / hist_count

    plt.figure(figsize=(6,6))
    plt.imshow(
        Hmean.T*100,
        origin='lower',
        extent=[-180, 180, -180, 180],
        aspect='equal',
        cmap='viridis'
    )
   # plt.title(f"Mean histogram across {hist_count} flies (axes flipped)")
    plt.xlabel("heading (deg)",fontsize=15)
    plt.ylabel("FSB-EB phase (deg)",fontsize=15)
    plt.colorbar(label="Mean count")
    plt.plot([-180,180],[fed,fed],color='r')
    plt.plot([-180,180],[-fed,-fed],color='r')
    
    plt.plot([-180,180],[fed2,fed2],color='r')
    plt.plot([-180,180],[-fed2,-fed2],color='r')
    plt.plot([0,0],[-180,180],color='w')
    plt.plot([-180,180],[0,0],color='w')
    plt.xticks(np.arange(-180,200,90),fontsize=15)
    plt.yticks(np.arange(-180,200,90),fontsize=15)
    plt.subplots_adjust(left=.2)
    plt.savefig(os.path.join(savedir,'HistogramOffset_'+datlab+ '_.png'))
    plt.savefig(os.path.join(savedir,'HistogramOffset_'+datlab+'_.pdf'))
#%% Phase progression upon plume exit
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
fc2_colour = np.array([81,156,205])/255
hdeltaj_colour = np.array([81,61,204])/255
hdeltac_colour = np.array([49,99,125])/255

plt.close('all')
minlen = 50
plt.figure()
for ib in range(3):
    if ib==0:
        data = all_flies_fc2
        colour = fc2_colour
    elif ib==1:
        data = all_flies_hdc
        colour = hdeltac_colour
    elif ib==2:
        data = all_flies_hdj
        colour = hdeltaj_colour
    meandat = np.zeros((600,len(data)))
    x = np.arange(0,600)/10
    for io in data:
       # plt.figure()
        cxa = data[io]
        amp = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
        jumps = cxa.get_jumps()
        amp = amp*-cxa.side
        jumps = cxa.get_entries_exits_like_jumps()
        jumplen = jumps[:,2]-jumps[:,1]
        jumps = jumps[jumplen>minlen]
        all_js = np.zeros((600,len(jumps)))
        for i,j in enumerate(jumps):
            dx = np.arange(j[1],j[2])
           # plt.plot(amp[dx],color='k',alpha=0.3)
            tamp = amp[dx]
            amplen  = np.min([600,len(tamp)])
            all_js[:amplen,i] = tamp[:amplen]
            
        all_js[all_js==0] = np.nan    
        meandat[:,int(io)] =stats.circmean(all_js,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    pmean = stats.circmean(meandat,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    
    for io in data:
        plt.scatter(x,meandat[:,int(io)],color=colour,alpha=0.3,s=5)
    plt.plot(x,pmean,color=colour,alpha=1,linewidth=3)
plt.xlim([0,20])
plt.ylim([-np.pi,np.pi])

plt.xlabel('time from plume exit (s)',fontsize=15)
plt.ylabel('allocentric phase (deg)',fontsize=15)
plt.title('phase',fontsize=15)
plt.yticks(np.arange(-np.pi,np.pi+.1,np.pi/2),fontsize=15,labels=np.arange(-180,181,90))
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.plot([0,20],[0,0],color='k',linestyle='--')
plt.plot([0,20],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
plt.savefig(os.path.join(savedir,'hDeltaC_hDeltaJ_phase_exit_comp.png'))

#%% Phase progression pulses
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
fc2_colour = np.array([81,156,205])/255
hdeltaj_colour = np.array([81,61,204])/255
hdeltac_colour = np.array([49,99,125])/255

plt.close('all')
minlen = 50
plt.figure()
for ib in range(3):
   
    if ib==0:
        data = all_flies_hdc_p
        dlen = int(len(all_flies_hdc_pulses)/2)
        odour = 'ACV'
        colour = hdeltac_colour
    elif ib==1:
        data = all_flies_hdj_p
        dlen = int(len(all_flies_hdj_pulses)/2)
        odour = 'ACV'
        colour = hdeltaj_colour
    meandat = np.zeros((600,dlen))
    x = np.arange(0,600)/10
    for io in range(dlen):
        print(io)
       # plt.figure()
        cxa = data[str(io)+'_'+odour]
        amp = ug.circ_subtract(cxa.pdat['offset_fsb_upper_phase'].to_numpy(),cxa.pdat['offset_eb_phase'].to_numpy())
        #jumps = cxa.get_jumps()
        #amp = amp*-cxa.side
        jumps = cxa.get_entries_exits_like_jumps()
        #jumplen = jumps[:,2]-jumps[:,1]
        #jumps = jumps[jumplen>minlen]
        all_js = np.zeros((600,len(jumps)))
        for i,j in enumerate(jumps):
            dx = np.arange(j[1],j[2])
           # plt.plot(amp[dx],color='k',alpha=0.3)
            tamp = amp[dx]
            amplen  = np.min([600,len(tamp)])
            all_js[:amplen,i] = tamp[:amplen]
            
        all_js[all_js==0] = np.nan    
        meandat[:,io] =stats.circmean(all_js,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    pmean = stats.circmean(meandat,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    
    for io in range(dlen):
        plt.scatter(x,meandat[:,io],color=colour,alpha=0.3,s=5)
    plt.plot(x,pmean,color=colour,alpha=1,linewidth=3)
plt.xlim([0,20])
plt.ylim([-np.pi,np.pi])

plt.xlabel('time from plume exit (s)',fontsize=15)
plt.ylabel('allocentric phase (deg)',fontsize=15)
plt.title('phase',fontsize=15)
plt.yticks(np.arange(-np.pi,np.pi+.1,np.pi/2),fontsize=15,labels=np.arange(-180,181,90))
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.plot([0,20],[0,0],color='k',linestyle='--')
plt.plot([0,20],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
#plt.savefig(os.path.join(savedir,'hDeltaC_hDeltaJ_phase_exit_comp.png'))


#%%
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
fc2_colour = np.array([81,156,205])/255
hdeltaj_colour = np.array([81,61,204])/255
hdeltac_colour = np.array([49,99,125])/255

plt.close('all')
minlen = 50
plt.figure()
for ib in range(3):
    if ib==0:
        data = all_flies_fc2
        colour = fc2_colour
    elif ib==1:
        data = all_flies_hdc
        colour = hdeltac_colour
    elif ib==2:
        data = all_flies_hdj
        colour = hdeltaj_colour
    meandat = np.zeros((600,len(data)))
    x = np.arange(0,600)/10
    for io in data:
       # plt.figure()
        cxa = data[io]
        amp = ug.circ_subtract(cxa.pdat['offset_fsb_upper_phase'].to_numpy(),cxa.pdat['offset_eb_phase'].to_numpy())
        jumps = cxa.get_jumps()
        amp = amp*-cxa.side
        jumps = cxa.get_entries_exits_like_jumps()
        jumplen = jumps[:,2]-jumps[:,1]
        jumps = jumps[jumplen>minlen]
        all_js = np.zeros((600,len(jumps)))
        for i,j in enumerate(jumps):
            dx = np.arange(j[1],j[2])
           # plt.plot(amp[dx],color='k',alpha=0.3)
            tamp = amp[dx]
            amplen  = np.min([600,len(tamp)])
            all_js[:amplen,i] = tamp[:amplen]
            
        all_js[all_js==0] = np.nan    
        meandat[:,int(io)] =stats.circmean(all_js,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    pmean = stats.circmean(meandat,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    
    for io in data:
        plt.scatter(x,meandat[:,int(io)],color=colour,alpha=0.3,s=5)
    plt.plot(x,pmean,color=colour,alpha=1,linewidth=3)
plt.xlim([0,20])
plt.ylim([-np.pi,np.pi])

plt.xlabel('time from plume exit (s)',fontsize=15)
plt.ylabel('fsb-eb phase (deg)',fontsize=15)
plt.title('phase',fontsize=15)
plt.yticks(np.arange(-np.pi,np.pi+.1,np.pi/2),fontsize=15,labels=np.arange(-180,181,90))
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.plot([0,20],[0,0],color='k',linestyle='--')
plt.plot([0,20],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
plt.savefig(os.path.join(savedir,'hDeltaC_hDeltaJ_fsb_eb_phase_exit_comp.png'))
#%% phase -progression good and bad trackers
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
fc2_colour = np.array([81,156,205])/255
hdeltaj_colour = np.array([81,61,204])/255
hdeltac_colour = np.array([49,99,125])/255

plt.close('all')
minlen = 50
plt.figure()
for ib in range(2):
    if ib==0:
        data = all_flies_hdj
        colour = hdeltaj_colour
    elif ib==1:
        data = all_flies_hdj_b
        colour = 'k'
    
    meandat = np.zeros((600,len(data)))
    x = np.arange(0,600)/10
    for io in data:
       # plt.figure()
        cxa = data[io]
        amp = ug.circ_subtract(cxa.pdat['offset_fsb_upper_phase'].to_numpy(),cxa.pdat['offset_eb_phase'].to_numpy())
       # amp = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
        try:   
            jumps = cxa.get_jumps()
        
            amp = amp*-cxa.side
        except:
            jumps = cxa.get_entries_exits_like_jumps()
            xf = cxa.ft2['ft_posx'].to_numpy()
            xj = xf[jumps[:,0]]
            xj = xj-xj[0]
            amp = amp*np.sign(xj[1])*-1
        jumps = cxa.get_entries_exits_like_jumps()
        jumplen = jumps[:,2]-jumps[:,1]
        jumps = jumps[jumplen>minlen]
        all_js = np.zeros((600,len(jumps)))
        for i,j in enumerate(jumps):
            dx = np.arange(j[1],j[2])
           # plt.plot(amp[dx],color='k',alpha=0.3)
            tamp = amp[dx]
            amplen  = np.min([600,len(tamp)])
            all_js[:amplen,i] = tamp[:amplen]
            
        all_js[all_js==0] = np.nan    
        meandat[:,int(io)] =stats.circmean(all_js,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    pmean = stats.circmean(meandat,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    
    for io in data:
        plt.scatter(x,meandat[:,int(io)],color=colour,alpha=0.3,s=5)
    plt.plot(x,pmean,color=colour,alpha=1,linewidth=3)
plt.xlim([0,20])
plt.ylim([-np.pi,np.pi])

plt.xlabel('time from plume exit (s)',fontsize=15)
plt.ylabel('fsb-eb phase (deg)',fontsize=15)
plt.title('phase',fontsize=15)
plt.yticks(np.arange(-np.pi,np.pi+.1,np.pi/2),fontsize=15,labels=np.arange(-180,181,90))
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.plot([0,20],[0,0],color='k',linestyle='--')
plt.plot([0,20],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
#%%
minlen = 50
meandat = np.zeros((600,len(all_flies_hdc)))
x = np.arange(0,600)/10
for io in all_flies_hdc:
   # plt.figure()
    cxa = all_flies_hdc[io]
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    all_js = np.zeros((600,len(jumps)))
    for i,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        all_js[:amplen,i] = tamp[:amplen]
        
    all_js[all_js==0] = np.nan    
    meandat[:,int(io)] = np.nanmean(all_js,axis=1)

plt.plot(x,meandat,color=hdeltac_colour,alpha=0.3)
plt.plot(x,np.mean(meandat,1),color=hdeltac_colour,alpha=1,linewidth=3)

meandat = np.zeros((600,len(all_flies_hdj)))
for io in all_flies_hdj:
   # plt.figure()
    cxa = all_flies_hdj[io]
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    all_js = np.zeros((600,len(jumps)))
    for i,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        all_js[:amplen,i] = tamp[:amplen]
        
    all_js[all_js==0] = np.nan    
    meandat[:,int(io)] = np.nanmean(all_js,axis=1)
    
plt.plot(x,meandat,color=hdeltaj_colour,alpha=0.3)
plt.plot(x,np.mean(meandat,axis=1),color=hdeltaj_colour,alpha=1,linewidth=3)
plt.xlim([0,20])
plt.ylim([0,0.8])
plt.xlabel('time from plume exit (s)',fontsize=15)
plt.ylabel('dF/F0',fontsize=15)
plt.title('mean fluorescence',fontsize=15)
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(0,1,.2),fontsize=15)
#%% Fluoresence recovery post plume exit comparison
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
fc2_colour = np.array([81,156,205])/255
hdeltaj_colour = np.array([81,61,204])/255
hdeltac_colour = np.array([49,99,125])/255

plt.close('all')
minlen = 50
meandat = np.zeros((600,len(all_flies_fc2)))
x = np.arange(0,600)/10
for io in all_flies_fc2:
   # plt.figure()
    cxa = all_flies_fc2[io]
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    all_js = np.zeros((600,len(jumps)))
    for i,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        all_js[:amplen,i] = tamp[:amplen]
        
    all_js[all_js==0] = np.nan    
    meandat[:,int(io)] = np.nanmean(all_js,axis=1)

plt.figure()
plt.plot(x,meandat,color=fc2_colour,alpha=0.3)
plt.plot(x,np.mean(meandat,1),color=fc2_colour,alpha=1,linewidth=3)


minlen = 50
meandat = np.zeros((600,len(all_flies_hdc)))
x = np.arange(0,600)/10
for io in all_flies_hdc:
   # plt.figure()
    cxa = all_flies_hdc[io]
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    all_js = np.zeros((600,len(jumps)))
    for i,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        all_js[:amplen,i] = tamp[:amplen]
        
    all_js[all_js==0] = np.nan    
    meandat[:,int(io)] = np.nanmean(all_js,axis=1)

plt.plot(x,meandat,color=hdeltac_colour,alpha=0.3)
plt.plot(x,np.mean(meandat,1),color=hdeltac_colour,alpha=1,linewidth=3)

meandat = np.zeros((600,len(all_flies_hdj)))
for io in all_flies_hdj:
   # plt.figure()
    cxa = all_flies_hdj[io]
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    all_js = np.zeros((600,len(jumps)))
    for i,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        all_js[:amplen,i] = tamp[:amplen]
        
    all_js[all_js==0] = np.nan    
    meandat[:,int(io)] = np.nanmean(all_js,axis=1)
    
plt.plot(x,meandat,color=hdeltaj_colour,alpha=0.3)
plt.plot(x,np.mean(meandat,axis=1),color=hdeltaj_colour,alpha=1,linewidth=3)
plt.xlim([0,20])
plt.ylim([0,0.8])
plt.xlabel('time from plume exit (s)',fontsize=15)
plt.ylabel('dF/F0',fontsize=15)
plt.title('mean fluorescence',fontsize=15)
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(0,1,.2),fontsize=15)
plt.savefig(os.path.join(savedir,'hDeltaC_hDeltaJ_meanFluor_comp.png'))
#%% Fluorescence recovery pulses
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
fc2_colour = np.array([81,156,205])/255
hdeltaj_colour = np.array([81,61,204])/255
hdeltac_colour = np.array([49,99,125])/255

plt.close('all')
minlen = 50
plt.figure()
odour = 'ACV'
for ib in range(3):
   
    if ib==0:
        data = all_flies_hdc_p
        dlen = int(len(all_flies_hdc_pulses)/2)
        
        colour = hdeltac_colour
    elif ib==1:
        data = all_flies_hdj_p
        dlen = int(len(all_flies_hdj_pulses)/2)
        
        colour = hdeltaj_colour
    meandat = np.zeros((600,dlen))
    x = np.arange(0,600)/10
    for io in range(dlen):
        print(io)
       # plt.figure()
        cxa = data[str(io)+'_'+odour]
        amp=np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        #amp = ug.circ_subtract(cxa.pdat['offset_fsb_upper_phase'].to_numpy(),cxa.pdat['offset_eb_phase'].to_numpy())
        #jumps = cxa.get_jumps()
        #amp = amp*-cxa.side
        jumps = cxa.get_entries_exits_like_jumps(odour=odour)
        #jumplen = jumps[:,2]-jumps[:,1]
        #jumps = jumps[jumplen>minlen]
        all_js = np.zeros((600,len(jumps)))
        for i,j in enumerate(jumps):
            dx = np.arange(j[1],j[2])
           # plt.plot(amp[dx],color='k',alpha=0.3)
            tamp = amp[dx]
            amplen  = np.min([600,len(tamp)])
            all_js[:amplen,i] = tamp[:amplen]
            
        all_js[all_js==0] = np.nan    
        meandat[:,io] =stats.circmean(all_js,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    pmean = stats.circmean(meandat,axis=1,nan_policy='omit',low=-np.pi,high=np.pi)
    
    for io in range(dlen):
        plt.plot(x,meandat[:,io],color=colour,alpha=0.3)
    plt.plot(x,pmean,color=colour,alpha=1,linewidth=3)
plt.xlim([0,20])
#plt.ylim([-np.pi,np.pi])

plt.xlim([0,20])
plt.ylim([0,0.8])
plt.xlabel('time from pulse end (s)',fontsize=15)
plt.ylabel('dF/F0',fontsize=15)
plt.title('mean fluorescence',fontsize=15)
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(0,1,.2),fontsize=15)
#plt.plot([0,20],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
plt.savefig(os.path.join(savedir,'hDeltaC_hDeltaJ_meanFluor_pulses_'+ odour+ '.png'))
#%%
#%% 
hdeltaj_colour = np.array([81,61,204])/255
hdeltac_colour = np.array([49,99,125])/255
plt.close('all')

plt.close('all')
minlen = 50
meandat = np.zeros((600,len(all_flies_fc2)))
x = np.arange(0,600)/10
for io in all_flies_fc2:
   # plt.figure()
    cxa = all_flies_fc2[io]
    amp = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    all_js = np.zeros((600,len(jumps)))
    for i,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        all_js[:amplen,i] = tamp[:amplen]
        
    all_js[all_js==0] = np.nan    
    meandat[:,int(io)] = np.nanmean(all_js,axis=1)

plt.figure()
plt.plot(x,meandat,color=fc2_colour,alpha=0.3)
plt.plot(x,np.mean(meandat,1),color=fc2_colour,alpha=1,linewidth=3)



minlen = 50
meandat = np.zeros((600,len(all_flies_hdc)))
x = np.arange(0,600)/10
for io in all_flies_hdc:
   # plt.figure()
    cxa = all_flies_hdc[io]
    #amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    amp = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    all_js = np.zeros((600,len(jumps)))
    for i,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        all_js[:amplen,i] = tamp[:amplen]
        
    all_js[all_js==0] = np.nan    
    meandat[:,int(io)] = np.nanmean(all_js,axis=1)

plt.plot(x,meandat,color=hdeltac_colour,alpha=0.5)
plt.plot(x,np.mean(meandat,1),color=hdeltac_colour,alpha=1,linewidth=2)

meandat = np.zeros((600,len(all_flies_hdj)))
for io in all_flies_hdj:
   # plt.figure()
    cxa = all_flies_hdj[io]
    #amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    amp = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    jumps = cxa.get_entries_exits_like_jumps()
    jumplen = jumps[:,2]-jumps[:,1]
    jumps = jumps[jumplen>minlen]
    all_js = np.zeros((600,len(jumps)))
    for i,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
       # plt.plot(amp[dx],color='k',alpha=0.3)
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        all_js[:amplen,i] = tamp[:amplen]
        
    all_js[all_js==0] = np.nan    
    meandat[:,int(io)] = np.nanmean(all_js,axis=1)
    
plt.plot(x,meandat,color=hdeltaj_colour,alpha=0.5)
plt.plot(x,np.mean(meandat,axis=1),color=hdeltaj_colour,alpha=1,linewidth=2)
plt.xlim([0,20])
plt.xlabel('time from plume exit (s)',fontsize=15)
plt.ylabel('pva',fontsize=15)
plt.title('pva',fontsize=15)

plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(0,.16,.05),fontsize=15)
plt.ylim([0.0,0.16])
plt.savefig(os.path.join(savedir,'hDeltaC_hDeltaJ_pva_comp.png'))

#%% Good and bad tracker example section
fc2_colour = np.array([81,156,205])/255
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
eg_fly =0
etp = etp_hdj_b[str(eg_fly)]
plt.close('all')
#Full trajectory

tstart = 60*5
tend = 60*10
tstart = 0
tend = 1000000
etp.example_trajectory(tstart,tend,boxes=([[60*5,60*15],[tstart,tend]]))
plt.ylim([-500,2500])
#plt.savefig(os.path.join(savedir,'Eg_traj_Fly_bad'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
phase = etp.cxa.pdat['offset_fsb_upper_phase'].to_numpy()
phase_eb = etp.cxa.pdat['offset_eb_phase'].to_numpy()

phase = etp.cxa.pdat['phase_fsb_upper']
phase_eb = etp.cxa.pdat['phase_eb']
heading = etp.cxa.ft2['ft_heading'].to_numpy()
ins = etp.cxa.ft2['instrip'].to_numpy().astype('float')
offset = ug.circ_subtract(phase_eb,heading)

#phase = ug.circ_subtract(phase,offset)
#phase_eb = ug.circ_subtract(phase_eb,offset)


tt = etp.cxa.pv2['relative_time'].to_numpy()
tdx = np.logical_and( tt>tstart,tt<tend)
t = tt[tdx]
tp = phase[tdx]
teb = phase_eb[tdx]
tins = ins[tdx]
tinsplot = tins*.5
offset_mn = stats.circmean(offset[tdx])
teb = ug.circ_subtract(teb,offset_mn)
tp = ug.circ_subtract(tp,offset_mn)


plt.figure()
plt.scatter(t,teb,color='k',s=5)
plt.scatter(t,ug.circ_subtract(teb,np.pi),color=[0.5,0.5,0.5],s=5)
plt.scatter(t,tp,color=fc2_colour,s=5)
plt.fill_between(t,np.zeros_like(tins)-np.pi-.5,tinsplot-np.pi-.5,color='r',linewidth=0,zorder=-1)


etp.example_heatmaps(tstart,tend)

#%%
eg_fly =0
etp = etp_hdj[str(eg_fly)]
#Full trajectory with boxes
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
tstart = 0
tend = 1000000
etp.example_trajectory(tstart,tend)
plt.ylim([-500,2500])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly_good'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))   
  
    
#%% Good and bad trackers
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
eg_fly =1
etp = etp_hdj_b[str(eg_fly)]
plt.close('all')
#Full trajectory
tstart = 0
tend = 1000000
etp.example_trajectory(tstart,tend)
plt.ylim([-500,2500])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly_bad'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
    


eg_fly =2
etp = etp_hdj[str(eg_fly)]
#Full trajectory with boxes
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
tstart = 0
tend = 1000000
etp.example_trajectory(tstart,tend)
plt.ylim([-500,2500])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly_good'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

#%% Example hdeltaC trajectories and heat
eg_fly = 5
etp = etp_hdc[str(eg_fly)]
plt.close('all')
#Full trajectory with boxes
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdc'
tstart = 0
tend = 1000000
etp.example_trajectory(tstart,tend,boxes=np.array([[60,110],[60*20.5,60*23]]))
plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Amenotaxis trajectory
plt.close('all')
tstart = 20
tend = 80
etp.example_trajectory(tstart,tend)
plt.xlim([-250,-100])
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
tstart = 60*20.5
tend = 60*23

etp.example_trajectory(tstart,tend)
plt.xlim([-75,75])
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

    
#%% Example hdeltaJ trajectories and heat
eg_fly =0
etp = etp_hdj[str(eg_fly)]
plt.close('all')
#Full trajectory with boxes
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
tstart = 0
tend = 1000000
etp.example_trajectory(tstart,tend,boxes=np.array([[60,110],[60*13.65,60*16.15]]))
plt.savefig(os.path.join(savedir,'Eg_traj_Fly'+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Amenotaxis trajectory
plt.close('all')
tstart = 20
tend = 80
etp.example_trajectory(tstart,tend)
plt.xlim([-80,0])

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
tstart = 60*13.65
tend = 60*16.15

etp.example_trajectory(tstart,tend)
plt.xlim([0,80])
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
    
    
    
    
    
    
    
    
    
    
    
#%% mean activity traces for jumps of good trackers
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
for i in all_flies_hdj:
    plt.figure()
    cxa = all_flies_hdj[i]
    pmean = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    pmnorm = (pmean-np.mean(pmean))/np.std(pmean)
    #cxa.mean_jump_heat(regions=['eb','fsb_upper'])
    fc = fci_regmodel(pmnorm,cxa.ft2,cxa.pv2)
    fc.mean_traj_heat_jump(pmnorm,set_cmx=True,cmx=1,set_cmn=True,cmn=-1)
    fc.mean_traj_jump_indi(pmnorm,-1,1)
    plt.ylim([-100,100])
    plt.savefig(os.path.join(savedir,'MeanAvFluor_'+cxa.name+'.pdf'))
    #fc.example_trajectory_jump(pmnorm,cxa.ft,cmin=0,cmax=1,xcent= 0,pw=5,jsize=3,cmap='coolwarm',selection=[])
#%%
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdc'
for i in all_flies_hdc:
    plt.figure()
    cxa = all_flies_hdc[i]
    pmean = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    pmnorm = (pmean-np.mean(pmean))/np.std(pmean)
    #cxa.mean_jump_heat(regions=['eb','fsb_upper'])
    fc = fci_regmodel(pmnorm,cxa.ft2,cxa.pv2)
    fc.mean_traj_heat_jump(pmnorm,set_cmx=True,cmx=1,set_cmn=True,cmn=-1)
    fc.mean_traj_jump_indi(pmnorm,-1,1)
    plt.ylim([-100,100])
    plt.savefig(os.path.join(savedir,'MeanAvFluor_'+cxa.name+'.pdf'))
    
#%% Plot phase when bump 'comes back' versus prior entry angle/epg phase
# EPG at 'bump off' versus fsb and bump on
plt.close('all')
ang = np.linspace(-np.pi,np.pi,num=100)
angs = np.sin(ang)*2
angc = np.cos(ang)*2
for i in all_flies_hdc:
    cxa = all_flies_hdc[i]
    pva = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    pva = sg.savgol_filter(pva,20,3)
    pva = (pva)/np.std(pva)
    pva = pva-np.min(pva)
    
    jumps = cxa.get_jumps(time_threshold=180)
    phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*-cxa.side
    phase = ug.savgol_circ(phase,20,3)
    phase_eb = cxa.pdat['offset_eb_phase'].to_numpy()*-cxa.side
    #phase = ug.circ_subtract(phase,phase_eb)
    offset = 0
    gridnum = np.ceil(np.sqrt(len(jumps))).astype(int)
    xoff = np.tile(np.arange(0,gridnum),gridnum)*10
    yoff =- xoff[np.argsort(xoff)]
    for i2,j in enumerate(jumps):
        plt.figure(int(i))
        dxpre = np.arange(j[0],j[1])
        dxpost = np.arange(j[1],j[2])
        x1 = np.arange(0,len(dxpre))
        x2 = np.arange(len(dxpre),len(dxpre)+len(dxpost))

        
        tpva = pva[dxpost]
        tphase = phase[dxpost]
        bstart,bsize = ug.find_blocks(tpva>.5)
        bstart = bstart[bsize>=2]
        thead = stats.circmean(phase_eb[dxpre[:5]],low=-np.pi,high=np.pi)
        head_ex = stats.circmean(phase_eb[dxpost[:5]],low=-np.pi,high=np.pi)
        
        


        xp = np.sin(phase[dxpost])*tpva+xoff[i2]
        yp = np.cos(phase[dxpost])*tpva+yoff[i2]
        plt.plot(angs+xoff[i2],angc+yoff[i2],color='k',zorder=-1,linewidth=0.5)
        plt.scatter(xp,yp,c=x2,s=10,alpha=0.5)
        plt.plot(np.array([0,0])+xoff[i2],np.array([-3,3])+yoff[i2],color='k',linestyle='--',zorder=-1)
        plt.plot(np.array([-3,3])+xoff[i2],np.array([0,0])+yoff[i2],color='k',linestyle='--',zorder=-1)
        xp = np.sin(thead)*3+xoff[i2] 
        yp = np.cos(thead)*3+yoff[i2]
        
        plt.plot([xoff[i2],xp],[yoff[i2],yp],color='r',linestyle='-')
        
        xp = np.sin(head_ex)*3+xoff[i2] 
        yp = np.cos(head_ex)*3+yoff[i2]
        
        plt.plot([xoff[i2],xp],[yoff[i2],yp],color=[0,0.5,1],linestyle='-')

    g=plt.gca()
    g.set_aspect('equal')
#%% 
plt.close('all')
ang = np.linspace(-np.pi,np.pi,num=100)
angs = np.sin(ang)*2
angc = np.cos(ang)*2
u = ug()
alldat = np.array([])
for i in all_flies_hdc:
    cxa = all_flies_hdc[i]
    pva = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    pva = sg.savgol_filter(pva,20,3)
    pva = (pva)/np.std(pva)
    pva = pva-np.min(pva)
    
    mnfluor = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    mnfluor = sg.savgol_filter(mnfluor,20,3)
    mnfluor =mnfluor/np.std(mnfluor)
    mnfluor = mnfluor-np.min(mnfluor)
    
    jumps = cxa.get_jumps(time_threshold=180) 
    #jumps = cxa.get_entries_exits_like_jumps()
    # jumplen = jumps[:,2]-jumps[:,1]
    # jumps = jumps[jumplen>600]
    phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*-cxa.side
    phase = ug.savgol_circ(phase,20,3)
    phase_eb = cxa.pdat['offset_eb_phase'].to_numpy()*-cxa.side
    _,_,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
    xf = cxa.ft2['ft_posx'].to_numpy()
    yf = cxa.ft2['ft_posy'].to_numpy()
    xf,yf = cxa.fictrac_repair(xf,yf)
    xf = xf*-cxa.side
    vd = np.append(0,vd)
    vd = vd/np.std(vd)
    #phase = ug.circ_subtract(phase,phase_eb)
    offset = 0
    gridnum = np.ceil(np.sqrt(len(jumps))).astype(int)
    xoff = np.tile(np.arange(0,gridnum),gridnum)*10
    yoff =- xoff[np.argsort(xoff)]
    for i2,j in enumerate(jumps):
        
        dxpre = np.arange(j[0],j[1])
        dxpost = np.arange(j[1],j[2])
        x1 = np.arange(0,len(dxpre))
        x2 = np.arange(len(dxpre),len(dxpre)+len(dxpost))

        tpva = pva[dxpost]
        tmf = mnfluor[dxpost]
        tphase = phase[dxpost]
        bstart,bsize = ug.find_blocks(tpva>.5)
        bstart = bstart[bsize>=2]
        thead = stats.circmean(phase_eb[dxpre[:5]],low=-np.pi,high=np.pi)
        head_ex = stats.circmean(phase_eb[dxpost[:5]],low=-np.pi,high=np.pi)
        
        peaks = sg.find_peaks(tpva,height=2,prominence=.25)
        pdx =peaks[0]
        
        peaksfl = sg.find_peaks(tmf,height=4,prominence=1)
        pdxf = peaksfl[0]
        # if len(pdx)==0:
        #     pdx = np.array(np.argmax(tpva))
        #     plen=1
        # else:
        #     plen = len(pdx)
            
        
        plen = len(pdx)
        txf = xf[dxpost]
        tyf = yf[dxpost]
        tyfp = yf[dxpre]
        txfp = xf[dxpre]
        tyfp = tyfp-tyf[0]
        txfp = txfp-txf[0]
        txf=txf-txf[0]
        tyf = tyf-tyf[0]
        # plt.figure(int(i))
        psmall = tphase[pdx]
        psmallf = tphase[pdxf]
        # x2small = x2[pdx]
        pvsmall = tpva[pdx]
        mfsmall =tmf[pdx]
        # psin = np.sin(psmall)*pvsmall
        # pcos = np.cos(psmall)*pvsmall
        # plt.plot(angs+xoff[i2],angc+yoff[i2],color='k',zorder=-1,linewidth=0.5)
        # plt.scatter(psin+xoff[i2],pcos+yoff[i2],c=x2small,vmin=0,vmax=600,s=4)
        
        # xp = np.sin(thead)*3+xoff[i2] 
        # yp = np.cos(thead)*3+yoff[i2]
        
        # plt.plot([xoff[i2],xp],[yoff[i2],yp],color='r',linestyle='-')
        
        # xp = np.sin(head_ex)*3+xoff[i2] 
        # yp = np.cos(head_ex)*3+yoff[i2]
        
        # plt.plot([xoff[i2],xp],[yoff[i2],yp],color=[0,0.5,1],linestyle='-')
        
        # plt.scatter(psmall,pvsmall,color='k')#c=np.arange(0,plen),vmin=0,vmax=9)
        # # for p,ip in enumerate(pvsmall):
        # #     plt.figure(p)
        # #     plt.scatter(psmall[p],pvsmall[p],color='k')
        # dat = np.zeros((plen,3))
        
        # dat[:,0] =psmall
        # dat[:,1] = pvsmall
        # dat[:,2] = np.arange(0,plen)
        # if len(alldat)==0:
        #     alldat = dat.copy()
        # else:
        #     alldat = np.append(alldat,dat,axis=0)
        
        
        
        plt.figure()
        plt.subplot(1,2,1)
        plt.plot(x2,tpva)
        plt.plot(x2,vd[dxpost],color=[0.5,0.5,0.5])
        plt.plot(x2,tmf,color='g')
        plt.scatter(x2[pdx],tpva[pdx],color='r')
        
        plt.scatter(x2[pdxf],tmf[pdxf],color='g')
        plt.title(i+ '_'+str(i2))
        plt.subplot(1,2,2)
        plt.plot(txf,tyf,color='k')
        g=plt.gca()
        g.set_aspect('equal')
        scale=5
        
        for ip, p in enumerate(psmall):
            plt.subplot(1,2,2)
            psin = np.sin(p)*scale*pvsmall[ip]
            pcos = np.cos(p)*scale*pvsmall[ip]
            xp = [txf[pdx[ip]],txf[pdx[ip]]+psin]
            yp = [tyf[pdx[ip]],tyf[pdx[ip]]+pcos]
            plt.plot(xp,yp,color='b')
            
            plt.plot([txfp[0],xp[0]],[tyfp[0],yp[0]],color='r')
            # if ip>0:
            #     plt.plot( [txf[pdx[ip]],txf[pdx[ip-1]]],[tyf[pdx[ip]],tyf[pdx[ip-1]]],color='r')
            
            
            plt.text(xp[0],yp[0],str(ip))
            plt.subplot(1,2,1)
            plt.text(x2[pdx[ip]],tpva[pdx[ip]],str(ip))
           
        for ip,p in enumerate(psmallf): 
            plt.subplot(1,2,2)
            psin = np.sin(p)*scale*pvsmall[ip]
            pcos = np.cos(p)*scale*pvsmall[ip]
            xp = [txf[pdxf[ip]],txf[pdxf[ip]]+psin]
            yp = [tyf[pdxf[ip]],tyf[pdxf[ip]]+pcos]
            plt.plot(xp,yp,color='g')
# cu = np.unique(alldat[:,2])
# for c in cu:
#     dx= alldat[:,2]==c
#     cmn = stats.circmean(alldat[dx,0],low=-np.pi,high=np.pi)
#     pmn =np.mean(alldat[dx,1])
#     plt.scatter(cmn,pmn,c=c,vmin=0,vmax=5,s=50,marker='+')
#%%
offset = 0
plt.close('all')
plt.figure(1,figsize=(5,10))
suffix = 'jumps'
timethresh = 10000000
for i in all_flies_fc2:
    cxa = all_flies_fc2[i]
    pva = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    pva = sg.savgol_filter(pva,20,3)
    pva = (pva)/np.std(pva)
    pva = pva-np.min(pva)
    if suffix=='jumps':
        jumps = cxa.get_jumps(time_threshold=180)
    else:
        jumps = cxa.get_jumps(time_threshold=180)
        jumps = cxa.get_entries_exits_like_jumps()
    phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*-cxa.side
    phase = ug.savgol_circ(phase,20,3)
    phase_eb = cxa.pdat['offset_eb_phase'].to_numpy()*-cxa.side
    #phase = ug.circ_subtract(phase,phase_eb)
    
    gridnum = np.ceil(np.sqrt(len(jumps))).astype(int)
    xoff = np.tile(np.arange(0,gridnum),gridnum)*10
    yoff =- xoff[np.argsort(xoff)]
    for i2,j in enumerate(jumps):
        #plt.figure(int(i))
        
        dxpre = np.arange(j[0],j[1])
        dxpost = np.arange(j[1],j[2])
        x1 = np.arange(0,len(dxpre))
        x2 = np.arange(len(dxpre),len(dxpre)+len(dxpost))

        
        tpva = pva[dxpost]
        tphase = phase[dxpost]
        bstart,bsize = ug.find_blocks(tpva>.5)
        bstart = bstart[bsize>=2]
        thead = stats.circmean(phase_eb[dxpre[:5]],low=-np.pi,high=np.pi)
        head_ex = stats.circmean(phase_eb[dxpost[:5]],low=-np.pi,high=np.pi)
        #plt.scatter(tphase,tpva,c=x2,s=2,vmin=0,vmax=600)
        plt.figure(1)
        #plt.plot([-np.pi,np.pi],[offset,offset],color='k',linewidth=0.5)
        plt.scatter(tphase[x2<timethresh],tpva[x2<timethresh]+offset,color='k',alpha=.1,s=3)
        #plt.scatter(tphase,tpva+offset,c=x2,s=2,vmin=0,vmax=300)
        #plt.figure(2)
        #plt.scatter(tphase[x2<50],tpva[x2<50]+offset,color='k',alpha=.1,s=3)
    offset +=-6
plt.plot([0,0],[offset+6,5],color='k',linestyle='--')
plt.plot([-np.pi/2,-np.pi/2],[offset+6,5],color='r',linestyle='--')
plt.plot([np.pi/2,np.pi/2],[offset+6,5],color=[0.5,0.5,0.5],linestyle='--')
plt.xlim([-np.pi,np.pi])
plt.xticks([-np.pi,0,np.pi],labels=[-180,0,180])
plt.xlabel('phase')
plt.ylabel('bump amp (norm pva)')

savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
plt.ylim([-37,8])
#plt.savefig(os.path.join(savedir,'ReturnPhaseAmp' +suffix+ str(timethresh)+'.png'))
#plt.savefig(os.path.join(savedir,'ReturnPhaseAmp' +suffix+str(timethresh)+'.pdf'))
#%%

#%% Plot phase when bump 'comes back' versus prior entry angle/epg phase
# EPG at 'bump off' versus fsb and bump on
plt.close('all')
ang = np.linspace(-np.pi,np.pi,num=100)
angs = np.sin(ang)*6
angc = np.cos(ang)*6
choice_text = np.array(['entry angle','exit angle','entry phase'])
for i in all_flies_fc2:
    # if i=='0':
    #     continue
    cxa = all_flies_fc2[i]
    pva = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    pva = sg.savgol_filter(pva,20,3)
    pva = (pva)/np.std(pva)
    pva = pva-np.min(pva)
    
    jumps = cxa.get_jumps(time_threshold=180)
    phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*-cxa.side
    phase = ug.savgol_circ(phase,20,3)
    phase_eb = cxa.pdat['offset_eb_phase'].to_numpy()*-cxa.side
    #phase = ug.circ_subtract(phase,phase_eb)
    offset = 0
    gridnum = np.ceil(np.sqrt(len(jumps))).astype(int)
    xoff = np.tile(np.arange(0,gridnum),gridnum)*15.5
    yoff =- xoff[np.argsort(xoff)]
    offset = 0
    for i2,j in enumerate(jumps):
        plt.figure(int(i))
        
        dxpre = np.arange(j[0],j[1])
        dxpost = np.arange(j[1],j[2])
        x1 = np.arange(0,len(dxpre))
        x2 = np.arange(len(dxpre),len(dxpre)+len(dxpost))

        
        tpva = pva[dxpost]
        tphase = phase[dxpost]
        bstart,bsize = ug.find_blocks(tpva>.5)
        bstart = bstart[bsize>=2]
        thead = stats.circmean(phase_eb[dxpre[:10]],low=-np.pi,high=np.pi)
        head_ex = stats.circmean(phase_eb[dxpost[:10]],low=-np.pi,high=np.pi)
        ent_phase = stats.circmean(phase[dxpre[:10]],low=-np.pi,high=np.pi)
        r = np.random.rand(len(dxpost))
        r = 1.5*x2/np.max(x2)
        hs = np.sin(phase_eb[dxpost])*(6+r)+xoff[i2]
        hc = np.cos(phase_eb[dxpost])*(6+r)+yoff[i2]
        plt.scatter(hs,hc,c=x2,alpha=.4,s=10)


        xp = np.sin(phase[dxpost])*tpva+xoff[i2]
        yp = np.cos(phase[dxpost])*tpva+yoff[i2]
        
        xy = np.append((np.sin(phase[dxpost])*tpva)[:,np.newaxis],(np.cos(phase[dxpost])*tpva)[:,np.newaxis],axis=1)
        
        plt.plot(angs+xoff[i2],angc+yoff[i2],color='k',zorder=-1,linewidth=0.5)
        plt.scatter(xp[tpva<2],yp[tpva<2],c=x2[tpva<2],s=10,alpha=.1,vmin=np.min(x2),vmax=np.max(x2))
        plt.scatter(xp[tpva>2],yp[tpva>2],c=x2[tpva>2],s=10,alpha=1,vmin=np.min(x2),vmax=np.max(x2))
        plt.plot(np.array([0,0])+xoff[i2],np.array([-3,3])+yoff[i2],color='k',linestyle='--',zorder=-1)
        plt.plot(np.array([-3,3])+xoff[i2],np.array([0,0])+yoff[i2],color='k',linestyle='--',zorder=-1)
        
        xp = np.sin(thead)*3+xoff[i2] 
        yp = np.cos(thead)*3+yoff[i2]
        vprior = np.array([np.sin(thead),np.cos(thead)])
        vp = ug.round_to_sig_figs(np.mean(np.matmul(xy,vprior)),2)
        
        plt.plot([xoff[i2],xp],[yoff[i2],yp],color='r',linestyle='-')
        
        xp = np.sin(head_ex)*3+xoff[i2] 
        yp = np.cos(head_ex)*3+yoff[i2]
        vex= np.array([np.sin(head_ex),np.cos(head_ex)])
        ve = ug.round_to_sig_figs(np.mean(np.matmul(xy,vex)),2)
        
        plt.plot([xoff[i2],xp],[yoff[i2],yp],color=[0,0.5,1],linestyle='-')
        
        xp = np.sin(ent_phase)*3+xoff[i2] 
        yp = np.cos(ent_phase)*3+yoff[i2]
        
        vent= np.array([np.sin(ent_phase),np.cos(ent_phase)])
        ven = ug.round_to_sig_figs(np.mean(np.matmul(xy,vent)),2)
        
        plt.plot([xoff[i2],xp],[yoff[i2],yp],color=[1,0.5,0],linestyle='-')
       # plt.text(xoff[i2],yoff[i2]-2,'Entry '+str(vp) +' Exit ' +str(ve) +' Phase '+ str(ven),horizontalalignment='center',fontsize=8)
        
        plt.figure(100+int(i))
        
        plt.plot(x2+offset,np.matmul(xy,vprior),color='r')
        plt.plot(x2+offset,-np.matmul(xy,vex),color=[0,0.5,1])
        plt.plot(x2+offset,np.matmul(xy,vent),color=[1,0.5,0])
        
        plt.figure(int(i))
        cat1 =np.matmul(xy,vprior)
        cat2 = -np.matmul(xy,vex)
        cat3 = np.matmul(xy,vent)
        
        c1b,c1bs = ug.find_blocks(cat1>2)
        c2b,c2bs = ug.find_blocks(cat2>2)
        c3b,c3bs = ug.find_blocks(cat3>2)
        
        c1b = c1b[c1bs>=10]
        c1bs = c1bs[c1bs>=10]
        c2b = c2b[c2bs>=10]
        c2bs = c2bs[c2bs>=10]
        c3b = c3b[c3bs>=10]
        c3bs = c3bs[c3bs>=10]
        bs_array = np.zeros(3,dtype=int)-1
        bs_array2 = np.zeros(3,dtype=int)-1
        bsize_array = np.zeros(3)
        if len(c1b)>0:
            bs_array[0] = c1b[-1]
            bs_array2[0] = c1b[-1]+c1bs[-1]
            bsize_array[0] = np.sum(cat1[bs_array[0]:bs_array[0]+c1bs[-1]])
        if len(c2b)>0:
            bs_array[1] = c2b[-1]
            bs_array2[1] = c2b[-1]+c2bs[-1]
            bsize_array[1] = np.sum(cat2[bs_array[1]:bs_array[1]+c2bs[-1]])
        if len(c3b)>0:
            bs_array[2] = c3b[-1]
            bs_array2[2] = c3b[-1]+c3bs[-1]
            bsize_array[2] = np.sum(cat3[bs_array[2]:bs_array[2]+c3bs[-1]])
            
        tm =np.max(bs_array2)
        if tm>0:
            tmm = np.max(tm-10,0)
            t_indices = np.where(bs_array2>tmm)[0]
            tbs = bsize_array[t_indices]
            
            tchoice = t_indices[np.where(tbs==np.max(tbs))[0]]
            print(len(tchoice))
            plt.text(xoff[i2],yoff[i2]+3,choice_text[tchoice[:]][:],horizontalalignment='center')
        else:
            plt.text(xoff[i2],yoff[i2]+3,'other',horizontalalignment='center')
        
        
        offset += np.max(x2)+10
        plt.figure(1000)
        plt.scatter(vp,ve,color='k')
    plt.figure(100+int(i))
    plt.plot([0,offset],[0,0],color='k',zorder=-1)
    plt.plot([0,offset],[2,2],color='k',zorder=-1)
    plt.plot([0,offset],[-2,-2],color='k',zorder=-1)
    plt.figure(int(i))
    g=plt.gca()
    g.set_aspect('equal')
plt.figure(1000)
plt.plot([-3,3],[-3,3],color='r')
#%% Look at first 25 entries
plt.close('all')
ang = np.linspace(-np.pi,np.pi,num=100)
angs = np.sin(ang)*6
angc = np.cos(ang)*6
for i in all_flies_hdc:
    cxa = all_flies_hdc[i]
    pva = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    pva = sg.savgol_filter(pva,20,3)
    pva = (pva)/np.std(pva)
    pva = pva-np.min(pva)
    jumps = cxa.get_jumps()
    jumps = cxa.get_cross_overs()
    #jumps = cxa.get_entries_exits_like_jumps()
   # jumps = jumps[:,:]
    phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*-cxa.side
    phase = ug.savgol_circ(phase,20,3)
    phase_eb = cxa.pdat['offset_eb_phase'].to_numpy()*-cxa.side
    #phase = ug.circ_subtract(phase,phase_eb)
    offset = 0
    gridnum = np.ceil(np.sqrt(len(jumps))).astype(int)
    xoff = np.tile(np.arange(0,gridnum),gridnum)*14
    yoff =- xoff[np.argsort(xoff)]
    for i2,j in enumerate(jumps):
        plt.figure(int(i))
        dxpre = np.arange(j[0],j[1])
        dxpost = np.arange(j[1],j[2])
        x1 = np.arange(0,len(dxpre))
        x2 = np.arange(len(dxpre),len(dxpre)+len(dxpost))

        
        tpva = pva[dxpost]
        tphase = phase[dxpost]
        bstart,bsize = ug.find_blocks(tpva>.5)
        bstart = bstart[bsize>=2]
        thead = stats.circmean(phase_eb[dxpre[:5]],low=-np.pi,high=np.pi)
        head_ex = stats.circmean(phase_eb[dxpost[:5]],low=-np.pi,high=np.pi)
        ent_phase = stats.circmean(phase[dxpre[:5]],low=-np.pi,high=np.pi)
        r = np.random.rand(len(dxpost))
        
        hs = np.sin(phase_eb[dxpost])*(6+r)+xoff[i2]
        hc = np.cos(phase_eb[dxpost])*(6+r)+yoff[i2]
        plt.scatter(hs,hc,c=x2,alpha=.4,s=10)


        xp = np.sin(phase[dxpost])*tpva+xoff[i2]
        yp = np.cos(phase[dxpost])*tpva+yoff[i2]
        plt.plot(angs+xoff[i2],angc+yoff[i2],color='k',zorder=-1,linewidth=0.5)
        plt.scatter(xp,yp,c=x2,s=10,alpha=0.5)
        plt.plot(np.array([0,0])+xoff[i2],np.array([-3,3])+yoff[i2],color='k',linestyle='--',zorder=-1)
        plt.plot(np.array([-3,3])+xoff[i2],np.array([0,0])+yoff[i2],color='k',linestyle='--',zorder=-1)
        
        xp = np.sin(thead)*3+xoff[i2] 
        yp = np.cos(thead)*3+yoff[i2]
        plt.plot([xoff[i2],xp],[yoff[i2],yp],color='r',linestyle='-')
        
        xp = np.sin(head_ex)*3+xoff[i2] 
        yp = np.cos(head_ex)*3+yoff[i2]
        plt.plot([xoff[i2],xp],[yoff[i2],yp],color=[0,0.5,1],linestyle='-')
        
        xp = np.sin(ent_phase)*3+xoff[i2] 
        yp = np.cos(ent_phase)*3+yoff[i2]
        plt.plot([xoff[i2],xp],[yoff[i2],yp],color=[1,0.5,0],linestyle='-')

    g=plt.gca()
    g.set_aspect('equal')

#%% plot details to check
plt.close('all')
for i in all_flies_hdc:
    cxa = all_flies_hdc[i]
    cxa.jump_return_details(flynum=int(i),jumplen=180)

#%% Plot details to check
plt.close('all')
for i in all_flies_hdj:
    cxa = all_flies_hdj[i]
    cxa.jump_return_details(flynum=int(i),jumplen=180)
#%% RDP algorithm for straight line walking
# Get phase during this
import rdp
x = cxa.ft2['ft_posx'].to_numpy()
y = cxa.ft2['ft_posy'].to_numpy()
ins = cxa.ft2['instrip'].to_numpy()
xy = np.append(x[:,np.newaxis],y[:,np.newaxis],axis=1)

o = rdp.rdp(xy,epsilon=5)
o_index = rdp.rdp(xy,epsilon=5,return_mask=True)
ow = np.where(o_index)[0]
phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
phase_eb = cxa.pdat['offset_eb_phase']
for i,io in enumerate(ow[1:]):
    dx = np.arange(ow[i],io)
    tins = ins[dx]
    if np.sum(tins)>0:
        negtins = tins==0
        bstart,bsize = ug.find_blocks(negtins)
        bm = np.argmax
    
    runlen = len(dx)/10
    tphase = phase[dx]
    tmn = stats.circmean(tphase,low=-np.pi,high=np.pi)
    if np.sum(tins)==0:
        plt.scatter(runlen,tmn,s=5,color='k')
    else:
        plt.scatter(runlen,tmn,s=5,color='r')
    
    #plt.figure()
    #plt.plot(tphase)

#%%
xy2 = xy.copy()
xy2[:,0] = xy2[:,0]*2
o = rdp.rdp(xy,epsilon=5)
plt.plot(x,y,color='k')
plt.plot(o[:,0],o[:,1],color='r')
g = plt.gca()
g.set_aspect('equal')


#%% Scatter fsb upper phase and eb phase
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
plt.close('all')
ticks =np.array( [-np.pi,-np.pi/2,0,np.pi/2,np.pi])
ticklabels = np.array([-180,-90,0,90,180])
for i in all_flies_hdj:
    cxa = all_flies_hdj[i]
    x = cxa.pdat['phase_eb']
    y = cxa.pdat['phase_fsb_upper']
    
    plt.figure()
    plt.scatter(x,y,s=2,alpha = 0.2,color='k')
    
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.xticks(ticks,labels=ticklabels,fontsize=15)
    plt.yticks(ticks,labels=ticklabels,fontsize=15)
    plt.xlabel('EB phase (deg)',fontsize=15)
    plt.ylabel('FSB phase (deg)',fontsize=15)
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.plot([-np.pi,0],[0,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.plot([0,np.pi],[-np.pi,0],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.subplots_adjust(bottom=.2,left=.2)
    savename = 'EB_vs_FSB_phase_good_tracker' + cxa.name+ '.png'
    plt.savefig(os.path.join(savedir,savename))
    savename = 'EB_vs_FSB_phase_good_tracker' + cxa.name+ '.pdf'
    plt.savefig(os.path.join(savedir,savename))
    
    
for i in all_flies_hdj_b:
    cxa = all_flies_hdj_b[i]
    x = cxa.pdat['phase_eb']
    y = cxa.pdat['phase_fsb_upper']
    #cxa.plot_traj_arrow(cxa.pdat['offset_fsb_upper_phase'],np.mean(cxa.pdat['wedges_fsb_upper']/4,axis=1),a_sep= 4)
    
    plt.figure()
    plt.scatter(x,y,s=2,alpha = 0.2,color='k')
    
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.xticks(ticks,labels=ticklabels,fontsize=15)
    plt.yticks(ticks,labels=ticklabels,fontsize=15)
    plt.xlabel('EB phase (deg)',fontsize=15)
    plt.ylabel('FSB phase(deg)',fontsize=15)
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.plot([-np.pi,0],[0,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.plot([0,np.pi],[-np.pi,0],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.subplots_adjust(bottom=.2,left=.2)
    
    savename = 'EB_vs_FSB_phase_bad_tracker' + cxa.name+ '.png'
    plt.savefig(os.path.join(savedir,savename))
    savename = 'EB_vs_FSB_phase_bad_tracker' + cxa.name+ '.pdf'
    plt.savefig(os.path.join(savedir,savename))
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdc'
for i in all_flies_hdc:
    cxa = all_flies_hdc[i]
    x = cxa.pdat['phase_eb']
    y = cxa.pdat['phase_fsb_upper']
    #cxa.plot_traj_arrow(cxa.pdat['offset_fsb_upper_phase'],np.mean(cxa.pdat['wedges_fsb_upper']/4,axis=1),a_sep= 4)
    
    plt.figure()
    plt.scatter(x,y,s=2,alpha = 0.2,color='k')
    
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.xticks(ticks,labels=ticklabels,fontsize=15)
    plt.yticks(ticks,labels=ticklabels,fontsize=15)
    plt.xlabel('EB phase (deg)',fontsize=15)
    plt.ylabel('FSB phase(deg)',fontsize=15)
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.plot([-np.pi,0],[0,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.plot([0,np.pi],[-np.pi,0],color=[0.0,.9,.6],linestyle='-',linewidth=2)
    plt.subplots_adjust(bottom=.2,left=.2)
    
    savename = 'EB_vs_FSB_phase_bad_tracker' + cxa.name+ '.png'
    plt.savefig(os.path.join(savedir,savename))
    savename = 'EB_vs_FSB_phase_bad_tracker' + cxa.name+ '.pdf'
    plt.savefig(os.path.join(savedir,savename))
#%% Scatter fsb upper phase and eb phase with inside and out
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj'
plt.close('all')
ticks =np.array( [-np.pi,-np.pi/2,0,np.pi/2,np.pi])
ticklabels = np.array([-180,-90,0,90,180])
suffix = ['out','in','pre']
for i in all_flies_hdj:
    cxa = all_flies_hdj[i]
    x = cxa.pdat['phase_eb']
    y = cxa.pdat['phase_fsb_upper']
    ins = cxa.ft2['instrip'].to_numpy()
    first = np.where(ins>0)[0][0]
    ins[:first] = 2
    for i in range(3):
        plt.figure()
        plt.scatter(x[ins==i],y[ins==i],s=2,alpha = 0.2,color='k')
        #plt.scatter(x[ins==1],y[ins==1],s=2,alpha = 0.2,color='r')
        plt.xlim([-np.pi,np.pi])
        plt.ylim([-np.pi,np.pi])
        plt.xticks(ticks,labels=ticklabels,fontsize=15)
        plt.yticks(ticks,labels=ticklabels,fontsize=15)
        plt.xlabel('EB phase (deg)',fontsize=15)
        plt.ylabel('FSB phase (deg)',fontsize=15)
        plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
        plt.plot([-np.pi,0],[0,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
        plt.plot([0,np.pi],[-np.pi,0],color=[0.0,.9,.6],linestyle='-',linewidth=2)
        plt.subplots_adjust(bottom=.2,left=.2)
        plt.title(suffix[i],fontsize=20)
    
    
        savename = 'EB_vs_FSB_phase_good_tracker'  + cxa.name + '_'+suffix[i]+ '.png'
        plt.savefig(os.path.join(savedir,savename))
        savename = 'EB_vs_FSB_phase_good_tracker_' + cxa.name + '_'+suffix[i]+ '.pdf'
        plt.savefig(os.path.join(savedir,savename))
    
    
for i in all_flies_hdj_b:
    cxa = all_flies_hdj_b[i]
    x = cxa.pdat['phase_eb']
    y = cxa.pdat['phase_fsb_upper']
    ins = cxa.ft2['instrip'].to_numpy()

    #cxa.plot_traj_arrow(cxa.pdat['offset_fsb_upper_phase'],np.mean(cxa.pdat['wedges_fsb_upper']/4,axis=1),a_sep= 4)
    first = np.where(ins>0)[0][0]
    ins[:first] = 2
    for i in range(3):
        plt.figure()
        plt.scatter(x[ins==i],y[ins==i],s=2,alpha = 0.2,color='k')
        #plt.scatter(x[ins==1],y[ins==1],s=2,alpha = 0.2,color='r')
        plt.xlim([-np.pi,np.pi])
        plt.ylim([-np.pi,np.pi])
        plt.xticks(ticks,labels=ticklabels,fontsize=15)
        plt.yticks(ticks,labels=ticklabels,fontsize=15)
        plt.xlabel('EB phase (deg)',fontsize=15)
        plt.ylabel('FSB phase (deg)',fontsize=15)
        plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
        plt.plot([-np.pi,0],[0,np.pi],color=[0.0,.9,.6],linestyle='-',linewidth=2)
        plt.plot([0,np.pi],[-np.pi,0],color=[0.0,.9,.6],linestyle='-',linewidth=2)
        plt.subplots_adjust(bottom=.2,left=.2)
        plt.title(suffix[i],fontsize=20)
        
        savename = 'EB_vs_FSB_phase_bad_tracker'  + cxa.name + '_'+suffix[i]+ '.png'
        plt.savefig(os.path.join(savedir,savename))
        savename = 'EB_vs_FSB_phase_bad_tracker_' + cxa.name + '_'+suffix[i]+ '.pdf'
        plt.savefig(os.path.join(savedir,savename))
#%%
from scipy.optimize import minimize

class model3_2:
    def __init__(self):
        self.version=1
    def predict_phase(self,phase,phase_eb,offset_pre,offset_in,offset_out,ins):
        phase_pred = np.zeros_like(phase)
        phase_pred[ins>0] = ug.circ_subtract(phase_eb[ins>0].copy(),offset_in)
        phase_pred[ins<1] = ug.circ_subtract(phase_eb[ins<1].copy(),offset_out)
        first = np.where(ins>0)[0][0]
        phase_pred[:first] = 0
        phase_pred[:first] = ug.circ_subtract(phase_eb[:first].copy(),offset_pre)
        return phase_pred
    
    def fit_model(self,phase,phase_eb,ins,delay_range):
        offsets = [0,0,0]
        bounds = [(-np.pi,np.pi),(-np.pi,np.pi),(-np.pi,np.pi)]
        mses = np.zeros(delay_range-1)
        models = {}
        for i in np.arange(1,delay_range,dtype=int):
            print('Delay: ',i)
            def objective(offsets):
                phase_pred = self.predict_phase(phase[i:],phase_eb[:-i].copy(),offsets[0],offsets[1],offsets[2],ins[:-i])
                return ug.circ_mse(phase[i:],phase_pred)
            models.update({str(i): minimize(objective,offsets,bounds=bounds,method='L-BFGS-B')})
            mses[i-1] = models[str(i)].fun
            
        best = np.argmin(mses)+1
        self.mses = mses
        self.res = models[str(best)]
        self.delay = best
#%%        
hdj_models = {}
for i in all_flies_hdj:
    cxa = all_flies_hdj[i]
    phase = cxa.pdat['phase_fsb_upper']
    phase_eb = cxa.pdat['phase_eb']    
    ins = cxa.ft2['instrip'].to_numpy()
    mdl3_2 = model3_2()
    mdl3_2.fit_model(phase,phase_eb,ins,30)
    hdj_models.update({i:mdl3_2})
    
hdj_models_b = {}
for i in all_flies_hdj:
    cxa = all_flies_hdj_b[i]
    phase = cxa.pdat['phase_fsb_upper']
    phase_eb = cxa.pdat['phase_eb']    
    ins = cxa.ft2['instrip'].to_numpy()
    mdl3_2 = model3_2()
    mdl3_2.fit_model(phase,phase_eb,ins,30)
    hdj_models_b.update({i:mdl3_2}) 
    
    
#%% Simple scatter of predictions
plt.close('all')
savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj\model'
for i in all_flies_hdj:
    cxa = all_flies_hdj[i]
    mdl = hdj_models[i]
    phase = cxa.pdat['phase_fsb_upper']
    phase_eb = cxa.pdat['phase_eb']    
    ins = cxa.ft2['instrip'].to_numpy().copy()
    
    phase_pred = mdl.predict_phase(phase[mdl.delay:],phase_eb[:-mdl.delay]
                                   ,mdl.res.x[0],mdl.res.x[1],mdl.res.x[2],ins[:-mdl.delay])
    plt.figure()
    plt.scatter(phase_pred,phase[mdl.delay:],color='k',s=2,alpha=0.2)
    uplt.standard_phase_plt()
    plt.xlabel('Predicted phase (deg)',fontsize=15)
    plt.ylabel('FSB phase (deg)',fontsize=15)
    savename = 'Pred_vs_FSB_phase_good_tracker' + cxa.name+ '.png'
    plt.savefig(os.path.join(savedir,savename))
    savename = 'Pred_vs_FSB_phase_good_tracker' + cxa.name+ '.pdf'
    plt.savefig(os.path.join(savedir,savename))
    
for i in all_flies_hdj_b:
    cxa = all_flies_hdj_b[i]
    mdl = hdj_models_b[i]
    phase = cxa.pdat['phase_fsb_upper']
    phase_eb = cxa.pdat['phase_eb']    
    ins = cxa.ft2['instrip'].to_numpy().copy()
    
    phase_pred = mdl.predict_phase(phase[mdl.delay:],phase_eb[:-mdl.delay]
                                   ,mdl.res.x[0],mdl.res.x[1],mdl.res.x[2],ins[:-mdl.delay])
    plt.figure()
    plt.scatter(phase_pred,phase[mdl.delay:],color='k',s=2,alpha=0.2)
    uplt.standard_phase_plt()
    plt.xlabel('Predicted phase (deg)',fontsize=15)
    plt.ylabel('FSB phase (deg)',fontsize=15)
    savename = 'Pred_vs_FSB_phase_bad_tracker' + cxa.name+ '.png'
    plt.savefig(os.path.join(savedir,savename))
    savename = 'Pred_vs_FSB_phase_bad_tracker' + cxa.name+ '.pdf'
    plt.savefig(os.path.join(savedir,savename))
#%% Plot MSE and parameters

plt.close('all')

savedir = r'Y:\Presentations\2025\11_LabMeeting\hdj\model'
x = np.arange(0,3)
x2 = np.arange(0,4)
for i in all_flies_hdj:
    cxa = all_flies_hdj[i]
    mdl = hdj_models[i]
    
    offsets = 180*mdl.res.x/np.pi
    offsets = offsets*np.sign(offsets)
    
    # Model offsets
    plt.figure(1,figsize=(5,5))
    plt.scatter(x,offsets,color=[0.4,0.6,.8])
    
    
    # Model errors
    plt.figure(2,figsize=(5,5))
    error= 180* np.arccos(1-mdl.res.fun)/np.pi
   # plt.scatter(1,error,color=[0.4,0.6,0.8])
    
    errors = np.zeros(4)
    errors[0] = error
    
    phase = cxa.pdat['phase_fsb_upper']
    phase_eb = cxa.pdat['phase_eb']    
    ins = cxa.ft2['instrip'].to_numpy().copy()
    phase2 = phase[mdl.delay:]
    
    phase_pred = mdl.predict_phase(phase[mdl.delay:],phase_eb[:-mdl.delay]
                                   ,mdl.res.x[0],mdl.res.x[1],mdl.res.x[2],ins[:-mdl.delay])
    
    ins2 = ins[mdl.delay:]
    first = np.where(ins2>0)[0][0]
    errors[1] = 180*np.arccos(1-ug.circ_mse(phase2[:first],phase_pred[:first]))/np.pi
    errors[2] = 180*np.arccos(1-ug.circ_mse(phase2[ins2>0],phase_pred[ins2>0]))/np.pi
    out = np.where(ins2<1)[0]
    out = out[out>first]
    errors[3] = 180* np.arccos(1-ug.circ_mse(phase2[out],phase_pred[out]))/np.pi
    plt.scatter(x2,errors,color=[.4,.6,.8])
    #Model delays
    plt.figure(3,figsize=(2,5))
    plt.scatter(0,mdl.delay/10,color=[.4,.6,.8])
    

for i in all_flies_hdj_b:
    cxa = all_flies_hdj_b[i]
    mdl = hdj_models_b[i]
    
    offsets = 180*mdl.res.x/np.pi
    offsets = offsets*np.sign(offsets)
    plt.figure(1)
    plt.scatter(x+.25,offsets,color=[.8,.2,.2])
    plt.figure(2)
    error=180* np.arccos(1-mdl.res.fun)/np.pi
    errors = np.zeros(4)
    errors[0] = error
    
    phase = cxa.pdat['phase_fsb_upper']
    phase_eb = cxa.pdat['phase_eb']    
    ins = cxa.ft2['instrip'].to_numpy().copy()
    phase2 = phase[mdl.delay:]
    
    phase_pred = mdl.predict_phase(phase[mdl.delay:],phase_eb[:-mdl.delay]
                                   ,mdl.res.x[0],mdl.res.x[1],mdl.res.x[2],ins[:-mdl.delay])
    
    ins2 = ins[mdl.delay:]
    first = np.where(ins2>0)[0][0]
    errors[1] = 180*np.arccos(1-ug.circ_mse(phase2[:first],phase_pred[:first]))/np.pi
    errors[2] = 180*np.arccos(1-ug.circ_mse(phase2[ins2>0],phase_pred[ins2>0]))/np.pi
    out = np.where(ins2<1)[0]
    out = out[out>first]
    errors[3] = 180* np.arccos(1-ug.circ_mse(phase2[out],phase_pred[out]))/np.pi
    plt.scatter(x2+.25,errors,color=[.8,.2,.2])


    plt.figure(3,figsize=(2,5))
    plt.scatter(0.25,mdl.delay/10,color=[.8,.2,.2])
plt.figure(1)
plt.ylim([0,180])
plt.xticks(x+0.125,labels =['pre','in','out'],fontsize=15)
plt.ylabel('offsets (deg)',fontsize=15)
plt.yticks(fontsize=15)
plt.subplots_adjust(0.3)

plt.figure(2)
plt.ylim([0,90])
plt.xticks(x2,labels=['total','pre','in','out'],fontsize=15)
plt.ylabel('angular error (deg)',fontsize=15)
#plt.xlim([0.5,2])
plt.yticks(fontsize=15)
plt.subplots_adjust(0.3)

plt.figure(3)
plt.ylim([0,2])
plt.ylabel('delay (s)',fontsize=15)
plt.yticks(fontsize=15)
plt.xticks([])
plt.subplots_adjust(0.4)
plt.xlim([-0.5,0.5])

savenames = ['offsets','errors','delay']
for i,s in enumerate(savenames):
    plt.figure(i+1)
    plt.savefig(os.path.join(savedir,s +'.png'))
    plt.savefig(os.path.join(savedir,s +'.pdf'))
    
#%% Example model fit
plt.close('all')
choice = str(1)
cxa = all_flies_hdj[choice]
mdl = hdj_models[choice]

phase = cxa.pdat['phase_fsb_upper']
phase_eb = cxa.pdat['phase_eb']    
ins = cxa.ft2['instrip'].to_numpy().copy()
phase2 = phase[mdl.delay:]
ins2 = ins[mdl.delay:].astype(float)
x= np.arange(0,len(phase2))/10


phase_pred = mdl.predict_phase(phase[mdl.delay:],phase_eb[:-mdl.delay]
                               ,mdl.res.x[0],mdl.res.x[1],mdl.res.x[2],ins[:-mdl.delay])
plt.fill_between(x,ins2*0.0-np.pi,2*np.pi*ins2-np.pi,color='r',linewidth=0)
s1 = plt.scatter(x,phase2,color='k',s=2, clip_on=False)
s1.set_clip_path(None)
s2 = plt.scatter(x,phase_pred,color=[.4,.6,.8],s=2, clip_on=False)
s2.set_clip_path(None)
error = np.abs(ug.circ_subtract(phase2,phase_pred))/np.pi

u =ug()
vx,vy,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'])
headings = ug.savgol_circ(cxa.ft2['ft_heading'].to_numpy(),50,3)
va = ug.get_ang_velocity(headings,cxa.pv2['relative_time'].to_numpy())
va = np.append(0,va)
va = va/np.std(va)
vd = np.append(0,vd)
vd = vd/np.std(vd)
vdn = vd[mdl.delay:]*.5-8
van = np.abs(va[mdl.delay:])*.5-11
plt.plot(x,vdn,color=[0.5,0.5,0.5])
plt.plot(x,error*2-6,color=[0,.4,.4])
plt.xlim([0,400])
plt.plot(x,van,color=[0.3,0.3,0.3])
plt.xlabel('time (s)')
plt.yticks([-11,-8,-6,-4,-np.pi,0,np.pi])
for artist in plt.gca().get_children():
    artist.set_clip_on(False)
#plt.savefig(os.path.join(savedir,'ExampleFitError.pdf'))

#plt.savefig(os.path.join(savedir,'ExampleFitError.svg'))

#%% Example model fit
plt.close('all')
choice = str(2)
cxa = all_flies_hdj_b[choice]
mdl = hdj_models_b[choice]

phase = cxa.pdat['phase_fsb_upper']
phase_eb = cxa.pdat['phase_eb']    
ins = cxa.ft2['instrip'].to_numpy().copy()
phase2 = phase[mdl.delay:]
ins2 = ins[mdl.delay:].astype(float)
x= np.arange(0,len(phase2))/10


phase_pred = mdl.predict_phase(phase[mdl.delay:],phase_eb[:-mdl.delay]
                               ,mdl.res.x[0],mdl.res.x[1],mdl.res.x[2],ins[:-mdl.delay])
plt.fill_between(x,ins2*0.0-np.pi,2*np.pi*ins2-np.pi,color='r',linewidth=0)
s1 = plt.scatter(x,phase2,color='k',s=2, clip_on=False)
s1.set_clip_path(None)
s2 = plt.scatter(x,phase_pred,color=[.4,.6,.8],s=2, clip_on=False)
s2.set_clip_path(None)
error = np.abs(ug.circ_subtract(phase2,phase_pred))/np.pi

u =ug()
vx,vy,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'])
headings = ug.savgol_circ(cxa.ft2['ft_heading'].to_numpy(),50,3)
va = ug.get_ang_velocity(headings,cxa.pv2['relative_time'].to_numpy())
va = np.append(0,va)
va = va/np.std(va)
vd = np.append(0,vd)
vd = vd/np.std(vd)
vdn = vd[mdl.delay:]*.5-8
van = np.abs(va[mdl.delay:])*.5-11
plt.plot(x,vdn,color=[0.5,0.5,0.5])
plt.plot(x,error*2-6,color=[0,.4,.4])
plt.xlim([0,400])
plt.plot(x,van,color=[0.3,0.3,0.3])
plt.xlabel('time (s)')
plt.yticks([-11,-8,-6,-4,-np.pi,0,np.pi])

#plt.savefig(os.path.join(savedir,'ExampleFitError_bad.pdf'))

#plt.savefig(os.path.join(savedir,'ExampleFitError_bad.svg'))
#%% HMM
cxa = all_flies_hdj['1']
from analysis_funs.CX_HMM import CX_HMM
hmm = CX_HMM(3) 
phase = cxa.pdat['phase_fsb_upper']
phase_eb = cxa.pdat['phase_eb']
hmm.fit(phase,phase_eb,verbose=True)
#%%
states = hmm.viterbi(phase,phase_eb)
x = np.arange(0,len(phase))
ustates = np.unique(states)
plt.scatter(x,phase,color='k',s=3)
#plt.scatter(x,phase_eb,color='r',s=3)
ypred = np.zeros_like(phase)
for u in ustates:
    dx = states==u
    plt.scatter(x[dx],hmm.wrap_ang(phase_eb[dx]+hmm.angles[u]),s=2)
    ypred[dx] = phase_eb[dx]+hmm.angles[u]
ins = cxa.ft2['instrip'].to_numpy()
plt.plot(x,ins*np.pi*2 -np.pi ,color='r')

print(ug.circ_mse(phase,ypred))

u = ug()
vx,vy,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'])
vd = vd/np.std(vd)
plt.plot(x[1:],vd/2 -6,color=[0.3,0.3,0.3])

















