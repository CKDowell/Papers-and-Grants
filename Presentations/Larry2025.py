# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 15:21:00 2025

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
    
    "hDeltaH_ACV_pulses": [r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f1\Trial4",
                           r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f2\Trial5",
                           r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f1\Trial5",
                           r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f2\Trial4"],
    "hDeltaH_Oct_pulses":[
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f1\Trial5",
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\251204\f2\Trial4",
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f1\Trial4",
        r"Y:\Data\FCI\Hedwig\hDeltaH_SS92512\260109\f2\Trial5"]
}
all_flies = {}

for group_name, datadirs in groups.items():
    all_flies[group_name] = {}
    
    for i, datadir in enumerate(datadirs):
        print(f"[{group_name}] {datadir}")
        etp = ET_paper(datadir)
        all_flies[group_name][i] = etp
#%% (1) Example raw data from first plume encounter
plt.close('all')
neurons = ['hDeltaC','hDeltaJ','hDeltaH']
tdat = [1,2,4]
savedir = r'Y:\Presentations\2026\01_Larry'
for i,n in enumerate(neurons):
    cxa = all_flies[n][tdat[i]].cxa
    cxa.simple_raw_plot(plotphase=False,regions = ['eb','fsb_upper'],yk='eb')
    explen = len(cxa.pdat['phase_eb'])
    ins = cxa.ft2['instrip'].to_numpy()
    pon = np.where(ins>0)[0][0]
    ylmin = pon-200 
    ylmax = pon+1200
    plt.ylim([ylmin,ylmax])
    plt.savefig(os.path.join(savedir,'Raw_first_odour_'+n+str(i)+'.pdf'))
#%% Default mode estimation
from scipy.optimize import minimize
plt.close('all')
# Model is a simple offset relative to heading
class model_hdj:
    def __init__(self):
        self.version=1
    def predict_phase(self,phase,phase_eb,offset):
        phase_pred = ug.circ_subtract(phase_eb,offset)
        return phase_pred
    
    def fit_model(self,phase,phase_eb,delay_range):
        offset = 0
        bounds = [(-np.pi,np.pi)]
        mses = np.zeros(delay_range-1)
        models = {}
        for i in np.arange(1,delay_range,dtype=int):
            
            def objective(offset):
                phase_pred = self.predict_phase(phase[i:],phase_eb[:-i],offset)
                return ug.circ_mse(phase[i:],phase_pred)
            models.update({str(i): minimize(objective,offset,bounds=bounds,method='L-BFGS-B')})
            mses[i-1] = models[str(i)].fun
            
        best = np.argmin(mses)+1
        self.mses = mses
        self.res = models[str(best)]
        self.delay = best

# Model is a single offset whose phase varies with animal being left or right of heading
class model_hdc:
    def __init__(self):
        self.version =1 
    def predict_phase(self,phase,phase_eb,heading,heading_off,offset):
        heading = ug.circ_subtract(heading,heading_off)
        lheading = heading<0 
        rheading = heading>0 
        phase_pred = np.zeros_like(phase_eb)
        phase_pred[lheading] = ug.circ_subtract(phase_eb[lheading],offset)
        phase_pred[rheading] = ug.circ_subtract(phase_eb[rheading],-offset)
        return phase_pred
    
    def fit_model(self,phase,phase_eb,heading,delay_range):
        offsets = [0.01,0]
        bounds = [(0,np.pi),(-np.pi,np.pi)]
        
        mses = np.zeros(delay_range-1)
        models = {}
        for i in np.arange(1,delay_range,dtype=int):
        
            def objective(offsets):
                phase_pred = self.predict_phase(phase[i:],phase_eb[:-i],heading[i:],offsets[0],offsets[1])
                return ug.circ_mse(phase[i:],phase_pred)
            models.update({str(i): minimize(objective,offsets,bounds=bounds,method='L-BFGS-B')})
            mses[i-1] = models[str(i)].fun
        best = np.argmin(mses)+1
        self.mses = mses
        self.res = models[str(best)]
        self.delay = best
groups = ['hDeltaC','hDeltaJ']
index_type = 'long_ret'
for g in groups:
    data = all_flies[g]
    for d in data:
        
        # Data
        cxa = data[d].cxa
        jumps = cxa.get_entries_exits_like_jumps()
        jlen = jumps[:,2]-jumps[:,1]
        jdx = jlen>600 # greater than a minute
        fitdx = np.array([],dtype=int)
        tjumps = jumps[jdx,:]
        for j in tjumps:
            tdx = np.arange(j[1],j[2])
            fitdx = np.append(fitdx,tdx)
        
        
        
        phase = cxa.pdat['phase_fsb_upper']
        phase = ug.savgol_circ(phase,20,3)
        eb = cxa.pdat['phase_eb']
        eb = ug.savgol_circ(eb,20,3)
        heading = cxa.ft2['ft_heading'].to_numpy()
        mh = ug.motor2heading(cxa.ft2['motor_step_command'].to_numpy())
        
        
        
        ins = cxa.ft2['instrip'].to_numpy()
        iw = np.where(ins>0)[0]
        if index_type=='post':
            dx = dx = np.arange(iw[-1],len(eb))
        elif index_type=='pre':
            dx = np.arange(0,iw[0])
        elif index_type=='all':
            dx = np.arange(0,len(eb))
        elif index_type == 'long_ret':
            fitdx = np.append(fitdx,np.arange(iw[-1],len(eb)))
            dx = np.unique(fitdx)
        # phase = phase[dx]
        # eb = eb[dx]
        # heading = heading[dx]
        # hDeltaJ model
        mdl = model_hdj()
        mdl.fit_model(phase[dx],eb[dx],50)
        
        pred = mdl.predict_phase(phase[mdl.delay:],eb[:-mdl.delay],mdl.res.x[0])
        errors = np.zeros(2)
        angles = np.zeros(2)
        errors[0] =mdl.res.fun
        angles[0] = 180*mdl.res.x[0]/np.pi
        ang_error1 = np.abs(ug.circ_subtract(phase[mdl.delay:],pred))
        
        # plot it
        plt.figure()
        
        x  = np.arange(0,len(eb))/10
        plt.scatter(x,phase,color='k',s=2)
        plt.scatter(x,eb,color=[0.5,0.5,0.5],s=2)
        #plt.scatter(x,ug.circ_subtract(eb,np.pi/4),color=[0.5,0.5,0.5],s=2)
       # plt.scatter(x,ug.circ_subtract(eb,-np.pi/4),color=[0.5,0.5,0.5],s=2)
        #plt.fill_between(x,ug.circ_subtract(eb,np.pi/4),ug.circ_subtract(eb,-np.pi/4),color=[0.5,0.5,0.5],zorder=-1)
        simple_pred = ug.circ_subtract(eb,np.pi*(135/180))
        rheading = heading>0 
        simple_pred[rheading] = ug.circ_subtract(eb[rheading],-np.pi*(135/180))
        #plt.scatter(x,simple_pred,color=[0.5,0.5,0.5],s=2)
        #plt.scatter(x[mdl.delay:],pred,color='g',s=2) # hdeltaJ model
        #plt.plot(x[mdl.delay:],ang_error1-2*np.pi,color='g')
        
        # hDeltaC model
        mdl = model_hdc()
        mdl.fit_model(phase[dx],eb[dx],mh[dx],50)
        pred = mdl.predict_phase(phase[mdl.delay:],eb[:-mdl.delay],mh[mdl.delay:],mdl.res.x[0],mdl.res.x[1])
        preds = ug.savgol_circ(pred,20,3)
        ang_error1 = np.abs(ug.circ_subtract(phase[mdl.delay:],preds))
        # plot it
        plt.scatter(x[mdl.delay:],preds,color='b',s=2)
        
        pred = mdl.predict_phase(phase[mdl.delay:],eb[:-mdl.delay],mh[mdl.delay:],mdl.res.x[0],-mdl.res.x[1])
        preds = ug.savgol_circ(pred,20,3)
        ang_error1 = np.abs(ug.circ_subtract(phase[mdl.delay:],preds))
        # plot it
        plt.scatter(x[mdl.delay:],preds,color=[0.5,0,1],s=2)
        
        plt.plot(x,ins*2*np.pi-np.pi,color='r')
        #plt.plot(x[mdl.delay:],ang_error1-2*np.pi,color='b')

        errors[1] =mdl.res.fun
        angles[1] = 180*mdl.res.x[1]/np.pi
        errors = np.round(errors,decimals=2)
        plt.plot(x,heading-np.pi*2,color='k')
        plt.plot(x,mh-np.pi*2,color=[0.5,0,0])
        plt.plot(x,0*mh-np.pi*2,color='r')
        plt.title(g+ ': ' +str(d) +' error hdj: ' +str(errors[0]) + ' ang: '+str(angles[0])+ ' error hdc: ' +str(errors[1])+ ' ang: '+str(angles[1]))
        print(g)
        print('hdj ang',angles[0],' hdc ang',angles[1])
        print(' hdc centre',mdl.res.x[0])
#%% HMM model with two modes
from analysis_funs.CX_HMM import CX_HMM
cxa =  all_flies['hDeltaC'][5].cxa
hmm = CX_HMM(3) 
phase = cxa.pdat['phase_fsb_upper']
phase_eb = cxa.pdat['phase_eb']
hmm.fit_hDeltaC(phase,phase_eb,verbose=True,max_iters=200)

#%%
ins = cxa.ft2['instrip'].to_numpy()
states,logp = hmm.viterbi(phase,phase_eb)
p = np.exp(logp)
ustates = np.unique(states)
#plt.scatter(x,phase_eb,color='r',s=3)
ypred = np.zeros_like(phase)
phase_pred = phase_eb
x= np.arange(0,len(phase_pred))
plt.scatter(x,phase,color='k',s=1)
colours = np.array([[0,0.3,1],[0,0.3,0.5],[0,0.8,0.3]])
for u in ustates:
    dx = states==u
    #plt.scatter(x[dx],hmm.wrap_ang(phase_pred[dx]+hmm.angles[u]),s=2)
    ypred[dx] = hmm.wrap_ang(phase_pred[dx]+hmm.angles[u])
ypreds = ug.savgol_circ(ypred,20,3)
for u in ustates:
    dx = states==u
    plt.scatter(x[dx],ypreds[dx],s=2,label=np.round(180*(hmm.angles[u])/np.pi),color=colours[u,:])
plt.legend()
plt.show()
plt.plot(x,ins*2*np.pi-np.pi,color='r')
#%% 
states,logp = hmm.viterbi(phase,phase_eb)
jumps = cxa.get_jumps()
pn = np.exp(logp)
p = np.zeros((len(pn),2))
p[:,0] = pn[:,0]+pn[:,1]
p[:,1] = pn[:,2]
data =np.zeros((100,p.shape[1],len(jumps)))
newtime = np.arange(0,100)
for i,j in enumerate(jumps):
    dx = np.arange(j[1],j[2])
    tdat = p[dx]
    oldtime = np.linspace(0,99,len(dx))
    for ir in range(p.shape[1]):
        data[:,ir,i] = np.interp(newtime,oldtime,tdat[:,ir])
pmean = np.mean(data,axis=2)
plt.plot(pmean)
#%%
ug.circ_mse(phase_short,ypred)
u = ug()
vx,vy,vd = u.get_velocity(cxa.ft2['ft_posx'].to_numpy(),cxa.ft2['ft_posy'].to_numpy(),cxa.pv2['relative_time'])
vd = vd/np.std(vd)
vs = np.append(0,vd)
plt.plot(x[1:],vd[mdl.winlen:]/2 -6,color=[0.3,0.3,0.3])
#%% Iterate through data running model
groups = ['FC2','hDeltaC','hDeltaJ']
out_dict = {}
for ig,g in enumerate(groups):
    print(g)
    data =  all_flies[g]
    angles = np.zeros((len(data),3))
    error = np.zeros(len(data))
    models = {}
    for i,d in enumerate(data):
        print(d)
        cxa=  data[d].cxa
        hmm = CX_HMM(3) 
        phase = cxa.pdat['phase_fsb_upper']
        phase_eb = cxa.pdat['phase_eb']
        hmm.fit_hDeltaC(phase,phase_eb,verbose=False,max_iters=200)
        angles[i,:] = hmm.angles
        
        states,logp = hmm.viterbi(phase,phase_eb)
        p = np.exp(logp)
        ustates = np.unique(states)
        #plt.scatter(x,phase_eb,color='r',s=3)
        ypred = np.zeros_like(phase)
        phase_pred = phase_eb
        x= np.arange(0,len(phase_pred))
        for u in ustates:
            dx = states==u
            ypred[dx] = hmm.wrap_ang(phase_pred[dx]+hmm.angles[u])
        error[i] = stats.circmean(np.abs(ug.circ_subtract(phase,ypred)),low=-np.pi,high=np.pi)
        models.update({str(i):hmm})
    out_dict.update({g:{'angles': angles,'error':error,'models':models}})
ug.save_pick(out_dict,os.path.join(r'Y:\Data\FCI\ConsolidatedData\hmm','hhm_models.pkl'))
#%% Load pickle
out_dict = ug.load_pick(os.path.join(r'Y:\Data\FCI\ConsolidatedData\hmm','hhm_models.pkl'))
#%% Plot angles
colours = np.array([[81,156,205],[49,99,125],[81,61,204]])/255
plt.close('all')
offset=0
for ig,g in enumerate(groups):
    plt.figure(1)
    angles =180*out_dict[g]['angles']/np.pi
    for i in [0,2]:
        
        plt.scatter(np.ones(len(angles))*i+offset,angles[:,i],color = colours[ig,:])
    offset+=.1
    plt.figure(0)
    plt.scatter(angles[:,0],angles[:,2],color=colours[ig,:])
    plt.figure(2)
    error= 180*out_dict[g]['error']/np.pi
    plt.scatter(np.ones(len(angles))*ig,error,color=colours[ig,:])
    angles =out_dict[g]['angles']
    amean = stats.circmean(angles,high=np.pi,low=-np.pi,axis=0)
    ameanT = np.round(180*amean/np.pi)
    ameans = np.sin(amean)
    ameanc = np.cos(amean)
    angs = np.sin(angles)
    angc = np.cos(angles)
    plt.figure(3)
    for i in range(3):
        for i2 in range(len(angles)):
            plt.plot(np.array([0,angs[i2,i]])+ig*2,[0,angc[i2,i]],color=colours[ig,:],alpha=0.3)
        plt.plot(np.array([0,ameans[i]])+ig*2,[0,ameanc[i]],color=colours[ig,:])
        plt.scatter(ameans[i]+ig*2,ameanc[i],color=colours[ig,:])
        plt.text(ameans[i]+ig*2,ameanc[i],str(ameanT[i]),horizontalalignment='center')
        plt.text(ig*2,1.5,g,horizontalalignment='center')
    g= plt.gca()
    g.set_aspect('equal')
    
#%% Plot states on phase
#%% Plot states on trajectory
colours = np.array([[1,0.3,0.3],[0,0.3,0.8],[0,0.8,0.3]])
plt.close('all')
g = 'FC2'
i =5
cxa = all_flies[g][i].cxa
x = cxa.ft2['ft_posx'].to_numpy()
y = cxa.ft2['ft_posy'].to_numpy()
model = out_dict[g]['models'][str(i)]
x,y = ug.fictrac_repair(x,y)

ins = cxa.ft2['instrip'].to_numpy()
plt.scatter(x[ins>0],y[ins>0],color=[0.5,0.5,0.5])
phase = cxa.pdat['phase_fsb_upper']
phase_eb = cxa.pdat['phase_eb']
states,logp = model.viterbi(phase,phase_eb)
ypred = np.zeros_like(phase)
plt.figure(1)
for i in range(3):
    dx = states==i
    px = x.copy()
    py = y.copy()
    px[~dx] =np.nan
    py[~dx] = np.nan
    plt.plot(px,py,color=colours[i,:])
    ypred[dx] = model.wrap_ang(phase_eb[dx]+model.angles[i])
g = plt.gca()
#uplt.coloured_line_simple(x,y,states,colours,0,2)
g.set_aspect('equal')
plt.figure(2)
t = np.arange(0,len(phase))/10
plt.scatter(t,phase,color='k',s=3)
plt.scatter(t,phase_eb,color=[0.5,0.5,0.5],s=3)
for i in range(3):
    dx = states==i
    plt.scatter(t[dx],ypred[dx],color=colours[i,:],s=3)
plt.plot(t,ins*np.pi*2-np.pi,color = 'r')
jumps = cxa.get_entries_exits_like_jumps()
for i,j in enumerate(jumps):
    plt.figure(1)
    plt.text(x[j[0]],y[j[0]],str(i))
    
    plt.figure(2)
    plt.text(t[j[0]],np.pi+.5,str(i))
amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
amp = amp/np.std(amp)
plt.plot(t,amp/3-np.pi*2,color='k')
for i in range(3):
    dx = states==i
    plt.scatter(t[dx],amp[dx]/3-np.pi*2,color=colours[i,:],s=5)
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
#%% State transitions
#1 Characterise
#2 document when - are they before the animal heads back to the plume
# - can they solely be explained by hypothetical PFNa inputs


#%% histogram of phase offsets
savedir = r'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'
plt.close('all')
groups = ['FC2','hDeltaC','hDeltaJ']
colours = np.array([[81,156,205],[49,99,125],[81,61,204]])/255
bins = np.linspace(-np.pi,np.pi,30)
plotbins = bins[1:]-np.mean(np.diff(bins))/2
outdata = {}

for ig,g in enumerate(groups):
    print(g)
    data =  all_flies[g]
    angles = np.zeros((len(data),3))
    error = np.zeros(len(data))
    
    tdata = np.zeros((len(bins)-1,len(data),6))
    for i,d in enumerate(data):
        print(d)
        cxa=  data[d].cxa
        jumps = cxa.get_jumps()
        jdx = np.array([],dtype=int)
        for j in jumps:
            jdx = np.append(jdx,np.arange(j[1],j[2]))
        pfsb = cxa.pdat['phase_fsb_upper']*-cxa.side
        peb = cxa.pdat['phase_eb']*-cxa.side
        diffdata = ug.circ_subtract(pfsb,peb)
        tdata[:,i,0],_ = np.histogram(diffdata,bins,density=True)
        ins = cxa.ft2['instrip'].to_numpy()
        befdx = np.arange(0,np.where(ins>0)[0][0])
        tdata[:,i,1],_ =np.histogram(diffdata[ins>0],bins,density=True)
        tdata[:,i,2],_ =np.histogram(diffdata[ins<1],bins,density=True)
        tdata[:,i,3],_ =np.histogram(diffdata[jdx],bins,density=True)
        tdata[:,i,4],_ = np.histogram(diffdata[befdx],bins,density=True)
        tdata[:,i,5],_ = np.histogram(diffdata[befdx[-300:]],bins,density=True)
    outdata.update({g:tdata})


titles= ['all','in odour','out of odour','jump returns','pre air','pre air last min']
savename = ['all','in_odour','out_odour','jump_return','pre_air','pre_air_last_min']
for ir in range(6):
    plt.figure()
    plt.title(titles[ir])
    for ig,g in enumerate(groups):
        tdata = outdata[g][:,:,ir]
        tmean = np.mean(tdata,axis=1)
        plt.plot(plotbins,tmean,color=colours[ig,:],label=g)
        tse = np.std(tdata,axis=1)/np.sqrt(tdata.shape[1])
        plt.fill_between(plotbins,tmean+tse,tmean-tse,alpha=0.5,color=colours[ig,:])
    plt.legend()
    
    plt.ylabel('probability')
    plt.xlabel('fsb - eb phase (deg)')
    plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
    plt.xlim([-np.pi,np.pi])
    
    plt.savefig(os.path.join(savedir,'phase_histogam_' +savename[ir]+'.png'))
#%% State transitions
g = 'hDeltaC'
i =2
cxa = all_flies[g][i].cxa
x = cxa.ft2['ft_posx'].to_numpy()
y = cxa.ft2['ft_posy'].to_numpy()
model = out_dict[g]['models'][str(i)]
x,y = ug.fictrac_repair(x,y)
ins = cxa.ft2['instrip'].to_numpy()

ins = cxa.ft2['instrip'].to_numpy()
st = np.where(ins>0)[0][0]
#x = x-x[0]
#y = y-y[0]
phase = cxa.pdat['phase_fsb_upper']
phase_eb = cxa.pdat['phase_eb']
states,logp = model.viterbi(phase,phase_eb)
state_diff = np.diff(states)
prior_state = np.where(np.abs(state_diff)>0)[0]
next_state = prior_state+1
t = np.arange(0,len(phase))
# plt.plot(t,states,color='k')
# plt.scatter(t[prior_state],states[prior_state],color='r')
# plt.scatter(t[next_state],states[next_state],color='g')


ypred = np.zeros_like(phase)
plt.figure(1)
for i in range(3):
    dx = states==i
    px = x.copy()
    py = y.copy()
    px[~dx] =np.nan
    py[~dx] = np.nan
   # plt.plot(px,py,color=colours[i,:])
    ypred[dx] = model.wrap_ang(phase_eb[dx]+model.angles[i])
colours = np.array([[1,0.3,0.3],[0,0.3,0.8],[0,0.8,0.3]])

t = np.arange(0,len(phase))/10
plt.scatter(t,phase,color='k',s=3)
plt.scatter(t,phase_eb,color=[0.5,0.5,0.5],s=3)
for i in range(3):
    dx = states==i
    plt.scatter(t[dx],ypred[dx],color=colours[i,:],s=3)

plt.plot(t,ins*np.pi*2-np.pi,color = 'r')
fsb_trans = np.array([],dtype=int)
eb_trans = np.array([],dtype=int)
for i,e in enumerate(prior_state):
    pdx = np.arange(e-10,e)
    adx = np.arange(e+1,np.min([e+11,len(phase)]))
    ebp = stats.circmean(phase_eb[pdx],high=np.pi,low=-np.pi)
    eba = stats.circmean(phase_eb[adx],high=np.pi,low=-np.pi)
    
    pp = stats.circmean(phase[pdx],high=np.pi,low=-np.pi)
    pa = stats.circmean(phase[adx],high=np.pi,low=-np.pi)
    
    e_diff = ug.circ_subtract(ebp,eba)
    p_diff = ug.circ_subtract(pp,pa)
    
    if np.abs(e_diff)>np.abs(p_diff):
        plt.scatter(t[e],ypred[e],marker='*',color='m')
        eb_trans = np.append(eb_trans,e+1)
    else:
        plt.scatter(t[e],ypred[e],marker='D',color='m')
        fsb_trans = np.append(fsb_trans,e+1)
    
eb_offset = cxa.pdat['offset_eb_phase'].to_numpy()
p_offset = cxa.pdat['offset_fsb_upper_phase'].to_numpy()

plt.figure()
cxa.plot_traj_arrow(cxa.pdat['offset_fsb_upper_phase'].to_numpy(),np.mean(cxa.pdat['wedges_fsb_upper']/2,axis=1),a_sep= 4)
#plt.plot(x,y,color='k')
#plt.scatter(x[ins>0],y[ins>0],color=[0.7,.7,.7])
plt.scatter(x[fsb_trans],y[fsb_trans],color='r',zorder=1)
# g = plt.gca()
# g.set_aspect('equal')

# for f in fsb_trans:
#     fdx = np.arange(f,f+3)
#     eb_mn = stats.circmean(eb_offset[fdx],high=np.pi,low=-np.pi)
#     p_mn = stats.circmean(p_offset[fdx],high=np.pi,low=-np.pi)
    
#     xs = np.sin(eb_mn)*5
#     yc =np.cos(eb_mn)*5
#     plt.plot(np.array([0,xs])+x[f],np.array([0,yc])+y[f],color='k')    
    
#     xs = np.sin(p_mn)*5
#     yc =np.cos(p_mn)*5
#     plt.plot(np.array([0,xs])+x[f],np.array([0,yc])+y[f],color='b') 

#%% Plot phase during jumps by percentile recovery

phase_percentiles = [20,40,60,80,100,120,140,160,180,200]
phase_percentiles = [100]
#cxa = all_flies['hDeltaJ'][3].cxa
plt.close('all')
groupdat = all_flies['hDeltaJ']
groups = ['hDeltaJ','hDeltaC']
colours = np.array([[81,156,205],[49,99,125],[81,61,204]])/255

for ig,g in enumerate(groups):
    groupdat = all_flies[g]
    for f in groupdat:
        cxa = groupdat[f].cxa
        jumps = cxa.get_jumps_w_prior()
        amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*-cxa.side
        phase_off = ug.circ_subtract(cxa.pdat['phase_fsb_upper'],cxa.pdat['phase_eb'])*-cxa.side
        rec_phase_all = np.zeros((1000,len(jumps),len(phase_percentiles)))*np.nan
        rec_phase_all_t = np.zeros((1000,len(jumps),len(phase_percentiles)))*np.nan
        for ij,j in enumerate(jumps):
            dxbef = np.arange(j[0],j[1])
            if len(dxbef)>5:
                dxbef = dxbef[-5:]
            befmean = np.mean(amp[dxbef])
            dxret = np.arange(j[2],j[3])
            ampret = amp[dxret]
            
            ampret_perc = 100*ampret/befmean
            dxall = np.arange(j[1],j[3])
            ampall_perc = 100*amp[dxall]/befmean
            #plt.plot(dxall-j[2],ampall_perc)
            #plt.plot(ampret_perc)
            
            ap_round = np.round(ampret_perc/20)*20
            phase_ret = phase[dxret]
            phaseo_ret = phase_off[dxret]
            plt_perc =np.array([])
            plt_t = np.array([])
            plt_phase = np.array([])
            plt_phaseo = np.array([])
            for ip,p in enumerate(phase_percentiles):
                pdx = ap_round==p
                if sum(pdx)==0:
                    continue
                pf = np.where(pdx)[0][0]
                #plt.scatter(ampret_perc[pf],phase_ret[pf],color='k',s=2)
                plt_perc = np.append(plt_perc,ampret_perc[pf])
                plt_phase = np.append(plt_phase,phase_ret[pf])
                plt_phaseo = np.append(plt_phaseo,phaseo_ret[pf])
                plt_t = np.append(plt_t,(dxret[pf]-j[2])/10)
                
                # plt.figure(3)
                # plt.scatter((dxret[pf:]-dxret[pf])/10,phase_ret[pf:],color=colours[ig,:],alpha=.5,s=2)
                rec_phase_all[dxret[pf:]-dxret[pf],ij,ip] = phase_ret[pf:]
                rec_phase_all_t[dxret[pf:]-j[2],ij,ip] = phase_ret[pf:]
            
            plt.figure(1)
            plt.subplot(1,2,1)
            plt.plot(plt_perc,plt_phase,color=colours[ig,:],marker='.',alpha=.5)
            plt.subplot(1,2,2)
            plt.plot(plt_perc,plt_phaseo,color=colours[ig,:],marker='.',alpha=.5)
    
            plt.figure(2)
            plt.subplot(1,2,1)
            plt.scatter(plt_t,plt_phase,color=colours[ig,:],marker='.',alpha=.5)
            plt.subplot(1,2,2)
            plt.scatter(plt_t,plt_phaseo,color=colours[ig,:],marker='.',alpha=.5)

        plt.figure(3)
        plt.subplot(1,2,1)
        pmean = stats.circmean(rec_phase_all,low=-np.pi,high=np.pi,axis=1,nan_policy='omit')
        pmean = pmean[~np.isnan(pmean)]
        t = np.arange(0,len(pmean))/10
        plt.scatter(t,pmean,color=colours[ig,:],s=4)
        plt.xlim([0,20])
        
        plt.subplot(1,2,2)
        pmean = stats.circmean(rec_phase_all_t,low=-np.pi,high=np.pi,axis=1,nan_policy='omit')
        pmean = pmean[~np.isnan(pmean)]
        t = np.arange(0,len(pmean))/10
        plt.scatter(t,pmean,color=colours[ig,:],s=4)
        plt.xlim([0,20])
#%% Plot mean phase and amp on trajectory
dchoices = ['hDeltaJ','hDeltaC','FC2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
numpoints = 12  
plt.close('all')
grandphase = np.zeros((200,len(dchoices)))
grandtraj = np.zeros((200,2))
grandamp = np.zeros((200,2,len(dchoices)))
for i,d1 in enumerate(dchoices):
    data=all_flies[d1]
        
    phases = np.zeros((200,len(data)))
    trajs = np.zeros((200,2,len(data)))
    amps = np.zeros((200,2,len(data)))
    for d in data:
        print(d1)
        etp = data[d]
        phase,traj,amp = etp.trajectory_mean(time_threshold2=0,regions=['fsb_upper'])
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
pdx = np.linspace(0,199,numpoints,dtype= 'int')

fx = np.array([-10,0,0,-10])
fy = np.array([-50,-50,0,0])
plt.fill(fx,fy,color=[0.7,0.7,.7])
plt.fill(fx-3,fy+50,color=[.7,.7,.7])
achoice=1
scale=15
for p in pdx:
    for i in range(len(dchoices)):
    
        tphase = grandphase[p,i]
        x = [grandtraj[p,0]+np.sin(tphase)*scale*grandamp[p,achoice,i],grandtraj[p,0]]
        y = [grandtraj[p,1]+np.cos(tphase)*scale*grandamp[p,achoice,i],grandtraj[p,1]]
        
        plt.plot(x,y,color=colours[i,:])
        
    vsx = np.sum(np.sin(grandphase[p,:2])*grandamp[p,achoice,:2])
    vsy = np.sum(np.cos(grandphase[p,:2])*grandamp[p,achoice,:2])
    
    ta = np.sqrt(vsx**2+vsy**2)
    vsx = grandamp[p,achoice,2]*(vsx/ta)*scale
    vsy = grandamp[p,achoice,2]*(vsy/ta)*scale
    
    x = np.array([vsx,0])+grandtraj[p,0]
    y = np.array([vsy,0])+grandtraj[p,1]    
    plt.plot(x,y,color=colours[i,:],linestyle='--')
g = plt.gca()
g.set_aspect('equal')
plt.ylim([-50,50])
savedir = r'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'
plt.savefig(os.path.join(savedir,'hDelta_FC2_vector_on_traj.pdf'))
#%% Pseudotime phase plots
achoice=1
dchoices = ['hDeltaJ','hDeltaC','FC2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
numpoints = 12  
plt.close('all')
grandphase = np.zeros((200,2,len(dchoices)))
grandtraj = np.zeros((200,2))
grandamp = np.zeros((200,2,2,len(dchoices)))
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
        for ip in range(2):
            tphase = grandphase[p,ip,i]
            x = np.array([np.sin(tphase)*grandamp[p,achoice,ip,i],0])+offset
            y = np.array([np.cos(tphase)*grandamp[p,achoice,ip,i],0])
            if ip>0:
                
                plt.plot(x,y,color=colours[i,:])
            else:
                plt.plot(x,y,color=[0.5,0.5,0.5])
            
                
                
g = plt.gca()
g.set_aspect('equal')
plt.fill([-.5,5.5,5.5,-.5],[-.5,-.5,.5,.5],color=[0.8,0.8,0.8],zorder=-1)
plt.savefig(os.path.join(savedir,'Jump_phases_pseudo_time.pdf'))
        
        
#%% Mean phase of jumps over time from plume exit

bins = np.arange(0,10.1,1)
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
grandmean_amp = np.zeros((len(bins)-1,3))
grandmean_phase = np.zeros((len(bins)-1,3))
for i,d in enumerate(dchoices):
    data = all_flies[d]
    tcolours = np.append(np.array([0.5,0.5,0.5])[:,np.newaxis],colours[i,:][:,np.newaxis],axis=1).T
    plt.figure(i+1)
    yoffset = 0
    flymeanphase = np.zeros((len(bins)-1,2,len(data)))
    flymeanamp = np.zeros((len(bins)-1,2,len(data)))
    for d in data:
        etp = data[d]
        phases,amps = etp.phase_time_leave_amp(int(np.max(bins)))
        
        offset = 0
        
        for b in range(0,len(bins)-1):
            bdx = np.arange(bins[b]*10,bins[b+1]*10,dtype=int)
            tp = stats.circmean(phases[bdx,:,:],low=-np.pi,high=np.pi,axis=0,nan_policy='omit')
            ta = np.nanmean(amps[bdx,:,:],axis=0)
            tpmean = stats.circmean(tp,low=-np.pi,high=np.pi,axis=1,nan_policy='omit')
            tamean = np.nanmean(ta,axis=1)
            
            x = np.sin(tp)*ta
            y = np.cos(tp)*ta
            
            xm = np.sin(tpmean)*tamean
            ym = np.cos(tpmean)*tamean
            flymeanphase[b,:,d] = tpmean
            flymeanamp[b,:,d]  = tamean
            for r in range(0,2):
                tx = np.append(x[r,:][:,np.newaxis],x[r,:][:,np.newaxis]*0,axis=1).T
                ty = np.append(y[r,:][:,np.newaxis],y[r,:][:,np.newaxis]*0,axis=1).T
                plt.plot(tx+offset,ty+yoffset,color=tcolours[r,:],alpha=0.25)
            
                tx = np.array([0,xm[r]])
                ty = np.array([0,ym[r]])
                plt.plot(tx+offset,ty+yoffset,color=tcolours[r,:],linewidth =2,zorder=5)
            offset+=1
        g = plt.gca()
        g.set_aspect('equal')            
        plt.xticks(np.linspace(0,offset-1,len(bins)-1),bins[:-1])
        yoffset = yoffset-1
        
    plt.figure(101)
    pmean = stats.circmean(flymeanphase,axis=2,low=-np.pi,high=np.pi)
    amean = np.nanmean(flymeanamp,axis=2)
    grandmean_amp[:,i] = amean[:,1]
    grandmean_phase[:,i] = pmean[:,1]
    offset = 0
    for b in range(0,len(bins)-1):
        for r in range(0,2):
            tx = np.array([0,np.sin(pmean[b,r])*amean[b,r]])
            ty = np.array([0,np.cos(pmean[b,r])*amean[b,r]])
            plt.plot(tx+offset,ty,color=tcolours[r,:])
        offset = offset+1
        
        
plt.figure(101)    
g= plt.gca()
g.set_aspect('equal')
plt.xticks(np.linspace(0,offset-1,len(bins)-1),bins[:-1])

xsum = np.sum(np.sin(grandmean_phase[:,:2])*grandmean_amp[:,:2],axis=1)
ysum = np.sum(np.cos(grandmean_phase[:,:2])*grandmean_amp[:,:2],axis=1)

offset = 0
for b in range(0,len(bins)-1):
    
    tx = np.array([0,xsum[b]])
    ty = np.array([0,ysum[b]])
    
    norm = grandmean_amp[b,2]/np.sqrt(tx[1]**2+ty[1]**2)
    plt.plot(tx*norm+offset,ty*norm,color='r',linestyle='--')
    offset = offset+1
savedir = r'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'
plt.savefig(os.path.join(savedir,'Mean_vector_sum.pdf'))

for i in range(1,4):
    print(i)
    plt.figure(i)
    plt.ylim([-7,1])
    plt.xlim([-1,10])
    plt.savefig(os.path.join(savedir,'Mean_vectors_indiflies_'+dchoices[i-1]+'.pdf'))
#%% Mean phase of jumps over time to plume entry

bins = np.arange(0,10.1,1)
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
grandmean_amp = np.zeros((len(bins)-1,3))
grandmean_phase = np.zeros((len(bins)-1,3))
for i,d in enumerate(dchoices):
    data = all_flies[d]
    tcolours = np.append(np.array([0.5,0.5,0.5])[:,np.newaxis],colours[i,:][:,np.newaxis],axis=1).T
    plt.figure(i+1)
    yoffset = 0
    flymeanphase = np.zeros((len(bins)-1,2,len(data)))
    flymeanamp = np.zeros((len(bins)-1,2,len(data)))
    for d in data:
        etp = data[d]
        phases,amps = etp.phase_time_amp(int(np.max(bins)))
        
        offset = 0
        
        for b in range(0,len(bins)-1):
            bdx = np.arange(bins[b]*10,bins[b+1]*10,dtype=int)
            tp = stats.circmean(phases[bdx,:,:],low=-np.pi,high=np.pi,axis=0,nan_policy='omit')
            ta = np.nanmean(amps[bdx,:,:],axis=0)
            tpmean = stats.circmean(tp,low=-np.pi,high=np.pi,axis=1,nan_policy='omit')
            tamean = np.nanmean(ta,axis=1)
            
            x = np.sin(tp)*ta
            y = np.cos(tp)*ta
            
            xm = np.sin(tpmean)*tamean
            ym = np.cos(tpmean)*tamean
            flymeanphase[b,:,d] = tpmean
            flymeanamp[b,:,d]  = tamean
            for r in range(0,2):
                tx = np.append(x[r,:][:,np.newaxis],x[r,:][:,np.newaxis]*0,axis=1).T
                ty = np.append(y[r,:][:,np.newaxis],y[r,:][:,np.newaxis]*0,axis=1).T
                plt.plot(tx+offset,ty+yoffset,color=tcolours[r,:],alpha=0.25)
            
                tx = np.array([0,xm[r]])
                ty = np.array([0,ym[r]])
                plt.plot(tx+offset,ty+yoffset,color=tcolours[r,:],linewidth =2,zorder=5)
            offset+=-1
        g = plt.gca()
        g.set_aspect('equal')            
        plt.xticks(np.linspace(offset+1,0,len(bins)-1),bins[:-1])
        yoffset = yoffset-1
        
    plt.figure(101)
    pmean = stats.circmean(flymeanphase,axis=2,low=-np.pi,high=np.pi)
    amean = np.nanmean(flymeanamp,axis=2)
    grandmean_amp[:,i] = amean[:,1]
    grandmean_phase[:,i] = pmean[:,1]
    offset = 0
    for b in range(0,len(bins)-1):
        for r in range(0,2):
            tx = np.array([0,np.sin(pmean[b,r])*amean[b,r]])
            ty = np.array([0,np.cos(pmean[b,r])*amean[b,r]])
            plt.plot(tx+offset,ty,color=tcolours[r,:])
        offset = offset-1
        
        
plt.figure(101)    
g= plt.gca()
g.set_aspect('equal')
plt.xticks(np.linspace(offset+1,0,len(bins)-1),bins[:-1])

xsum = np.sum(np.sin(grandmean_phase[:,:2])*grandmean_amp[:,:2],axis=1)
ysum = np.sum(np.cos(grandmean_phase[:,:2])*grandmean_amp[:,:2],axis=1)

offset = 0
for b in range(0,len(bins)-1):
    
    tx = np.array([0,xsum[b]])
    ty = np.array([0,ysum[b]])
    
    norm = grandmean_amp[b,2]/np.sqrt(tx[1]**2+ty[1]**2)
    plt.plot(tx*norm+offset,ty*norm,color=colours[2,:],linestyle='--')
    offset = offset-1
savedir = r'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'
plt.savefig(os.path.join(savedir,'Mean_vector_sum_to_entry.pdf'))

for i in range(1,4):
    print(i)
    plt.figure(i)
    plt.ylim([-7,1])
    plt.xlim([-10,1])
    plt.savefig(os.path.join(savedir,'Mean_vectors_to_entry_indiflies_'+dchoices[i-1]+'.pdf'))


#%% Analysis of phase before first turn back to plume
from scipy import ndimage as nd
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
from rdp import rdp
angthresh = .1
d = 'FC2'
f = 5
cxa = all_flies[d][f].cxa
jumps = cxa.get_jumps()
x = cxa.ft2['ft_posx'].to_numpy()
y = cxa.ft2['ft_posy'].to_numpy()
angvel = ug.get_ang_velocity(cxa.ft2['ft_heading'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
angvelsmth = nd.gaussian_filter(angvel,3)**2
phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
eb = cxa.pdat['offset_eb_phase'].to_numpy()
heading = cxa.ft2['ft_heading'].to_numpy()
u = ug()
vx,vy,vd = u.get_velocity(x,y,cxa.pv2['relative_time'].to_numpy())
minssize = 1
for ju in jumps:
    dx = np.arange(ju[1],ju[2])
    tx = x[dx]
    ty = y[dx]
    tx = tx-tx[0]
    ty = ty-ty[0]
    tp = phase[dx]
    te = eb[dx]
    tvd = vd[dx]
    th = heading[dx]
    av = angvel[dx]
    avs = angvelsmth[dx]
    
    xy = np.append(tx[:,np.newaxis],ty[:,np.newaxis],1)
    r = rdp(xy,epsilon=2)
    rdx = np.where(np.in1d(tx,r[:,0]))[0]
    still = np.abs(av)<.1
    
    sb,sbs = ug.find_blocks(still,mergeblocks=True,merg_threshold=2)
    sb = sb[sbs>=minssize]
    sbs = sbs[sbs>=minssize]
    
    if len(rdx)>2:
        sdx = np.arange(max(rdx[1]-5,0),rdx[1])
        pmn = stats.circmean(tp[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
        hmn = stats.circmean(th[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
        rvec = r[1,:]-r[0,:]
        rcom = rvec[0]*1j+rvec[1]
        ran = np.angle(rcom)*cxa.side*-1
        plt.figure(101)
        plt.scatter(hmn,pmn,color='k')
        plt.scatter(ran,pmn,color='r')
    
    # if len(sb)>0:
    #     fstdx = np.arange(sb[0],sb[0]+sbs[0])
    # # if len(rdx)>2:
    # #     fstill = np.argmin(np.abs(sb-rdx[1]))
    # #     fstdx =np.arange(sb[fstill],sb[fstill]+sbs[fstill])
    
    #     if len(fstdx)>5:
    #         fstdx2 = fstdx[-5:]
    #     else:
    #         fstdx2 = fstdx
            
    #     pmn = stats.circmean(tp[fstdx2],high=np.pi,low=-np.pi)*cxa.side*-1
    #     hmn = stats.circmean(th[fstdx2],high=np.pi,low=-np.pi)*cxa.side*-1
    #     plt.figure(101)
    #     plt.scatter(hmn,pmn,color='k')
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(tx,ty,color='k')
    plt.plot(r[:,0],r[:,1],color='r')
    
    for r1 in rdx:
        xp = np.array([0,np.sin(tp[r1])])+tx[r1]
        yp = np.array([0,np.cos(tp[r1])])+ty[r1]
        plt.plot(xp,yp,color='b')
    
    
    plt.gca().set_aspect('equal')
    plt.subplot(1,2,2)
    plt.plot(-8+av/2,color='k')
    plt.plot(-8+avs/2,color='r')
    plt.scatter(rdx,rdx*0-4,color='b')
    t= np.arange(0,len(tp))
    plt.scatter(t,tp,color='b',s=2)
    plt.scatter(t,th,color='k',s=2)
    plt.scatter(t[still],th[still],color='g',s=2)
    # if len(sb)>0:
    #     plt.scatter(t[fstdx],th[fstdx],color='r',s=5,zorder=-1)
    plt.plot(t,t*0,color='k')
    
    
    plt.plot(t,-8+tvd/5,color=[.5,.5,.5])


plt.figure(101)
plt.xlim([-np.pi,np.pi])
plt.ylim([-np.pi,np.pi])
plt.plot([-np.pi,np.pi],[0,0],color='k')
plt.plot([0,0],[-np.pi,np.pi],color='k')
plt.plot([-np.pi/2,-np.pi/2],[-np.pi,np.pi],color='r',linestyle='--')
plt.plot([-np.pi,np.pi],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
#%% plot rdp example
from scipy import ndimage as nd
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
from rdp import rdp
angthresh = .1
d = 'FC2'
f = 5
cxa = all_flies[d][f].cxa
jumps = cxa.get_jumps()
x = cxa.ft2['ft_posx'].to_numpy()*-cxa.side
y = cxa.ft2['ft_posy'].to_numpy()
angvel = ug.get_ang_velocity(cxa.ft2['ft_heading'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
angvelsmth = nd.gaussian_filter(angvel,3)**2
phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
eb = cxa.pdat['offset_eb_phase'].to_numpy()
heading = cxa.ft2['ft_heading'].to_numpy()
u = ug()
vx,vy,vd = u.get_velocity(x,y,cxa.pv2['relative_time'].to_numpy())
minssize = 1
for ju in jumps:
    dx = np.arange(ju[1],ju[2])
    
    tx = x[dx]
    ty = y[dx]
    tx = tx-tx[0]
    ty = ty-ty[0]
    xy = np.append(tx[:,np.newaxis],ty[:,np.newaxis],1)
    r = rdp(xy,epsilon=2)
    
    dx2 = np.arange(ju[0],ju[2])
    xplt = x[dx2]-x[ju[1]]
    yplt = y[dx2]-y[ju[1]]
    
    plt.figure()
    plt.plot(xplt,yplt,color='k')
    ymin = np.min(yplt)
    ymax = np.max(yplt)
    xvec = np.array([0,0,-10,-10])
    plt.fill(xvec,np.array([0,ymin,ymin,0]),color=[0.7,0.7,0.7])
    plt.fill(xvec-3,np.array([0,ymax,ymax,0]),color=[0.7,0.7,0.7])
    plt.plot(r[:,0],r[:,1],color='r')
    plt.scatter(r[:,0],r[:,1],color='r',zorder=2,s=30)
    plt.gca().set_aspect('equal')
#%% Heading and phase before first RDP turning point
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']

minssize = 1
xlabs = ['heading (deg)','rdp vector (deg)']
for di,d in enumerate(dchoices):
    plt.figure(di)
    for s in range(1,3):
        plt.subplot(1,2,s)
        plt.xlim([-np.pi,np.pi])
        plt.ylim([-np.pi,np.pi])
        plt.plot([-np.pi,np.pi],[0,0],color='k')
        plt.plot([0,0],[-np.pi,np.pi],color='k')
        plt.plot([-np.pi/2,-np.pi/2],[-np.pi,np.pi],color='r',linestyle='--')
        plt.plot([-np.pi,np.pi],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
        plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,-0,90,180])
        plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,-0,90,180])
        plt.ylabel(d+' phase (deg)')
        plt.xlabel(xlabs[s-1])
        plt.gca().set_aspect('equal')
        
    data = all_flies[d]
    for tf in data:
        cxa = data[tf].cxa
        jumps = cxa.get_jumps()
        x = cxa.ft2['ft_posx'].to_numpy()
        y = cxa.ft2['ft_posy'].to_numpy()
        angvel = ug.get_ang_velocity(cxa.ft2['ft_heading'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
        angvelsmth = nd.gaussian_filter(angvel,3)**2
        phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
        eb = cxa.pdat['offset_eb_phase'].to_numpy()
        heading = cxa.ft2['ft_heading'].to_numpy()
        u = ug()
        vx,vy,vd = u.get_velocity(x,y,cxa.pv2['relative_time'].to_numpy())
        for ju in jumps:
            dx = np.arange(ju[1],ju[2])
            tx = x[dx]
            ty = y[dx]
            tx = tx-tx[0]
            ty = ty-ty[0]
            tp = phase[dx]
            te = eb[dx]
            tvd = vd[dx]
            th = heading[dx]
            av = angvel[dx]
            avs = angvelsmth[dx]
            
            xy = np.append(tx[:,np.newaxis],ty[:,np.newaxis],1)
            r = rdp(xy,epsilon=2)
            rdx = np.where(np.in1d(tx,r[:,0]))[0]
            # still = np.abs(av)<.1
            
            # sb,sbs = ug.find_blocks(still,mergeblocks=True,merg_threshold=2)
            # sb = sb[sbs>=minssize]
            # sbs = sbs[sbs>=minssize]
            
            if len(rdx)>2:
                sdx = np.arange(max(rdx[1]-3,0),rdx[1])
                pmn = stats.circmean(tp[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
                hmn = stats.circmean(th[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
                rvec = r[1,:]-r[0,:]
                rcom = rvec[0]*1j+rvec[1]
                ran = np.angle(rcom)*cxa.side*-1
                
                plt.subplot(1,2,1)
                plt.scatter(hmn,pmn,color=colours[di,:],s=10)
                plt.subplot(1,2,2)
                plt.scatter(ran,pmn,color=colours[di,:],s=10)
for di,d in enumerate(dchoices):
    plt.figure(di)
    plt.savefig(os.path.join(savedir,'RDP_first_'+d+'.pdf'))
#%% Amplitude and timing of first RDP turning point
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
savedir = r'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'

minssize = 1
xlabs = ['heading (deg)','rdp vector (deg)']

for di,d in enumerate(dchoices):
    data = all_flies[d]
    dmean = np.zeros(len(data))
    tmean = np.zeros(len(data))
    inplume_amp = np.zeros(len(data))
    for tf in data:
        cxa = data[tf].cxa
        jumps = cxa.get_jumps()
        x = cxa.ft2['ft_posx'].to_numpy()
        y = cxa.ft2['ft_posy'].to_numpy()
        angvel = ug.get_ang_velocity(cxa.ft2['ft_heading'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
        angvelsmth = nd.gaussian_filter(angvel,3)**2
        amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        #phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
        #eb = cxa.pdat['offset_eb_phase'].to_numpy()
        #heading = cxa.ft2['ft_heading'].to_numpy()
        u = ug()
        vx,vy,vd = u.get_velocity(x,y,cxa.pv2['relative_time'].to_numpy())
        amps = np.array([])
        times = np.array([])
        plamps = np.array([])
        for ju in jumps:
            dx = np.arange(ju[1],ju[2])
            tx = x[dx]
            ty = y[dx]
            tx = tx-tx[0]
            ty = ty-ty[0]
            #tp = phase[dx]
            #te = eb[dx]
            tvd = vd[dx]
            
            ta = amp[dx]
            
            
            xy = np.append(tx[:,np.newaxis],ty[:,np.newaxis],1)
            r = rdp(xy,epsilon=2)
            rdx = np.where(np.in1d(tx,r[:,0]))[0]
            # still = np.abs(av)<.1
            
            # sb,sbs = ug.find_blocks(still,mergeblocks=True,merg_threshold=2)
            # sb = sb[sbs>=minssize]
            # sbs = sbs[sbs>=minssize]
            
            if len(rdx)>2:
                sdx = np.arange(max(rdx[1]-3,0),rdx[1])
                pmn = np.mean(ta[sdx])
                plt.figure(1)
                plt.scatter(rdx[1]*.1,pmn,color=colours[di,:],s=20)
                amps = np.append(amps,pmn)
                times = np.append(times,rdx[1]*.1)
                plamps = np.append(plamps,np.mean(amp[np.arange(ju[1]-3,ju[1])]))
                # pmn = stats.circmean(tp[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
                # hmn = stats.circmean(th[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
                # rvec = r[1,:]-r[0,:]
                # rcom = rvec[0]*1j+rvec[1]
                # ran = np.angle(rcom)*cxa.side*-1
                
                # plt.subplot(1,2,1)
                # plt.scatter(hmn,pmn,color=colours[di,:],s=10)
                # plt.subplot(1,2,2)
                # plt.scatter(ran,pmn,color=colours[di,:],s=10)
        dmean[tf] = np.median(amps)
        tmean[tf] = np.median(times)
        inplume_amp[tf] = np.median(plamps)
    plt.figure(2) 
    
    
    #plt.scatter(tmean,dmean,color=colours[di,:],s=20,label=d,alpha=.5)
    #plt.scatter(np.median(tmean),np.median(dmean),color=colours[di,:],s=100)
    pltx = np.zeros((len(tmean),2))
    pltx[:,0] = tmean
    plty = np.append(dmean[:,np.newaxis],inplume_amp[:,np.newaxis],axis=1)
    plt.plot(pltx.T,plty.T,color=colours[di,:],alpha=.25)
    plt.scatter(np.median(pltx,axis=0),np.median(plty,axis=0),100,color=colours[di,:])
    
    #plt.scatter(0,np.median(inplume_amp),color=colours[di,:],s=100)
    plt.plot([0,np.median(tmean)],[np.median(inplume_amp),np.median(dmean)],color=colours[di,:],linewidth=5)
plt.legend()
plt.xlabel('time of first turn (s)')
plt.ylabel('dF/F')
plt.ylim([0,.65])
plt.savefig(os.path.join(savedir,'RDP_first_Amp.pdf'))
#%% Heading and phase after last RDP turning point
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']

minssize = 1
xlabs = ['heading (rad)','rdp vector (rad)']
min_num = 3
for di,d in enumerate(dchoices):
    plt.figure(di)
    for s in range(1,3):
        plt.subplot(1,2,s)
        plt.xlim([-np.pi,np.pi])
        plt.ylim([-np.pi,np.pi])
        plt.plot([-np.pi,np.pi],[0,0],color='k')
        plt.plot([0,0],[-np.pi,np.pi],color='k')
        plt.plot([-np.pi/2,-np.pi/2],[-np.pi,np.pi],color='r',linestyle='--')
        plt.plot([-np.pi,np.pi],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
        plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,-0,90,180])
        plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,-0,90,180])
        plt.ylabel(d+' phase (deg)')
        plt.xlabel(xlabs[s-1])
        plt.gca().set_aspect('equal')
    data = all_flies[d]
    for tf in data:
        cxa = data[tf].cxa
        jumps = cxa.get_jumps()
        x = cxa.ft2['ft_posx'].to_numpy()
        y = cxa.ft2['ft_posy'].to_numpy()
        angvel = ug.get_ang_velocity(cxa.ft2['ft_heading'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
        angvelsmth = nd.gaussian_filter(angvel,3)**2
        phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
        eb = cxa.pdat['offset_eb_phase'].to_numpy()
        heading = cxa.ft2['ft_heading'].to_numpy()
        u = ug()
        vx,vy,vd = u.get_velocity(x,y,cxa.pv2['relative_time'].to_numpy())
        for ju in jumps:
            dx = np.arange(ju[1],ju[2])
            tx = x[dx]
            ty = y[dx]
            tx = tx-tx[0]
            ty = ty-ty[0]
            tp = phase[dx]
            te = eb[dx]
            tvd = vd[dx]
            th = heading[dx]
            av = angvel[dx]
            avs = angvelsmth[dx]
            
            xy = np.append(tx[:,np.newaxis],ty[:,np.newaxis],1)
            r = rdp(xy,epsilon=2)
            rdx = np.where(np.in1d(tx,r[:,0]))[0]
            # still = np.abs(av)<.1
            
            # sb,sbs = ug.find_blocks(still,mergeblocks=True,merg_threshold=2)
            # sb = sb[sbs>=minssize]
            # sbs = sbs[sbs>=minssize]
            
            if len(rdx)>min_num:
                sdx = np.arange(rdx[-2],min(rdx[-2]+3,len(tp)))
                pmn = stats.circmean(tp[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
                hmn = stats.circmean(th[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
                rvec = r[-1,:]-r[-2,:]
                rcom = rvec[0]*1j+rvec[1]
                ran = np.angle(rcom)*cxa.side*-1
                
                plt.subplot(1,2,1)
                plt.scatter(hmn,pmn,color=colours[di,:],s=10)
                plt.subplot(1,2,2)
                plt.scatter(ran,pmn,color=colours[di,:],s=10)
for di,d in enumerate(dchoices):
    plt.figure(di)
    plt.savefig(os.path.join(savedir,'RDP_last_'+d+str(min_num)+'.pdf'))

#%% Heading and phase for other points
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']

minssize = 1
xlabs = ['heading (deg)','rdp vector (deg)']
for di,d in enumerate(dchoices):
    plt.figure(di)
    for s in range(1,3):
        plt.subplot(1,2,s)
        plt.xlim([-np.pi,np.pi])
        plt.ylim([-np.pi,np.pi])
        plt.plot([-np.pi,np.pi],[0,0],color='k')
        plt.plot([0,0],[-np.pi,np.pi],color='k')
        plt.plot([-np.pi/2,-np.pi/2],[-np.pi,np.pi],color='r',linestyle='--')
        plt.plot([-np.pi,np.pi],[-np.pi/2,-np.pi/2],color='r',linestyle='--')
        plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,-0,90,180])
        plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,-0,90,180])
        plt.ylabel(d+' phase (deg)')
        plt.xlabel(xlabs[s-1])
        plt.gca().set_aspect('equal')
    data = all_flies[d]
    for tf in data:
        cxa = data[tf].cxa
        jumps = cxa.get_jumps()
        x = cxa.ft2['ft_posx'].to_numpy()
        y = cxa.ft2['ft_posy'].to_numpy()
        angvel = ug.get_ang_velocity(cxa.ft2['ft_heading'].to_numpy(),cxa.pv2['relative_time'].to_numpy())
        angvelsmth = nd.gaussian_filter(angvel,3)**2
        phase = cxa.pdat['offset_fsb_upper_phase'].to_numpy()
        eb = cxa.pdat['offset_eb_phase'].to_numpy()
        heading = cxa.ft2['ft_heading'].to_numpy()
        u = ug()
        vx,vy,vd = u.get_velocity(x,y,cxa.pv2['relative_time'].to_numpy())
        for ju in jumps:
            dx = np.arange(ju[1],ju[2])
            tx = x[dx]
            ty = y[dx]
            tx = tx-tx[0]
            ty = ty-ty[0]
            tp = phase[dx]
            te = eb[dx]
            tvd = vd[dx]
            th = heading[dx]
            av = angvel[dx]
            avs = angvelsmth[dx]
            
            xy = np.append(tx[:,np.newaxis],ty[:,np.newaxis],1)
            r = rdp(xy,epsilon=2)
            xdf = np.where(np.diff(tx)==0)[0]+1
            tx[xdf]= -1
            rdx = np.where(np.logical_and(np.in1d(tx,r[:,0]),np.in1d(ty,r[:,1])))[0]
            # still = np.abs(av)<.1
            
            # sb,sbs = ug.find_blocks(still,mergeblocks=True,merg_threshold=2)
            # sb = sb[sbs>=minssize]
            # sbs = sbs[sbs>=minssize]
            
            if len(rdx)>4:
                rdx2 = rdx[1:-1]
                r2 = r[1:-1,:]
                for ir in range(len(rdx2)-2):
                    sdx = np.arange(rdx2[ir+1],min(rdx2[ir+1]+3,rdx2[ir+2]))
                    sdx = np.arange(max(rdx2[ir+1]-3,rdx[ir]),rdx2[ir+1])
                    pmn = stats.circmean(tp[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
                    hmn = stats.circmean(th[sdx],high=np.pi,low=-np.pi)*cxa.side*-1
                    rvec = r2[ir+1,:]-r2[ir,:]
                    rcom = rvec[0]*1j+rvec[1]
                    ran = np.angle(rcom)*cxa.side*-1
                
                    plt.subplot(1,2,1)
                    plt.scatter(hmn,pmn,color=colours[di,:],s=10)
                    plt.subplot(1,2,2)
                    plt.scatter(ran,pmn,color=colours[di,:],s=10)


for di,d in enumerate(dchoices):
    plt.figure(di)
    plt.savefig(os.path.join(savedir,'RDP_inter_'+d+'.pdf'))
    
#%% fluorescence recovery timecourse
savedir = r'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'

colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255
plt.close('all')
minlen = 50
x = np.arange(0,600)/10
dchoices = ['hDeltaJ','hDeltaC','FC2']

for di,d in enumerate(dchoices):
    data = all_flies[d]
    meandat = np.zeros((600,len(data)))
    for f in data:
        cxa = data[f].cxa
        
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
        meandat[:,int(f)] = np.nanmean(all_js,axis=1)

    
    plt.plot(x[:200],meandat[:200,:],color=colours[di,:],alpha=0.3)
    plt.plot(x[:200],np.nanmean(meandat[:200,:],1),color=colours[di,:],alpha=1,linewidth=3)
plt.ylabel('mean dF/F',fontsize=15)
plt.xlabel('time from plume exit (s)',fontsize=15)
plt.xlim([0,20])
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(0,1,.2),fontsize=15)
plt.savefig(os.path.join(savedir,'FluorRecoveryAllEntries.png'))
plt.savefig(os.path.join(savedir,'FluorRecoveryAllEntries.pdf'))
#%% fluorescence recovery timecourse ACV pulses

colours = np.array([[81,61,204],#hDeltaJ
                    [49,99,125],#hDeltaC
                    [60, 181, 136],#hDeltaH
                    [81,156,205] # FC2
                    ])/255
plt.close('all')
minlen = 50
x = np.arange(0,600)/10
dchoices = ['hDeltaJ_ACV_pulses','hDeltaC_ACV_pulses','hDeltaH_ACV_pulses']

for di,d in enumerate(dchoices):
    data = all_flies[d]
    meandat = np.zeros((600,len(data)))
    for f in data:
        cxa = data[f].cxa
        
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
        meandat[:,int(f)] = np.nanmean(all_js,axis=1)

    
    plt.plot(x[:200],meandat[:200,:],color=colours[di,:],alpha=0.3)
    plt.plot(x[:200],np.nanmean(meandat[:200,:],1),color=colours[di,:],alpha=1,linewidth=3)
plt.ylabel('mean dF/F',fontsize=15)
plt.xlabel('time from plume exit (s)',fontsize=15)
plt.xlim([0,20])
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(0,1,.2),fontsize=15)
plt.savefig(os.path.join(savedir,'FluorRecoveryACV_pulses.png'))
plt.savefig(os.path.join(savedir,'FluorRecoveryACV_pulses.pdf'))
#%% oct pulses

colours = np.array([[81,61,204],#hDeltaJ
                    [49,99,125],#hDeltaC
                    [60, 181, 136],#hDeltaH
                    [81,156,205] # FC2
                    ])/255
plt.close('all')
minlen = 50
x = np.arange(0,600)/10
dchoices = ['hDeltaJ_Oct_pulses','hDeltaC_Oct_pulses','hDeltaH_Oct_pulses']

for di,d in enumerate(dchoices):
    data = all_flies[d]
    meandat = np.zeros((600,len(data)))
    for f in data:
        cxa = data[f].cxa
        
        amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        jumps = cxa.get_entries_exits_like_jumps(odour='Oct')
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
        meandat[:,int(f)] = np.nanmean(all_js,axis=1)

    
    plt.plot(x[:200],meandat[:200,:],color=colours[di,:],alpha=0.3)
    plt.plot(x[:200],np.nanmean(meandat[:200,:],1),color=colours[di,:],alpha=1,linewidth=3)
plt.ylabel('mean dF/F',fontsize=15)
plt.xlabel('time from plume exit (s)',fontsize=15)
plt.xlim([0,20])
plt.xticks(np.arange(0,21,5),fontsize=15)
plt.yticks(np.arange(0,1,.2),fontsize=15)
plt.savefig(os.path.join(savedir,'FluorRecoveryOct_pulses.png'))
plt.savefig(os.path.join(savedir,'FluorRecoveryOct_pulses.pdf'))

#%% Fluorescence recovery, first plume entry
colours = np.array([[81,61,204],#hDeltaJ
                    [49,99,125],#hDeltaC
                    [60, 181, 136],#hDeltaH
                    [81,156,205] # FC2
                    ])/255
plt.close('all')
minlen = 50
x = np.arange(0,600)/10
dchoices = ['hDeltaJ','hDeltaC','hDeltaH','FC2']
dtimebase = np.arange(0,600,5)/10
for di,d in enumerate(dchoices):
    data = all_flies[d]
    meandat = np.zeros((600,len(data)))
    for f in data:
        cxa = data[f].cxa
        
        amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        jumps = cxa.get_entries_exits_like_jumps()
        
        
        j = jumps[0,:]
        dx = np.arange(j[1],j[2])
        tamp = amp[dx]
        amplen  = np.min([600,len(tamp)])
        meandat[:amplen,f] = tamp[:amplen]
        
        # for i,j in enumerate(jumps):
        #     dx = np.arange(j[1],j[2])
        #    # plt.plot(amp[dx],color='k',alpha=0.3)
        #     tamp = amp[dx]
        #     amplen  = np.min([600,len(tamp)])
        #     all_js[:amplen,i] = tamp[:amplen]
            
        # all_js[all_js==0] = np.nan    
        # meandat[:,int(f)] = np.nanmean(all_js,axis=1)

    meandat[meandat[:]==0] = np.nan
    meandat_d = np.zeros((len(dtimebase),meandat.shape[1]))
    for m in range(meandat.shape[1]):
        meandat_d[:,m] = np.interp(dtimebase,x,meandat[:,m])
    
    plt.figure(1)
    plt.plot(x[:200],meandat[:200,:],color=colours[di,:],alpha=0.3)
    plt.plot(x[:200],np.nanmean(meandat[:200,:],1),color=colours[di,:],alpha=1,linewidth=3)
    
    plt.figure(2)
    plt.plot(dtimebase,meandat_d,color=colours[di,:],alpha=0.3)
    plt.plot(dtimebase,np.nanmean(meandat_d,1),color=colours[di,:],alpha=1,linewidth=3)
    
for i in range(1,3):
    plt.figure(i)
    plt.ylabel('mean dF/F',fontsize=15)
    plt.xlabel('time from plume exit (s)',fontsize=15)
    plt.xlim([0,20])
    plt.xticks(np.arange(0,21,5),fontsize=15)
    plt.yticks(np.arange(0,1,.2),fontsize=15)
plt.figure(1)
plt.savefig(os.path.join(savedir,'FluorRecoveryFirstEntry.png'))
plt.savefig(os.path.join(savedir,'FluorRecoveryFirstEntry.pdf'))
#%% Fluorescence recovery, first plume with entry
colours = np.array([[81,61,204],#hDeltaJ
                    [49,99,125],#hDeltaC
                    [60, 181, 136],#hDeltaH
                    [81,156,205] # FC2
                    ])/255
plt.close('all')
minlen = 50
x = np.arange(-300,600)/10
dchoices = ['hDeltaJ','hDeltaC','hDeltaH','FC2']
dtimebase = np.arange(-300,600,5)/10
for di,d in enumerate(dchoices):
    data = all_flies[d]
    meandat = np.zeros((len(x),len(data)))
    for f in data:
        cxa = data[f].cxa
        
        amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
        jumps = cxa.get_entries_exits_like_jumps()
        
        
        j = jumps[0,:]
        dx = np.arange(np.max([j[0],j[1]-300]),np.min([j[2],j[1]+600]))
        dxt = (dx-j[1])/10
        tamp = amp[dx]
        
        
        
        
        
        
        # amplen  = np.min([len(x),len(tamp)])
        # tamp = tamp[:amplen]
        # dxt = dxt[:amplen]
        tdx = np.in1d(x,dxt)
        meandat[tdx,f] = tamp
        
        # for i,j in enumerate(jumps):
        #     dx = np.arange(j[1],j[2])
        #    # plt.plot(amp[dx],color='k',alpha=0.3)
        #     tamp = amp[dx]
        #     amplen  = np.min([600,len(tamp)])
        #     all_js[:amplen,i] = tamp[:amplen]
            
        # all_js[all_js==0] = np.nan    
        # meandat[:,int(f)] = np.nanmean(all_js,axis=1)

    meandat[meandat[:]==0] = np.nan
    meandat_d = np.zeros((len(dtimebase),meandat.shape[1]))
    for m in range(meandat.shape[1]):
        meandat_d[:,m] = np.interp(dtimebase,x,meandat[:,m])
    
    plt.figure(1)
    plt.plot(x,meandat,color=colours[di,:],alpha=0.3)
    plt.plot(x,np.nanmean(meandat,1),color=colours[di,:],alpha=1,linewidth=3)
    
    plt.figure(2)
    plt.plot(dtimebase,meandat_d,color=colours[di,:],alpha=0.3)
    plt.plot(dtimebase,np.nanmean(meandat_d,1),color=colours[di,:],alpha=1,linewidth=3)
    
for i in range(1,3):
    plt.figure(i)
    plt.ylabel('mean dF/F',fontsize=15)
    plt.xlabel('time from plume exit (s)',fontsize=15)
    plt.xlim([-2,20])
    plt.xticks(np.arange(0,21,5),fontsize=15)
    plt.yticks(np.arange(0,1,.2),fontsize=15)
#%% Average return duration
plt.close('all')
dchoices = ['hDeltaJ','hDeltaC','FC2']
all_jumps = np.array([])
for d in dchoices:
    data = all_flies[d]
    for f in data:
        cxa = data[f].cxa
        jumps = cxa.get_jumps()
        jlen = (jumps[:,2]-jumps[:,1])/10
        all_jumps = np.append(all_jumps,jlen)
        
plt.hist(all_jumps,np.arange(0,60,1),color='k',density=True)
medjump = np.median(all_jumps)
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
#%% hDeltaC error from expected phase


f =5
data = all_flies['hDeltaC']
cxa = data[f].cxa
epg = cxa.pdat['phase_eb']
fsb = cxa.pdat['phase_fsb_upper']
epg = ug.savgol_circ(epg,20,3)
fsb = ug.savgol_circ(fsb,20,3)
#pred1 = ug.circ_subtract(epg,(3/4)*np.pi)
#pred2 = ug.circ_subtract(epg,-(3/4)*np.pi)
pred1 = ug.circ_subtract(epg,(3/4)*np.pi)
pred2 = ug.circ_subtract(epg,-(3/4)*np.pi)
t = np.arange(0,len(epg))/10
ins = cxa.ft2['instrip'].to_numpy()
plt.figure()
plt.scatter(t,pred1,color='g',s=3)
plt.scatter(t,pred2,color='b',s=3)
plt.scatter(t,epg,color=[0.5,0.5,0.5],s=3)
plt.scatter(t,fsb,color='k',s=3)

plt.figure()
err1 = np.abs(ug.circ_subtract(fsb,epg))
err2 = np.abs(ug.circ_subtract(fsb,pred1))
err3 = np.abs(ug.circ_subtract(fsb,pred2))

# err1 = ug.circ_subtract(fsb,epg)
# err2 = ug.circ_subtract(fsb,pred1)
# err3 = ug.circ_subtract(fsb,pred2)
total_error = np.min(np.append(err2[:,np.newaxis],err3[:,np.newaxis],axis=1),axis=1)

plt.scatter(t,err1,color='k',s=3)
plt.scatter(t,err2,color='g',s=3)
plt.scatter(t,err3,color='b',s=3)
plt.scatter(t,total_error,color='m',s=3)
plt.plot(t,ins*np.pi,color='r')
plt.yticks([0,np.pi/4,np.pi/2,3*np.pi/4,np.pi],labels =[0,45,90,135,180])


x = cxa.ft2['ft_posx'].to_numpy()
y = cxa.ft2['ft_posy'].to_numpy()
x,y = ug.fictrac_repair(x,y)





colour  =total_error.copy()
cmin = 0 
cmax= 2
colour[colour<cmin] = cmin
colour[colour>cmax] = cmax 
colour[0] = cmin # Hacky fix to get over scaling issues
colour[-1] = cmax # Hacky fix to get over scaling issues

plt.figure()
ax = plt.gca()
uplt.coloured_line(x,y,colour,ax,cmap='coolwarm')
plt.show()
plt.gca().set_aspect('equal')
plt.scatter(x[ins>0],y[ins>0],color=[.7,.7,.7])

cxa.plot_traj_arrow_heat(['fsb_upper'],total_error,a_sep=20,traindat=False,fulldat=True,colormap='coolwarm',cmin=0,cmax=2)
#%% Default mode error during jump returns
savedir = 'Y:\Data\FCI\FCI_summaries\hDelta_FC2_Comparison'
plt.close('all')
data = all_flies['hDeltaC']
inttime = 50
plotdata = np.zeros((inttime,len(data)))
plotamp = np.zeros((inttime,len(data)))
trajdata = np.zeros((inttime,2,len(data)))
for f in data:
    cxa = data[f].cxa
    etp = data[f]
    epg = cxa.pdat['phase_eb']
    fsb = cxa.pdat['phase_fsb_upper']
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    pred1 = ug.circ_subtract(epg,(3/4)*np.pi)
    pred2 = ug.circ_subtract(epg,-(3/4)*np.pi)
    
    # pred1 = ug.circ_subtract(epg,(2/3)*np.pi)
    # pred2 = ug.circ_subtract(epg,-(2/3)*np.pi)
    
    err2 = np.abs(ug.circ_subtract(fsb,pred1))
    err3 = np.abs(ug.circ_subtract(fsb,pred2))
    total_error = np.min(np.append(err2[:,np.newaxis],err3[:,np.newaxis],axis=1),axis=1)
    total_error = 180*total_error/np.pi
    jumps = cxa.get_jumps()
    dmat = np.zeros((inttime,len(jumps)))
    dmata = np.zeros((inttime,len(jumps)))
    
    phase,traj,_ = etp.trajectory_mean(time_threshold2=0,regions=['eb','fsb_upper'],bins=inttime)

    trajdata[:,:,f] = np.mean(traj[inttime:],axis=2)
    
    for ij,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
        oldtime = dx-dx[0]
        newtime = np.linspace(oldtime[0],oldtime[-1],inttime)
        terr = total_error[dx]
        dmat[:,ij] = np.interp(newtime,oldtime,terr)
        dmata[:,ij] = np.interp(newtime,oldtime,amp[dx])
        
        
    plotdata[:,f] = np.mean(dmat,axis=1)
    plotamp[:,f] = np.mean(dmata,axis=1)
    
fig,ax1 = plt.subplots()
ax1.plot(plotdata,color='k',alpha=.2)
ax1.plot(np.mean(plotdata,axis=1),color='k')
ax1.set_xlabel('Pseudo time',fontsize=15)
ax1.set_ylabel('Error from default (deg)',fontsize=15)
ax1.set_xticks([0,25,49],labels = [0,.5,1],fontsize=15)
ax1.set_yticks(np.arange(30,91,10),labels=np.arange(30,91,10),fontsize=15)
ax1.set_xlim([0,49])
ax2 = ax1.twinx()
ax2.plot(plotamp,color='r',alpha=.2)
ax2.plot(np.mean(plotamp,axis=1),color='r')
ax2.set_ylabel('mean dF/F',fontsize=15)
ax2.set_ylim([0,0.6])
ax2.set_yticks([0,.25,.5],labels = [0,.25,.5],fontsize=15)
trajmean = np.mean(trajdata,axis=2)



plt.savefig(os.path.join(savedir,'hDeltaC_default_over_return.pdf'))

plt.figure()
uplt.coloured_line_simple(trajmean[:,0],trajmean[:,1],np.mean(plotdata,axis=1),cmap='coolwarm',cmin=45,cmax=90)
plt.plot([-3,-3],[0,35],color=[.7,.7,.7])
uplt.coloured_line_simple(trajmean[:,0]+7,trajmean[:,1],np.mean(plotamp,axis=1),cmap='coolwarm',cmin=0,cmax=.5)    
plt.plot(np.array([-3,-3])+7,[0,35],color=[.7,.7,.7])
plt.gca().set_aspect('equal')
plt.savefig(os.path.join(savedir,'hDeltaC_default_over_return_traj.pdf'))

        


#%%

amp = (amp-np.mean(amp))/np.std(amp)
cxa.plot_traj_arrow_heat(['fsb_upper'],amp,a_sep=20,traindat=False,fulldat=True,colormap='coolwarm',cmin=-2,cmax=2)

#%% 

data = all_flies['hDeltaC']
bins = np.arange(-3,3,.2)
plotbins = bins[1:]-np.mean(np.diff(bins))/1

plotdata = np.zeros((len(bins)-1,len(data)))
for f in data:
    cxa = data[f].cxa
    
    epg = cxa.pdat['phase_eb']
    fsb = cxa.pdat['phase_fsb_upper']
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    amp = (amp-np.mean(amp))/np.std(amp)
    pred1 = ug.circ_subtract(epg,(3/4)*np.pi)
    pred2 = ug.circ_subtract(epg,-(3/4)*np.pi)
    err1 = np.abs(ug.circ_subtract(fsb,epg))
    err2 = np.abs(ug.circ_subtract(fsb,pred1))
    err3 = np.abs(ug.circ_subtract(fsb,pred2))

    # err1 = ug.circ_subtract(fsb,epg)
    # err2 = ug.circ_subtract(fsb,pred1)
    # err3 = ug.circ_subtract(fsb,pred2)
    total_error = np.min(np.append(err2[:,np.newaxis],err3[:,np.newaxis],axis=1),axis=1)
    ind = np.digitize(amp,bins)
    for i in range(len(bins)-1):
        plotdata[i,f] = np.mean(total_error[ind==i])
    
    
    
plotdata= 180*plotdata/np.pi
plt.plot(plotbins,plotdata,color='k',alpha=.25)
plt.plot(plotbins,np.nanmean(plotdata,axis=1),color='k')    
plt.xlabel('mean dF/F (zscore')
plt.ylabel('mean error (deg)')
#%%
colours = np.array([[81,61,204],#hDeltaJ
                    [49,99,125],#hDeltaC
                    #[60, 181, 136],#hDeltaH
                    [81,156,205] # FC2
                    ])/255
dchoices = ['hDeltaJ','hDeltaC']
amp = ['dFF','pva']
for a in amp:
    plt.figure()
    for di, d in enumerate(dchoices):
        data = all_flies[d]
        bins = np.arange(-3,3,.25)
        plotbins = bins[1:]-np.mean(np.diff(bins))/2
        
        plotdata = np.zeros((len(bins)-1,len(data)))
        for f in data:
            cxa = data[f].cxa
            
            epg = cxa.pdat['phase_eb']
            fsb = cxa.pdat['phase_fsb_upper']
            if a =='dFF':
                amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
            elif a=='pva':
                amp = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
            amp = (amp-np.mean(amp))/np.std(amp)
            pred1 = ug.circ_subtract(epg,(3/4)*np.pi)
            pred2 = ug.circ_subtract(epg,-(3/4)*np.pi)
            err1 = np.abs(ug.circ_subtract(fsb,epg))
            err2 = np.abs(ug.circ_subtract(fsb,pred1))
            err3 = np.abs(ug.circ_subtract(fsb,pred2))
        
            # err1 = ug.circ_subtract(fsb,epg)
            # err2 = ug.circ_subtract(fsb,pred1)
            # err3 = ug.circ_subtract(fsb,pred2)
            total_error = np.abs(err1)
            ind = np.digitize(amp,bins)
            for i in range(len(bins)-1):
                plotdata[i,f] = np.median(total_error[ind==i])
                plotdata[i,f] = stats.mode(np.round(total_error[ind==i],decimals=1))[0]
            
            
        plotdata= 180*plotdata/np.pi
        plt.plot(plotbins,plotdata,color=colours[di],alpha=.2)
        plt.plot(plotbins,np.nanmean(plotdata,axis=1),color=colours[di,:],linewidth=2)    
        
        
    if a =='dFF':
        plt.xlabel('mean dF/F (zscore)',fontsize=15)
    elif a=='pva':
        plt.xlabel('pva (zscore)',fontsize=15)
    plt.ylabel('modal FSB-EPG phase (deg)',fontsize=15)
    
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim([-3,3])
    plt.ylim([0,180])
    plt.savefig(os.path.join(savedir,  'ModalPhase_w_'+a+'.png' ))
#%% Default error across jumps for hDeltaC
a = 'dFF'
plt.close('all')
for f in data:
    cxa = data[f].cxa
    
    epg = cxa.pdat['phase_eb']
    fsb = cxa.pdat['phase_fsb_upper']
    if a =='dFF':
        amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    elif a=='pva':
        amp = ug.get_pvas(cxa.pdat['wedges_fsb_upper'])
    #amp = (amp-np.mean(amp))/np.std(amp)
    pred1 = ug.circ_subtract(epg,(3/4)*np.pi)
    pred2 = ug.circ_subtract(epg,-(3/4)*np.pi)
    err1 = np.abs(ug.circ_subtract(fsb,epg))
    total_error = np.abs(err1)
    jumps = cxa.get_jumps()
    plt.figure()
    for ij,j in enumerate(jumps):
        dx = np.arange(j[0],j[2])
        terror = total_error[dx]
        plt.scatter(amp[dx],terror,s=1,color='k')


#%% hDeltaJ recovery swing direction

dchoice = 'hDeltaJ'
data = all_flies[dchoice]
cxa = data[1].cxa
jumps = cxa.get_jumps()
phase = cxa.pdat['phase_fsb_upper']
eb = cxa.pdat['phase_eb']
plt.close('all')
for ij,j in enumerate(jumps):
    dx = np.arange(j[1],j[2])
    tp = phase[dx]
    teb = eb[dx]
    twed = np.flipud(cxa.pdat['wedges_fsb_upper'][dx,:].T)
    plt.figure()
    plt.subplot(2,1,1)
    plt.scatter(dx,tp,color='b',s=3)
    plt.scatter(dx,teb,color='k',s=3)
    plt.subplot(2,1,2)
    plt.imshow(twed,interpolation='None',aspect='auto',vmin=0,vmax=1)
#%% hDeltaJ default lag

dchoice = 'hDeltaJ'
data = all_flies[dchoice]
for f in data:
    cxa = data[f].cxa
    phase = cxa.pdat['phase_fsb_upper']
    eb = cxa.pdat['phase_eb']
    ins = cxa.ft2['instrip'].to_numpy()
    st = np.where(ins>0)[0][0]
    dx = np.arange(0,len(eb))
    lags,corr = ug.circular_xcorr_lags(phase[dx],eb[dx])
    plt.plot(lags/10,np.abs(corr),color='k')
    ldx = np.where(np.abs(lags/10)<10)[0]
    pks = np.argmax(np.abs(corr[ldx]))
    print(lags[ldx[pks]]/10)
plt.plot([0,0],[0,1],color='r')
plt.xlim([-20,20])
#%% Phase relationships upon first stop
plt.close('all')
for f in data:
    cxa = data[f].cxa
    u = ug()
    x = cxa.ft2['ft_posx'].to_numpy()
    y = cxa.ft2['ft_posy'].to_numpy()
    t = cxa.pv2['relative_time'].to_numpy()
    dx,dy,dd = u.get_velocity(x,y,t)
    dd = np.append(0,dd)
    jumps=cxa.get_jumps()
    phase= cxa.pdat['phase_fsb_upper']*-cxa.side
    eb = cxa.pdat['phase_eb']*-cxa.side
    
    phase= cxa.pdat['offset_fsb_upper_phase'].to_numpy()*-cxa.side
    eb = cxa.pdat['offset_eb_phase'].to_numpy()*-cxa.side
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    
    print(cxa.side)
    jumps = cxa.get_entries_exits_like_jumps()
    for j in jumps:
        #j = jumps[12,:]
        
        dx  = np.arange(j[1],j[2])
        #plt.plot(amp[dx])
        #plt.plot(dd[dx]/10)
        stops = dd[dx]<.5        
        
        # Compute phase at stops
        bstart,bsize = ug.find_blocks(stops,mergeblocks=True,merg_threshold=2)
        tphase = phase[dx]
        teb = eb[dx]
        minsize = 5 
        bstart = bstart[bsize>=minsize]
        bsize = bsize[bsize>=minsize]
        if len( bsize)==0:
            continue
        sphases = np.zeros((len(bsize),2))
        for ib,b in enumerate(bstart):
            #plt.plot([b,b+bsize[ib]],[1,1],color='k')
            dx2 = np.arange(b,b+bsize[ib])
            sphases[ib,0] = stats.circmean(tphase[dx2],high=np.pi,low=-np.pi)
            sphases[ib,1] = stats.circmean(teb[dx2],high=np.pi,low=-np.pi)
        
        # Compute phase/eb on entry
        dxe = np.arange(j[0]-2,j[0]+2)
        pent = stats.circmean(phase[dxe],high=np.pi,low=-np.pi)
        eent = stats.circmean(eb[dxe],high=np.pi,low=-np.pi)
        plt.figure(0)
        plt.scatter(ug.circ_subtract(pent,np.pi),sphases[0,0],color='b',s=5)
        plt.figure(1)
        plt.scatter(eent,sphases[0,0],color='k',s=5)
        plt.figure(2)
        plt.scatter(ug.circ_subtract(sphases[0,1],np.pi),sphases[0,0],color='r',s=5)
        
        
xlabs = ['minus fsb entry phase','eb entry phase','minus phase']    
for i in range(3):
    plt.figure(i)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.ylabel('Stopping phase')
    plt.xlabel(xlabs[i])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r',linestyle='--')
#%% Phase relationships upon recovery
plt.close('all')
flycolours =[[0.5,0.5,0.5],'r','g','b','m','k']
to_plume = np.array([],dtype=float)
to_plume_ent = np.array([],dtype=float)
fcol = np.array([],dtype=float)
aw_plume = np.array([],dtype=float)
aw_plume_ent = np.array([],dtype=float)
fcolaw = np.array([],dtype=float)
timethresh = 30
dchoice = 'hDeltaJ'
data = all_flies[dchoice]
for f in data:
    cxa = data[f].cxa
    u = ug()
    x = cxa.ft2['ft_posx'].to_numpy()
    y = cxa.ft2['ft_posy'].to_numpy()
    t = cxa.pv2['relative_time'].to_numpy()
    dx,dy,dd = u.get_velocity(x,y,t)
    dd = np.append(0,dd)
    jumps=cxa.get_jumps()
    jumps = cxa.get_entries_exits_like_jumps()
    jdx = (jumps[:,2]-jumps[:,1])>=timethresh
    jumps = jumps[jdx,:]
    phase= cxa.pdat['phase_fsb_upper']*-cxa.side
    eb = cxa.pdat['phase_eb']*-cxa.side
    
    phase= cxa.pdat['offset_fsb_upper_phase'].to_numpy()*-cxa.side
    eb = cxa.pdat['offset_eb_phase'].to_numpy()*-cxa.side
    eb = cxa.ft2['ft_heading'].to_numpy()*-cxa.side
    amp = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    
    print(cxa.side)
    #jumps = cxa.get_entries_exits_like_jumps()
    
    for j in jumps:
        #j = jumps[12,:]
        
        dx  = np.arange(j[1],j[2])
        tamp = amp[dx]
        tphase = phase[dx]
        teb = eb[dx]
        if len(tamp)<30:
            continue
        #plt.figure()
        #plt.plot(amp[dx])
        #plt.plot(dd[dx]/10)
        stops = dd[dx]<.5        
        
        amin = np.argmin(tamp[:50])
        am = np.argmax(tamp[amin:50])+amin
        #plt.scatter(am,tamp[am],color='k')
        
        mdx = np.arange(am-2,min(am+2,len(tphase)))
        mphase = stats.circmean(tphase[mdx],low=-np.pi,high=np.pi)
        meb  = stats.circmean(teb[mdx],low=-np.pi,high=np.pi)
        
        # # Compute phase at stops
        # bstart,bsize = ug.find_blocks(stops,mergeblocks=True,merg_threshold=2)
        # tphase = phase[dx]
        # 
        # minsize = 5 
        # bstart = bstart[bsize>=minsize]
        # bsize = bsize[bsize>=minsize]
        # if len( bsize)==0:
        #     continue
        # sphases = np.zeros((len(bsize),2))
        # for ib,b in enumerate(bstart):
        #     #plt.plot([b,b+bsize[ib]],[1,1],color='k')
        #     dx2 = np.arange(b,b+bsize[ib])
        #     sphases[ib,0] = stats.circmean(tphase[dx2],high=np.pi,low=-np.pi)
        #     sphases[ib,1] = stats.circmean(teb[dx2],high=np.pi,low=-np.pi)
        
        # # Compute phase/eb on entry
        dxe = np.arange(j[0]-2,j[0]+2)
        pent = stats.circmean(phase[dxe],high=np.pi,low=-np.pi)
        eent = stats.circmean(eb[dxe],high=np.pi,low=-np.pi)
        if mphase<0:
           # print(mphase)
            to_plume = np.append(to_plume,mphase)
            to_plume_ent = np.append(to_plume_ent,eent)
            fcol = np.append(fcol,f)
        elif mphase>0:
            #print(mphase)
            aw_plume = np.append(aw_plume,mphase)
            
            print(mphase)
            print(aw_plume[-1])
            aw_plume_ent = np.append(aw_plume_ent,eent)
            fcolaw = np.append(fcolaw,f)
        #plt.figure(0)
        #plt.scatter(pent,mphase,color=flycolours[f],s=20)
        plt.figure(f)
        plt.scatter(eent,mphase,color='b',s=30)
        plt.scatter(ug.circ_subtract(meb,np.pi),mphase,color='k',s=30)
        plt.plot([eent,ug.circ_subtract(meb,np.pi) ],[mphase,mphase],color=[0.5,0.5,0.5],zorder=-1)
        #plt.figure(2)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r',linestyle='--')
        #plt.scatter(ug.circ_subtract(meb,np.pi),mphase,color=flycolours[f],s=20)
    
plt.figure()
plt.scatter(to_plume_ent,to_plume,c=fcol,cmap='tab10')
plt.scatter(aw_plume_ent,aw_plume,color='r')
plt.scatter(aw_plume_ent,aw_plume-2*np.pi,color='r')
plt.xlim([-np.pi,np.pi])
plt.ylim([-np.pi,np.pi])
plt.ylabel('recovery phase')
plt.xlabel('eb entry phase')
#plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r',linestyle='--')
#plt.plot([-3/4*np.pi,np.pi],[-np.pi,3/4*np.pi])
#plt.plot([-2/3*np.pi,np.pi],[-np.pi,2/3*np.pi])
# xlabs = ['fsb entry phase','eb entry phase','minus phase']    
# for i in range(3):
#     plt.figure(i)
#     plt.xlim([-np.pi,np.pi])
#     plt.ylim([-np.pi,np.pi])
#     plt.ylabel('recovery phase')
#     plt.xlabel(xlabs[i])
#     plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r',linestyle='--')
#%% Evidence for two bumps for hDeltaJ
plt.close('all')
dchoice = 'hDeltaJ'
f =0
cxa = all_flies[dchoice][f].cxa
jumps = cxa.get_entries_exits_like_jumps()
jumps = cxa.get_jumps()
binsize = .25
binsize = int(binsize*10)
offset=.5
wedges = cxa.pdat['wedges_fsb_upper']
wedges_eb = cxa.pdat['wedges_eb']
for ij,j in enumerate(jumps):
    dxe  = np.arange(j[0]-5,j[2])
    
    twed = wedges[dxe,:]
    twede = wedges_eb[dxe,:]
    n_bins = twed.shape[0]//binsize
    onbin = (j[1]-j[0]+5)//binsize
    wed_trimmed =twed[:n_bins*binsize,:]
    wed_binned = wed_trimmed.reshape(n_bins,binsize,16,-1).mean(axis=1)
    
    wed_trimmede =twede[:n_bins*binsize,:]
    wed_binnede = wed_trimmede.reshape(n_bins,binsize,16,-1).mean(axis=1)
    
    
    dxi  = np.arange(j[0],j[1])
    
        # plt.figure()
        # plt.subplot(1,2,1)
        # plt.imshow(twed)
        # plt.subplot(1,2,2)
        # plt.imshow(wed_binned)
    plt.figure()
    if n_bins>3:
        for w in range(len(wed_binned)):
            tw = wed_binned[w,:]
            twe = wed_binnede[w,:]
            x = np.arange(w*16,(w+1)*16)
            plt.plot(x,twe,color='k')
            plt.plot(x,tw+offset,color='b')
            if w==1:
                plt.plot([x[-1]+.5,x[-1]+.5],[-.1,1.2],color='r')
            if w==onbin:
                plt.plot([x[-1]+.5,x[-1]+.5],[-.1,1.2],color='r')
            if np.mod(w,2)==0:
                plt.fill_between(x,x*0-.1,x*0+1.2,color=[0.8,0.8,0.8],zorder=-1)
                plt.text(x[7],-.2,str((w-1)*binsize/10))
        plt.title(str(ij))
        
cxa.jump_return_details()
#%% EPG bump offset and bias
dchoices = ['hDeltaJ','hDeltaC']
#plt.close('all')
cnt=0
for d in dchoices:
    data= all_flies[d]
    for f in data:
        cxa = data[f].cxa
        ins = cxa.ft2['instrip'].to_numpy()
        ins = np.arange(0,np.where(ins>0)[0][0])
        cxa.get_jumps()
        epg = cxa.pdat['phase_eb'][ins]
        heading = cxa.ft2['ft_heading'].to_numpy()[ins]
        offset= ug.circ_subtract(epg,heading)
        #plt.figure()
        #plt.hist(offset,bins=np.linspace(-np.pi,np.pi,100))
        omn = stats.circmean(offset,high=-np.pi,low=np.pi)
        plt.figure(101)
        plt.scatter(-cxa.side,omn,color='k')
        # plt.figure(cnt)
        # x = cxa.ft2['ft_posx'].to_numpy()
        # y = cxa.ft2['ft_posy'].to_numpy()
        # x,y = cxa.fictrac_repair(x,y)
        # plt.plot(x[ins],y[ins],color='k')
        # plt.gca().set_aspect('equal')
        # plt.title('Offset: ' +str(180*omn/np.pi)+' Side: ' +str(-cxa.side))
        
        cnt+=1
plt.figure(101)
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180])
plt.ylim([-np.pi,np.pi])
plt.ylabel('Mean heading EPG offset (degrees)')
plt.xticks([-1,1],labels=['left tracker','right tracker'])
plt.plot([-1,1],[0,0],color='r')

