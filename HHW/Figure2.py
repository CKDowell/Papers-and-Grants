# -*- coding: utf-8 -*-
"""
Created on Mon May 20 13:50:44 2024

@author: dowel

This script will output Figure 2 for the HHW application

Panels:
    Example

"""

#%%
from analysis_funs.regression import fci_regmodel
from analysis_funs.optogenetics import opto 
import os
import matplotlib.pyplot as plt 
from src.utilities import imaging as im
from skimage import io, data, registration, filters, measure
from scipy import signal as sg
from analysis_funs.CX_imaging import CX
import numpy as np
import pandas as pd
plt.rcParams['pdf.fonttype'] = 42 
#%% 
datadirs = ["Y:\Data\FCI\Hedwig\\TH_GC7f\\240306\\f1\\Trial2",
            "Y:\Data\FCI\Hedwig\\TH_GC7f\\240314\\f2\\Trial2",
            "Y:\\Data\\FCI\\Hedwig\\TH_GC7f\\240529\\f3\\Trial1"]
savedir = "Y:\\Applications\\HHW\\Figures\\Figure2"
#%% Plot example fluorescence with velocity
dchoice = 0
colours = np.array([[79,0,148],[212,0,149],[255,170,239]])/255
datadir = datadirs[dchoice]
d = datadir.split("\\")
name = d[-3] + '_' + d[-2] + '_' + d[-1]
cx = CX(name,['fsb_layer'],datadir)
pv2, ft, ft2, ix = cx.load_postprocessing()


regchoice = ['odour onset', 'odour offset', 'in odour', 
                            'cos heading pos','cos heading neg', 'sin heading pos', 'sin heading neg',
                        'angular velocity pos','angular velocity neg',
                        'translational vel',
                        'ramp down since exit','ramp to entry']
plt.figure(figsize=(18,8))
ilist = [1,2,0]

for i in range(3):
    tstr = str(ilist[i]) + '_fsb_layer'
    y = pv2[tstr].to_numpy()
    fc = fci_regmodel(y,ft2,pv2)
    fc.rebaseline(span=500,plotfig=False)
    fc.run(regchoice)
    if i==0:
       istr = ft2['instrip'].astype(float).to_numpy()
       ts = fc.ts.to_numpy()
       plt.fill_between(ts,istr*3.5,istr*0,color= [.7,0.7,0.7],linewidth=0) 
    
    if ilist[i]>0:
        plt.plot(fc.ts,fc.ca+ilist[i],color = colours[ilist[i],:])
    else:
        fc.plot_flur_w_regressors(['translational vel'],cacol=colours[ilist[i],:])
plt.xlim([740,860])
plt.savefig(os.path.join(savedir,'EgFluor.pdf'))
#%% Plot example trajectory of animal
fc.simple_trajectory(dx=np.logical_and(ts>740,ts<860 ))
plt.savefig(os.path.join(savedir,'EgTraj.pdf'))
#%% Pearson correlation for each animal
regchoice = ['odour onset', 'odour offset', 'in odour', 
                        'angular velocity abs','translational vel']

rhos = np.zeros((len(regchoice),3,len(datadirs)))
for ids,datadir in enumerate(datadirs):
    
    d = datadir.split("\\")
    name = d[-3] + '_' + d[-2] + '_' + d[-1]
    cx = CX(name,['fsb_layer'],datadir)
    pv2, ft, ft2, ix = cx.load_postprocessing()
    for i in range(3):
        tstr = str(i) + '_fsb_layer'
        y = pv2[tstr].to_numpy()
        fc = fci_regmodel(y,ft2,pv2)
        fc.rebaseline(span=500,plotfig=False)
        
        fc.run_pearson(regchoice)
        x = np.arange(len(regchoice))
        #plt.plot(x,fc.pearson_rho,color=colours[i,:])
        #plt.scatter(x,fc.pearson_rho,color=colours[i,:],alpha=1,s=10)
        rhos[:,i,ids] = fc.pearson_rho
rmean = np.mean(rhos,axis=2)
plt.plot([0, len(regchoice)-1],[0,0],color='k',linestyle='--')
for ids,datadir in enumerate(datadirs):
    for i in range(3):
        plt.scatter(x,rhos[:,i,ids],color=colours[i,:],s=7.5)
for i in range(3):
    plt.plot(x,rmean[:,i],color=colours[i,:],linewidth=2)
plt.xticks(x,labels=regchoice,rotation=45)
plt.ylim([-0.5,1])
plt.yticks(np.arange(-0.5,1.25,0.25))
plt.subplots_adjust(bottom=0.3)
plt.subplots_adjust(left=0.3)
plt.ylabel("Pearson's rho")
plt.xlabel('Regressor')
plt.box('off')
plt.savefig(os.path.join(savedir,'PearsonCorr.pdf'))
# %% Delta R2 models
regchoice = ['odour onset', 'odour offset', 'in odour', 
                        'angular velocity abs','translational vel']
dR2_w = np.zeros((3,len(regchoice),len(datadirs)),dtype = float)
coeffs = np.zeros_like(dR2_w)
r2 = np.zeros((3,len(datadirs)))
corr = np.zeros((3,3,len(datadirs)))
for ix, d in enumerate(datadirs):
    cx = CX(name,['fsb_layer'],d)
    pv2, ft, ft2, i1 = cx.load_postprocessing()
    y = pv2[['0_fsb_layer','1_fsb_layer','2_fsb_layer']]
    for i in range(3): 
        fc = fci_regmodel(pv2[str(i) + '_fsb_layer'].to_numpy().flatten(),ft2,pv2)
        fc.rebaseline(span=500)# sorts out drift of data
        y[str(i) + '_fsb_layer'] = fc.ca
        fc.run(regchoice)
        r2[i,ix] = fc.r2
        fc.run_dR2(20,fc.xft)
        dR2_w[i,:,ix] =  fc.dR2_mean
        coeffs[i,:,ix] = fc.coeff_cv[:-1]
        mn,t = fc.plot_mean_flur('odour_onset',taf=10,output=True,plotting=False)
        if i==0:
            mean_trace = np.zeros((len(mn),3))
        mean_trace[:,i] = mn
    if ix==0:
        trace_dict = {'mn_trace_' +str(ix) :  mean_trace,
                  'ts_'+str(ix): t}
    else :
        trace_dict.update({'mn_trace_' +str(ix) :  mean_trace,
                           'ts_'+str(ix): t})
    corr[:,:,ix] = y.corr()
    
    
#%%
plt.figure()
sign_dR2 = -dR2_w*np.sign(coeffs)
for ix in range(len(datadirs)):
    for i in range(3):
        plt.plot(np.linspace(0,len(regchoice)-1,len(regchoice)),sign_dR2[i,:,ix],color=colours[i,:])
plt.xticks(np.linspace(0,len(regchoice)-1,len(regchoice)),labels=regchoice,rotation=90)
plt.subplots_adjust(bottom=0.4)
plt.legend(['Lower', 'Middle','Upper'])
plt.xlabel('Regressors')
plt.ylabel('delta R2 * sign(coeffs)')
