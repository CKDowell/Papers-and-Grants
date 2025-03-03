# -*- coding: utf-8 -*-
"""
Created on Mon May 20 17:43:02 2024

@author: dowel

Script makes panels for Figure 4

Panels include:

"""

import numpy as np
import pandas as pd
import src.utilities.funcs as fc
from analysis_funs.optogenetics import opto 
import os
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 
#%%
datadirs = ["Y:\Data\Optogenetics\Gr64-f\Gr64af_Horizontal\\240118\\f2\Trial1",
    "Y:\Data\Optogenetics\Gr64-f\Gr64af_Horizontal\\240121\\f3\Trial1",
"Y:\Data\Optogenetics\Gr64-f\Gr64af_Horizontal\\240124\\f1\Trial1"]
meta_data = {'stim_type': 'plume',
              'act_inhib':'act',
    'ledOny': -float(50)/2,
              'ledOffy':float(50)/2,
              'ledOnx': 10,
              'ledOffx': 30,
              'LEDoutplume': False,
              'LEDinplume': True,
              'PlumeWidth': float(50),
              'PlumeAngle': 90,
              'RepeatInterval':250
              }
#%%
dchoice= 1
datadir = datadirs[dchoice]
lname = os.listdir(datadir)
savepath = os.path.join(datadir,lname[0])
df = fc.read_log(savepath)


inst = df['instrip'].to_numpy()
exst = np.where(inst>0)[0][0]

op = opto()
x = df['ft_posx'].to_numpy()
y = df['ft_posy'].to_numpy()

x,y = op.fictrac_repair(x,y)
x = x[exst:]
y = y[exst:]
x = x-x[0]
y= y -y[0]


led_stim = df['led1_stpt'][exst:].to_numpy()
lsd = np.diff(led_stim)
ledon = np.where(lsd<0)[0]
ledoff = np.where(lsd>0)[0]
#%%
def gradient_plot(x,y,ledon,ledoff,t_real):
    tconst = 120
    for i,lon in enumerate(ledon):
        plt.plot(x[lon:ledoff[i]],y[lon:ledoff[i]],color='r',linewidth=0.5)
        if i<(len(ledoff)-1):
            intarr = np.arange(ledoff[i],ledon[i+1])
        else:
            intarr = np.arange(ledoff[i],len(t_real)-2)
        t = t_real[intarr]
        t = t-t[0]
        exp_dec = np.exp(-t/tconst)
        xseg = x[intarr]
        yseg = y[intarr]
        for r in range(len(xseg)-1):
            plt.plot(xseg[r:r+2],yseg[r:r+2],color=[exp_dec[r],0,0],linewidth=0.5)
    

#%% plot data
savedir = 'Y:\\Applications\\HHW\\Figures\\Figure4'
t_real = op.get_time(df)
t_real = t_real[exst:]
t_real = t_real-t_real[0]

x_on = meta_data['ledOnx']
x_off = meta_data['ledOffx']
y_on = meta_data['ledOny']
y_off= meta_data['ledOffy']
y_stm = [y_on, y_off, y_off,y_on]
rep_int = meta_data['RepeatInterval']


offset = 0
plt.figure()
        
ax =  plt.gca()
ax.set_aspect('equal',adjustable='box')

plt.fill([min(x),min(x),max(x),max(x)],[-25-offset,25-offset,25-offset,-25-offset],color=[0.7,0.7,0.7])
plt.plot(x[:ledon[0]],y[:ledon[0]]-offset,color='k',linewidth=0.5)

xrange = [min(x), max(x)]
y_stm = [y_on-offset, y_off-offset, y_off-offset,y_on-offset]
if xrange[0]<0:
    xr = -np.arange(0,np.abs(xrange[0]),rep_int)
    for i in xr:
        x_stm = [i-x_on, i-x_on,i-x_off,i-x_off]
        plt.fill(x_stm,y_stm,color=[1,0.5,0.5])

tranges = np.array([[ ledon[0],np.where(t_real>805)[0][0]], 
                    [np.where(t_real>805)[0][0],np.where(t_real>1007.5)[0][0]],
                    [np.where(t_real>1007.5)[0][0],len(t_real)-1]])

for t in tranges:
    print(t)
    trange = np.arange(t[0],t[1])
    offset = offset+125
    plt.fill([min(x),min(x),max(x),max(x)],[-25-offset,25-offset,25-offset,-25-offset],color=[0.7,0.7,0.7])
    tled = led_stim[trange]
    
    plt.scatter(x[trange][tled<1],y[trange][tled<1]-offset,color='r',linewidth=0.5,zorder=2)
    plt.plot(x[trange],y[trange]-offset,color='k',linewidth=0.5)
    #gradient_plot(x[trange],y[trange]-offset,np.where(lsd[trange[:-1]]<0)[0], np.where(lsd[trange[:-1]]>0)[0],t_real[trange])
    xrange = [min(x), max(x)]
    y_stm = [y_on-offset, y_off-offset, y_off-offset,y_on-offset]
    if xrange[0]<0:
        xr = -np.arange(0,np.abs(xrange[0]),rep_int)
        for i in xr:
            x_stm = [i-x_on, i-x_on,i-x_off,i-x_off]
            plt.fill(x_stm,y_stm,color=[1,0.5,0.5])
    

#
plt.figure()
t_realp = -t_real
tconst = 120
#plt.fill([min(x),min(x),max(x),max(x)],[-25,25,25,-25],color=[0.7,0.7,0.7])
plt.plot(x[:ledon[0]],t_realp[:ledon[0]],color='k',linewidth=0.5)
xrange = [min(x), max(x)]
if xrange[0]<0:
    xr = -np.arange(0,np.abs(xrange[0]),rep_int)
    for i in xr:
        x_stm = [i-x_on, i-x_on,i-x_off,i-x_off]
        #plt.fill(x_stm,y_stm,color=[1,0.5,0.5])




for i,lon in enumerate(ledon):
    plt.scatter(x[lon:ledoff[i]],t_realp[lon:ledoff[i]],color='r',linewidth=0.5)
    if i<(len(ledoff)-1):
        intarr = np.arange(ledoff[i],ledon[i+1])
    else:
        intarr = np.arange(ledoff[i],len(t_real)-2)
    t = t_real[intarr]
    t = t-t[0]
    exp_dec = np.exp(-t/tconst)
    xseg = x[intarr]
    yseg = t_realp[intarr]
    plt.plot(xseg,yseg,color='k',linewidth=0.5)
    # for r in range(len(xseg)-1):
    #     plt.plot(xseg[r:r+2],yseg[r:r+2],color=[exp_dec[r],0,0],linewidth=0.5)




xrange = [min(x), max(x)]
if xrange[0]<0:
    xr = -np.arange(0,np.abs(xrange[0]),rep_int)
    for i in xr:
        x_stm = [i-x_on, i-x_on,i-x_off,i-x_off]
        plt.plot([i-x_on,i-x_on],[min(t_realp),0],linestyle='--',color='r')
        plt.plot([i-x_off,i-x_off],[min(t_realp),0],linestyle='--',color='r')


tranges = np.array([[0, ledon[0]],
    [ ledon[0],np.where(t_real>805)[0][0]], 
                    [np.where(t_real>805)[0][0],np.where(t_real>1007.5)[0][0]],
                    [np.where(t_real>1007.5)[0][0],len(t_real)-1]])
offset = 10
for t in tranges:
    plt.plot([min(x)-offset,min(x)-offset],t_realp[t],color='k')
    offset = offset+10
plt.show()
plt.xlim([min(x)-offset,max(x)])
plt.savefig(os.path.join(savedir,'X_thru_time.pdf'))

f = plt.figure(1)
ax = f.gca()
ax.set_xlim([min(x)-offset,max(x)])
plt.savefig(os.path.join(savedir,'Trajs.pdf'))

