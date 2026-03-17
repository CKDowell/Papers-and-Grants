# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 16:39:42 2025

@author: dowel
"""

import numpy as np
import os

import analysis_funs.utilities.funcs as fc
from analysis_funs.optogenetics import opto 
import matplotlib.pyplot as plt
import pandas as pd
from Utilities.utils_general import utils_general as ug
plt.rcParams['pdf.fonttype'] = 42 
plt.rcParams['font.sans-serif'] = 'Arial'
#%%
rootdirs = [
    r"Y:\Data\Behaviour\AndyPaper\0_degrees",
            
            r"Y:\Data\Behaviour\AndyPaper\45_degrees",
            
            r"Y:\Data\Behaviour\AndyPaper\90_degrees",
            r"Y:\Data\Behaviour\AndyPaper\jumps",
            ]


#%% Determine the efficiency of each return to plume Gather data
plt.close('all')
tilts = [0,45,90,0]
data_dict = {}
for ir, r in enumerate(rootdirs):
    rparts = r.split('\\')
    files = os.listdir(r)
    data_array = np.zeros((300,len(files)))
    for i,f in enumerate(files):
        print(r,'___',f)
        # try:
        df = fc.read_log(os.path.join(r,f))
        # except:
            
        #     df = pd.read_table(os.path.join(r,f), delimiter=',', engine='python', on_bad_lines='skip')
        op = opto()
        effs = op.get_ret_efficiency(df,tilt=tilts[ir])
        data_array[:len(effs),i] = effs
    dsum = np.sum(data_array==0,axis=1)
    zero = np.where(dsum==len(files))[0][0]   
    data_array = data_array[:zero,:]
    data_dict.update({rparts[-1]:data_array})
#%% Plot data normalised
savedir = r"Y:\Papers, Review, Theses\RutaLabPapers\ET_2024\MyContribution\SupplementalPanels\ExtFig10"
parts = ['0_degrees','90_degrees','45_degrees','jumps']
plt.close('all')
nmin= 5
for p in parts:
    tdata = data_dict[p]
    tsum = np.sum(tdata>0,axis=1)
    t_thresh = np.where(tsum>=5)[0][-1]
    tdata = tdata[:t_thresh,:]
    tdata[tdata==0] = np.nan
    
    tdata = tdata/tdata[0,:]
    x = np.arange(0,len(tdata))[:,np.newaxis]
    x = np.tile(x,(1,tdata.shape[1]))
    plt.figure()
    plt.plot([0,t_thresh],[1,1],color='r',linestyle='--')
    plt.scatter(x[:],tdata[:],color='k',s=5,alpha=0.2)
    pmean = np.nanmean(tdata,axis=1)
    plt.plot(x[:,0],pmean,color='k',linewidth=2)
    plt.title(p)
    plt.yscale('log')
    plt.ylabel('normalised return efficiency')
    plt.xlabel('return number')
    plt.ylim([.5*10**-1,11])
    plt.savefig(os.path.join(savedir,p+'norm_return_eff.png'))
    plt.savefig(os.path.join(savedir,p+'norm_return_eff.pdf'))
#%% Plot data 
savedir = r"Y:\Papers, Review, Theses\RutaLabPapers\ET_2024\MyContribution\SupplementalPanels\ExtFig10"
parts = ['0_degrees','90_degrees','45_degrees','jumps']
plt.close('all')
nmin= 5
for p in parts:
    tdata = data_dict[p]
    tsum = np.sum(tdata>0,axis=1)
    t_thresh = np.where(tsum>=5)[0][-1]
    tdata = tdata[:t_thresh,:]
    tdata[tdata==0] = np.nan
    
    #tdata = tdata/tdata[0,:]
    x = np.arange(0,len(tdata))[:,np.newaxis]
    x = np.tile(x,(1,tdata.shape[1]))
    plt.figure()
    #plt.plot([0,t_thresh],[1,1],color='r',linestyle='--')
    plt.scatter(x[:],tdata[:],color='k',s=5,alpha=0.2)
    pmean = np.nanmean(tdata,axis=1)
    plt.plot(x[:,0],pmean,color='k',linewidth=2)
    plt.title(p)
    #plt.yscale('log')
    plt.ylabel('return efficiency')
    plt.xlabel('return number')
    plt.ylim([0,1])
    plt.savefig(os.path.join(savedir,p+'return_eff.png'))
    plt.savefig(os.path.join(savedir,p+'return_eff.pdf'))
#%% Get return path length
data_dict = {}
for ir, r in enumerate(rootdirs):
    rparts = r.split('\\')
    files = os.listdir(r)
    data_array = np.zeros((300,len(files)))
    for i,f in enumerate(files):
        print(r,'___',f)
        # try:
        df = fc.read_log(os.path.join(r,f))
        # except:
            
        #     df = pd.read_table(os.path.join(r,f), delimiter=',', engine='python', on_bad_lines='skip')
        op = opto()
        path_len = op.get_ret_path_len(df,min_duration=1)
        data_array[:len(path_len),i] = path_len
    dsum = np.sum(data_array==0,axis=1)
    zero = np.where(dsum==len(files))[0][0]   
    data_array = data_array[:zero,:]
    data_dict.update({rparts[-1]:data_array})
#%%
savedir = r"Y:\Papers, Review, Theses\RutaLabPapers\ET_2024\MyContribution\SupplementalPanels\ExtFig10"
parts = ['0_degrees','90_degrees','45_degrees','jumps']
plt.close('all')
nmin= 5
for p in parts:
    tdata = data_dict[p]
    tsum = np.sum(tdata>0,axis=1)
    t_thresh = np.where(tsum>=5)[0][-1]
    tdata = tdata[:t_thresh,:]
    tdata[tdata==0] = np.nan
    
    #tdata = tdata/tdata[0,:]
    x = np.arange(0,len(tdata))[:,np.newaxis]
    x = np.tile(x,(1,tdata.shape[1]))
    plt.figure()
    #plt.plot([0,t_thresh],[1,1],color='r',linestyle='--')
    plt.scatter(x[:],tdata[:],color='k',s=5,alpha=0.2)
    pmean = np.nanmean(tdata,axis=1)
    plt.plot(x[:,0],pmean,color='k',linewidth=2)
    plt.title(p)
    #plt.yscale('log')
    #plt.ylabel('return efficiency')
    plt.xlabel('return number')
    #plt.ylim([0,1])
    plt.yscale('log')
    plt.xlim([-.5,4.5])
#%% Plot overlaid jump trajectories
plt.close('all')
colours = np.array([[43,134,167],[219,0,34]])/255
r = r"Y:\Data\Behaviour\AndyPaper\jumps"
files = os.listdir(r)


normsamp = 500
timepoints = [0,1,4]
inarray = np.zeros((normsamp,2,len(files),len(timepoints)))
outarray = np.zeros((normsamp,2,len(files),len(timepoints)))
newtime = np.arange(0,normsamp)
for i,f in enumerate(files):
    print(r,'___',f)
    # try:
    df = fc.read_log(os.path.join(r,f))
    op = opto()
    e_e = op.get_entries_exits(df,ent_duration=1)
    x = df['ft_posx'].to_numpy()
    y = df['ft_posy'].to_numpy()
    x,y = op.fictrac_repair(x,y)
    dx = np.arange(e_e[0,0],e_e[0,2])
    tx = x[dx]
    ty = y[dx]
    tx = tx-x[e_e[0,1]]
    ty = ty-y[e_e[0,1]]
    side = -np.sign(tx[-1])
    tx= tx*side
    xoffset= 0
    for it,t in enumerate(timepoints):
        
        dxi = np.arange(e_e[t,0],e_e[t,1])
        tx = x[dxi]
        ty = y[dxi]
        tx = tx-x[e_e[t,1]]
        
        ty = ty-y[e_e[t,1]]
        oldtime = np.linspace(0,normsamp-1,len(tx))
        tx= tx*side+xoffset
        inarray[:,0,i,it] = np.interp(newtime,oldtime,tx)
        inarray[:,1,i,it] = np.interp(newtime,oldtime,ty)
        
        
        plt.plot(tx,ty,color=colours[1,:],alpha=0.1,zorder=1)
        
        dxo = np.arange(e_e[t,1],e_e[t,2])
        
        tx = x[dxo]
        ty = y[dxo]
        tx = tx-x[e_e[t,1]]
        ty = ty-y[e_e[t,1]]
        tx= tx*side+xoffset
        oldtime = np.linspace(0,normsamp-1,len(tx))
        outarray[:,0,i,it] = np.interp(newtime,oldtime,tx)
        outarray[:,1,i,it] = np.interp(newtime,oldtime,ty)
        
        plt.plot(tx,ty,color=colours[0,:],alpha=0.1,zorder=1)
    
        xoffset +=250
    

    

in_mean = np.mean(inarray,axis=2)
out_mean= np.mean(outarray,axis=2)
xoffset = 0
for i in range(len(timepoints)):
    
    
    

    plt.fill(np.array([-50,0,0,-50])+xoffset,[-1000,-1000,0,0],color=[0.7,0.7,0.7],alpha=0.25,zorder=-2)
    plt.fill(np.array([-70,-20,-20,-70])+xoffset,[1000,1000,-1000,-1000],color=[0.7,0.7,0.7],alpha=0.25,zorder=-1)
    
    plt.plot(in_mean[:,0,i],in_mean[:,1,i],color=colours[1,:])
    plt.plot(out_mean[:,0,i],out_mean[:,1,i],color=colours[0,:])
    xoffset+=250


g = plt.gca()
g.set_aspect('equal')
plt.ylim([-400,500])
plt.savefig(os.path.join(savedir,'Jump_trajs.pdf'))
#%% 45 degree
plt.close('all')
colours = np.array([[43,134,167],[219,0,34]])/255
r = r"Y:\Data\Behaviour\AndyPaper\45_degrees"
files = os.listdir(r)


normsamp = 500
timepoints = [0,1,4]
inarray = np.zeros((normsamp,2,len(files),len(timepoints)))
outarray = np.zeros((normsamp,2,len(files),len(timepoints)))
newtime = np.arange(0,normsamp)
for i,f in enumerate(files):
    print(r,'___',f)
    # try:
    df = fc.read_log(os.path.join(r,f))
    op = opto()
    e_e = op.get_entries_exits(df,ent_duration=1)
    x = df['ft_posx'].to_numpy()
    y = df['ft_posy'].to_numpy()
    x,y = op.fictrac_repair(x,y)
    dx = np.arange(e_e[-1,0],e_e[-1,2])
    tx = x[dx]
    ty = y[dx]
    tx = tx-x[e_e[0,1]]
    ty = ty-y[e_e[0,1]]
    side = -(np.sign(tx[-1])*np.sign(ty[-1]))
    tx= tx*side
    xoffset= 0
    #plt.figure()
    #plt.title(f)
    for it,t in enumerate(timepoints):
        
        dxi = np.arange(e_e[t,0],e_e[t,1])
        tx = x[dxi]
        ty = y[dxi]
        tx = tx-x[e_e[t,1]]
        
        ty = ty-y[e_e[t,1]]
        oldtime = np.linspace(0,normsamp-1,len(tx))
        tx= tx*side+xoffset
        inarray[:,0,i,it] = np.interp(newtime,oldtime,tx)
        inarray[:,1,i,it] = np.interp(newtime,oldtime,ty)
        
        
        plt.plot(tx,ty,color=colours[1,:],alpha=0.1,zorder=1)
        
        dxo = np.arange(e_e[t,1],e_e[t,2])
        
        tx = x[dxo]
        ty = y[dxo]
        tx = tx-x[e_e[t,1]]
        ty = ty-y[e_e[t,1]]
        tx= tx*side+xoffset
        oldtime = np.linspace(0,normsamp-1,len(tx))
        outarray[:,0,i,it] = np.interp(newtime,oldtime,tx)
        outarray[:,1,i,it] = np.interp(newtime,oldtime,ty)
        
        plt.plot(tx,ty,color=colours[0,:],alpha=0.1,zorder=1)
    
        xoffset +=250
    g = plt.gca()
    g.set_aspect('equal')

   
   
in_mean = np.mean(inarray,axis=2)
out_mean= np.mean(outarray,axis=2)
xoffset = 0
mult = np.sqrt(2)/2
for i in range(len(timepoints)):
    
    
    plt.fill(np.array([-200*mult,200*mult,(200-50)*mult,(-200-50)*mult])+xoffset,[200*mult,-200*mult,-(200+50)*mult,(200-50)*mult],color=[0.7,0.7,0.7],alpha=0.25,zorder=-2)

    #plt.fill(np.array([-50*mult,0,0,-50*mult])+xoffset,[-1000,-1000,0,0],color=[0.7,0.7,0.7],alpha=0.25,zorder=-2)
    #plt.fill(np.array([-70,-20,-20,-70])+xoffset,[1000,1000,-1000,-1000],color=[0.7,0.7,0.7],alpha=0.25,zorder=-1)
    
    plt.plot(in_mean[:,0,i],in_mean[:,1,i],color=colours[1,:])
    plt.plot(out_mean[:,0,i],out_mean[:,1,i],color=colours[0,:])
    xoffset+=250


g = plt.gca()
g.set_aspect('equal')
plt.ylim([-400,500])
plt.savefig(os.path.join(savedir,'degree45_trajs.pdf'))
#%% Plot overlaid horizontal
plt.close('all')
colours = np.array([[43,134,167],[219,0,34]])/255
r = r"Y:\Data\Behaviour\AndyPaper\90_degrees"
files = os.listdir(r)
xoffset= 250


plt.fill([-200,-200,500,500],[-50,0,0,-50],color=[0.7,0.7,0.7],alpha=0.25,zorder=1)

normsamp = 500
timepoints = [0,1,4]
inarray = np.zeros((normsamp,2,len(files),len(timepoints)))
outarray = np.zeros((normsamp,2,len(files),len(timepoints)))
newtime = np.arange(0,normsamp)

for i,f in enumerate(files):
    print(r,'___',f)
    # try:
    df = fc.read_log(os.path.join(r,f))
    op = opto()
    e_e = op.get_entries_exits(df,ent_duration=1)
    x = df['ft_posx'].to_numpy()
    y = df['ft_posy'].to_numpy()
    x,y = op.fictrac_repair(x,y)
    dx = np.arange(e_e[0,0],e_e[0,2])
    tx = x[dx]
    ty = y[dx]
    tx = tx-x[e_e[0,1]]
    ty = ty-y[e_e[0,1]]
    side = -np.sign(tx[-1])
    
    xoffset= 0
    #plt.figure()
    #plt.title(f)
    
    for it,t in enumerate(timepoints):
        print(e_e.shape)
        if e_e.shape[0]==3 and it==2:
            outarray[:,:,i,2] = np.nan
            inarray[:,:,i,2] = np.nan
            continue
        dxo = np.arange(e_e[t,1],e_e[t,2])
        dxi = np.arange(e_e[t,0],e_e[t,1])
        
        side = -np.sign(x[dxo][-1]-x[dxo][0])
        
        tx = x[dxi]
        ty = y[dxi]
        
        tx = tx-x[e_e[t,1]]
        
        ty = ty-y[e_e[t,1]]
        oldtime = np.linspace(0,normsamp-1,len(tx))
        tx= tx*side+xoffset
        inarray[:,0,i,it] = np.interp(newtime,oldtime,tx)
        inarray[:,1,i,it] = np.interp(newtime,oldtime,ty)
        
        
        plt.plot(tx,ty,color=colours[1,:],alpha=0.1,zorder=1)
        
        
        
        tx = x[dxo]
        
        ty = y[dxo]
        tx = tx-x[e_e[t,1]]
        ty = ty-y[e_e[t,1]]
        tx= tx*side+xoffset
        oldtime = np.linspace(0,normsamp-1,len(tx))
        outarray[:,0,i,it] = np.interp(newtime,oldtime,tx)
        outarray[:,1,i,it] = np.interp(newtime,oldtime,ty)
        
        plt.plot(tx,ty,color=colours[0,:],alpha=0.1,zorder=1)
        xoffset +=250
    
    
    
    # dxi = np.arange(e_e[0,0],e_e[0,1])
    # tx = x[dxi]
    # ty = y[dxi]
    # tx = tx-x[e_e[0,1]]
    # ty = ty-y[e_e[0,1]]
    # oldtime = np.linspace(0,normsamp-1,len(tx))
    # tx= tx*side
    # inarray[:,0,i,0] = np.interp(newtime,oldtime,tx)
    # inarray[:,1,i,0] = np.interp(newtime,oldtime,ty)
    
    
    # plt.plot(tx,ty,color=colours[1,:],alpha=0.1,zorder=1)
    
    # dxo = np.arange(e_e[0,1],e_e[0,2])
    
    # tx = x[dxo]
    # ty = y[dxo]
    # tx = tx-x[e_e[0,1]]
    # ty = ty-y[e_e[0,1]]
    # tx= tx*side
    # oldtime = np.linspace(0,normsamp-1,len(tx))
    # outarray[:,0,i,0] = np.interp(newtime,oldtime,tx)
    # outarray[:,1,i,0] = np.interp(newtime,oldtime,ty)
    
    # plt.plot(tx,ty,color=colours[0,:],alpha=0.1,zorder=1)
    
    
    
    
    # dx = np.arange(e_e[0,0],e_e[1,2])
    # tx = x[dx]
    # ty = y[dx]
    # tx = tx-x[e_e[1,1]]
    # ty = ty-y[e_e[1,1]]
    # side = -np.sign(tx[-1])
    # tx= tx*side+xoffset
    
    
    # dxi = np.arange(e_e[1,0],e_e[1,1])
    # tx = x[dxi]
    # ty = y[dxi]
    # tx = tx-x[e_e[1,1]]
    # ty = ty-y[e_e[1,1]]
    # tx= tx*side+xoffset
    
    # oldtime = np.linspace(0,normsamp-1,len(tx))
    # inarray[:,0,i,1] = np.interp(newtime,oldtime,tx)
    # inarray[:,1,i,1] = np.interp(newtime,oldtime,ty)
    
    # plt.plot(tx,ty,color=colours[1,:],alpha=0.1,zorder=1)
    
    # dxo = np.arange(e_e[1,1],e_e[1,2])
    
    # tx = x[dxo]
    # ty = y[dxo]
    # tx = tx-x[e_e[1,1]]
    # ty = ty-y[e_e[1,1]]
    # tx= tx*side+xoffset
    # oldtime = np.linspace(0,normsamp-1,len(tx))
    # outarray[:,0,i,1] = np.interp(newtime,oldtime,tx)
    # outarray[:,1,i,1] = np.interp(newtime,oldtime,ty)
    # plt.plot(tx,ty,color=colours[0,:],alpha=0.1,zorder=1)
       
in_mean = np.nanmean(inarray,axis=2)
out_mean= np.nanmean(outarray,axis=2)
for i in range(3):
    plt.plot(in_mean[:,0,i],in_mean[:,1,i],color=colours[1,:])
    plt.plot(out_mean[:,0,i],out_mean[:,1,i],color=colours[0,:])


g = plt.gca()
g.set_aspect('equal')
plt.ylim([-400,500])
plt.savefig(os.path.join(savedir,'Horizontal_trajs.pdf'))