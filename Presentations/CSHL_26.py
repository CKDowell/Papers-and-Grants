# -*- coding: utf-8 -*-

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
#%% Vector sum example 1
savedir = r'Y:\Presentations\2026\03_CSHL\Panels'

colours = np.array([
    [81,156,205],
    [0,170,105],
    [80,80,191]])/255
goal_angle =np.pi* -(110/180)
vector_1 = np.pi*-5/180
vector_2 = np.pi*-145/180
b = (np.sin(goal_angle)*np.cos(vector_1)-np.cos(goal_angle)*np.sin(vector_1))/(np.sin(vector_2)*np.cos(vector_1)-np.cos(vector_2)*np.sin(vector_1))
a = (np.cos(goal_angle)-b*np.cos(vector_2))/np.cos(vector_1)

x = np.linspace(-np.pi,np.pi,500)
g = np.cos(x-goal_angle)
plt.plot(x,g,color=colours[0,:],zorder=5,linewidth=4)

v1 = np.cos(x-vector_1)*a
plt.plot(x,v1,linestyle='--',color=colours[1,:],linewidth=2)

v2 = np.cos(x-vector_2)*b
plt.plot(x,v2,linestyle='--',color=colours[2,:],linewidth=2)
plt.xlim([-np.pi,np.pi])
plt.savefig(os.path.join(savedir,'VectorSumSin.pdf'))
#plt.plot(x,v2+v1)

#%% Vector sum example 2
savedir = r'Y:\Presentations\2026\03_CSHL\Panels'

colours = np.array([
    [81,156,205],
    [0,170,105],
    [80,80,191]])/255
goal_angle =np.pi* -(30/180)
vector_1 = np.pi*-5/180
vector_2 = np.pi*-145/180
b = (np.sin(goal_angle)*np.cos(vector_1)-np.cos(goal_angle)*np.sin(vector_1))/(np.sin(vector_2)*np.cos(vector_1)-np.cos(vector_2)*np.sin(vector_1))
a = (np.cos(goal_angle)-b*np.cos(vector_2))/np.cos(vector_1)

x = np.linspace(-np.pi,np.pi,500)
g = np.cos(x-goal_angle)
plt.plot(x,g,color=colours[0,:],zorder=5,linewidth=4)

v1 = np.cos(x-vector_1)*a
plt.plot(x,v1,linestyle='--',color=colours[1,:],linewidth=2)

v2 = np.cos(x-vector_2)*b
plt.plot(x,v2,linestyle='--',color=colours[2,:],linewidth=2)
plt.xlim([-np.pi,np.pi])
plt.savefig(os.path.join(savedir,'VectorSumSin2.pdf'))
#plt.plot(x,v2+v1)
#%% #%% Vector sum example 3
savedir = r'Y:\Presentations\2026\03_CSHL\Panels'

colours = np.array([
    [81,156,205],
    [0,170,105],
    [80,80,191]])/255
goal_angle =np.pi* -(110/180)
vector_1 = np.pi*-80/180
vector_2 = np.pi*-145/180
#b = (np.sin(goal_angle)*np.cos(vector_1)-np.cos(goal_angle)*np.sin(vector_1))/(np.sin(vector_2)*np.cos(vector_1)-np.cos(vector_2)*np.sin(vector_1))
#a = (np.cos(goal_angle)-b*np.cos(vector_2))/np.cos(vector_1)

x = np.linspace(-np.pi,np.pi,500)
g = np.cos(x-goal_angle)
#plt.plot(x,g,color=colours[0,:],zorder=5,linewidth=4)

v1 = np.cos(x-vector_1)*a
plt.plot(x,v1,linestyle='--',color=colours[1,:],linewidth=2)

v2 = np.cos(x-vector_2)*b
plt.plot(x,v2,linestyle='--',color=colours[2,:],linewidth=2)
plt.xlim([-np.pi,np.pi])

plt.plot(x,v1+v2,color=colours[0,:],zorder=5,linewidth=4)

plt.savefig(os.path.join(savedir,'VectorSumSin3.pdf'))
#%% Vector modulation movie
savedir = r'Y:\Presentations\2026\03_CSHL\Panels'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

colours = np.array([
    [81,156,205],
    [0,170,105],
    [80,80,191]])/255

vector_1 = np.pi * -5/180
vector_2 = np.pi * -145/180

x = np.linspace(-np.pi, np.pi, 500)

fig, ax = plt.subplots()

g_line, = ax.plot([], [], color=colours[0,:], zorder=5, linewidth=4)
v1_line, = ax.plot([], [], linestyle='--', color=colours[1,:], linewidth=2)
v2_line, = ax.plot([], [], linestyle='--', color=colours[2,:], linewidth=2)

ax.set_xlim(-np.pi, np.pi)
ax.set_ylim(-2, 2)
ax.axis('off')

goal_angles = np.linspace(vector_2, vector_1, 200)
goal_angles = np.append(goal_angles,np.flipud(goal_angles))
def update(frame):
    goal_angle = goal_angles[frame]

    b = (np.sin(goal_angle)*np.cos(vector_1) - np.cos(goal_angle)*np.sin(vector_1)) / \
        (np.sin(vector_2)*np.cos(vector_1) - np.cos(vector_2)*np.sin(vector_1))

    a = (np.cos(goal_angle) - b*np.cos(vector_2)) / np.cos(vector_1)

    g = np.cos(x - goal_angle)
    v1 = np.cos(x - vector_1) * a
    v2 = np.cos(x - vector_2) * b

    g_line.set_data(x, g)
    v1_line.set_data(x, v1)
    v2_line.set_data(x, v2)

    return g_line, v1_line, v2_line

ani = FuncAnimation(fig, update, frames=len(goal_angles), blit=True)
savepath = os.path.join(savedir,"goal_angle_movie.mp4")

ani.save(savepath, fps=30)


#%% Load fly data

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
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260425\f2\Trial3', # 3 jumps
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260401\f1\Trial2', # Some nice tracking 3 jumps
        r'Y:\Data\FCI\Hedwig\68A10_60D05_FC2_GC8s_RC3\260429\f2\Trial1',# 17 jumps
        ],
    
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
        if group_name=='hDeltaC_FC2':
            etp = ET_paper(datadir,regions=['fsb1_ch1','fsb2_ch2'])
        elif group_name=='hDeltaC_FC2_EPG':
            etp = ET_paper(datadir,regions=['eb_ch1','fsb1_ch1','fsb2_ch2'])
        else:
            etp = ET_paper(datadir)
        all_flies[group_name][i] = etp
        
        
        
        
#%% Plot example data from FC2
savedir = r'Y:\Presentations\2026\03_CSHL\ExampleDataPanels'
eg_fly = 3
cell_type = 'FC2'
etp = all_flies[cell_type][eg_fly] # 2

# Jump trajectory
plt.close('all')
#tstart =660
tstart = 60*4.2 + 10
tend = 60*4.2 + 60
etp.example_trajectory(tstart,tend)
#plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Jump heats
etp.example_heatmaps(tstart,tend,cmaps=['Greys','Blues'])
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly_'+cell_type +str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'Phase_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Plume_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))



colours = np.array([[0,0,0],[81,156,205]])/255
etp.example_trajectory_arrows(tstart,tend,colours,samp_period=2,boxes=[],regions=['eb','fsb_upper'])
plt.savefig(os.path.join(savedir,'Eg_traj_arrows_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

#%% Plot example data from FC2 - amenotaxis
savedir = r'Y:\Presentations\2026\03_CSHL\ExampleDataPanels'
eg_fly = 3
cell_type = 'FC2'
etp = all_flies[cell_type][eg_fly] # 2

# Jump trajectory
plt.close('all')
#tstart =660
tstart = 60
tend = 110
etp.example_trajectory(tstart,tend)
#plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Jump heats
etp.example_heatmaps(tstart,tend,cmaps=['Greys','Blues'])
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly_'+cell_type +str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'Phase_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Plume_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

etp.example_trajectory_arrows(tstart,tend,colours,samp_period=2,boxes=[],regions=['eb','fsb_upper'])
plt.savefig(os.path.join(savedir,'Eg_traj_arrows_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
#%% Plot example data from hDeltaJ
savedir = r'Y:\Presentations\2026\03_CSHL\ExampleDataPanels'
eg_fly = 6
cell_type = 'hDeltaJ'
etp = all_flies[cell_type][eg_fly] # 2
t_duration = 78
t_duration = 50
# Jump trajectory
plt.close('all')
#tstart =660
tstart = 895
tend = tstart+t_duration
etp.example_trajectory(tstart,tend)
#plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Jump heats
etp.example_heatmaps(tstart,tend,cmaps=['Greys','Purples'])
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly_'+cell_type +str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'Phase_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Plume_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

colours = colours = np.array([[0,0,0],[81,61,204]])/255

etp.example_trajectory_arrows(tstart,tend,colours,samp_period=2,boxes=[],regions=['eb','fsb_upper'])
plt.savefig(os.path.join(savedir,'Eg_traj_arrows_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
#%% 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def sigmoid_cmap(rgb, N=256, k=10,cent=0.5):
    """
    rgb : target colour in 0–255
    N   : number of colour steps
    k   : steepness of sigmoid
    """
    
    rgb = np.array(rgb) / 255
    white = np.array([1,1,1])

    x = np.linspace(0,1,N)

    # sigmoid centered at 0.5
    s = 1/(1 + np.exp(-k*(x-cent)))

    # normalize to 0–1
    s = (s - s.min())/(s.max()-s.min())

    colors = white + (rgb - white) * s[:,None]

    return LinearSegmentedColormap.from_list("sigmoid_cmap", colors)

# create colormap
cmap = sigmoid_cmap([49,99,125], k=5)
#%%
def exp_cmap(rgb, N=256, T=0.25):
    """
    rgb : target color [R,G,B] in 0–255
    N   : number of color levels
    T   : time constant controlling curvature
    """
    
    rgb = np.array(rgb)/255
    white = np.array([1,1,1])

    x = np.linspace(0,1,N)
    
    # exponential transform
    y = 1 - np.exp(-x/T)
    
    # normalize to 0–1
    y = y / y.max()

    colors = white + (rgb - white) * y[:,None]

    return LinearSegmentedColormap.from_list("exp_cmap", colors)
#%% Plot example data from hDeltaC
eg_fly = 4 
cell_type = 'hDeltaC'
etp = all_flies[cell_type][eg_fly] # 1
t_duration = 78
t_duration = 50
# Jump trajectory
plt.close('all')
#tstart =660
tstart = 858
tend = tstart+t_duration
etp.example_trajectory(tstart,tend)
#plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

# Jump heats

from matplotlib.colors import LinearSegmentedColormap

# target colour
rgb = np.array([49, 99, 125]) / 255

# create colormap
white_to_blue = LinearSegmentedColormap.from_list(
    "white_to_blue",
    [(1, 1, 1), tuple(rgb)]
)


target = np.array([49,99,125])/255

colors = np.linspace([1,1,1], target, 256)
cmap = LinearSegmentedColormap.from_list("custom", colors)
cmap = sigmoid_cmap([49,99,125], k=3,cent=.75)
#cmap = exp_cmap([49,99,125],T=.25)
#cmap = white_to_rgb_cmap([49,99,125])

etp.example_heatmaps(tstart,tend,cmaps=['Greys',cmap])
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly_'+cell_type +str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'Phase_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Plume_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))

colours = colours = np.array([[0,0,0],[49,99,125]])/255
etp.example_trajectory_arrows(tstart,tend,colours,samp_period=2,boxes=[],regions=['eb','fsb_upper'])
plt.savefig(os.path.join(savedir,'Eg_traj_arrows_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
#%%
cmap = sigmoid_cmap([49,99,125], k=3,cent=.75)

data = np.random.rand(20,20)
plt.figure()
plt.imshow(data, cmap=cmap)
plt.colorbar()
plt.show()
plt.savefig(os.path.join(savedir,'Colorbar_DarkBlue.pdf'))

plt.figure()
plt.imshow(data, cmap='Purples')
plt.colorbar()
plt.show()
plt.savefig(os.path.join(savedir,'Colorbar_Purple.pdf'))
#%%
from matplotlib import cm

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from colorspacious import cspace_convert

def white_to_rgb_cmap(rgb, N=256):
    rgb = np.array(rgb)/255

    # convert endpoints to LAB
    white_lab = cspace_convert([[1,1,1]], "sRGB1", "CIELab")[0]
    target_lab = cspace_convert([rgb], "sRGB1", "CIELab")[0]

    # interpolate in LAB
    t = np.linspace(0,1,N)[:,None]
    lab_ramp = white_lab + (target_lab-white_lab)*t

    # convert back to RGB
    rgb_ramp = np.clip(cspace_convert(lab_ramp, "CIELab", "sRGB1"),0,1)

    return LinearSegmentedColormap.from_list("custom_white_seq", rgb_ramp)

# example
cmap = white_to_rgb_cmap([49,99,125])
data = np.random.rand(20,20)
plt.imshow(data, cmap=cmap)
plt.colorbar()
plt.show()
#%% Pseudo time plots

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

#%%
savedir = r'Y:\Presentations\2026\03_CSHL\Panels'
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

#%% Dual colour data
colours = np.array([[81,61,204],[49,99,125],[81,156,205]])/255

# Plot jump return phase relationships 
dchoice = 'hDeltaC_FC2'
fchoice =0  

offset = 0
flies = [0,1,2]
L = 20
regions = ['fsb1_ch1','fsb2_ch2']
for f in flies:
    etp = all_flies[dchoice][f]
    cxa = etp.cxa
    jumps = cxa.get_jumps()
    fc2 = cxa.pdat['phase_fsb1_ch1']
    hdc = cxa.pdat['phase_fsb2_ch2']
    for ij,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
        tf = fc2[dx]
        th = hdc[dx]
        fsub =stats.circmean(tf[:5],low=-np.pi,high=np.pi)
        ft = ug.circ_subtract(tf,fsub)
        th = ug.circ_subtract(th,fsub)
        x = np.arange(0,len(ft))
        plt.fill([x[0],x[-1],x[-1],x[0]],np.array([-np.pi,-np.pi,np.pi,np.pi])+offset,color=[.8,.8,.8])
        plt.scatter(x,ft+offset,color=colours[2,:],s=5)
        # plt.scatter(x,th,color=colours[1,:],s=5)
        #plt.plot(x,ft,color=colours[2,:])
        plt.scatter(x,th+offset,color=colours[1,:],s=5)
        
        phase,traj,amp = etp.trajectory_mean(regions=regions)
        #Predict
        
#         zf = np.exp(1j*ft)
#         zh = np.exp(1j*th)
#         #zf = zf[:len(zf)-np.mod(len(zf),L)]
#         #zh = zh[:len(zh)-np.mod(len(zh),L)]
#         X = np.column_stack(
#     [zh[L-k:-k] if k > 0 else zh[L:] for k in range(L)]
# )
        

#         w = np.linalg.lstsq(X, zf[:-L], rcond=None)[0]
    
#         phihat = np.angle(X @ w)
    
#         plt.scatter(x[:len(phihat)],phihat+offset,color='r',s=5)



        offset = offset-2*np.pi-1
#%%
savedir = r'Y:\Presentations\2026\03_CSHL\Panels'
x = np.linspace(0,1,100)
x  = np.append(x,np.linspace(x[-1]+x[1],5+x[1],100))
cdx = [2,1]

allmean = np.zeros((200,2,len(flies)))
alltraj = np.zeros((200,2,len(flies)))
ampmean = np.zeros((200,2,len(flies)))
for f in flies:
    etp = all_flies[dchoice][f]
    for r in regions:
        etp.cxa.pdat['offset_'+r+'_phase'] =pd.Series( np.squeeze(etp.cxa.pdat['phase_'+r]))
        etp.cxa.pdat['wedges_offset_'+r] = etp.cxa.pdat['wedges_'+r]
    
    phase,traj,amp = etp.trajectory_mean(time_threshold=60,time_threshold2=0,regions=regions)
    psubtract = stats.circmean(phase[90:100,0,:],axis=0,low=-np.pi,high=np.pi)
    for ip in range(phase.shape[2]):
        phase[:,0,ip] = ug.circ_subtract(phase[:,0,ip],psubtract[ip])
        phase[:,1,ip] = ug.circ_subtract(phase[:,1,ip],psubtract[ip])
    pmean = stats.circmean(phase,high=np.pi,low=-np.pi,axis=2,nan_policy='omit')
    allmean[:,:,f] = pmean
    alltraj[:,:,f] = np.mean(traj,axis=2)
    ampmean[:,:,f] = np.mean(amp[:,1,:,:],axis=2)
    # for i2 in range(phase.shape[2]):
    #     psub = stats.circmean(phase[95:105,0,i2],high=np.pi,low=-np.pi)
    #     for i in range(2):
    #         tp = ug.circ_subtract(phase[:,i,i2],psub)
    #         plt.scatter(tp,x,color=colours[cdx[i],:],s=5)
cchoice =[2,1]
plt.figure()
plt.fill([0,1,1,0],[-np.pi,-np.pi,np.pi,np.pi],color=[.8,.8,.8],alpha=.5,linewidth=0)
plt.fill([1,5,5,1],[-np.pi,-np.pi,0,0],color=[.9,.1,.1],alpha=.2,linewidth=0)
plt.fill([1,5,5,1],[np.pi,np.pi,0,0],color=[.4,.8,.96],alpha=.2,linewidth=0)
plt.plot([0,max(x)],[0,0],color='k',linestyle='--')
#plt.plot([0,max(x)],[-np.pi/2,-np.pi/2],color='r',linestyle='--')

for i in range(2):
    pm = stats.circmean(allmean,low=-np.pi,high=np.pi,axis=2)
    for f in flies:
        plt.scatter(x,allmean[:,i,f],color=colours[cchoice[i],:],s=10,alpha=.25)
        
    plt.scatter(x,pm[:,i],color=colours[cchoice[i],:],zorder=5)

plt.xticks([])
plt.xlim([0,5])
plt.ylim([-np.pi,np.pi])
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],labels=[-180,-90,0,90,180],fontsize=15)
plt.ylabel('phase [deg]',fontsize=15)
plt.savefig(os.path.join(savedir,'DualColourPseudotimeProgression.pdf'))


#%%
#%% hDeltaC FC2 phase offset jumps from re-entry

plt.close('all')
regions2 = ['fsb1_ch1','fsb2_ch2']
plen = 20*10
x = np.arange(0,plen)/10
pltmean = np.zeros((plen,len(datadirs)))*np.nan
pltmeani = np.zeros((plen,2,len(datadirs)))*np.nan
flies = [0,1,2]
for i,d in enumerate(flies):
    cxa = all_flies[dchoice][d].cxa
    jumps = cxa.get_jumps(time_threshold=20)
    #jumps = cxa.get_jumps()
    #jumps = cxa.get_entries_exits_like_jumps()
    fc2 = cxa.pdat['phase_fsb1_ch1']*-cxa.side
    hdc = np.squeeze(cxa.pdat['phase_fsb2_ch2'])*-cxa.side
    offsets = ug.circ_subtract(hdc,fc2)
    plotmat = np.zeros((plen,len(jumps)))*np.nan
    plotmat_indi = np.zeros((plen,2,len(jumps)))*np.nan
    for ij,j in enumerate(jumps):
        dx = np.arange(j[1],j[2])
        offs = fc2[dx[0]]
        if len(dx)>plen:
            dx = dx[-plen:]
            
        plotmat[-len(dx):,ij] = offsets[dx]
        plotmat_indi[-len(dx):,0,ij] = ug.circ_subtract(fc2[dx],offs)
        plotmat_indi[-len(dx):,1,ij] = ug.circ_subtract(hdc[dx],offs)
    pltmean[:,i] = stats.circmean(plotmat,low=-np.pi,high=np.pi,axis=1,nan_policy='omit')
    pltmeani[:,:,i] = stats.circmean(plotmat_indi,low=-np.pi,high=np.pi,axis=2,nan_policy='omit')
    
plt.figure(101)
plt.plot(x,pltmean,color='k',alpha=0.5)
plt.plot(x,stats.circmean(pltmean,axis=1,low=-np.pi,high=np.pi,nan_policy='omit'),color='k')

pimean = stats.circmean(pltmeani,axis=2,low=-np.pi,high=np.pi,nan_policy='omit')


plt.figure(102,figsize=(3,6))
for i in range(len(flies)):
    plt.scatter(pltmeani[50:,0,i],x[50:]-20,color=colours[2,:],alpha=0.8,s=10)
    #plt.plot(x,pimean[:,0],color=colours[1,:])
    
    plt.scatter(pltmeani[50:,1,i],x[50:]-20,color=colours[1,:],alpha=0.8,s=10)


#plt.plot(x,pimean[:,1],color=colours[0,:])
plt.plot([0,0],[-15,0],color='k',linestyle='--')
plt.plot([-np.pi,np.pi],[0,0],color='k',linestyle='--')
plt.yticks([-15,-10,-5,0])
#plt.plot([-np.pi/2,-np.pi/2],[-20,0],color='r',linestyle='--')
plt.xticks(np.arange(-np.pi,np.pi+.1,np.pi/2),labels=np.arange(-180,181,90))
plt.xlabel('Phase (deg)')
plt.xlim([-np.pi,np.pi])
plt.ylabel('Time to return (s)')

plt.savefig(os.path.join(savedir, 'DualColour_Average.pdf'))


plt.figure(101)
plt.ylim([0,np.pi])
    

for i in range(len(datadirs)):
    plt.figure()
    plt.plot(x,pltmeani[:,0,i],color='b',alpha=0.5,marker='o',linestyle='None')

    plt.plot(x,pltmeani[:,1,i],color='r',alpha=0.5,marker='o',linestyle='None')



#%% Plot example data from dual colour
savedir = r'Y:\Presentations\2026\03_CSHL\ExampleDataPanels'

#eg_fly = 2
#tstart =1000
#tstart = 978
#tend = tstart+150
eg_fly = 2
cell_type = 'hDeltaC_FC2'
etp = all_flies[cell_type][eg_fly] # 2

# Jump trajectory
plt.close('all')
tstart =1000
#tstart = 978
tend = tstart+100

regions = ['fsb1_ch1','fsb2_ch2']
for r in regions:
    etp.cxa.pdat['offset_'+r+'_phase'] =pd.Series( np.squeeze(etp.cxa.pdat['phase_'+r]))
    etp.cxa.pdat['wedges_offset_'+r] = etp.cxa.pdat['wedges_'+r]
etp.example_trajectory(tstart,tend)
#plt.xlim([-50,15])
#plt.savefig(os.path.join(savedir,'Eg_traj_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
cmap = sigmoid_cmap([49,99,125], k=3,cent=.75)

# Jump heats
tcolours = np.zeros((4,3))
tcolours[2] = np.array([237,30,36])/255
tcolours[:2] = colours[[2,1],:]
etp.example_heatmaps(tstart,tend,cmaps=['Blues',cmap],regions=regions,colours=tcolours)
# plt.figure(2)
# plt.savefig(os.path.join(savedir,'EPG_heat_Fly_'+cell_type +str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
# plt.figure(3)
# plt.savefig(os.path.join(savedir,'FC2_heat_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
# plt.figure(4)
# plt.savefig(os.path.join(savedir,'Phase_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
# plt.figure(5)
# plt.savefig(os.path.join(savedir,'Plume_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))


#%% Plot example data from dual colour Zoom
savedir = r'Y:\Presentations\2026\03_CSHL\ExampleDataPanels'
eg_fly = 2
cell_type = 'hDeltaC_FC2'
etp = all_flies[cell_type][eg_fly] # 2

# Jump trajectory
plt.close('all')
tstart =1135.5

#tstart = 978
tend = 1150

regions = ['fsb1_ch1','fsb2_ch2']
for r in regions:
    etp.cxa.pdat['offset_'+r+'_phase'] =pd.Series( np.squeeze(etp.cxa.pdat['phase_'+r]))
    etp.cxa.pdat['wedges_offset_'+r] = etp.cxa.pdat['wedges_'+r]
etp.example_trajectory(tstart,tend)
#plt.xlim([-50,15])
plt.savefig(os.path.join(savedir,'Eg_traj_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
cmap = sigmoid_cmap([49,99,125], k=3,cent=.75)

# Jump heats
tcolours = np.zeros((4,3))
tcolours[2] = np.array([237,30,36])/255
tcolours[:2] = colours[[2,1],:]
etp.example_heatmaps(tstart,tend,cmaps=['Blues',cmap],regions=regions,colours=tcolours)
plt.figure(2)
plt.savefig(os.path.join(savedir,'EPG_heat_Fly_'+cell_type +str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(3)
plt.savefig(os.path.join(savedir,'FC2_heat_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(4)
plt.savefig(os.path.join(savedir,'Phase_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))
plt.figure(5)
plt.savefig(os.path.join(savedir,'Plume_Fly_'+cell_type+str(eg_fly)+'_'+str(tstart) +'_'+str(tend)+'.pdf'))


#%% 

plt.close('all')
data = all_flies['hDeltaC']
for f in data:
    cxa = data[f].cxa
    ins = cxa.ft2['instrip'].to_numpy()
    y = np.mean(cxa.pdat['wedges_fsb_upper'],axis=1)
    cxa.plot_traj_arrow_new(['eb','fsb_upper'],a_sep=10)
    plt.figure()
    plt.plot(y,color='k')
    plt.plot(ins,color='r')




#%% Dual colour plots - when is FC2 paying attention to hDC

# Metrics - 1. circular correlation
#           2. Offset in signal
#           3. Dot product

data = all_flies['hDeltaC_FC2'][2]
cxa = data.cxa
import pingouin as pg
fc2 = cxa.pdat['phase_fsb1_ch1']
hdc = cxa.pdat['phase_fsb2_ch2'].squeeze()
heading = cxa.ft2['ft_heading'].to_numpy()



C = pg.circ_corrcc(fc2,hdc)
C2 = pg.circ_corrcc(heading,fc2)

winlength = 50
step =1
n = len(fc2)
starts = range(0,n-winlength+1,step)
ins = cxa.ft2['instrip'].to_numpy()


def time_varying_circ_corr(data1,data2,winlength):
    indices = np.array([])
    correlations = np.array([])
    for s in starts:
        
        end = s+winlength
        r = pg.circ_corrcc(data1[s:end],data2[s:end])
        idx = (s + end - 1) / 2
        indices = np.append(indices,idx)
        correlations = np.append(correlations,r[0])

    return indices, correlations


t= np.arange(0,len(fc2))/10

i_hdc,c_hdc = time_varying_circ_corr(hdc,fc2,30)
plt.plot(i_hdc/10,c_hdc,color='k')
i_hd,c_hd = time_varying_circ_corr(heading,fc2,30)
plt.plot(i_hd/10,c_hd,color='b')
plt.plot(t,ins,color='r')
colours = uplt.columnar_colours()
plt.scatter(t,fc2+2*np.pi,color=colours[2,:],s=3)
plt.scatter(t,hdc+2*np.pi,color=colours[1,:],s=3)

#%%
savedir = r'Y:\Presentations\2026\04_LabMeeting_Frankfurt\DualColour_hDC_FC2'
plt.close('all')
for i in range(3):
    cxa = all_flies['hDeltaC_FC2'][i].cxa
    jumps = cxa.get_jumps()
    fc2 = cxa.pdat['phase_fsb1_ch1']
    hdc = cxa.pdat['phase_fsb2_ch2'].squeeze()
    amp = np.mean(cxa.pdat['wedges_fsb2_ch2'],axis=1)
    amp = (amp-np.mean(amp))/np.std(amp)
    heading = cxa.ft2['ft_heading'].to_numpy()
    t= np.arange(0,len(fc2))/10
    ins = cxa.ft2['instrip'].to_numpy()

    dot = np.cos(ug.circ_subtract(hdc,fc2))
    plt.figure()
    plt.hist(dot,bins=100)
    plt.figure()
    plt.scatter(t,fc2+2*np.pi,color=colours[2,:],s=3)
    plt.scatter(t,hdc+2*np.pi,color=colours[1,:],s=3)
    plt.scatter(t,heading+2*np.pi,color='k',s=3)
    plt.plot(t,dot,color='k')
    plt.plot(t,ins,color='r')
    plt.scatter(t[jumps[:,1]],np.ones(len(jumps)),color='m')
    plt.figure(figsize=(15,15))
    for ij,j in enumerate(jumps):
        
        
        dx = np.arange(j[0],j[2])
        tp = np.where(dx==j[1])[0]
        tfc2 = fc2[dx]
        thdc = hdc[dx]
        tamp = amp[dx]
        poffset = tfc2[tp]
        tfc2 = ug.circ_subtract(tfc2,poffset)
        thdc = ug.circ_subtract(thdc,poffset)
        tt = t[dx]
        tt = tt-t[j[1]]
        offset =- np.pi*2.5*ij 
        plt.fill(np.array([tt[0],0,0,tt[0]]),np.array([-np.pi,-np.pi,np.pi,np.pi])+offset,color=[1,.8,.8])
        plt.fill(np.array([tt[0],tt[-1],tt[-1],tt[0]]),np.array([-np.pi,-np.pi,np.pi,np.pi])+offset,color=[.8,.8,.8],zorder=-1)
        
        plt.scatter(tt,tfc2+offset,color=colours[2,:],s=5)
        plt.scatter(tt,thdc+offset,color=colours[1,:],s=5)
        plt.plot([tt[0],tt[-1]],np.array([0,0])+offset,color='k',linestyle='--')
        #plt.plot([tt[0],tt[-1]],np.array([-np.pi,-np.pi])+offset,color='k')
        #plt.plot([tt[0],tt[-1]],np.array([+np.pi,+np.pi])+offset,color='k')
        uplt.coloured_line_simple(tt,tt*0-1.2*np.pi+offset,tamp,'coolwarm',-2,2,linewidth=3)
    plt.ylim([-130,5])
    plt.xlim([-30,60])
    plt.savefig(os.path.join(savedir,'jump_lines_'+cxa.name+'.pdf'))
    plt.savefig(os.path.join(savedir,'jump_lines_'+cxa.name+'.png'))
    
for i in range(2):
    cxa = all_flies['hDeltaC_FC2_EPG'][i].cxa
    jumps = cxa.get_jumps()
    fc2 = cxa.pdat['phase_fsb1_ch1']
    hdc = cxa.pdat['phase_fsb2_ch2'].squeeze()
    epg = cxa.pdat['phase_eb_ch1']
    heading = cxa.ft2['ft_heading'].to_numpy()
    
    amp = np.mean(cxa.pdat['wedges_fsb2_ch2'],axis=1)
    amp = (amp-np.mean(amp))/np.std(amp)
    
    fc2_s = ug.savgol_circ(cxa.pdat['phase_fsb1_ch1'],30,3)
    hdc_s = ug.savgol_circ(cxa.pdat['phase_fsb2_ch2'].squeeze(),30,3)
    epg_s = ug.savgol_circ(cxa.pdat['phase_eb_ch1'],30,3)
    
    t= np.arange(0,len(fc2))/10
    ins = cxa.ft2['instrip'].to_numpy()

    dot = np.cos(ug.circ_subtract(hdc_s,fc2_s))
    dot2 = np.cos(ug.circ_subtract(epg_s,fc2_s))
    plt.figure()
    plt.subplot(1,2,1)
    plt.hist(dot,bins=100)
    plt.subplot(1,2,2)
    plt.hist(dot2,bins=100)
    plt.figure()
    plt.scatter(dot,dot2,s=1)
    
    plt.figure()
    plt.scatter(t,fc2+2*np.pi,color=colours[2,:],s=3)
    plt.scatter(t,hdc+2*np.pi,color=colours[1,:],s=3)
    plt.scatter(t,epg+2*np.pi,color='k',s=3)
    plt.plot(t,dot,color=colours[1,:])
    plt.plot(t,dot2,color='k')
    plt.plot(t,ins,color='r')
    plt.scatter(t[jumps[:,1]],np.ones(len(jumps)),color='m')
    
    
    
    plt.figure(figsize=(15,15))
    poffset = stats.circmean(ug.circ_subtract(epg,heading),low=-np.pi,high=np.pi)
    for ij,j in enumerate(jumps):
        
        
        dx = np.arange(j[0],j[2])
        tp = np.where(dx==j[1])[0]
        tfc2 = fc2[dx]
        thdc = hdc[dx]
        tepg = epg[dx]
        tamp = amp[dx]
        #poffset = tfc2[tp]
        tfc2 = ug.circ_subtract(tfc2,poffset)*-cxa.side
        thdc = ug.circ_subtract(thdc,poffset)*-cxa.side
        tepg = ug.circ_subtract(tepg,poffset)*-cxa.side
        tt = t[dx]
        tt = tt-t[j[1]]
        offset = -np.pi*2.5*ij 
        plt.fill(np.array([tt[0],0,0,tt[0]]),np.array([-np.pi,-np.pi,np.pi,np.pi])+offset,color=[1,.8,.8])
        plt.fill(np.array([tt[0],tt[-1],tt[-1],tt[0]]),np.array([-np.pi,-np.pi,np.pi,np.pi])+offset,color=[.8,.8,.8],zorder=-1)
        plt.scatter(tt,tfc2+offset,color=colours[2,:],s=5)
        plt.scatter(tt,thdc+offset,color=colours[1,:],s=5)
        plt.scatter(tt,tepg+offset,color='k',s=5)
        plt.plot([tt[0],tt[-1]],np.array([0,0])+offset,color='k',linestyle='--')
        #plt.plot([tt[0],tt[-1]],np.array([-np.pi,-np.pi])+offset,color='k')
        #plt.plot([tt[0],tt[-1]],np.array([+np.pi,+np.pi])+offset,color='k')
        uplt.coloured_line_simple(tt,tt*0-1.2*np.pi+offset,tamp,'coolwarm',-2,2,linewidth=3)
    plt.ylim([-130,5])
    plt.xlim([-30,60])
    plt.savefig(os.path.join(savedir,'jump_lines_'+cxa.name+'.pdf'))
    plt.savefig(os.path.join(savedir,'jump_lines_'+cxa.name+'.png'))
#%% Example FC2, EPG, hDeltaC
etp = all_flies['hDeltaC_FC2_EPG'][1]

tcolours = uplt.columnar_colours()
tcolours[0,:] = [0,0,0]
etp.plt_tmp(regions=['eb_ch1','fsb2_ch2','fsb1_ch1'],scalar=7,colours=tcolours,phase_num=5)
plt.savefig(os.path.join(savedir,'Average_jump'+etp.cxa.name+'.pdf'))
#%% Identify default mode of hDC and whether offset with FC2 aligns

plt.close('all')
plt.rcParams.update({'font.size': 15})
cxa=  etp.cxa
fc2 = cxa.pdat['phase_fsb1_ch1']
hdc = cxa.pdat['phase_fsb2_ch2'].squeeze()
epg = cxa.pdat['phase_eb_ch1']
heading = cxa.ft2['ft_heading'].to_numpy()
offset = stats.circmean(ug.circ_subtract(epg,heading),low=-np.pi,high=np.pi)
fc2 = ug.circ_subtract(fc2,offset)
epg = ug.circ_subtract(epg,offset)
hdc = ug.circ_subtract(hdc,offset)

ins = cxa.ft2['instrip'].to_numpy()
# dx = np.where(ins>0)[0][0]
# fc2= fc2[:dx]
# hdc = hdc[:dx]
# epg = epg[:dx]
plt.figure()
plt.scatter(hdc,fc2,s=2,color='k',alpha=.2)
plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='b')
plt.plot([-3*np.pi/4,-3*np.pi/4],[-np.pi,np.pi],color='r',linestyle='-')
plt.plot([3*np.pi/4,3*np.pi/4],[-np.pi,np.pi],color='r',linestyle='-')
plt.plot([-np.pi,np.pi],[-3*np.pi/4,-3*np.pi/4],color='r',linestyle='-')
plt.plot([-np.pi,np.pi],[3*np.pi/4,3*np.pi/4],color='r',linestyle='-')
plt.plot([-np.pi,np.pi],[np.pi/2,np.pi/2],color='g',linestyle='-',linewidth=2)
plt.xlabel('hdc phase (deg)')
plt.ylabel('fc2 phase (deg)')
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xlim([-np.pi,np.pi])
plt.ylim([-np.pi,np.pi])
plt.savefig(os.path.join(savedir,'hdc_fc2_scatter_'+cxa.name+'.png'))


plt.figure()
plt.scatter(epg,fc2,s=2,color='k',alpha=.2)
plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='b')
plt.plot([-np.pi,np.pi],[-3*np.pi/4,-3*np.pi/4],color='r',linestyle='-')
plt.plot([-np.pi,np.pi],[3*np.pi/4,3*np.pi/4],color='r',linestyle='-')
plt.plot([-np.pi,np.pi],[np.pi/2,np.pi/2],color='g',linestyle='-',linewidth=2)
plt.xlabel('epg phase (deg)')
plt.ylabel('fc2 phase (deg)')
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xlim([-np.pi,np.pi])
plt.ylim([-np.pi,np.pi])
plt.savefig(os.path.join(savedir,'epg_fc2_scatter_'+cxa.name+'.png'))


plt.figure()
plt.scatter(epg,hdc,s=2,color='k',alpha=.2)
plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='b')
plt.plot([-np.pi,np.pi],[-3*np.pi/4,-3*np.pi/4],color='r',linestyle='-')
plt.plot([-np.pi,np.pi],[3*np.pi/4,3*np.pi/4],color='r',linestyle='-')
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xlabel('epg phase (deg)')
plt.ylabel('hdc phase (deg)')
plt.xlim([-np.pi,np.pi])
plt.ylim([-np.pi,np.pi])
plt.savefig(os.path.join(savedir,'epg_hdc_scatter_'+cxa.name+'.png'))



plt.figure()
plt.title('hDeltaC')
plt.hist(hdc,100,density=True)
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xlim([-np.pi,np.pi])
plt.plot([0,0],[0,0.5],color='k')
plt.plot([np.pi/2,np.pi/2],[0,0.5],color='g')
plt.plot([3*np.pi/4,3*np.pi/4],[0,0.5],color='r')
plt.plot([-3*np.pi/4,-3*np.pi/4],[0,0.5],color='r')
plt.xlabel('phase (deg)')
plt.savefig(os.path.join(savedir,'hdc_hist_'+cxa.name+'.png'))


plt.figure()
plt.title('FC2')
plt.hist(fc2,100,density=True)
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xlim([-np.pi,np.pi])
plt.plot([0,0],[0,0.5],color='k')
plt.plot([np.pi/2,np.pi/2],[0,0.5],color='g')
plt.plot([3*np.pi/4,3*np.pi/4],[0,0.5],color='r')
plt.plot([-3*np.pi/4,-3*np.pi/4],[0,0.5],color='r')
plt.xlabel('phase (deg)')
plt.savefig(os.path.join(savedir,'fc2_hist_'+cxa.name+'.png'))



plt.figure()
plt.title('EPG')
plt.hist(epg,100,density=True)
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])
plt.xlim([-np.pi,np.pi])
plt.plot([0,0],[0,0.5],color='k')
plt.xlabel('phase (deg)')
plt.savefig(os.path.join(savedir,'epg_hist_'+cxa.name+'.png'))



#%% Example hDeltaC and FC2

cxa = all_flies['hDeltaC_FC2'][2].cxa
jumps = cxa.get_jumps()
fc2 = cxa.pdat['phase_fsb1_ch1']
hdc = cxa.pdat['phase_fsb2_ch2'].squeeze()
amp = np.mean(cxa.pdat['wedges_fsb2_ch2'],axis=1)
amp = (amp-np.mean(amp))/np.std(amp)
heading = cxa.ft2['ft_heading'].to_numpy()
t= np.arange(0,len(fc2))/10
ins = cxa.ft2['instrip'].to_numpy().astype(float)

dot = np.cos(ug.circ_subtract(hdc,fc2))
plt.figure(figsize=(16.88,2.69))
plt.scatter(t,fc2,color=colours[2,:],s=3)
plt.scatter(t,hdc,color=colours[1,:],s=3)

plt.fill_between(t,ins*0-4.5, ins*.5-4.5,color=[1,.2,.2],linewidth=0)
plt.scatter(t[jumps[:,1]],np.ones(len(jumps))*.5-4.5,color='k',marker='*',s=30)
plt.xlim([1100,1400])
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[-180,-90,0,90,180])

plt.savefig(os.path.join(savedir,'ExamplePhase_timecourse'+cxa.name+'.pdf'))



#%% Pre air offsets
plt.close('all')
dchoice = 'FC2'
data = all_flies[dchoice] 
for d in data:
    cxa = data[d].cxa
    cxa.get_jumps()
    tphase = cxa.pdat['phase_fsb_upper']
    heading = cxa.ft2['ft_heading'].to_numpy()
    teb = cxa.pdat['phase_eb']
    ins = cxa.ft2['instrip'].to_numpy()
    mfc1 = cxa.ft2['mfc1_stpt'].to_numpy()
    mstart = np.where(mfc1>0)[0][0]
    istart  = np.where(ins>0)[0][0]
    dx = np.arange(mstart,istart)
    offset = ug.circ_subtract(tphase,teb)
    
    x = np.arange(0,len(dx))
    h_offset = stats.circmean(ug.circ_subtract(teb,heading),low=-np.pi,high=np.pi)
    
    teb = ug.circ_subtract(teb,h_offset)
    tphase = ug.circ_subtract(tphase,h_offset)
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(teb[dx],tphase[dx],c=x,s=5)    
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    
    plt.subplot(1,2,2)
    plt.scatter(teb[dx[-1]:],tphase[dx[-1]:],color='k',s=1,alpha=.2)    
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,0],[-np.pi,np.pi],color='r')
    
    plt.title(cxa.side)
    # allo_downwind = ug.circ_subtract(teb,h_offset+np.pi)
   # cxa.plot_traj_arrow_new(['eb','fsb_upper'],a_sep=5,colours =np.array([[0,0,0],[0,0,1]]))
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(heading[dx],ug.circ_subtract(teb[dx],tphase[dx]),color='k',s=5)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')

    plt.ylabel('eb - fsb')
    plt.subplot(1,2,2)
    plt.scatter(heading[dx[-1]:],ug.circ_subtract(teb[dx[-1]:],tphase[dx[-1]:]),color='k',s=2,alpha=.3)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    plt.xlabel('heading')
    # plt.scatter(x,teb[dx],color='k',s=5)
    # plt.scatter(x,tphase[dx],color='b',s=5)
    # plt.ylim([-np.pi,np.pi])
    
    # x2 = np.arange(0,len(tphase))
    # plt.figure()
    # plt.scatter(x2[dx[-1]:],teb[dx[-1]:],color='k',s=5)
    # plt.scatter(x2[dx[-1]:],tphase[dx[-1]:],color='b',s=5)
    # plt.ylim([-np.pi,np.pi])
    # plt.scatter(x,offset[dx],color='k',s=5)
    # plt.scatter(x,allo_downwind[dx],color='b',s=5)
    
#%% Dual colour - with EPG
#plt.close('all')
dchoice = 'hDeltaC_FC2_EPG'
data = all_flies[dchoice] 
for d in data:
    cxa = data[d].cxa
    jumps = cxa.get_jumps()
    tphase = cxa.pdat['phase_fsb1_ch1']
    tphase2 = cxa.pdat['phase_fsb2_ch2']
    heading = cxa.ft2['ft_heading'].to_numpy()
    teb = cxa.pdat['phase_eb_ch1']
    ins = cxa.ft2['instrip'].to_numpy()
    mfc1 = cxa.ft2['mfc1_stpt'].to_numpy()
    mstart = np.where(mfc1>0)[0][0]
    istart  = np.where(ins>0)[0][0]
    dx = np.arange(mstart,istart)
    # dx = np.array([],dtype='int')
    # for j in jumps:
    #     dx = np.append(dx,np.arange(j[1],j[2]))
    offset = ug.circ_subtract(tphase,teb)
    
    x = np.arange(0,len(dx))
    h_offset = stats.circmean(ug.circ_subtract(teb,heading),low=-np.pi,high=np.pi)
    
    teb = ug.circ_subtract(teb,h_offset)
    tphase = ug.circ_subtract(tphase,h_offset)
    tphase2 = ug.circ_subtract(tphase2,h_offset)
    
    plt.figure()
    plt.subplot(2,2,1)
    plt.scatter(teb[dx],tphase[dx],color='k',s=5)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    plt.xlabel('epg')
    plt.ylabel('fc2')

    plt.subplot(2,2,2)
    plt.scatter(teb[dx],tphase2[dx],color='k',s=5)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    plt.xlabel('epg')
    plt.ylabel('hdc')
    
    plt.subplot(2,2,3)
    plt.scatter(tphase2[dx],tphase[dx],color='k',s=5)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    plt.xlabel('hdc')
    plt.ylabel('fc2')
    
    # plt.figure()
    # plt.subplot(2,2,1)
    # plt.scatter(teb[dx[-1]:],tphase[dx[-1]:],color='k',s=1,alpha=.2)
    # plt.xlim([-np.pi,np.pi])
    # plt.ylim([-np.pi,np.pi])
    # plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    # plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    # plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    # plt.xlabel('epg')
    # plt.ylabel('fc2')

    # plt.subplot(2,2,2)
    # plt.scatter(teb[dx[-1]:],tphase2[dx[-1]:],color='k',s=1,alpha=.2)
    # plt.xlim([-np.pi,np.pi])
    # plt.ylim([-np.pi,np.pi])
    # plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    # plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    # plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    # plt.xlabel('epg')
    # plt.ylabel('hdc')
    
    # plt.subplot(2,2,3)
    # plt.scatter(tphase2[dx[-1]:],tphase[dx[-1]:],color='k',s=1,alpha=.2)
    # plt.xlim([-np.pi,np.pi])
    # plt.ylim([-np.pi,np.pi])
    # plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    # plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    # plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    # plt.xlabel('hdc')
    # plt.ylabel('fc2')
    
    # plt.figure()
    # plt.subplot(1,2,1)
    # plt.scatter(heading,ug.circ_subtract(teb,tphase2),s=1,color='k')
    # plt.xlim([-np.pi,np.pi])
    # plt.ylim([-np.pi,np.pi])
    # plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    # plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    # plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    
    # plt.subplot(1,2,2)
    # plt.scatter(heading,ug.circ_subtract(teb,tphase),s=1,color='k')
    # plt.xlim([-np.pi,np.pi])
    # plt.ylim([-np.pi,np.pi])
    # plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    # plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    # plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
#%% Dual colour - without EPG
plt.close('all')
dchoice = 'hDeltaC_FC2'
data = all_flies[dchoice] 
for d in data:
    cxa = data[d].cxa
    cxa.get_jumps()
    tphase = cxa.pdat['phase_fsb1_ch1']
    tphase2 = cxa.pdat['phase_fsb2_ch2'].squeeze()
    heading = cxa.ft2['ft_heading'].to_numpy()
    #teb = cxa.pdat['phase_eb_ch1']
    ins = cxa.ft2['instrip'].to_numpy()
    mfc1 = cxa.ft2['mfc1_stpt'].to_numpy()
    mstart = np.where(mfc1>0)[0][0]
    istart  = np.where(ins>0)[0][0]
    dx = np.arange(mstart,istart)
    #offset = ug.circ_subtract(tphase,teb)
    
    x = np.arange(0,len(dx))
    #h_offset = stats.circmean(ug.circ_subtract(teb,heading),low=-np.pi,high=np.pi)
    
    #teb = ug.circ_subtract(teb,h_offset)
    #tphase = ug.circ_subtract(tphase,h_offset)
    #tphase2 = ug.circ_subtract(tphase2,h_offset)
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(tphase2[dx],tphase[dx],color='k',s=5)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    plt.xlabel('hdc')
    plt.ylabel('fc2')
    plt.subplot(1,2,2)
    plt.scatter(tphase2[dx[-1]:],tphase[dx[-1]:],color='k',s=1,alpha=.2)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    plt.xlabel('hdc')
    plt.ylabel('fc2')
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(heading[dx],ug.circ_subtract(tphase2[dx],tphase[dx]),color='k',s=5)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    plt.xlabel('hdc')
    plt.ylabel('fc2')
    plt.subplot(1,2,2)
    plt.scatter(heading[dx[-1]:],ug.circ_subtract(tphase2[dx[-1]:],tphase[dx[-1]:]),color='k',s=1,alpha=.2)
    plt.xlim([-np.pi,np.pi])
    plt.ylim([-np.pi,np.pi])
    plt.plot([-np.pi,np.pi],[-np.pi,np.pi],color='r')
    plt.plot([0,np.pi],[-np.pi,0],linestyle='--',color='r')
    plt.plot([-np.pi,0],[0,np.pi],linestyle='--',color='r')
    plt.xlabel('heading')
    plt.ylabel('hdc-fc2')