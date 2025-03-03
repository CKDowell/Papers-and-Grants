# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 14:26:23 2024

@author: dowel

Plots panels for figure 0
Main panel is an edge tracking trajectory with zoomed in portion showing 
leave/returns
"""


#%%
import numpy as np
import pandas as pd
import src.utilities.funcs as fc
from analysis_funs.optogenetics import opto 
import os
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 
savedir = "Y:\\Applications\\HHW\\Figures\\Figure0"
#%%
plt.close('all')
# Inhibition of PFL3 but there is no obvs phenotype
datadir = "Y:\Data\Optogenetics\PFL3\PFL3_inhibition_SS60291\\240304\\f2\Trial1"
lname = os.listdir(datadir)
savepath = os.path.join(datadir,lname[0])
df = fc.read_log(savepath)
x = df['ft_posx'].to_numpy()
y = df['ft_posy'].to_numpy()
ins = df['instrip'].to_numpy()
insd = np.diff(ins.astype(int))
pon = np.where(insd>0)[0]+1
poff = np.where(insd<0)[0]+1

start = np.where(ins>0)[0][0]
x = x-x[start]
y = y -y[start]
ym = np.max(y)
plt.fill([-25,25,25,-25],[0,0,ym,ym],color=[0.7,0.7,0.7])

for i,p in enumerate(pon):
    plt.plot(x[p:poff[i]+1],y[p:poff[i]+1],color='r')
    if i<(len(pon)-1):
        plt.plot(x[poff[i]:pon[i+1]],y[poff[i]:pon[i+1]],color='k')
    
plt.plot(x[poff[-1]:],y[poff[-1]:],color='k')
#plt.plot(x,y,color='k')
plt.plot([0,50],[430,430],color='k')
plt.plot([0,50],[465,465],color='k')
plt.plot([50,50],[430,465],color='k')
plt.plot([0,0],[430,465],color='k')
plt.gca().set_aspect('equal')
plt.show()
savename = os.path.join(savedir,'FullEgTraj.pdf')
plt.savefig(savename)

plt.figure()
plt.fill([-25,25,25,-25],[0,0,ym,ym],color=[0.7,0.7,0.7])
for i,p in enumerate(pon):
    plt.plot(x[p:poff[i]+1],y[p:poff[i]+1],color='r')
    if i<(len(pon)-1):
        plt.plot(x[poff[i]:pon[i+1]],y[poff[i]:pon[i+1]],color='k')
    
plt.plot(x[poff[-1]:],y[poff[-1]:],color='k')
plt.xlim([0,50])
plt.ylim([430,465])
plt.gca().set_aspect('equal')
plt.show()
savename = os.path.join(savedir,'TrajSegment.pdf')
plt.savefig(savename)
#%% Video of plume
testdir = r'Y:\Data\Optogenetics\37G12_PFR_a\37G12_PFR_a_inhibition\Test'
flies = [
     "241013\\f1\\Trial1",
    "241013\\f2\\Trial1",
    "241013\\f3\\Trial1",
    "241014\\f1\\Trial1",
    "241014\\f2\\Trial1",
    "241014\\f3\\Trial1",
    #"241014\\f4\\Trial1",
    "241014\\f5\\Trial1",
    "241015\\f4\\Trial2",
    "241015\\f5\\Trial2",
    "241015\\f6\\Trial2",
    
    "241023\\f1\\Trial1",#0 thresh
    "241023\\f2\\Trial1",# crap walker
    
    "241015\\f4\\Trial1",# led off
    "241015\\f5\\Trial1",# crap walker
    
    "241107\\f1\\Trial1",# 0 thresh
    "241107\\f3\\Trial1",# led off
    "241107\\f4\\Trial1",# -1000 mm threshold
    
    "241114\\f1\\Trial1",
    "241114\\f3\\Trial1",
    "241114\\f5\\Trial1",
    
    ]


plt.close('all')
# Inhibition of PFL3 but there is no obvs phenotype
for f in flies:
    #datadir = "Y:\Data\Optogenetics\PFL3\PFL3_inhibition_SS60291\\240304\\f2\Trial1"
    searchdir = os.path.join(testdir,f)
    indir = os.listdir(searchdir)
    datadir= os.path.join(searchdir,indir[0])
    files = os.listdir(datadir)
    savepath = os.path.join(datadir,files[0])
    df = fc.read_log(savepath)
    x = df['ft_posx'].to_numpy()
    y = df['ft_posy'].to_numpy()
    ins = df['instrip'].to_numpy()
    insd = np.diff(ins.astype(int))
    pon = np.where(insd>0)[0]+1
    poff = np.where(insd<0)[0]+1
    
    start = np.where(ins>0)[0][0]
    x = x-x[start]
    y = y -y[start]
    ym = np.max(y)
    plt.figure()
    plt.title(f)
    plt.fill([-25,25,25,-25],[0,0,ym,ym],color=[0.7,0.7,0.7])
    
    for i,p in enumerate(pon):
        plt.plot(x[p:poff[i]+1],y[p:poff[i]+1],color='r')
        if i<(len(pon)-1):
            plt.plot(x[poff[i]:pon[i+1]],y[poff[i]:pon[i+1]],color='k')
        
    plt.plot(x[poff[-1]:],y[poff[-1]:],color='k')
    #plt.plot(x,y,color='k')
    plt.plot([0,50],[430,430],color='k')
    plt.plot([0,50],[465,465],color='k')
    plt.plot([50,50],[430,465],color='k')
    plt.plot([0,0],[430,465],color='k')
    plt.gca().set_aspect('equal')
    plt.show()

#%%

searchdir = r'Y:\Data\Optogenetics\37G12_PFR_a\37G12_PFR_a_inhibition\Test\241015\\f4\\Trial2'
indir = os.listdir(searchdir)
datadir= os.path.join(searchdir,indir[0])
files = os.listdir(datadir)
savepath = os.path.join(datadir,files[0])
df = fc.read_log(savepath)

savedir = r'Y:\Presentations\2025\02_BSV'
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg'
import networkx as nx

#mpl.use("TkAgg") 
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.animation import ImageMagickFileWriter, ImageMagickWriter,FFMpegWriter
# Your specific x and y values
x = pd.Series.to_numpy(df['ft_posx'])
y = pd.Series.to_numpy(df['ft_posy'])

pon = pd.Series.to_numpy(df['instrip']>0)
pw = np.where(pon)

x = df['ft_posx'].to_numpy()
y = df['ft_posy'].to_numpy()
ins = df['instrip'].to_numpy()
insd = np.diff(ins.astype(int))
pon = np.where(insd>0)[0]+1
poff = np.where(insd<0)[0]+1

start = np.where(ins>0)[0][0]
x = x-x[start]
y = y -y[start]
ym = np.max(y)
x = x[start:]
y = y[start:]
ins = ins[start:].astype('float')
yrange = [min(y), max(y)]
xrange = [min(x), max(x)]


fig, ax = plt.subplots(figsize=(20,20))
plt.fill([-25,25,25,-25],[yrange[0],yrange[0],yrange[1],yrange[1]],color=[0.7,0.7,0.7])
sc = ax.scatter([],[],color='k',s=0.5)
sc2 = ax.scatter([],[],color=[0.9,0.1,0.1],s=0.5)
# line, = ax.plot([], [], lw=2,color=[1,0.2,0.2])
# line2, = ax.plot([],[],lw=2,color=[0.2,0.2,0.2])


# Set axis limits
#ax.set_xlim(xrange[0], xrange[1])
#ax.set_ylim(yrange[0], 1000)
ax.set_aspect('equal')
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
def update(frame):
    # Update the line plot with the current x and y values
    #line2.set_data(x[:frame], y[:frame])
    # if frame>100:
    #     line.set_data(x[frame-100:frame], y[frame-100:frame])
    # else:
    #     line.set_data(x[:frame], y[:frame])
    
    tx = x[:frame]
    ty = y[:frame]
    tins = ins[:frame]
    
    sc.set_offsets(np.column_stack((tx[tins<1],ty[tins<1])))
    sc2.set_offsets(np.column_stack((tx[tins>0],ty[tins>0])))
    ax.set_xlim([-150, 150])
    if frame>0:
        ax.set_ylim([ty[-1]-50,ty[-1]+50])
    #sc.set_color(ins[:frame])
    # if z[frame]<1:
    #     sc.set_offsets(np.column_stack((x[frame],y[frame])))
    
# Create animation
animation = mpl.animation.FuncAnimation(fig, update, frames=len(x), interval=0.01)
#savedir = "Y:\Data\Optogenetics\\Gr64-f\\Gr64af_Horizontal_Outside\\SummaryFigures"
path_to_convert = r'C:\\Program Files\\ImageMagick-7.1.1-Q16-HDRI\\convert.exe'
path_to_magick = r'C:\\Program Files\\ImageMagick-7.1.1-Q16-HDRI\\magick.exe'
#writer = ImageMagickWriter(fps=100, metadata=dict(artist='Me'), bitrate=1800)
writer = FFMpegWriter(fps=300)

animation.save(os.path.join(savedir,'EgStraightPlume.avi'), writer=writer)

plt.show()    
#%%











x = x-x[pw[0][0]]
y = y-y[pw[0][0]]
x = x[pw[0][0]:]
y = y[pw[0][0]:]
z = z[pw[0][0]:]
# Create initial line plot


line, = ax.plot([], [], lw=2,color=[0.2,0.4,1])  # Empty line plot with line width specified



sc = ax.scatter([],[],color='red')
plt.box('False')
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
pa = meta_data['PlumeWidth']
yrange = [min(y), max(y)]
xrange = [min(x), max(x)]
x_plm = [xrange[0], xrange[0], xrange[1],xrange[1]]
y_plm = [-pa/2, pa/2, pa/2, -pa/2]
plt.fill(x_plm,y_plm,color =[0.8,0.8,0.8])
x_on = meta_data['ledOnx']
x_off = meta_data['ledOffx']

y_on = meta_data['ledOny']
y_off= meta_data['ledOffy']
y_stm = [y_on, y_off, y_off,y_on]

rep_int = meta_data['RepeatInterval']

a_s = meta_data['act_inhib']
if a_s=='act':
    led_colour = [1,0.8,0.8]
elif a_s=='inhib':
    led_colour = [0.8, 1, 0.8]
if xrange[0]<0:
    
    xr = -np.arange(0,np.abs(xrange[0]),rep_int)
    print(xr)
    for i in xr:
        
        x_stm = [i-x_on, i-x_on,i-x_off,i-x_off]
        ax.fill(x_stm,y_stm,color=led_colour)
if xrange[1]>0:
    
    xr = np.arange(0,np.abs(xrange[1]),rep_int)
    print(xr)
    for i in xr:
        
        x_stm = [i+x_on, i+x_on,i+x_off,i+x_off]
        ax.fill(x_stm,y_stm,color=led_colour)

ax.set_aspect('equal')
# Set axis limits
ax.set_xlim(xrange[0], xrange[1])
ax.set_ylim(yrange[0], yrange[1])  # Adjust the y-axis limits based on your data

# Animation update function
def update(frame):
    # Update the line plot with the current x and y values
    line2.set_data(x[:frame], y[:frame])
    if frame>100:
        line.set_data(x[frame-100:frame], y[frame-100:frame])
    else:
        line.set_data(x[:frame], y[:frame])
    
    if z[frame]<1:
        sc.set_offsets(np.column_stack((x[frame],y[frame])))
    
# Create animation
animation = mpl.animation.FuncAnimation(fig, update, frames=len(x), interval=0.01)
savedir = "Y:\Data\Optogenetics\\Gr64-f\\Gr64af_Horizontal_Outside\\SummaryFigures"
path_to_convert = r'C:\\Program Files\\ImageMagick-7.1.1-Q16-HDRI\\convert.exe'
path_to_magick = r'C:\\Program Files\\ImageMagick-7.1.1-Q16-HDRI\\magick.exe'
#writer = ImageMagickWriter(fps=100, metadata=dict(artist='Me'), bitrate=1800)
writer = FFMpegWriter(fps=300)
#writer.program = path_to_magick

animation.save(os.path.join(savedir,'Outisde_Stim_fly1_better.avi'), writer=writer)

plt.show()         
