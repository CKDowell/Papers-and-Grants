# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:18:49 2026

@author: dowel

Script plots panels for the hDeltaC figure


"""

#%%
savedir = r'Y:\Papers, Review, Theses\EarlyDraft26\hDeltaC_Figure\panels'

#%% Panel A - connectivity
#%% Load up functions
from neuprint import Client
c = Client('neuprint.janelia.org',dataset='hemibrain:v1.2.1', token ='41275a6e3b110e7065d2d6fd22650fed7bddd1da2d2ddd395262ac3245c48a46')
c.fetch_version()
from neuprint import NeuronCriteria as NC, SynapseCriteria as SC, fetch_synapses, fetch_synapse_connections, fetch_neurons,fetch_adjacencies
from scipy import stats
from sklearn import mixture
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 
import matplotlib as mpl
mpl.rcParams['path.simplify'] = True
plt.rcParams['path.simplify'] = True
import os
#%% 
def get_den_ax_loc(cell,num_gauss=3):
    """
    Function determines axon and dendrite location from cell type
    Returns: 
        ax_term: axon terminal mean location
        den_term: mean locations of dendrite terminals
    How it works:
        Function takes the mean pre-synapse location to be the axon terminal
        Then makes a guassian mixture model on post synapse terminals
        Where the two don't overlap are where dendrites are
    Future extensions: 
        Use gaussian mixture modelling to get terminal location for
        cells with multiple axon terminal sites
    Limitations: 
        You should know in advance how the neurite tree looks to specify the
        number of gaussians
    
    """
    
    syndf = fetch_synapses(cell)
    pre_post = syndf['type']
    pre = pd.Series.to_numpy(pre_post=='pre')
    post = pd.Series.to_numpy(pre_post=='post')
    syn_locs = np.array([syndf['x'], syndf['y'], syndf['z']])
    
    cdx = pd.Series.to_numpy(syndf['confidence']>0.9) # Can relax this if needs be
    syn_locs = np.transpose(syn_locs)
    pre_locs = syn_locs[pre&cdx,:]
    post_locs = syn_locs[post&cdx,:]
    
    # get mean pre location
    mn_pre = np.mean(pre_locs,axis=0)
    sd_pre = np.std(pre_locs,axis=0)
    #
    gm = mixture.GaussianMixture(n_components=num_gauss, random_state=0).fit(post_locs)
    # Gaussian mixture model with two gaussians
    mn_post = gm.means_
    
    ax_term = mn_pre
    den_dx = np.sqrt(np.sum(np.square(mn_post-mn_pre),axis=1))
    dx = [True, True, True]
    ax_dx = np.argmin(den_dx)
    dx[ax_dx] = False
    den_term = mn_post[dx,:]
    
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.scatter(post_locs[:,0],post_locs[:,1],post_locs[:,2],color='k')
    # ax.scatter(pre_locs[:,0],pre_locs[:,1],pre_locs[:,2],color='r')
    # ax.scatter(mn_post[:,0],mn_post[:,1],mn_post[:,2],color='b')
    # plt.show()
    
    return ax_term, den_term
def get_syn_pair_loc(cell_pre,cell_post):
    # First get synapse information of post-synaptic neuron
    ax_term, den_term = get_den_ax_loc(cell_post,num_gauss=3)
    # Get synapse information about the pair
    df_syn = fetch_synapse_connections(cell_pre,cell_post)
    syn_post_big = np.array([df_syn['x_post'],df_syn['y_post'],df_syn['z_post']]).transpose()
    syn_loc= {}
    for i_n,n in enumerate(cell_pre):
        dx = pd.Series.to_numpy(df_syn['bodyId_pre']==n)
        syn_post = syn_post_big[dx,:]
        dshape = np.shape(den_term)
        syn_shape = np.shape(syn_post)
        ax_den = np.empty([syn_shape[0],dshape[0]+1],dtype=float)
        ax_diff = np.sqrt(np.sum(np.square(syn_post-ax_term),axis=1))
        ax_den[:,0] = ax_diff
        for i in range(dshape[0]):
            d_diff = np.sqrt(np.sum(np.square(syn_post-den_term[i,:]),axis=1)) 
            ax_den[:,i+1] = d_diff
            
        d_code = np.argmin(ax_den,axis=1)
        ax = d_code==0    
        syn_loc[i_n] = {'ax_den': ax_den, 'd_code': d_code, 'is_axon': ax}    
    return syn_loc,ax_term,den_term 

def ax_dendrite_hDelta(hdelta):
    """
        Parameters
        ----------
        hdelta : 
        Type of hDelta neuron
       
        Returns
        Axon vs denrite location of all pre synapses by type
        -------
        None.
       
    """   
    # Get all hdeltaC type ids
    criteria = NC(type=hdelta)
    n_df1,ndf2 = fetch_neurons(criteria)
    post_ids = n_df1['bodyId']
    # Get a summary of inputs to create the matrix
    neuron_df, conn_df, =  fetch_adjacencies(None,criteria)
    uni_type = np.unique(pd.Series.to_numpy(neuron_df['type'],'str'))
    
    ax_den_mat = np.zeros([len(uni_type),2],dtype=float)
    ax_den_mat_all = np.zeros([len(uni_type),3,len(post_ids)],dtype=float)
    ax_den_mat_adv = np.zeros([len(uni_type),3],dtype=float)
    
    # Iterate through each hdelta C neuron
    for i, n in enumerate(post_ids):
        print(str(i),  ' of ' ,str(len(post_ids)))
        dx = conn_df['bodyId_post']==n
        cell_pre = pd.Series.to_numpy(conn_df['bodyId_pre'][dx],dtype='int64')
        syn_loc,ax_term,den_term = get_syn_pair_loc(cell_pre,n)
        v = np.argmin(den_term[:,2])+1
        d = np.argmax(den_term[:,2])+1
        for i_s in range(len(cell_pre)):
            print(i_s)
            s = syn_loc[i_s]
            axnum = np.sum(s['d_code']==0)
            denum = np.sum(s['d_code']!=0)
            
            de_d = np.sum(s['d_code']==d)#dorsal or ventral dendrite
            de_v = np.sum(s['d_code']==v)
            
            cdx = neuron_df['bodyId']==cell_pre[i_s]
            tp = pd.Series.to_numpy(neuron_df['type'][cdx])
            tdx = uni_type ==tp
            ax_den_mat[tdx,0] = ax_den_mat[tdx,0]+axnum 
            ax_den_mat[tdx,1] = ax_den_mat[tdx,1]+denum 
            ax_den_mat_adv[tdx,0] = ax_den_mat_adv[tdx,0]+axnum 
            ax_den_mat_adv[tdx,1] = ax_den_mat_adv[tdx,1]+de_d 
            ax_den_mat_adv[tdx,2] = ax_den_mat_adv[tdx,2]+de_v 
            #print(s['d_code'])
            
            # if (de_d+de_v)!=denum:
            #     print('count off')
            #     print('d',d)
            #     print('v',v)
                
            #     return de_d
            
            ax_den_mat_all[tdx,0,i] = ax_den_mat_all[tdx,0,i]+axnum 
            ax_den_mat_all[tdx,1,i] = ax_den_mat_all[tdx,1,i]+de_d 
            ax_den_mat_all[tdx,2,i] = ax_den_mat_all[tdx,2,i]+de_v
            
    h_delta_ax_den = {'AxonDendrite': ax_den_mat,'AxonDendriteAdv': ax_den_mat_adv,
                      'AxonDendriteAll': ax_den_mat_all,'Post_types': uni_type,'hDelta': hdelta}
    return h_delta_ax_den

#%% Panel A main plot
h = 'hDeltaC'
tmin = 10

tans = ['FB','OA']



# Gather data
h_delta_ax_den = ax_dendrite_hDelta(h)
h_all = h_delta_ax_den['AxonDendriteAll']
h_s = np.shape(h_all)

types = h_delta_ax_den['Post_types']
tot_counts = np.sum(h_delta_ax_den['AxonDendrite'],axis=1)
I = np.argsort(-tot_counts)
tot_counts = tot_counts[I]/h_s[2]
types = types[I]


AD = h_delta_ax_den['AxonDendrite'][I]/h_s[2]
# Threshold min mean number 
types = types[tot_counts>=tmin]
AD = AD[tot_counts>=tmin]

#Separate tangentials and columnars
ordered_index = np.array([],dtype='int')
for c in tans:
    match = np.where(np.char.startswith(types,c))[0]
    ordered_index = np.append(ordered_index,match)

for r in range(2):
    plt.figure()
    if r==0:
        ttypes = types[ordered_index]
        tAD = AD[ordered_index,:]
        sav_suff = 'tan'
    else:
        ttypes = types.copy()
        ttypes = np.delete(ttypes,ordered_index)
        tAD = AD.copy()
        tAD = np.delete(tAD,ordered_index,axis=0)
        sav_suff = 'col'
    for i,n in enumerate(ttypes):
        plt.fill_between([0,tAD[i,1]],[-i-0.4,-i-0.4],[-i+0.4,-i+0.4],color='k')
        plt.fill_between([0,-tAD[i,0]],[-i-0.4,-i-0.4],[-i+0.4,-i+0.4],color=[0.5,0.5,0.5])
        # plt.plot([0,AD[i,1]],[-i,-i],color='k')
        # plt.plot([0,-AD[i,0]],[-i,-i],color=[0.5,0.5,0.5])
    

    plt.xlim([-100,100])
    plt.ylim([-28,1])
    plt.yticks(np.arange(-i,1),labels=np.flipud(ttypes))
    plt.xlabel('Mean syn count (-ax,+den)')

    savename = os.path.join(savedir,'DenAx_Count_'+ h +sav_suff+'.png')
    plt.savefig(savename)
    savename = os.path.join(savedir,'DenAx__Count'+ h +sav_suff+'.pdf')
    plt.savefig(savename)

#%% Load neural data
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

    "hDeltaC": [r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250411\f1\Trial1",# Phase recording is not the best - strong pointer
    #r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f1\Trial3",# Not many jumps, weak pointer
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250605\f2\Trial2", # Strong pointer
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250714\f1\Trial2",# Strong pointer, just like FC2
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f1\Trial1", # Points away ******
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250718\f2\Trial3", # Strong pointer
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250721\f1\Trial2",# Strong pointer
    r"Y:\Data\FCI\Hedwig\hDeltaC_SS02863\250907\f1\Trial2",# Good pointer
       
    ],

    
    
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

#%% Panel B - scatter plot of phases
plt.close('all')
dchoice = 'hDeltaC'
data = all_flies[dchoice]
# Scatter pre
epgall = np.array([])
hdcall = np.array([])
for f in data:
    cxa = data[f].cxa
    cxa.get_jumps()
    j = cxa.get_entries_exits_like_jumps()
    if dchoice =='hDeltaC':
        thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        epg = cxa.pdat['offset_eb_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
    elif dchoice =='hDeltaC_FC2_EPG':
        thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
    
    
    hdcall = np.append(hdcall,thdc)
   
    epgall = np.append(epgall,epg)
    
epgall = 180*epgall/np.pi
hdcall = 180*hdcall/np.pi
plotbins = np.linspace(-180,180,21)
plt.hist2d(epgall,hdcall,plotbins,cmap='grey_r',density=True,vmin=0,vmax=4e-5)
plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='--')
plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='--')
plt.yticks([-180,-90,0,90,180])
plt.xticks([-180,-90,0,90,180])
plt.xlabel('EPG phase (deg)')
plt.ylabel('hDeltaC phas (deg)')
plt.colorbar()
plt.title('Pre odour')
plt.savefig(os.path.join(savedir,'PreAir_hist.png'))
plt.savefig(os.path.join(savedir,'PreAir_hist.pdf'))

# post bins all not just jumps
pbins = [0,2,4,6,8,16,32]

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
        if dchoice =='hDeltaC':
            thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
            epg = cxa.pdat['offset_eb_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        elif dchoice =='hDeltaC_FC2_EPG':
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
    plt.hist2d(epgall,hdcall,plotbins,cmap='grey_r',density=True,vmin=0,vmax=4e-5)
    plt.yticks([-180,-90,0,90,180])
    plt.xticks([-180,-90,0,90,180])
    plt.xlabel('EPG phase (deg)')
    plt.ylabel('hDeltaC phas (deg)')
    plt.plot([-180,180],[-90,-90],color='r',linestyle='--')
    plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='--')
    plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='--')
    plt.colorbar()
    plt.title('Post odour '+str(b1/10) +' - '+str(b2/10))
    # post bins
    plt.savefig(os.path.join(savedir,'PostAirAll_hist_'+ str(b1/10) +' - '+str(b2/10) + '.png'))
    plt.savefig(os.path.join(savedir,'PostAirAll_hist_'+ str(b1/10) +' - '+str(b2/10)   + '.pdf'))

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
        if dchoice =='hDeltaC':
            thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
            epg = cxa.pdat['offset_eb_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        elif dchoice =='hDeltaC_FC2_EPG':
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
    plt.hist2d(epgall,hdcall,plotbins,cmap='grey_r',density=True,vmin=0,vmax=4e-5)
    plt.yticks([-180,-90,0,90,180])
    plt.xticks([-180,-90,0,90,180])
    plt.xlabel('EPG phase (deg)')
    plt.ylabel('hDeltaC phas (deg)')
    plt.plot([-180,180],[-90,-90],color='r',linestyle='--')
    plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='--')
    plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='--')
    plt.colorbar()
    plt.title('Post odour '+str(b1/10) +' - '+str(b2/10))
    # post bins
    plt.savefig(os.path.join(savedir,'PostAirJump_hist_'+ str(b1/10) +' - '+str(b2/10) + '.png'))
    plt.savefig(os.path.join(savedir,'PostAirJump_hist_'+ str(b1/10) +' - '+str(b2/10)   + '.pdf'))




#%% Panel B - 2D hist of phases with EPG normalization


plt.close('all')
dchoice = 'hDeltaC'
data = all_flies[dchoice]
# Scatter pre
epgall = np.array([])
hdcall = np.array([])
vmax = .15
for f in data:
    cxa = data[f].cxa
    cxa.get_jumps()
    j = cxa.get_entries_exits_like_jumps()
    if dchoice =='hDeltaC':
        thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        epg = cxa.pdat['offset_eb_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
    elif dchoice =='hDeltaC_FC2_EPG':
        thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
        epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
    
    hdcall = np.append(hdcall,thdc)
    
    epgall = np.append(epgall,epg)
    
epgall = 180*epgall/np.pi
hdcall = 180*hdcall/np.pi
plotbins = np.linspace(-180,180,21)
h,xb,yb = np.histogram2d(epgall,hdcall,plotbins)
hs = np.sum(h,axis=1)[:,np.newaxis]
hsT = np.tile(hs,(1,len(hs)))
hnorm = h/hsT
extent = [xb[0], xb[-1], yb[0], yb[-1]]
plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=vmax)


plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='--')
plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='--')
plt.yticks([-180,-90,0,90,180])
plt.xticks([-180,-90,0,90,180])
plt.xlabel('EPG phase (deg)')
plt.ylabel('hDeltaC phas (deg)')
plt.colorbar()
plt.title('Pre odour')



# post bins all not just jumps
pbins = [0,2,4,6,8,16,32,60]

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
        if dchoice =='hDeltaC':
            thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
            epg = cxa.pdat['offset_eb_phase'].to_numpy()*cxa.side*-1
        elif dchoice =='hDeltaC_FC2_EPG':
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
    plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=vmax)
    plt.yticks([-180,-90,0,90,180])
    plt.xticks([-180,-90,0,90,180])
    plt.xlabel('EPG phase (deg)')
    plt.ylabel('hDeltaC phas (deg)')
    plt.plot([-180,180],[-90,-90],color='r',linestyle='--')
    plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='--')
    plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='--')
    plt.colorbar()
    plt.title('Post odour '+str(b1/10) +' - '+str(b2/10))
    # post bins
    plt.savefig(os.path.join(savedir,'PostAirAll_hist_'+ str(b1/10) +' - '+str(b2/10) + '_colnorm.png'))
    plt.savefig(os.path.join(savedir,'PostAirAll_hist_'+ str(b1/10) +' - '+str(b2/10)   + '_colnorm.pdf'))

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
        if dchoice =='hDeltaC':
            thdc = cxa.pdat['offset_fsb_upper_phase'].to_numpy()*cxa.side*-1
            epg = cxa.pdat['offset_eb_phase'].to_numpy()*cxa.side*-1
        elif dchoice =='hDeltaC_FC2_EPG':
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
    plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=vmax)
    
    
    plt.yticks([-180,-90,0,90,180])
    plt.xticks([-180,-90,0,90,180])
    plt.xlabel('EPG phase (deg)')
    plt.ylabel('hDeltaC phas (deg)')
    plt.plot([-180,180],[-90,-90],color='r',linestyle='--')
    plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='--')
    plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='--')
    plt.colorbar()
    plt.title('Post odour '+str(b1/10) +' - '+str(b2/10))
    # post bins
    plt.savefig(os.path.join(savedir,'PostAirJump_hist_'+ str(b1/10) +' - '+str(b2/10) + '_colnorm.png'))
    plt.savefig(os.path.join(savedir,'PostAirJump_hist_'+ str(b1/10) +' - '+str(b2/10)   + '_colnorm.pdf'))
#%%  Panel B - 2D hist of FC2 and hDeltaC phases with EPG normalization
pbins = [0,2,8,16,32,60]

plt.close('all')
dchoice = 'hDeltaC_FC2_EPG'
data = all_flies[dchoice]
# Scatter pre
epgall = np.array([])
hdcall = np.array([])
fc2all = np.array([])
vmax = .25
bins = 16
for f in data:
    cxa = data[f].cxa
    cxa.get_jumps()
    j = cxa.get_entries_exits_like_jumps()
    
    
    thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
    epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
    tfc2 = cxa.pdat['offset_fsb1_ch1_phase'].to_numpy()[:j[0,0]]*cxa.side*-1
    hdcall = np.append(hdcall,thdc)
    fc2all = np.append(fc2all,tfc2)
    epgall = np.append(epgall,epg)
    
epgall = 180*epgall/np.pi
hdcall = 180*hdcall/np.pi
fc2all = 180*fc2all/np.pi
plotbins = np.linspace(-180,180,bins)
plt.figure()
for i in range(3):
    plt.subplot(2,2,i+1)
    if i==0:
        h,xb,yb = np.histogram2d(epgall,hdcall,plotbins)
        medians, edges, bin_numbers = stats.binned_statistic(
             epgall,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
        )        
    elif i ==1:
        h,xb,yb = np.histogram2d(epgall,fc2all,plotbins)
        medians, edges, bin_numbers = stats.binned_statistic(
             epgall,fc2all, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
        )
    elif i==2:
        h,xb,yb = np.histogram2d(fc2all,hdcall,plotbins)
        medians, edges, bin_numbers = stats.binned_statistic(
             fc2all,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
        )
    hs = np.sum(h,axis=1)[:,np.newaxis]
    hsT = np.tile(hs,(1,len(hs)))
    hnorm = h/hsT
    extent = [xb[0], xb[-1], yb[0], yb[-1]]
    
    
    plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=vmax)
    pedge = plotbins[1:]-np.mean(np.diff(plotbins))/2
    plt.scatter(pedge,medians,s=10,color='r')
    
    plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='-')
    plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='-')
    plt.plot([-180,180],[-180,180],color=[.2,.2,1],linestyle='--')
    plt.yticks([-180,-90,0,90,180])
    plt.xticks([-180,-90,0,90,180])
    plt.xlabel('EPG phase (deg)')
    if i==0:
        plt.ylabel('hDeltaC phase (deg)')
    elif i ==1:
        plt.ylabel('FC2 phase (deg)')
    elif i==2:
        plt.xlabel('FC2 phase (deg)')
        plt.ylabel('hDeltaC phase (deg)')
    plt.colorbar()
    plt.title('Pre odour')



# post bins all not just jumps
pbins = [0,2,4,6,8,16,32,60]

for ip,p in enumerate(pbins):
    epgall = np.array([])
    hdcall = np.array([])
    fc2all = np.array([])
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
        
        
        thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()*cxa.side*-1
        epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()*cxa.side*-1
        tfc2 = cxa.pdat['offset_fsb1_ch1_phase'].to_numpy()*cxa.side*-1
            
            
        dx = np.array([],dtype='int')
        
        for ij,j in enumerate(jumps):
            d1 = min(j[1]+b1,j[2])
            d2 = min(j[1]+b2,j[2])
            if d1==d2:
                continue
            dx = np.append(dx,np.arange(d1,d2))
            
        
        
        hdcall = np.append(hdcall,thdc[dx])
        epgall = np.append(epgall,epg[dx])
        fc2all = np.append(fc2all,tfc2[dx])

    epgall = 180*epgall/np.pi
    hdcall = 180*hdcall/np.pi
    fc2all = 180*fc2all/np.pi
    #plotbins = np.linspace(-180,180,20)
    plt.figure()
    for i in range(3):
        plt.subplot(2,2,i+1)
        if i==0:
            h,xb,yb = np.histogram2d(epgall,hdcall,plotbins)
            medians, edges, bin_numbers = stats.binned_statistic(
                 epgall,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
            )        
        elif i ==1:
            h,xb,yb = np.histogram2d(epgall,fc2all,plotbins)
            medians, edges, bin_numbers = stats.binned_statistic(
                 epgall,fc2all, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
            )
        elif i==2:
            h,xb,yb = np.histogram2d(fc2all,hdcall,plotbins)
            medians, edges, bin_numbers = stats.binned_statistic(
                 fc2all,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
            )
        hs = np.sum(h,axis=1)[:,np.newaxis]
        hsT = np.tile(hs,(1,len(hs)))
        hnorm = h/hsT
        extent = [xb[0], xb[-1], yb[0], yb[-1]]
        
        
        plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=vmax)
        pedge = plotbins[1:]-np.mean(np.diff(plotbins))/2
        plt.scatter(pedge,medians,s=10,color='r')
        
        plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='-')
        plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='-')
        plt.plot([-180,180],[-180,180],color=[.2,.2,1],linestyle='--')
        plt.yticks([-180,-90,0,90,180])
        plt.xticks([-180,-90,0,90,180])
        plt.xlabel('EPG phase (deg)')
        if i==0:
            plt.ylabel('hDeltaC phase (deg)')
        elif i ==1:
            plt.ylabel('FC2 phase (deg)')
        elif i==2:
            plt.xlabel('FC2 phase (deg)')
            plt.ylabel('hDeltaC phase (deg)')
        plt.colorbar()
        plt.title('Post odour '+str(b1/10) +' - '+str(b2/10))
    # post bins
    #plt.savefig(os.path.join(savedir,'PostAirAll_hist_'+ str(b1/10) +' - '+str(b2/10) + '_colnorm.png'))
   #plt.savefig(os.path.join(savedir,'PostAirAll_hist_'+ str(b1/10) +' - '+str(b2/10)   + '_colnorm.pdf'))

# for just jump returns

for ip,p in enumerate(pbins):
    epgall = np.array([])
    hdcall = np.array([])
    fc2all = np.array([])
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
        thdc = cxa.pdat['offset_fsb2_ch2_phase'].to_numpy()*cxa.side*-1
        epg = cxa.pdat['offset_eb_ch1_phase'].to_numpy()*cxa.side*-1
        tfc2 = cxa.pdat['offset_fsb1_ch1_phase'].to_numpy()*cxa.side*-1
        dx = np.array([],dtype='int')
        
        for ij,j in enumerate(jumps):
            d1 = min(j[1]+b1,j[2])
            d2 = min(j[1]+b2,j[2])
            if d1==d2:
                continue
            dx = np.append(dx,np.arange(d1,d2))
            
        
        
        hdcall = np.append(hdcall,thdc[dx])
        epgall = np.append(epgall,epg[dx])
        fc2all = np.append(fc2all,tfc2[dx])

    epgall = 180*epgall/np.pi
    hdcall = 180*hdcall/np.pi
    fc2all = 180*fc2all/np.pi
   # plotbins = np.linspace(-180,180,20)
    plt.figure()
    for i in range(3):
        plt.subplot(2,2,i+1)
        if i==0:
            h,xb,yb = np.histogram2d(epgall,hdcall,plotbins)
            medians, edges, bin_numbers = stats.binned_statistic(
                 epgall,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
            )        
        elif i ==1:
            h,xb,yb = np.histogram2d(epgall,fc2all,plotbins)
            medians, edges, bin_numbers = stats.binned_statistic(
                 epgall,fc2all, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
            )
        elif i==2:
            h,xb,yb = np.histogram2d(fc2all,hdcall,plotbins)
            medians, edges, bin_numbers = stats.binned_statistic(
                 fc2all,hdcall, statistic=lambda a: stats.circmean(a,high=180,low=-180), bins=plotbins
            )
        hs = np.sum(h,axis=1)[:,np.newaxis]
        hsT = np.tile(hs,(1,len(hs)))
        hnorm = h/hsT
        extent = [xb[0], xb[-1], yb[0], yb[-1]]
        
        
        plt.imshow(hnorm.T, interpolation='none',origin='lower',extent=extent,cmap='grey_r',vmin=0,vmax=vmax)
        pedge = plotbins[1:]-np.mean(np.diff(plotbins))/2
        plt.scatter(pedge,medians,s=10,color='r')
        
        plt.plot([0,0],[-180,180],color=[.2,.2,1],linestyle='-')
        plt.plot([-180,180],[0,0],color=[.2,.2,1],linestyle='-')
        plt.plot([-180,180],[-180,180],color=[.2,.2,1],linestyle='--')
        plt.yticks([-180,-90,0,90,180])
        plt.xticks([-180,-90,0,90,180])
        plt.xlabel('EPG phase (deg)')
        if i==0:
            plt.ylabel('hDeltaC phase (deg)')
        elif i ==1:
            plt.ylabel('FC2 phase (deg)')
        elif i==2:
            plt.xlabel('FC2 phase (deg)')
            plt.ylabel('hDeltaC phase (deg)')
        plt.colorbar()
        plt.title('Post odour '+str(b1/10) +' - '+str(b2/10))
    # post bins
    #plt.savefig(os.path.join(savedir,'PostAirJump_hist_'+ str(b1/10) +' - '+str(b2/10) + '_colnorm.png'))
    #plt.savefig(os.path.join(savedir,'PostAirJump_hist_'+ str(b1/10) +' - '+str(b2/10)   + '_colnorm.pdf'))
#%% 
