# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 14:35:23 2024

@author: dowel
"""


import numpy as np
import matplotlib.pyplot as plt 
from analysis_funs.CX_analysis_col import CX_a
import pandas as pd
from Utilities.utils_general import utils_general as ug
from scipy.stats import circmean, circstd
from Utilities.utils_plotting import uplt
from scipy import stats
from rdp import rdp
#%%
class ET_paper:
    def __init__(self,datadir,regions=['eb','fsb_upper','fsb_lower']):
        self.cxa = CX_a(datadir,regions= regions,denovo=False)
        self.colours =np.array([[81,156,204],[84,39,143],[0,0,0],[237,30,36],[6,149,207]])/255
        self.colours2 =np.array([[106,207,246],[237,30,36],[229,229,229],[6,149,207]])/255
    def example_trajectory(self,tstart,tend,boxes=[],heat=False,col=None,cmap='coolwarm',cmin=-1,cmax=1):
        colours = np.array([[106,207,246],[237,30,36],[204,204,204],[6,149,207]])/255
        ft2 = self.cxa.ft2
        x = ft2['ft_posx'].to_numpy()
        y = ft2['ft_posy'].to_numpy()
        x,y = self.cxa.fictrac_repair(x,y)
        ins = ft2['instrip'].to_numpy()
        jumps = ft2['jump'].to_numpy()
        tt = self.cxa.pv2['relative_time'].to_numpy()
        tdx = np.logical_and( tt>tstart,tt<tend)
        inplume = ins>0
        st  = np.where(inplume)[0][0]
        x = x-x[st-1]
        y = y-y[st-1]
        ends = np.where(inplume)[0][-1]
        yend = y[ends]
        
        x = x[tdx]
        y = y[tdx]
        ins = ins[tdx]
        jumps = jumps[tdx]
        inplume = inplume[tdx]
        
        jumps = jumps-np.mod(jumps,3)
        jd = np.diff(jumps)
        jn = np.where(np.abs(jd)>0)[0]+1
        print(jumps[jn])
        jkeep = np.where(np.diff(jn)>1)[0]
       
       #jn = jn[jkeep]
        print(jumps[jn])
        #x,y = self.fictrac_repair(x,y)
        
        
        xrange = np.max(x)-np.min(x)
        yrange = np.max(y)-np.min(y)
        
        mrange = np.max([xrange,yrange])+100
        y_med = yrange/2
        x_med = xrange/2
        ylims = [y_med-mrange/2, y_med+mrange/2]
   
        xlims = [x_med-mrange/2, x_med+mrange/2]

        plt.rcParams['ps.fonttype'] = 42 
        fig = plt.figure(figsize=(15,15))
        
        ax = fig.add_subplot(111)
        #ax.scatter(x[inplume],y[inplume],color=[0.5, 0.5, 0.5])
        
        yj = y[jn]
        yj = np.append(yj,np.min([y[-1],yend]))
        plt.fill()
        
        plume_cents = [-630,-420,-210,0,210,420,630]
        for p in plume_cents:
            tj = 0
            x1 = p+5+tj
            x2 = p-5+tj
            y1 = 0
            y2 = yj[0]
            try:
                xvec = np.array([x1,x2,x2,x1])+jumps[jn[0]-1]
            except:
                xvec = np.array([x1,x2,x2,x1])
            yvec = [y1,y1,y2,y2]
            
            
            cents = np.arange(0,np.max(np.abs(x)),210)
            cents = np.append(np.flipud(-cents[1:]),cents)
           # cents = [-1050,-840,-630,-420,-210, 0,210,420,630,840,1050]
            
            plt.fill(xvec,yvec,color=colours[2,:])
            for c in cents:
                plt.fill(xvec+c,yvec,color=colours[2,:])
            for i,j in enumerate(jn):
                
                tj = jumps[j]
                print(tj)
                x1 = p+5+tj
                x2 = p-5+tj
                y1 = yj[i]
                y2 = yj[i+1]
                xvec = np.array([x1,x2,x2,x1])
                yvec = [y1,y1,y2,y2]
                for c in cents:
                    plt.fill(xvec+c,yvec,color=colours[2,:])
                    plt.plot([xvec[1],xvec[1]-20],yvec[0:2],color=colours[3,:],linewidth=0.75)
                    plt.scatter(xvec[1]-20,yvec[0],color=colours[3,:],marker='<')

        

        # plot in and out plume
        
        if heat:
            
            ax = plt.gca()
            hplot = col.copy()
            hplot[hplot<cmin] = cmin
            hplot[hplot>cmax] = cmax
            hplot[0] = cmin
            hplot[-1] = cmax
            hplot = hplot[tdx]
            uplt.coloured_line(x,y,hplot,ax,cmap=cmap)
            
            
        else:
            plt.plot(x,y,color=colours[0,:],linewidth=0.5)
            pstart,psize = ug.find_blocks(inplume)
            for i,p in enumerate(pstart):
                plt.plot(x[p:p+psize[i]],y[p:p+psize[i]],color=colours[1,:],linewidth=0.75)

        plt.xlim(xlims)
        plt.ylim(ylims)
        plt.xlabel('x position (mm)')
        plt.ylabel('y position (mm)')
        x1 = np.min(x)-30
        x2 = np.max(x)+30
        plt.xlim([x1,x2])
        plt.ylim([min(y),max(y)])
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        for b in boxes:
            dx1 = ug.find_nearest(tt,b[0])
            print(dx1)
            dx2 = ug.find_nearest(tt,b[1])
            y1 = y[dx1]
            y2 = y[dx2]
            x1 = y[dx1]
            plt.plot(np.array([-30,30,30,-30,-30])+x1,[y1,y1,y2,y2,y1],color='k',linestyle='--',linewidth=0.75)
        
        
        plt.show()
    def example_trajectory_arrows(self,tstart,tend,colours,scale=10,samp_period=10,boxes=[],regions=['eb','fsb_upper'],heat=False,col=None,cmap='coolwarm',cmin=-1,cmax=1):
        self.example_trajectory(tstart,tend,heat=heat,col=col,cmap=cmap,cmin=cmin,cmax=cmax)
        ft2 = self.cxa.ft2
        x = ft2['ft_posx'].to_numpy()
        y = ft2['ft_posy'].to_numpy()
        x,y = self.cxa.fictrac_repair(x,y)
        ins = ft2['instrip'].to_numpy()
        tt = self.cxa.pv2['relative_time'].to_numpy()
        tdx = np.logical_and( tt>tstart,tt<tend)
        inplume = ins>0
        st  = np.where(inplume)[0][0]
        x = x-x[st-1]
        y = y-y[st-1]
          
        x = x[tdx]
        y = y[tdx]
        
        for ir,r in enumerate(regions):
            tphase = self.cxa.pdat['offset_'+r+'_phase'].to_numpy()[tdx]
            dx = np.arange(0,len(tphase),samp_period*10)
            for d in dx:
                tx = np.sin(tphase[d])*scale
                ty = np.cos(tphase[d])*scale
                px = np.array([0,tx])+x[d]
                py = np.array([0,ty])+y[d]
                plt.plot(px,py,color=colours[ir,:])
        
    def example_trajectory_pva(self,tstart,tend,boxes=[],region='fsb_upper',cmap='coolwarm',crange=[-1,1]):
        colours = np.array([[106,207,246],[237,30,36],[168,170,173],[6,149,207]])/255
        ft2 = self.cxa.ft2
        x = ft2['ft_posx'].to_numpy()
        y = ft2['ft_posy'].to_numpy()
        x,y = self.cxa.fictrac_repair(x,y)
        ins = ft2['instrip'].to_numpy()
        jumps = ft2['jump'].to_numpy()
        jumps2 = self.cxa.get_jumps()
        
        tt = self.cxa.pv2['relative_time'].to_numpy()
        wedges= self.cxa.pdat['wedges_'+region]
        pva = ug.get_pvas(wedges)
        
        tdx = np.logical_and( tt>tstart,tt<tend)
        tdxw = np.where(tdx)[0]
        jkdx = np.logical_and(jumps2[:,1]>tdxw[0],jumps2[:,1]<tdxw[-1])
        jumps2 = jumps2[jkdx,:]
        jumps2 = jumps2-tdxw[0]
        print(jumps2)
        
        inplume = ins>0
        st  = np.where(inplume)[0][0]
        x = x-x[st-1]
        y = y-y[st-1]
        ends = np.where(inplume)[0][-1]
        yend = y[ends]
        x = x[tdx]
        y = y[tdx]
        pva = pva[tdx]
        ins = ins[tdx]
        jumps = jumps[tdx]
        inplume = inplume[tdx]
        
        jumps = jumps-np.mod(jumps,3)
        jd = np.diff(jumps)
        jn = np.where(np.abs(jd)>0)[0]+1
        print(jumps[jn])
        jkeep = np.where(np.diff(jn)>1)[0]
       
       #jn = jn[jkeep]
        print(jumps[jn])
        #x,y = self.fictrac_repair(x,y)
        
        
        xrange = np.max(x)-np.min(x)
        yrange = np.max(y)-np.min(y)
        
        mrange = np.max([xrange,yrange])+100
        y_med = yrange/2
        x_med = xrange/2
        ylims = [y_med-mrange/2, y_med+mrange/2]
   
        xlims = [x_med-mrange/2, x_med+mrange/2]

        plt.rcParams['ps.fonttype'] = 42 
        fig = plt.figure(figsize=(15,15))
        
        ax = fig.add_subplot(111)
        #ax.scatter(x[inplume],y[inplume],color=[0.5, 0.5, 0.5])
        
        yj = y[jn]
        print(yend)
        print(y[-1])
        
        yj = np.append(yj,np.min(y[-1],yend))
        plt.fill()
        tj = 0
        x1 = 0+5+tj
        x2 = 0-5+tj
        y1 = 0
        y2 = yj[0]
        xvec = np.array([x1,x2,x2,x1])
        yvec = [y1,y1,y2,y2]
        
        cents = [-630,-420,-210, 0,210,420,630]
        
        plt.fill(xvec,yvec,color=colours[2,:])
        for c in cents:
            plt.fill(xvec+c,yvec,color=colours[2,:])
        for i,j in enumerate(jumps2):
            
            tj = jumps[j[1]]
            print(tj)
            x1 = 0+5+tj
            x2 = 0-5+tj
            y1 = yj[i]
            y2 = yj[i+1]
            xvec = np.array([x1,x2,x2,x1])
            yvec = [y1,y1,y2,y2]
            for c in cents:
                plt.fill(xvec+c,yvec,color=colours[2,:])
                plt.plot([xvec[1],xvec[1]-20],yvec[0:2],color=colours[3,:],linewidth=0.75)
                plt.scatter(xvec[1]-20,yvec[0],color=colours[3,:],marker='<')

            jbef = np.arange(j[1]-10,j[1])
            jaf = np.arange(j[2]-10,j[2])
            plt.plot(x[jbef],y[jbef],color='g',linewidth=5)
            plt.plot(x[jaf],y[jaf],color='g',linewidth=5)
        # plot in and out plume
        
        uplt.coloured_line_simple(x,y,pva,cmap,crange[0],crange[1])
        # plt.plot(x,y,color=colours[0,:],linewidth=0.5)
        pstart,psize = ug.find_blocks(inplume)
        # for i,p in enumerate(pstart):
        #     plt.plot(x[p:p+psize[i]],y[p:p+psize[i]],color=colours[1,:],linewidth=0.75)

        plt.xlim(xlims)
        plt.ylim(ylims)
        plt.xlabel('x position (mm)')
        plt.ylabel('y position (mm)')
        x1 = np.min(x)-30
        x2 = np.max(x)+30
        plt.xlim([x1,x2])
        plt.ylim([min(y),max(y)])
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        for b in boxes:
            dx1 = ug.find_nearest(tt,b[0])
            print(dx1)
            dx2 = ug.find_nearest(tt,b[1])
            y1 = y[dx1]
            y2 = y[dx2]
            plt.plot([-30,30,30,-30,-30],[y1,y1,y2,y2,y1],color='k',linestyle='--',linewidth=0.75)
        
        
        plt.show()
    def example_heatmaps(self,tstart,tend,regions=['eb','fsb_upper'],cmaps=['Greys','Blues'],yoking=True,
                         colours = np.array([[0,0,0],[81,156,204],[237,30,36],[6,149,207]])/255 ):
        pdat = self.cxa.pdat
        tt = self.cxa.pv2['relative_time'].to_numpy()
        tdx = np.logical_and( tt>tstart,tt<tend)
        for i,r in enumerate(regions):
            plt.figure(figsize=(1,20))
            if yoking:
                dat = pdat['wedges_offset_'+r][tdx,:]
            else:
                dat = pdat['wedges_'+r][tdx,:]
            dat = np.flipud(dat)
            
            plt.imshow(dat,cmap=cmaps[i],interpolation='None',aspect='auto',vmin=0,vmax=1)
        
        plt.figure(figsize=(1,20))
        
        heading = self.cxa.ft2['ft_heading'].to_numpy()
        heading = heading[tdx]*180/np.pi
        y = tt[tdx]
        y = y-y[0]
        #plt.plot(heading,y,color='k')
        for i,r in enumerate(regions):
            if yoking:
                x = pdat['offset_'+r+'_phase'].to_numpy()
            else:
                x = pdat['phase_'+r]
            x = x[tdx]
            x = 180*x/np.pi
            
            
            plt.scatter(x,y,color=colours[i,:],s=5)
        
        plt.xlim([-180,180])
        plt.ylim([min(y),max(y)])
        plt.xticks([-180,0,180])
        
        plt.figure(figsize=(1,20))
        x = self.cxa.ft2['instrip'][tdx].to_numpy()
        #plt.plot(x,y,color=colours[2,:])
        plt.fill_betweenx(y,x.astype('float'),np.zeros_like(x,dtype='float'),color=colours[2,:],step='mid',linewidth=0, edgecolor='none')
        jumps = self.cxa.get_jumps(time_threshold=10000)
        jumps = jumps[:,1]
        twhere = np.where(tdx)[0]
        jumps = jumps-twhere[0]
        jumps = jumps[jumps>0]
        jumps = jumps[jumps<len(y)]
        print(jumps)
        plt.scatter(x[jumps]*0+1.1,y[jumps],marker ='<',color=colours[3,:])
        plt.ylim([min(y),max(y)])
        for i,j in enumerate(jumps):
            plt.plot([1.1,1.3],[y[j],y[j]],color=colours[3,:])
    def mean_wedge(self,dataperiod,timespan,regions=['eb','fsb_upper']):
        ft2 = self.cxa.ft2
        pdat = self.cxa.pdat
        wedge_mean = np.zeros((16,len(regions)))
        
        tdx = self.get_dxs(dataperiod,timespan)

        if dataperiod=='jumps' or dataperiod =='pre jumps':
            if self.cxa.side>0:
                for i,r in enumerate(regions):
                    pdat['wedges_offset_'+r] = np.fliplr(pdat['wedges_offset_'+r])
            
                
                
        
        for i,r in enumerate(regions):
            
            wedge_mean[:,i] = np.mean(pdat['wedges_offset_'+r][tdx,:],axis=0)
        return wedge_mean
    def mean_wedge_rotate(self,dataperiod,timespan,regions=['eb','fsb_upper']):
        ft2 = self.cxa.ft2
        pdat = self.cxa.pdat
        
        
        
        wedge_mean = np.zeros((16,len(regions)))
        
        tdx = self.get_dxs(dataperiod,timespan)
    
        if dataperiod=='jumps' or dataperiod =='pre jumps':
            if self.cxa.side>0:
                for i,r in enumerate(regions):
                    pdat['wedges_'+r] = np.fliplr(pdat['wedges_offset_'+r])
                    pdat['phase_'+r] = -pdat['phase_'+r]
                
                
        
        for i,r in enumerate(regions):
            tw = pdat['wedges_'+r]
            tp = pdat['phase_'+r]
            rect_p = (tp+np.pi)/(2*np.pi)
            rpr = np.round(rect_p*15)
            rpr = rpr.astype('int')
            tw_new = np.zeros(tw.shape)
            for ip,p in enumerate(tp):
                
                o = rpr[ip]
                
                tw_new[ip,:] = np.append(tw[ip,o:],tw[ip,:o])
            tw_new = np.append(tw_new[:,9:],tw_new[:,:9],axis=1)
            # plt.figure()
            # plt.imshow(tw,interpolation='None',aspect='auto',vmin=0,vmax=1)
            # plt.figure()
            # plt.imshow(tw_new,interpolation='None',aspect='auto',vmin=0,vmax=1)
            wedge_mean[:,i] = np.mean(tw_new[tdx,:],axis=0)
        return wedge_mean
    def output_phases(self,dataperiod,timespan,flyid,regions=['eb','fsb_upper']):
        
        pdat = self.cxa.pdat
        tdx = self.get_dxs(dataperiod,timespan)
        for i,r in enumerate(regions):
            phase = pdat['offset_'+r+'_phase'].to_numpy()
            if dataperiod=='jumps':
                phase = -phase*self.cxa.side
            if i==0:
                phases = np.zeros((len(tdx),len(regions)+1))
            phases[:,i] =phase[tdx]
            
        phases[:,-1] = flyid
        return phases
    def phase_histogram(self,dataperiod,timespan,bins,regions=['eb','fsb_upper']):
        pdat = self.cxa.pdat
        tdx = self.get_dxs(dataperiod,timespan)
        
        
        pdiff = pdat['phase_'+regions[1]]-pdat['phase_'+regions[0]]
        
        pdiff = pdiff[tdx]
        if dataperiod=='jumps' or dataperiod=='pre jumps':
            pdiff=-pdiff*self.cxa.side
        pcos = np.cos(pdiff)
        psin = np.sin(pdiff)
        pdiffa = np.arctan2(psin,pcos)
        h = np.histogram(pdiffa,bins = bins,density=True)
        
        
        ph = h[0]
        hists = np.zeros((len(ph),3))
        hists[:,0] = ph/np.sum(ph)
        for i,r in enumerate(regions):
            phase = pdat['offset_'+r+'_phase'].to_numpy()
            phase = phase[tdx]
            if dataperiod=='jumps' or dataperiod == 'pre jumps':
                phase=-phase*self.cxa.side
            h = np.histogram(phase,bins=bins,density=True)
            ph = h[0]
            hists[:,i+1] = ph/np.sum(ph)
        return hists
    def get_dxs(self,dataperiod,timespan):
        ft2 = self.cxa.ft2
        if dataperiod=='pre air':
            ins = ft2['instrip'].to_numpy()
            st = np.where(ins>0)[0][0]
            tdx =  np.arange(0,st,dtype='int')
        elif dataperiod =='pre air anemotaxis':
            ins = ft2['instrip'].to_numpy()
            from rdp import rdp
            from Utilities.utils_general import utils_general as utg
            x = ft2['ft_posx'].to_numpy()
            y = ft2['ft_posy'].to_numpy()
            
            x,y = self.cxa.fictrac_repair(x, y)
            
            st = np.where(ins>0)[0][0]
            y = y-y[st]
            x = x[0:st]
            y = y[0:st]
            x = x[:,np.newaxis]
            y = y[:,np.newaxis]
            arr = np.append(x,y,axis=1)
            arr2 = rdp(arr,epsilon=10)
            ardx  = np.array([utg.find_nearest(i,y) for i in arr2[:,1]])
            ad = np.diff(arr2,axis=0)
            
            ad_abs = np.sqrt(np.sum(np.square(ad),axis=1))
            
            adx = np.where(ad_abs>50)[0].astype('int') # from Mussels-Pires et al 2024
            
            tlist = [np.arange(ardx[i],ardx[i+1]) for i in adx]
            
            tdx = np.array([])
            for t in tlist:
                tdx = np.append(tdx,t)
            tdx = tdx.astype('int')
            #tdx = np.unique(tdx,dtype='int')
            
            
        elif dataperiod=='jumps':
            jumps = self.cxa.get_jumps()
            
            tdx = np.array([],dtype='int')
            tt = self.cxa.pv2['relative_time']
            td = np.mean(np.diff(tt))
            onesec = np.round(1/td)
            bins = int(onesec*timespan)
            for j in jumps:
                jrange = np.arange(j[1],j[2],dtype='int')
                if len(jrange)>bins:
                    jrange = jrange[-bins:]
                
                tdx = np.append(tdx,jrange)
        elif dataperiod=='pre jumps':
            jumps = self.cxa.get_jumps()
            
            tdx = np.array([],dtype='int')
            tt = self.cxa.pv2['relative_time']
            td = np.mean(np.diff(tt))
            onesec = np.round(1/td)
            bins = int(onesec*timespan)
            for j in jumps:
                jrange = np.arange(j[0],j[1],dtype='int')
                if len(jrange)>bins:
                    jrange = jrange[-bins:]
                
                tdx = np.append(tdx,jrange)
                
        return tdx
        
    def phase_time(self,t_length,regions = ['eb','fsb_upper'],rebase=False):
        jumps = self.cxa.get_jumps()
        tbins = t_length*10
        phases = np.empty((tbins,4,len(jumps)))
        for i,j in enumerate(jumps):
           # print(j)
            jend = j[2]
            jstart = j[1]
            jrange = np.arange(jstart,jend,dtype=int)
            if len(jrange)>tbins:
                
                #print('Too long', len(jrange), tbins)
                jrange = jrange[-tbins:]
                #print(jrange)
            for ir, r in enumerate(regions):
                tphase = -self.cxa.pdat['offset_'+r+'_phase'].to_numpy()*self.cxa.side
                
                
                if rebase:
                    rdx = np.arange(j[1]-1,j[1]+5)
                    toff = -self.cxa.pdat['offset_'+r+'_phase'].to_numpy()[rdx]*self.cxa.side
                    toff = circmean(toff,low=-np.pi,high=np.pi)
                    phases[-len(jrange):,ir,i] = ug.circ_subtract(tphase[jrange],toff)
                else: 
                    phases[-len(jrange):,ir,i] = tphase[jrange]
                
            tphase = -self.cxa.ft2['ft_heading'].to_numpy()*self.cxa.side
            phases[0:len(jrange),2,i] = tphase[jrange]
            
        phases[:,3,:] = phases[:,1,:]-phases[:,0,:]
        #phases = np.flipud(phases)
        return phases
    
    def phase_time_amp(self,t_length,regions = ['eb','fsb_upper']):
        # from plume entry
        jumps = self.cxa.get_jumps()
        tbins = t_length*10
        phases = np.empty((tbins,len(regions),len(jumps)))*np.nan
        amps = np.empty((tbins,len(regions),len(jumps)))*np.nan
        for ir,r in enumerate(regions):
            phase = self.cxa.pdat['offset_'+r+'_phase'].to_numpy()*-1*self.cxa.side
            amp = np.mean(self.cxa.pdat['wedges_'+r],axis=1)
            for i,j in enumerate(jumps):
                dx = np.arange(max(j[2]-tbins,j[1]),j[2])
               
                dx2 = dx-np.min(dx)
                phases[dx2,ir,i] = np.flipud(phase[dx])
                amps[dx2,ir,i] = np.flipud(amp[dx])
                
        return phases,amps
    
    def phase_time_leave_amp(self,t_length,regions=['eb','fsb_upper']):
        # from plume exit
        jumps = self.cxa.get_jumps()
        tbins = t_length*10
        phases = np.empty((tbins,len(regions),len(jumps)))*np.nan
        amps = np.empty((tbins,len(regions),len(jumps)))*np.nan
        for ir,r in enumerate(regions):
            phase = self.cxa.pdat['offset_'+r+'_phase'].to_numpy()*-1*self.cxa.side
            amp = np.mean(self.cxa.pdat['wedges_'+r],axis=1)
            for i,j in enumerate(jumps):
                dx = np.arange(j[1],min(j[2],j[1]+tbins))
               
                dx2 = dx-j[1]
                phases[dx2,ir,i] = phase[dx]
                amps[dx2,ir,i] = amp[dx]
                
        return phases,amps
            
        
        
    def phase_time_inside(self,t_length,regions=['eb','fsb_upper']):
        jumps = self.cxa.get_jumps()
        tbins = t_length*10
        phases = np.empty((tbins,4,len(jumps)))
        inplume = self.cxa.ft2['instrip'].to_numpy()
        ipdiff = np.diff(inplume)
        exits = np.where(ipdiff<0)[0]
        #print(exits)
        for i,j in enumerate(jumps):
            jend = j[2]
            #print(i)
            next_exit = np.min(exits[exits>jend])
            jrange = np.arange(jend,next_exit)
            #print(jrange)
            if len(jrange)>tbins:
                jrange = jrange[:tbins]
            for ir,r in enumerate(regions):
                tphase = -self.cxa.pdat['offset_'+r+'_phase'].to_numpy()*self.cxa.side
                phases[:len(jrange),ir,i] = tphase[jrange]
            tphase = -self.cxa.ft2['ft_heading'].to_numpy()*self.cxa.side
            phases[0:len(jrange),2,i] = tphase[jrange]
            
        phases[:,3,:] = phases[:,1,:]-phases[:,0,:]
        #phases = np.flipud(phases)
        return phases
                
    def amp_time(self,t_length,regions = ['eb','fsb_upper']):
        jumps = self.cxa.get_jumps()
        tbins = t_length*10
        amps = np.empty((tbins,len(regions),3,len(jumps)))
        for i,j in enumerate(jumps):
            # print(j)
            jend = j[2]
            jstart = j[1]
            jrange = np.arange(jstart,jend,dtype=int)
            if len(jrange)>tbins:
                
                #print('Too long', len(jrange), tbins)
                jrange = jrange[-tbins:]
                #print(jrange
            for ir, r in enumerate(regions):
                ampsum =np.sum(self.cxa.pdat['wedges_' + r],axis=1)
                ampmean = np.mean(self.cxa.pdat['wedges_' + r],axis=1)
                amppva = self.cxa.pdat['amp_'+r]
                amps[-len(jrange):,ir,0,i] = ampsum[jrange]
                amps[-len(jrange):,ir,1,i] = ampmean[jrange]
                amps[-len(jrange):,ir,2,i] = amppva[jrange]
            
        return amps
        
    def trajectory_mean_variable(self,variable,side_mult=False, bins=100,time_threshold=60,time_threshold2=0,min_inplume=5):    
        jumps = self.cxa.get_jumps(time_threshold=time_threshold,time_threshold2=time_threshold2)
        inplume = jumps[:,1]-jumps[:,0]
        jdx = inplume>=min_inplume
        jumps = jumps[jdx,:]
        if len(variable.shape)==1:
            variable  = variable[:,np.newaxis]
            
            
        retvar = np.zeros((bins,variable.shape[1],len(jumps)))
        levvar = np.zeros((bins,variable.shape[1],len(jumps)))
        varnum = variable.shape[1]
        if side_mult:
            variable = variable*-self.cxa.side
        retxy = np.zeros((bins,2,len(jumps)))
        levxy = np.zeros((bins,2,len(jumps)))
        
        
        newtime = np.linspace(0,1,bins)
        
        ins = self.cxa.ft2['instrip']
        inst = np.where(ins)[0][0]
        x = self.cxa.ft2['ft_posx'].to_numpy()
        y = self.cxa.ft2['ft_posy'].to_numpy()
        #x,y = self.cxa.fictrac_repair(x,y)
        x = x-x[inst]
        y = y-y[inst]
        x = -x*self.cxa.side
                
        for i,j in enumerate(jumps):
            bdx = np.arange(j[0],j[1])
            adx = np.arange(j[1],j[2])
            
            tx = x-x[j[1]-1]
            ty = y-y[j[1]-1]
            
            rx = tx[adx]
            ry = ty[adx]
            
            lx = tx[bdx]
            ly = ty[bdx]
            
            oldtime = np.linspace(0,1,len(rx))
            retxy[:,0,i] = np.interp(newtime,oldtime,rx)
            retxy[:,1,i] = np.interp(newtime,oldtime,ry)
            
            oldtime = np.linspace(0,1,len(lx))
            levxy[:,0,i] = np.interp(newtime,oldtime,lx)
            levxy[:,1,i] = np.interp(newtime,oldtime,ly)
            
            for iv in range(varnum):
                rv = variable[adx,iv]
                lv = variable[bdx,iv]
        
                
                oldtime = np.linspace(0,1,len(rv))
                retvar[:,iv,i] = np.interp(newtime,oldtime,rv)

                oldtime = np.linspace(0,1,len(lv))
                levvar[:,iv,i] = np.interp(newtime,oldtime,lv)
           
                
            
            
        trajs = np.append(levxy,retxy,axis=0)
        outvar = np.append(levvar,retvar,axis=0)
        return trajs,outvar
        
        
    def trajectory_mean(self,bins = 100,regions=['eb','fsb_upper','ft_heading'],time_threshold=60,time_threshold2=0,min_inplume=5):
        
        print(regions)
        jumps = self.cxa.get_jumps(time_threshold=time_threshold,time_threshold2=time_threshold2)
        inplume = jumps[:,1]-jumps[:,0]
        jdx = inplume>=min_inplume
        jumps = jumps[jdx,:]
        
        retphase = np.zeros((bins,len(regions),len(jumps)))
        levphase = np.zeros((bins,len(regions),len(jumps)))
        retxy = np.zeros((bins,2,len(jumps)))
        levxy = np.zeros((bins,2,len(jumps)))
        retamp = np.zeros( (bins,2,len(regions),len(jumps)))
        levamp = np.zeros((bins,2,len(regions),len(jumps)))
        
        newtime = np.linspace(0,1,bins)
        jseries = self.cxa.ft2['jump'].to_numpy()
        
        ins = self.cxa.ft2['instrip']
        inst = np.where(ins)[0][0]
        x = self.cxa.ft2['ft_posx'].to_numpy()
        y = self.cxa.ft2['ft_posy'].to_numpy()
        #x,y = self.cxa.fictrac_repair(x,y)
        x = x-x[inst]
        y = y-y[inst]
        x = -x*self.cxa.side
        pva = np.zeros((len(x),len(regions)))
        mnflr = np.zeros((len(x),len(regions)))
        
        for ir,r in enumerate(regions):
            pva[:,ir] = ug.get_pvas(self.cxa.pdat['wedges_'+r])
            mnflr[:,ir] = np.mean(self.cxa.pdat['wedges_'+r],axis=1)
        
        for i,j in enumerate(jumps):
            bdx = np.arange(j[0],j[1])
            adx = np.arange(j[1],j[2])
            tj = jseries[j[2]]
            
            tx = x-x[j[1]-1]
            ty = y-y[j[1]-1]
            
            rx = tx[adx]
            ry = ty[adx]
            
            lx = tx[bdx]
            ly = ty[bdx]
            
            oldtime = np.linspace(0,1,len(rx))
            retxy[:,0,i] = np.interp(newtime,oldtime,rx)
            retxy[:,1,i] = np.interp(newtime,oldtime,ry)
            
            oldtime = np.linspace(0,1,len(lx))
            levxy[:,0,i] = np.interp(newtime,oldtime,lx)
            levxy[:,1,i] = np.interp(newtime,oldtime,ly)
            for ip,p in enumerate(regions):
                if p=='ft_heading':
                    phase = self.cxa.ft2[p].to_numpy()
                else:    
                    phase = self.cxa.pdat['offset_'+p+'_phase'].to_numpy()
                
                    
                    
                    
                phase = -phase*self.cxa.side
                rp = phase[adx]
                lp = phase[bdx] 
                
                rpva = pva[adx,ip]
                lpva = pva[bdx,ip]
                
                rmnflr = mnflr[adx,ip]
                lmnflr = mnflr[bdx,ip]
                
                
                oldtime = np.linspace(0,1,len(rp))
                retphase[:,ip,i] = ug.circ_interp(newtime,oldtime,rp)
                retamp[:,0,ip,i] = np.interp(newtime,oldtime,rpva)
                retamp[:,1,ip,i] = np.interp(newtime,oldtime,rmnflr)
                
                oldtime = np.linspace(0,1,len(lp))
                levphase[:,ip,i] = ug.circ_interp(newtime,oldtime,lp)
                levamp[:,0,ip,i] = np.interp(newtime,oldtime,lpva)
                levamp[:,1,ip,i] = np.interp(newtime,oldtime,lmnflr)
                
            
            
        phases = np.append(levphase,retphase,axis=0)
        trajs = np.append(levxy,retxy,axis=0)
        amps = np.append(levamp,retamp,axis=0)
        return phases,trajs,amps
    
    def plt_tmp(self,scalar=1,bins=100,phase_num=10,regions=['eb','fsb_upper'],colours=np.array([[0,0,0],[81,156,205]])/255):
        colours2 = self.colours2
        phases,trajs,amps = self.trajectory_mean(bins=bins,regions=regions)
        tmean =np.mean(trajs,axis=2)
        pmean = circmean(phases,low=-np.pi,high=np.pi,axis=2)
        #trajs = trajs-trajs[bins,:]
        
        pedge = tmean[bins,0]
        ymn = np.min(np.min(trajs[0,1,:]))
        ymx = np.max(trajs[-1,1,:])
        
        print(ymx)
        fx = np.array([-5,5,5,-5])-(5-pedge)
        fy = np.array([ymn,ymn,0,0])
        plt.fill(fx,fy,color=colours2[2,:])
        fy2 = np.array([0,0,ymx,ymx])
        plt.fill(fx-3,fy2,color=colours2[2,:])
        
        plt.plot(fx[:2]-3,[0,0],color=self.colours[4,:])
        plt.plot(trajs[:100,0,:],trajs[:100,1,:],color = colours2[1,:],alpha = 0.1)
        plt.plot(trajs[100:,0,:],trajs[100:,1,:],color = colours2[0,:],alpha = 0.1)
        
        plt.plot(tmean[:100,0],tmean[:100,1],color = colours2[1,:],alpha = 1)
        plt.plot(tmean[99:,0],tmean[99:,1],color = colours2[0,:],alpha = 1)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.show()
        
        c = self.colours
        c[1,:] = colours2[3,:]
        c[0,:] = [0,0,0]
        p_points = np.linspace(0,bins*2-1,phase_num*2,dtype='int')
        for i,p in enumerate(p_points):
            for r in range(len(regions)):
                tp = pmean[p,r]
                tx =tmean[p,0]
                ty = tmean[p,1]
                tx2 = np.sin(tp)*scalar+tx
                ty2 = np.cos(tp)*scalar+ty
                plt.plot([tx,tx2],[ty,ty2],color=colours[r,:])
                
    def plt_tmp_all(self,scalar =1,bins=100,phase_num=10):
        colours2 = self.colours2
        phases,trajs = self.trajectory_mean(bins=bins)
        offset = 0
        
        c = self.colours
        c[1,:] = colours2[3,:]
        c[0,:] = [0,0,0]
        for i in range(trajs.shape[2]):
            
            x  =trajs[:,0,i]
            y = trajs[:,1,i]
            tphase = phases[:,:,i]
            pedge = x[49]
            ymn = np.min(y)
            ymx = np.max(y)
            xmx = np.max(x)
            fx = np.array([-5,5,5,-5])-5#-(5-pedge)
            fy = np.array([ymn,ymn,0,0])
            plt.fill(fx+offset,fy,color=colours2[2,:])
            fy2 = np.array([0,0,ymx,ymx])
            plt.fill(fx-3+offset,fy2,color=colours2[2,:])
            plt.plot(x[:100]+offset,y[:100],color = colours2[1,:])
            plt.plot(x[99:]+offset,y[99:],color = colours2[0,:])
            
            p_points = np.linspace(0,bins*2-1,phase_num*2,dtype='int')
            for i,p in enumerate(p_points):
                for r in range(2):
                    tp = tphase[p,r]
                    tx =x[p]+offset
                    ty = y[p]
                    tx2 = np.sin(tp)*scalar+tx
                    ty2 = np.cos(tp)*scalar+ty
                    plt.plot([tx,tx2],[ty,ty2],color=self.colours[r,:])
                    
            offset += xmx+ 20
            
        
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        plt.show()
       
        
        
        
       
        
        # c = self.colours
        # c[1,:] = colours2[3,:]
        # c[0,:] = [0,0,0]
        # p_points = np.linspace(0,bins*2-1,phase_num*2,dtype='int')
        # for i,p in enumerate(p_points):
        #     for r in range(2):
        #         tp = pmean[p,r]
        #         tx =tmean[p,0]
        #         ty = tmean[p,1]
        #         tx2 = np.sin(tp)*scalar+tx
        #         ty2 = np.cos(tp)*scalar+ty
        #         plt.plot([tx,tx2],[ty,ty2],color=self.colours[r,:])

    def turn2plume_metrics(self,regions=['eb','fsb_upper'],epsilon=1.25,run_len=5,plotting=False,to_plume=True,jumps='jumps'):
        """
            Function decomposes returns into RDP defined linear segments. Gets the point where the turn
            to plume happens and extracts metrics in a .5 s window prior to that event
            Metrics are the phase, bump amplitude, pva 
            
            Strategy:
                1. Use RDP algo to decompose traj into turns
                2. Identify turns towards the plume followed by runs above a threshold length
                3. If animal is still then advance turn point to the next time they move
                4. Extract parameters:
                    - Phases
                    - PVA
                    - Amp
                    - Heading
        
        """
        
        
        cxa = self.cxa
        # Get entries and exits
        if jumps =='jumps':
            jumps = cxa.get_jumps()
        elif jumps=='all':
            jumps = cxa.get_entries_exits_like_jumps()
            
        # Get metrics
        heading = cxa.ft2['ft_heading'].to_numpy()
        x = cxa.ft2['ft_posx'].to_numpy()
        y = cxa.ft2['ft_posy'].to_numpy()
        x,y = ug.fictrac_repair(x, y)
        xy = np.append(x[:,np.newaxis],y[:,np.newaxis],axis=1)
        xy[:,0] = xy[:,0]*-cxa.side
        phases = np.zeros((len(heading),len(regions)))
        pvas = np.zeros_like(phases)
        amps  = np.zeros_like(phases)
        for i,r in enumerate(regions):
            phases[:,i] = cxa.pdat['offset_'+r+'_phase'].to_numpy()*-cxa.side
            pvas[:,i] = ug.get_pvas(cxa.pdat['wedges_'+r])
            amps[:,i] = np.mean(cxa.pdat['wedges_'+r],axis=1)
            
       # pvas = (pvas-np.mean(pvas,axis=0))/np.std(pvas,axis=0) # zscore pvas
        u =ug()
        _,_,dd = u.get_velocity(x,y,cxa.pv2['relative_time'].to_numpy())
        still = ug.get_still(1,np.abs(dd))
        first_data = True
        
        if to_plume:
            rmult = 1 
        else:
            rmult = -1
        
        for i,j in enumerate(jumps):
            dx = np.arange(j[1],j[2])
            
            tp = phases[dx,:]
            tpv = pvas[dx,:]
            ta = amps[dx,:]
            th = heading[dx]
            txy = xy[dx]
            txy = txy-txy[0,:]
            ts = still[dx]
            tb,tbs = ug.find_blocks(ts)
            tbe = tb+tbs-1
            
            rdx = rdp(txy,epsilon,return_mask=True)
            rdxdx = dx[rdx]
            rxy = txy[rdx]
            rxyd = np.diff(rxy,axis=0)
            rdxw = np.where(rdx)[0]
            rdist = np.sqrt(np.sum(rxyd**2,axis=1))
            
            rang = np.angle(rxyd[:,1]+rxyd[:,0]*1j)
            
            if plotting:
                plt.figure()
                plt.subplot(1,2,1)
                plt.plot(txy[:,0],txy[:,1],color='k')
                plt.plot(rxy[:,0],rxy[:,1],color='r')
                plt.scatter(txy[ts>0,0],txy[ts>0,1],color='b')
                plt.gca().set_aspect('equal')
            turncount = -1
            for ri,r in enumerate(rang):
                if plotting:
                    plt.subplot(1,2,1)
                    plt.text(rxy[ri,0],rxy[ri,1],str(np.round(rdist[ri],decimals=1))+'  '+str(np.round(180*r/np.pi)))
                rd = rdist[ri]
                
                # Selection criteria 
                # negative r - animal needs to turn to the plume
                # run length needs to exceed a threshold
                # cannot be the first run - this is special since it is the exit trajectory
                
                if r*rmult<0 and rd>=run_len and ri>0: # minimum run length and number of runs, multiply by to plume or not
                    turncount +=1
                    # refine ri
                    
                    rrxy = rxy[ri,:]
                    tebdiff = tbe-rdxw[ri]
                    tebdiff = tebdiff[tebdiff<=0]

                    if ts[rdxw[ri]]>0: # check if animal is stil
                        instill = True
                        tebdiff = tbe-rdxw[ri] # difference between rdp timepoint and end of still block
                        tbe2 = tbe[tebdiff>=0] # f
                        ri2 = tbe2[0]
                    elif len(tebdiff)==0:
                        instill = False
                    else:
                    
                        last_still = tbe[np.argmax(tebdiff)]
                        dt = np.max(tebdiff)
                        lsxy = txy[last_still]
                        
                        still_d = np.sqrt(np.sum((rrxy-lsxy)**2))
                        if dt>-5 and still_d<1:
                            ri2 = last_still
                            instill= True
        
                        else:
                            instill=False
                    
                    rtime = rdxw[ri]
                    if plotting:
                        plt.plot(rxy[ri:ri+2,0],rxy[ri:ri+2,1],color='g')
                        plt.subplot(1,2,2)
                        plt.plot(ts,color='k')
                        plt.plot(tp[:,1],color='b')
                        plt.plot([rtime,rtime],[-np.pi,np.pi],color='r')
                    
                    if instill:
                        if plotting:
                            plt.plot([ri2,ri2],[-np.pi,np.pi],color='g')
                            #plt.scatter(txy[ri2,0],txy[ri2,1],color='g')
                        
                        dx_r = np.arange(ri2-5,ri2)
                    else:
                        dx_r = np.arange(rtime-5,rtime)
                        
                   # print(i,ri,instill)
                    
                    # Now gather data
                    if first_data:
                        ophase = stats.circmean(tp[dx_r,:],low=-np.pi,high=np.pi,axis=0)[np.newaxis,:]
                        opva = np.mean(tpv[dx_r,:],axis=0)[np.newaxis,:]
                        oamp = np.mean(ta[dx_r,:],axis=0)[np.newaxis,:]
                        oheading = stats.circmean(th[dx_r],low=-np.pi,high=np.pi)
                        index = np.array([i,ri,turncount,rdxdx[ri]])[np.newaxis,:]
                        rdpang = np.array(r)
                        first_data =False
                    else:
                        ophase = np.append(ophase,
                                  stats.circmean(tp[dx_r,:],low=-np.pi,high=np.pi,axis=0)[np.newaxis,:],axis=0)
                        opva = np.append(opva,
                                         np.mean(tpv[dx_r,:],axis=0)[np.newaxis,:],axis=0)
                        oamp = np.append(oamp,
                                        np.mean(ta[dx_r,:],axis=0)[np.newaxis,:],axis=0)
                        oheading = np.append(oheading,stats.circmean(th[dx_r],low=-np.pi,high=np.pi))
                        index = np.append(index,np.array([i,ri,turncount,rdxdx[ri]])[np.newaxis,:],axis=0)
                        rdpang = np.append(rdpang,r)
                        
        if first_data:
            ophase = np.array([])
            opva= np.array([])
            oamp =  np.array([])
            oheading = np.array([])
            index = np.array([])
            rdpang = np.array([])
        output_dictionary = {
            'index': index,
            'phase': ophase,
            'pva': opva,
            'mean_amp': oamp,
            'heading': oheading,
            'rdpang': rdpang
                             }
            
        return output_dictionary
            
    def distance_plume_boundary(self):
        # Iterate through all exits and get max distance from plume boundary
        # Remember to exclude crossings to additional plumes. Take into 
        # account the jumps.
        ft2 = self.cxa.ft2
        ins = ft2['instrip'].to_numpy()
        x = ft2['ft_posx'].to_numpy()
        y = ft2['ft_posy'].to_numpy()
        
        
        
        insdiff = np.diff(ins)
        tt = self.cxa.pv2['relative_time']
        son = np.where(insdiff>0)[0]
        soff = np.where(insdiff<0)[0]
        max_dist = np.array([])
        for i,s in enumerate(soff[:-1]):
            son2 = son[i+1]
            slen = tt[son2]-tt[s]
            if slen<1:
                print('Exit too short')
                continue
            sdx = np.arange(s,son2,dtype='int')
            tx = x[sdx]
            tx_diff = np.abs(tx[-1]-tx[0])
            
            if tx_diff>190:
                print('Got to another plume')
                continue
            
            dtx = tx-tx[0]
            tmx = np.max(np.abs(dtx))
            max_dist = np.append(max_dist,tmx)
            
        max_dist_off = np.array([])
        for i,s in enumerate(son):
            if i==0:
                continue
        
            sf = soff[i]
            slen = tt[sf]-tt[s]
            if slen<1:
                print('In plume too short')
                continue
            
            
            sdx = np.arange(s,sf,dtype='int')
            
            tx = x[sdx]
            tx_diff = np.abs(tx[-1]-tx[0])
            
            dtx = tx-tx[0]
            tmx = np.max(np.abs(dtx))
            max_dist_off = np.append(max_dist_off,tmx)
        
        return max_dist,max_dist_off
    
    
    