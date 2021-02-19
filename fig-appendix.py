#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 10:24:25 2020

"""
import pickle,gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from class_averageStats import averageStats
from class_SIRmodel import SIRmodel
from copy import deepcopy

prefix = 'nc5-'
imagedir = 'output/images/'+prefix+'app-'
outputdir = 'output/'+prefix

benchkwargs = {
    "q_seed" : 2443,
    "p_proby"     : [0.0],
    "p_probr"     : [1/6.5],
    "p_probc"     : [[0, 0.054, 0, 0, 0]], #prob of contagion by type and state
    "p_probd"     : [0],
    "p_infradius" : 0.013016,
    "p_avgstep"   : 0.034,
    "p_stdstep"   : 0,
    "p_cluster"   : 'cluster',
    "q_popsize"   : 25600,
    "q_days"      : 300,
    "q_printOption" : 0.6,
    'g_maxinf'    : 1,
} 
import matplotlib
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]})

figsize = 3.4
#%% Appendix: different number of initial clusters

vals = [1,2,3,5,6,10,20]
values = []
peakactive = []
peakday = []
maxinf = []
ncont = []
lastday = []

for idx,ncl in enumerate(vals):
    try: 
        if ncl==1:
            file=gzip.open(outputdir+'basesim.pickle.gz','rb')
        elif ncl == 20 :
            file=gzip.open(outputdir+'random.pickle.gz','rb')
        else :
            file=gzip.open(outputdir+'allclusters.pickle.gz','rb')
    except:
        print('[',idx,']',ncl,' not found')
        continue
    values = np.append(values,ncl)
    m = pickle.load(file)
    file.close()
    for mod in m:
        if len(mod.q_firstinf) != ncl:
            m.remove(mod)
    

    bavg = averageStats(m)
    peakactive = np.append(peakactive, max(bavg.prinf))
    peakday = np.append(peakday, np.where(bavg.prinf==max(bavg.prinf))[0][0])
    maxinf = np.append(maxinf,bavg.ninstate[4])
    lastday = np.append(lastday, bavg.maxdays)

    print('-----------------------------------------\nSTATS\n')
    print('___'+str(np.round(ncl,2))+'_ density ______')
    print('Pct S I Y R D: '+str(np.round(bavg.ninstate,3)))
    print('Max % I+Y : '+str(max(bavg.prinf)))
    print('Time of max',np.where(bavg.prinf==max(bavg.prinf))[0])
#    print('Average # of contacts: ',ncont[idx])
    
#%%    
    
fig,((ax,bx)) = plt.subplots(1,2,figsize=(figsize*2,figsize))

dens = np.array(values)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylim(0,1.02)

# ax.set_xlim(0.4,4.5)
ax.set_xlabel('Number of initial clusters')
ax.set_ylabel('Percent infected')
ax.plot(dens,maxinf,c='chocolate',label='Total cases')
ax.plot(dens,peakactive,'--',c='olive',label='Peak active cases')
ax.legend(loc='center right')
#ax.set_xticks(np.arange(0,21,4))

bx2 = bx.twinx()
bx2.set_ylim(0,33)
ax.set_xlim(0,21)
bx.set_xlim(0,21)
bx.spines['top'].set_visible(False)
bx2.spines['top'].set_visible(False)
bx.set_ylabel('Days')

bx.set_xlabel('Number of initial clusters')
bx.plot(dens,lastday,c='chocolate',label='Days until steady state')
bx.plot(dens,peakday,'--',c='olive',label='Day of peak (right scale)')
#bx.plot([],[],'--',c='chocolate',label='Number of contacts (right scale)')
bx.legend(loc='center right')

plt.savefig(imagedir+'mclusters.pdf')
plt.show()

#%% Appendix: goodness of fit

file = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file)
file.close()

from class_spatialModels import spatialSAYDR
mod3 = spatialSAYDR()

import pandas as pd
df = pd.DataFrame()
for idx,mod in enumerate(b):
    df[idx] = mod.igrowth[1:]

ma = df.rolling(3).mean()

g_plotdays = 55
malomdr = pd.DataFrame(b[0].mldr[0:g_plotdays]) #this is already a moving avg
malomdr3 = pd.DataFrame(mod3.mldr[0:g_plotdays]) #this is already a moving avg
#malomdr = malomdr.rolling(4,win_type='triang').mean()
malomdr = np.array(malomdr[:])
meanBoots = ma[ma.columns[0:1]].mean(axis=1)
  
fig,ax = plt.subplots(figsize=(figsize*4/3,figsize))
minday = 0

### draw legend with empty data
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks(np.arange(minday,g_plotdays,10))
#ax.set_ylim(0.0,0.44)
#ax.set_yticks(np.arange(0,0.8*m[0].q_popsize,10000))
#ax.set_ylim(4,g_maxinf)
ax.set_xlabel('Days ')
ax.set_ylabel('Growth rate')

# draw the actual chart. Rescale position so that it fits the square properly
ax.plot(np.arange(g_plotdays),meanBoots[:g_plotdays],c='DarkOrange',label='Infections growth, Baseline')
#ax.plot(np.arange(1,g_plotdays),malomdr[1:g_plotdays][:,0],c='olive',label='Fatalities growth, Lombardy')
ax.plot(np.arange(1,g_plotdays),malomdr3[1:g_plotdays][0],c='olive',label='Fatalities growth, Lombardy')
ax.legend()
plt.savefig(imagedir+'baseline-lombardy-match.pdf')  


#%% Appendix: Rescaling a city

from class_spatialModels import spatialSAYDR

# baseline model
kwargs = deepcopy(benchkwargs)
b = spatialSAYDR(**kwargs)

# model with twice the city side (4* the population)
sizemult = 2
kwargs2 = deepcopy(benchkwargs)
kwargs2['q_popsize'] = np.int(kwargs['q_popsize'] * sizemult**2)
kwargs2['q_citysize'] = sizemult
fromedge = sizemult-0.25
kwargs2['q_firstinf'] = [[0.25,0.25],[0.25,fromedge],[fromedge,0.25],[fromedge,fromedge]]
d = spatialSAYDR(**kwargs2)

printdays = [10, 20, 30, 40, 50, 70]

# Create 2x2 sub plots
def genfig(day=0):
    fig = plt.figure(figsize=(figsize*2,figsize))
    axd = fig.add_subplot(1,2,2)
    axb = fig.add_subplot(2,4,6)
    axr = fig.add_subplot(2,4,2)
    axb.set_xlim(0,1)
    axb.set_ylim(0,1)
    axd.set_ylim(0,sizemult)
    axd.set_xlim(0,sizemult)
    axr.set_xlim(0,120)
    axr.set_ylim(0,1.5)
    
    axb.set_yticks([])
    axd.set_yticks([])
    axb.set_xticks([])
    axd.set_xticks([])
    axr.set_ylabel('% Population')
    axr.set_yticks([0,0.5,1])
    axr.plot([],[],c='SaddleBrown',label='small city')
    axr.plot([],[],c='olive',label='big city')
    axr.plot([],[],'--',c='black',label='active')
    axd.scatter([],[],c='darkorange',label='Infected')
    axd.scatter([],[],c='olive',label='Recovered')
    axr.legend(loc='upper left',ncol=1,prop={'size': 3},title='Day'+str(day))
    axd.legend(loc='upper left',ncol=1)
    axr.spines['right'].set_visible(False)
    axr.spines['top'].set_visible(False)
    plt.rcParams['legend.title_fontsize'] = 'small'
    return fig,axb,axd,axr

fig = genfig()
from celluloid import Camera

fig,axb,axd,axr = genfig('')
camera = Camera(fig)



# start up the simulation
b.computeStats(0)
d.computeStats(0)
for day in range(85):
    
    #uncomment to save individual figures
    #comment to save movie
    fig,axb,axd,axr = genfig('')
    
    d.computeStats(day,1)
    b.computeStats(day,1)
    if np.mod(day,1)==0 :
        axr.legend(loc='upper left',ncol=1,prop={'size': 8},title='Day'+str(day))
#        axb.scatter(b.pos[:,0][b.state==b.S],b.pos[:,1][b.state==b.S],c='lightgrey',s=.2)
        axb.scatter(b.pos[:,0][b.state==b.I],b.pos[:,1][b.state==b.I],c='darkorange',s=1)
        axb.scatter(b.pos[:,0][b.state==b.R],b.pos[:,1][b.state==b.R],c='olive',s=.2)
#        axd.scatter(d.pos[:,0][d.state==d.S],d.pos[:,1][d.state==d.S],c='lightgrey',s=.2)
        axd.scatter(d.pos[:,0][d.state==d.I],d.pos[:,1][d.state==d.I],c='darkorange',s=1)
        axd.scatter(d.pos[:,0][d.state==d.R],d.pos[:,1][d.state==d.R],c='olive',s=.2)
        axr.plot(np.arange(day),b.prtinf[0:day],c='SaddleBrown')
        axr.plot(np.arange(day),d.prtinf[0:day],c='olive')
        axr.plot(np.arange(day),b.prinf[0:day],'--',c='SaddleBrown')
        axr.plot(np.arange(day),d.prinf[0:day],'--',c='olive')
        
        if day in printdays:
            plt.savefig(imagedir+'sizecomp-day'+str(day)+'.pdf')
            plt.show()
        camera.snap()
    if np.sum(d.state==d.I)+np.sum(d.state==d.Y) <= 15:
        break

    b.aDayInTheLife(day) 
    d.aDayInTheLife(day) 
      
anim = camera.animate(blit=True, interval=150, repeat_delay=500)
anim.save(imagedir+'animation-sizecomp.gif')

#%% Appendix: city size

from class_averageStats import avgDist,averageStats
file = gzip.open(outputdir+'allsizes.pickle.gz','rb')
allmodels = pickle.load(file)
file.close()

file = gzip.open(outputdir+'basesize2.pickle.gz', 'rb')
mod2 = pickle.load(file)
file.close()

for model in mod2:
    allmodels.append(model)
    
#values = np.append(np.arange(10000,40000,5000),np.arange(40000,150000,10000))
values = []
for model in allmodels:
    if model.q_popsize not in values:
        values.append(model.q_popsize)
        

peakactive = np.zeros(len(values))
peakday = np.zeros(len(values))
maxinf = np.zeros(len(values))
ncont = np.zeros(len(values))
lastday = np.zeros(len(values))

    
for idx,size in enumerate(values):
    print(size)
    modsize = list()
    for model in allmodels:
        if model.q_popsize == size:
            modsize.append(model)
    
    
    bavg = averageStats(modsize)
    peakactive[idx] = bavg.maxCases
    peakday[idx] = max(np.where(bavg.prinf==max(bavg.prinf)))[0]
    maxinf[idx] = max(bavg.prtinf)
    lastday[idx] =bavg.minlastday
    
    # I didn't save positions in the big pickle so need to recreated them
    # to compute # of contacts
    pos = np.random.rand(modsize[0].q_popsize,2)*modsize[0].q_citysize
    ravg = avgDist(pos=pos,radius=modsize[0].p_infradius,maxc=2000)[1]
    ncont[idx] = ravg
    
  
fig,((ax,bx)) = plt.subplots(1,2,figsize=(figsize*2,figsize))

dens = np.dot(values,1/25600)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylim(0,1.2)
ax.set_xlim(0.3,4.2)

# ax.set_xlim(0.4,4.5)
ax.set_xlabel('City size relative to baseline ')
ax.set_ylabel('Percent infected')
ax.set_xticks(np.arange(0.5,4.2,0.5))
ax.plot(dens,peakactive,c='olive',label='Peak active cases')
ax.plot(dens,maxinf,'--',c='saddlebrown',label='Total cases')
ax.legend(loc='center left')

bx2 = bx.twinx()
bx2.set_ylim(0,30)
bx.set_xticks(np.arange(0.5,4.2,0.5))
bx.set_xlim(0.3,4.2)
bx.spines['top'].set_visible(False)
bx2.spines['top'].set_visible(False)
bx.set_xlabel('City size relative to baseline ')
bx.plot(dens,lastday,c='saddlebrown',label='Days until steady state')
bx.plot([],[],':',linewidth=2,c='darkorange',label='\# of contacts (right scale)')
bx.plot(dens,peakday,'--',c='olive',label='Day of peak')
bx2.plot(dens,ncont,':',c='darkorange',linewidth=2 )
bx.legend(loc='upper left')

plt.savefig(imagedir+'sizes.pdf')

#%% Densities

from class_averageStats import avgDist
from class_spatialModels import spatialSAYDR

file1 = gzip.open(outputdir+'app-density.pickle.gz','rb')
allmods=  pickle.load(file1)
file1.close()

values = []
for model in allmods:
    if model.q_citysize not in values:
        values.append(model.q_citysize)
values=np.array(values)        

nv = len(values)
peakactive = np.zeros(nv)
peakday = np.zeros(nv)
maxinf = np.zeros(nv)
ncont = np.zeros(nv)
lastday = np.zeros(nv)

 
for idx,citysize in enumerate(values) :
    if citysize > 1.7:
        continue
    print(idx,citysize)
    modh = list()
    for mod in allmods:
        if mod.q_citysize == citysize :
            modh.append(mod)
    values[idx] = citysize
    modav = averageStats(modh)
    peakactive[idx] = modav.maxCases
    peakday[idx] = np.max(np.where(modav.prinf==np.max(modav.prinf)))
    maxinf[idx] = np.max(modav.prtinf)
    lastday[idx] = modav.minlastday
    
    kwargs = benchkwargs
    kwargs['q_citysize'] = citysize
    mm = spatialSAYDR(**kwargs)
    ncont[idx] = avgDist(radius=mm.p_infradius,pos=mm.pos,maxc=2000)[1]
    

     #%%
dens = 1/(values**2)
    
fig,((ax,bx)) = plt.subplots(1,2,figsize=(figsize*2,figsize))

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)


#ax.set_xlim(0.2,3.5)
ax.set_xlabel('City density relative to baseline ')
ax.set_ylabel('Percent infected')
ax.plot(dens,maxinf,c='saddlebrown',label='Total cases')
ax.plot(dens,peakactive,'--',c='darkorange',label='Peak active cases')
ax.legend(loc='lower right')


bx2 = bx.twinx()
bx2.set_ylim(0,53)
#bx2.set_xlim(0.2,3.5)

bx.spines['top'].set_visible(False)
bx2.spines['top'].set_visible(False)
bx.set_xlabel('City density relative to baseline ')
bx.plot(dens[:-3],peakday[:-3],'--',c='olive',label='Day of peak')
bx.plot(dens[:-3],lastday[:-3],c='darkorange',label='Days until steady state')
bx.plot([],[],':',c='saddlebrown', label='\# of contacts (right scale)')
bx2.plot(dens[:-3],ncont[:-3], ':', linewidth = 1.9, c='saddlebrown' )
bx.legend(loc='upper right')


fig.tight_layout()
plt.savefig(imagedir+'densities.pdf')

