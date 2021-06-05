#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 19:36:26 2020

Code in this file replicates all figures in

Alberto Bisin and Andrea Moro,
"Learning Epidemiology by Doing: The Empirical Implications of
a Spatial SIR Model with Behavioral Responses,"
NBER Working Paper 27590, June 2020
"""

import pickle,gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from class_averageStats import averageStats
from class_SIRmodel import SIRmodel
from copy import deepcopy
import matplotlib
import matplotlib.gridspec as gridspec

# import labellines
# this can be commented out if your system does not have or cannot
# read a LaTeX distribution.
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]})

prefix= 'nc5-'
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
    "p_ininf"     : 30,
    "q_popsize"   : 25600,
    "q_days"      : 300,
    "q_printOption" : 0.6,
    'g_maxinf'    : 1,
}

imagedir = 'output/images/'+prefix+'b'

cluster = benchkwargs['p_cluster']
popSize = benchkwargs['q_popsize']
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]

imagedir = 'output/images/'+prefix

fsize = 3.3

#%% Spatial progression of infections, benchmarkmodel

from class_spatialModels import spatialSAYDR

printdays = [3, 10, 20, 30]
kwargs = deepcopy(benchkwargs)
kwargs['q_days'] = np.max(printdays) + 1

m = spatialSAYDR(**kwargs)

for day in range(1,m.q_days):

    print(m.computeStats(day,1))
    if day in printdays:
        m.drawPositions(day,savefig='baseline',S=True,folder=imagedir,ext='.png')
        pass
    m.aDayInTheLife(day)

    if np.sum(m.state==m.I)+np.sum(m.state==m.Y) <=15:
        break

#%% Spatial progression of infections, random initial locations

printdays = [3, 10, 20, 30]
kwargs = deepcopy(benchkwargs)
kwargs['p_cluster']   = 'random'
kwargs['q_days'] = np.max(printdays) + 1

m = spatialSAYDR(**kwargs)

for day in range(1,m.q_days):

    print(m.computeStats(day,1))
    if day in printdays:
        m.drawPositions(day,savefig='randomcluster',S=True,folder=imagedir,ext='.png')
        pass
    m.aDayInTheLife(day)

    if np.sum(m.state==m.I)+np.sum(m.state==m.Y) <=15:
        break

#%% Spatial progression of infections, no movements

printdays = [10,20,30,50,150,250]

kwargs = deepcopy(benchkwargs)
kwargs['p_avgstep'] = 0
kwargs['q_days'] = np.max(printdays) + 1

m = spatialSAYDR(**kwargs)

for day in range(1,m.q_days):

    print(m.computeStats(day,1))
    if day in printdays:
        m.drawPositions(day,savefig='nomove',S=True,folder=imagedir,ext='.png')
        pass
    m.aDayInTheLife(day)

    if np.sum(m.state==m.I)+np.sum(m.state==m.Y) <=15:
        break

#%% Initial location in space, heterogeneous density model

from class_spatialModels import spSAYDR_hetDensity

mod1 = spSAYDR_hetDensity(q_lambda=1,**benchkwargs)
mod1.drawPositions(0,{'savefig':'hetdens','folder':imagedir,'S': True,
                      'colors':['grey','orange','orange','red','green'], 'ext':'.png'}, )


#%% Statistics reported in paper, section Local Herd Immunity ( \label{sec:SIRcomp})

# Benchmark model
file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
bavg = averageStats(b)
file1.close()

print('R0, average of '+str(len(b))+' replications: ',np.round(bavg.R0,2))
print('R0, std of '+str(len(b))+' replications: ',np.round(bavg.R0std,2))

avgmaxinf = np.average(list(map(lambda x: np.max(x.prtinf),b)))
stdmaxinf = np.std(list(map(lambda x: np.max(x.prtinf),b)))

print('Avg of steady state infected ', np.round(avgmaxinf,2))
print('Std of steady state infected ', np.round(stdmaxinf,3))

#%% SIR and Spatial-SIR comparison of infection dynamics

plt.rcParams['text.usetex'] = True #Let TeX do the typsetting
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath'] #Force sans-serif math mode (for axes labels)
plt.rcParams['font.family'] = 'serif' # ... for regular text
plt.rcParams['font.sans-serif'] = 'Computer Modern Roman' # Choose a nice font here


# Benchmark model
file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
bavg = averageStats(b)
file1.close()

dfb = pd.DataFrame()
logsb = pd.DataFrame()
for idx,m in enumerate(b) :
    dfb[idx] = np.sum(m.nstates[:,1:3],axis=1)
    dfb[str(idx)+'L1'] = dfb[idx].shift()
    logsb[str(idx)] = np.log(dfb[idx]) - np.log(dfb[str(idx)+'L1'])

# Random locations
file = gzip.open(outputdir+'allrandmodels.pickle.gz','rb')
rand = pickle.load(file)
file.close()
ravg = averageStats(rand)

df = pd.DataFrame()
logs = pd.DataFrame()
for idx,m in enumerate(rand) :
    df[idx] = np.sum(m.nstates[:,1:3],axis=1)
    df[str(idx)+'L1'] = df[idx].shift()
    logs[str(idx)] = np.log(df[idx]) - np.log(df[str(idx)+'L1'])

# SIR model
aa = SIRmodel(contRate*13.5,recRate,q_popsize=popSize,q_init=10)
dfSIR = pd.DataFrame(aa.SIR,index=np.arange(aa.day+1),columns=['S','I','R','L'])
dfSIR['ILag'] = dfSIR['I'].shift()
dfSIR['ld'] = np.log(dfSIR['I']) - np.log(dfSIR['ILag'])
dfSIRld = np.array(dfSIR['ld'])

days=aa.SIR[:,0].size

fig = plt.figure(figsize=(2*fsize,fsize*1))
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)

ax = fig.add_subplot(spec[0,0])
ax2 = fig.add_subplot(spec[0,1])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Days ')
ax.set_ylim(-0.2,0.5)
ax.set_xlim(-1,95)
l1=ax.plot(np.arange(len(dfSIRld)),dfSIRld,':', linewidth=2, color='saddlebrown',label ='SIR')[0]
l2=ax.plot(np.arange(days),logs.mean(axis=1)[:days],'--', color='darkorange',
        label ='\parbox{10em}{Spatial-SIR\\newline with random positions}')[0]
l3=ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Spatial-SIR')[0]

days=aa.SIR[:,0].size
sstateSIR = str(np.round(1-aa.SIR[days-1,0]/25600,2))
sstater = str(np.round(np.max(ravg.prtinf),2))
sstateb = str(np.round(np.max(bavg.prtinf),2))

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Days ')
ax2.set_xlim(-1,92)
ax2.set_ylim(0,0.52)
ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,':', linewidth=2, color='saddlebrown',label = sstateSIR)
ax2.plot(np.arange(days),ravg.prinf[:days],'--', color='darkorange', label =sstater)
ax2.plot(np.arange(days),bavg.prinf[:days],color='olive',label = sstateb)

ax2.legend(title='$\\frac{R_\\infty}{N}$', bbox_to_anchor=(.52,.3), title_fontsize=15)

line_labels=['SIR',
             '\parbox{10em}{Spatial-SIR\\newline with random positions}',
             'Spatial-SIR']

fig.legend(bbox_to_anchor=(.62,.82), labels=line_labels, fontsize=12,framealpha=1)
ax.set_title('Growth rate')
ax2.set_title('Infected')

fig.tight_layout()
plt.savefig(imagedir+'density_contagion2.pdf')

plt.close()
#%% Random locations

# Benchmark model
file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
bavg = averageStats(b)
file1.close()

dfb = pd.DataFrame()
logsb = pd.DataFrame()
for idx,m in enumerate(b) :
    dfb[idx] = np.sum(m.nstates[:,1:3],axis=1)
    dfb[str(idx)+'L1'] = dfb[idx].shift()
    logsb[str(idx)] = np.log(dfb[idx]) - np.log(dfb[str(idx)+'L1'])

# Random locations
file = gzip.open(outputdir+'random.pickle.gz','rb')
rand = pickle.load(file)
file.close()
ravg = averageStats(rand)

df = pd.DataFrame()
logs = pd.DataFrame()
for idx,m in enumerate(rand) :
    df[idx] = np.sum(m.nstates[:,1:3],axis=1)
    df[str(idx)+'L1'] = df[idx].shift()
    logs[str(idx)] = np.log(df[idx]) - np.log(df[str(idx)+'L1'])


fig = plt.figure(figsize=(2*fsize,fsize*1))
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)

ax = fig.add_subplot(spec[0,0])
ax2 = fig.add_subplot(spec[0,1])

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Days ')
ax.set_ylim(-0.2,0.5)
ax.set_xlim(-1,95)
l0= ax.plot(np.arange(days),logs.mean(axis=1)[:days],'--', color='darkorange',
        label ='\parbox{10em}{Spatial-SIR \newline with random outbreaks}')[0]
l1=ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Spatial-SIR')[0]

days=aa.SIR[:,0].size
sstateSIR = str(np.round(1-aa.SIR[days-1,0]/25600,2))
sstater = str(np.round(np.max(ravg.prtinf),2))
sstateb = str(np.round(np.max(bavg.prtinf),2))

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Days ')
ax2.set_xlim(-1,92)
ax2.set_ylim(0,0.3)
ax2.set_yticks(np.arange(0,0.31,.1))
ax2.plot(np.arange(days),ravg.prinf[:days],'--', color='darkorange', label =sstater)
ax2.plot(np.arange(days),bavg.prinf[:days],color='olive',label = sstateb)

ax2.legend(title='$\\frac{R_\\infty}{N}$', bbox_to_anchor=(.52,.42), title_fontsize=15)

line_labels=['Random outbreaks',
             'Baseline Spatial-SIR']

fig.legend([l0,l1],line_labels, bbox_to_anchor=(.62,.88), fontsize= 12, borderaxespad=0.5,framealpha=1)
ax.set_title('Growth rate')
ax2.set_title('Infected')


fig.tight_layout()

plt.savefig(imagedir+'short-random-rates.pdf')
plt.close()

print('Peak active, baseline',np.round(max(bavg.prinf),2)
      ,'day ',np.where(max(bavg.prinf)==bavg.prinf)[0])
print('Peak active, random',np.round(max(ravg.prinf),2)
      ,'day ',np.where(max(ravg.prinf)==ravg.prinf)[0])
avgmaxinf = np.average(list(map(lambda x: np.max(x.prtinf),b)))
print('Steady state infected baseline', np.round(avgmaxinf,4))
avgmaxinf2 = np.average(list(map(lambda x: np.max(x.prtinf),rand)))
print('Steady state infected random', np.round(avgmaxinf2,4))

#%%  SIR model city size

file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
file2 = gzip.open(outputdir+'basesize2.pickle.gz','rb')
fileq = gzip.open(outputdir+'quartersim.pickle.gz','rb')

b = pickle.load(file1)
s_bavg = averageStats(b)
file1.close()

b2 = pickle.load(file2)
s_b2avg = averageStats(b2)
file2.close()

q = pickle.load(fileq)
s_qavg = averageStats(q)
fileq.close()

import pandas as pd

popSize = benchkwargs['q_popsize']
cluster = 30
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]
bavg = SIRmodel(contRate*13.5, recRate, q_popsize=popSize, q_init=cluster)
b2avg = SIRmodel(contRate*13.5, recRate, q_popsize=popSize*4, q_init=cluster)
sqavg = SIRmodel(contRate*13.5, recRate, q_popsize=popSize/4, q_init=cluster)

maxday=bavg.day
day2= b2avg.day
dayq= sqavg.day

fig,(ax12,ax2) = plt.subplots(1,2,figsize=(2*fsize,fsize))
g_plotdays=150
g_maxinf=0.5

ax12.spines['right'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax12.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax12.set_ylim(0,g_maxinf)
ax2.set_ylim(0,g_maxinf)
ax2.set_yticks(np.arange(0,g_maxinf+.01,.1))
ax12.set_yticks(np.arange(0,g_maxinf+.01,.1))


ax12.set_xlabel('Days ')
ax2.set_xlabel('Days ')

s_qavg.prtinf[s_qavg.minlastday:] = s_qavg.prtinf[s_qavg.minlastday]
s_bavg.prtinf[s_bavg.minlastday:] = s_bavg.prtinf[s_bavg.minlastday]

#spatial SIR

m2 = str(np.round(np.max(s_b2avg.prtinf),2))
ma = str(np.round(np.max(s_bavg.prtinf),2))
mq = str(np.round(np.max(s_qavg.prtinf),2))
ll3 = ax12.plot(np.arange(g_plotdays),s_qavg.prinf[0:g_plotdays],'--',c='darkorange',label=mq)[0]
ll2 = ax12.plot(np.arange(g_plotdays),s_bavg.prinf[0:g_plotdays],c='olive',label=ma)[0]
ll1 = ax12.plot(np.arange(g_plotdays),s_b2avg.prinf[0:g_plotdays],':',c='saddlebrown',linewidth=1.9, label=m2)[0]

ax12.legend(title='${R_\\infty}/{N}$', bbox_to_anchor=(.23,.25), title_fontsize=11)
#ax12.set_title(title='Infected')
g_plotdays = 70

ax2.set_xticks(np.arange(0,201,20))



ms2 = str(np.round(np.max(b2avg.prtinf),2))
msa = str(np.round(np.max(bavg.prtinf),2))
msq = str(np.round(np.max(sqavg.prtinf),2))
ax2.plot(np.arange(g_plotdays),sqavg.prinf[0:g_plotdays],'--',c='darkorange',label=msq)
ax2.plot(np.arange(g_plotdays),bavg.prinf[0:g_plotdays],c='olive',label=msa)
ax2.plot(np.arange(g_plotdays),b2avg.prinf[0:g_plotdays],':',c='saddlebrown',linewidth=1.9,label= ms2)
ax2.legend(title='${R_\\infty}/{N}$', bbox_to_anchor=(.48,.22), title_fontsize=11)
line_labels=[             '$N=1/4*$Baseline',
             'Baseline',
    '$N=4*$Baseline',]

fig.legend([ll3,ll2,ll1],line_labels, bbox_to_anchor=(.61,.85), fontsize= 12, borderaxespad=0.5,framealpha=1)
ax12.set_title('Spatial-SIR')
ax2.set_title('SIR')

fig.tight_layout()
plt.savefig(imagedir+'SIR-citysize-rates.pdf')
plt.close()

print('Peak active, baseline',np.round(max(s_bavg.prinf),2)
      ,'day ',np.where(max(s_bavg.prinf)==s_bavg.prinf)[0])
print('Peak active, 2*',np.round(max(s_b2avg.prinf),2)
      ,'day ',np.where(max(s_b2avg.prinf)==s_b2avg.prinf)[0])
print('Peak active, 1/4',np.round(max(s_qavg.prinf),2)
      ,'day ',np.where(max(s_qavg.prinf)==s_qavg.prinf)[0])
avgmaxinf = np.average(list(map(lambda x: np.max(x.prtinf),b)))
print('Steady state infected baseline', np.round(avgmaxinf,4))
avgmaxinf2 = np.average(list(map(lambda x: np.max(x.prtinf),b2)))
print('Steady state infected 2*', np.round(avgmaxinf2,4))
avgmaxinf3 = np.average(list(map(lambda x: np.max(x.prtinf),q)))
print('Steady state infected 1.4', np.round(avgmaxinf3,4))

print('Sir DAYS to peak',np.where(max(bavg.prinf)==bavg.prinf)[0]
      ,np.where(max(b2avg.prinf)==b2avg.prinf)[0],
      np.where(max(sqavg.prinf)==sqavg.prinf)[0])

#%% Baseline vs 6*contagious, 1/6*density


file = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file)
file.close()

file = gzip.open(outputdir+'benchmark_6x_cont.pickle.gz','rb')
benchplushuge  = pickle.load(file)
file.close()

# benchmark plus (Y contagious)
dfb = pd.DataFrame()
logsb = pd.DataFrame()
for idx,m in enumerate(b) :
    dfb[idx] = np.sum(m.nstates[:,1:3],axis=1)
    dfb[str(idx)+'L1'] = dfb[idx].shift()
    logsb[str(idx)] = np.log(dfb[idx]) - np.log(dfb[str(idx)+'L1'])

# benchmark plus huge
df6 = pd.DataFrame()
logs6 = pd.DataFrame()
for idx,m in enumerate(benchplushuge) :
    df6[idx] = np.sum(m.nstates[:,1:3],axis=1)
    df6[str(idx)+'L1'] = df6[idx].shift()
    logs6[str(idx)] = np.log(df6[idx]) - np.log(df6[str(idx)+'L1'])

bavg = averageStats(b)
bavg.prtinf[bavg.minlastday:] = bavg.prtinf[bavg.minlastday]
b6 = averageStats(benchplushuge)
b6.prtinf[b6.minlastday:] = b6.prtinf[b6.minlastday]

# first plot: add bench plus with less density
days = 213

fig,(ax,ax2) = plt.subplots(1,2,figsize=(2*fsize,fsize))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

ax.set_xlabel('Days')
ax2.set_xlabel('Days')
ll2=ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Baseline')[0]
ll1=ax.plot(np.arange(days),logs6.mean(axis=1)[:days],'--', color='darkorange',label ='1/6*density, 6*contagious')[0]

days = 213

qqq= str(np.round(np.max(b[0].prtinf),2))
ppp= str(np.round(np.max(benchplushuge[0].prtinf),2))
ax.set_ylim(-0.2,0.4)
ax2.set_ylim(0,0.17)
ax2.plot(np.arange(days),bavg.prinf[:days],color='olive',label =qqq)
ax2.plot(np.arange(days),b6.prinf[:days],'--', color='darkorange',label =ppp)

ax2.set_yticks(np.arange(0,0.17,0.04))

plt.setp(ax.get_title(), multialignment='center')

ax.set_title('Growth rate')
ax2.set_title('Infected')
ax2.legend(title='$\\frac{R_\\infty}{N}$', bbox_to_anchor=(.58,.35), title_fontsize=15)
line_labels=['Baseline Spatial-SIR',
    '\\parbox{6.5em}{1/6*density,\\newline 6*contagious}',]

fig.legend([ll2,ll1],line_labels, bbox_to_anchor=(.63,.84), fontsize= 12, borderaxespad=0.5,framealpha=1)


fig.tight_layout()
plt.savefig(imagedir+'short-density_contagion1.pdf')
plt.close()

#%% Different city density (different size same population)

file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
s_bavg = averageStats(b)
file1.close()

file1 = gzip.open(outputdir+'density.pickle.gz','rb')
allmods=  pickle.load(file1)
file1.close()

modh = list()
mod2 = list()
for mod in allmods:
    if np.round(mod.q_citysize,1) == 0.7 :
        mod2.append(mod)
    if np.round(mod.q_citysize,1) == 1.4 :
        modh.append(mod)

s_b2avg = averageStats(mod2)
s_qavg = averageStats(modh)

popSize = benchkwargs['q_popsize']
cluster = 30
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]
bavg = SIRmodel(contRate*13.5, recRate, q_popsize=popSize, q_init=cluster)
b2avg = SIRmodel(contRate*13.5*2, recRate, q_popsize=popSize*4, q_init=cluster)
sqavg = SIRmodel(contRate*13.5/2, recRate, q_popsize=popSize/4, q_init=cluster)


print('Sir DAYS to peak',np.where(max(bavg.prinf)==bavg.prinf)[0]
      ,np.where(max(b2avg.prinf)==b2avg.prinf)[0],
      np.where(max(sqavg.prinf)==sqavg.prinf)[0])
print('Peak active, 1*',np.round(max(bavg.prinf),2)
      ,'day ',np.where(max(bavg.prinf)==bavg.prinf)[0])
print('Peak active, 2*',np.round(max(b2avg.prinf),2)
      ,'day ',np.where(max(b2avg.prinf)==b2avg.prinf)[0])
print('Peak active, 1/2*',np.round(max(sqavg.prinf),2)
      ,'day ',np.where(max(sqavg.prinf)==sqavg.prinf)[0])

avgmaxinf = np.round(np.max(b[0].prtinf),2)
avgmaxinf2 = np.round(np.max(mod2[0].prtinf),2)
avgmaxinf3 = np.round(np.max(modh[0].prtinf),2)
print('Steady state infected 1*', np.round(avgmaxinf,4))
print('Steady state infected 1/2*', np.round(avgmaxinf3,4))
print('Steady state infected 2*', np.round(avgmaxinf2,4))

sirmaxa = np.round(np.max(bavg.SIR[:,2])/(25600),2)
sirmax2 = np.round(np.max(b2avg.SIR[:,2])/(25600*4),2)
sirmaxh = np.round(np.max(sqavg.SIR[:,2])/(25600/4),2)

print('Sir DAYS to peak',np.where(max(bavg.prinf)==bavg.prinf)[0]
      ,np.where(max(b2avg.prinf)==b2avg.prinf)[0],
      np.where(max(sqavg.prinf)==sqavg.prinf)[0])

print('Peak active, 1*  ',np.round(max(bavg.prinf),2)
      ,'day ',np.where(max(bavg.prinf)==bavg.prinf)[0])
print('Peak active, 2*  ',np.round(max(b2avg.prinf),2)
      ,'day ',np.where(max(b2avg.prinf)==b2avg.prinf)[0])
print('Peak active, 1/2*',np.round(max(sqavg.prinf),2)
      ,'day ',np.where(max(sqavg.prinf)==sqavg.prinf)[0])


print('Steady state infected 1*  ', np.round(sirmaxa,4))
print('Steady state infected 2*  ', np.round(sirmax2,4))
print('Steady state infected 1/2*', np.round(sirmaxh,4))


maxday=bavg.day
day2= b2avg.day
dayq= sqavg.day

fig,(ax12,ax2) = plt.subplots(1,2,figsize=(2*fsize,fsize))
g_plotdays=150
g_maxinf=0.81

ax12.spines['right'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax12.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax12.set_ylim(0,g_maxinf)
ax2.set_ylim(0,g_maxinf)
ax2.set_yticks(np.arange(0,g_maxinf+.01,.2))
ax12.set_yticks(np.arange(0,g_maxinf+.01,.2))


ax12.set_xlabel('Days ')
ax2.set_xlabel('Days ')

s_qavg.prtinf[s_qavg.minlastday:] = s_qavg.prtinf[s_qavg.minlastday]
s_bavg.prtinf[s_bavg.minlastday:] = s_bavg.prtinf[s_bavg.minlastday]

#spatial SIR

ll3 = ax12.plot(np.arange(g_plotdays),s_b2avg.prinf[0:g_plotdays],'--',c='darkorange',label=avgmaxinf2)[0]
ll2 = ax12.plot(np.arange(g_plotdays),s_bavg.prinf[0:g_plotdays],c='olive',label=avgmaxinf)[0]
ll1 = ax12.plot(np.arange(g_plotdays),s_qavg.prinf[0:g_plotdays],':',c='saddlebrown',linewidth=1.9, label=avgmaxinf3)[0]

#ax12.set_title(title='Infected')
g_plotdays = 70

ax2.set_xticks(np.arange(0,201,20))
ax2.plot(np.arange(g_plotdays),b2avg.prinf[0:g_plotdays],'--',c='darkorange',label=sirmax2)
ax2.plot(np.arange(g_plotdays),bavg.prinf[0:g_plotdays],c='olive',label=sirmaxa)
ax2.plot(np.arange(g_plotdays),sqavg.prinf[0:g_plotdays],':',c='saddlebrown',linewidth=1.9,label= sirmaxh)


ax12.legend(title='${R_\\infty}/{N}$', bbox_to_anchor=(.23,.18), title_fontsize=11)
ax2.legend(title='${R_\\infty}/{N}$', bbox_to_anchor=(.48,.22), title_fontsize=11)
line_labels=['$2*$ density',
             'Baseline',
    '$1/2*$ density',]

fig.legend([ll3,ll2,ll1],line_labels, bbox_to_anchor=(.6,.88), fontsize= 12, borderaxespad=0.5,framealpha=1)
ax12.set_title('Spatial-SIR')
ax2.set_title('SIR')

fig.tight_layout()
plt.savefig(imagedir+'short-3densities.pdf')
plt.close()

#%% Heterogeneous density

file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
file1.close()
bavg = averageStats(b)

fileq = outputdir+'hetdens-lambda-1.pickle.gz'
file = gzip.open(fileq,'rb')
mod1 = pickle.load(file)
file.close()
b1avg = averageStats(mod1)

fig,ax2 = plt.subplots(figsize=(fsize,fsize))
ax2.set_ylim(0,0.23)
ax2.set_yticks(np.arange(0,0.21,0.05))
g_plotdays = b1avg.minlastday
bavg.prtinf[bavg.minlastday:] = bavg.prtinf[bavg.minlastday]

### draw legend with empty data
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Days ')
ax2.set_title('Infected')

msa = '$\\frac{R_\\infty}{N}=$'+str(np.round(np.max(bavg.prtinf),2))
ms1 = '$\\frac{R_\\infty}{N}=$'+str(np.round(np.max(b1avg.prtinf),2))

ax2.plot(np.arange(g_plotdays),bavg.prinf[0:g_plotdays],'',c='olive',linewidth=1.3, label='Baseline Spatial-SIR, '+msa)
ax2.plot(np.arange(g_plotdays),b1avg.prinf[0:g_plotdays],'--',c='darkorange',linewidth=1.7, label="Heterogeneous density, "+ms1)
ax2.legend(loc='upper right')

fig.tight_layout()
plt.savefig(imagedir+'hetdens1.pdf')
plt.close()

#%% Different movement speed

# Baseline model
file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
bavg = averageStats(b)
file1.close()

file = gzip.open(outputdir+'allmovements.pickle.gz','rb')
models = pickle.load(file)
file.close()

df = pd.DataFrame()

df['speed'] =list(map(lambda x: x.p_avgstep, models))
df['seed'] = list(map(lambda x: x.q_seed,models))
df['lastday'] = list(map(lambda x: x.lastday,models))
df['maxinf'] = list(map(lambda x: max(x.prtinf),models))
df['peakinf'] = list(map(lambda x: max(x.prinf),models))

days = 601
prtinf = pd.DataFrame(list(map(lambda x: x.prtinf[0:days],models)))
prinf = pd.DataFrame(list(map(lambda x: x.prinf[0:days],models)))

df1 = pd.concat([df,prtinf],axis=1)
df2 = pd.concat([df,prinf],axis=1)

dfav = df1.groupby('speed').mean()
df2av = df2.groupby('speed').mean()

dfn = dfav.to_numpy()[:,5:]
df2n = df2av.to_numpy()[:,5:]

b0 = df2n[0,:] #0 speed averages
b20 = df2n[1,:] #20% speed

fig,(ax2) = plt.subplots(1,1,figsize=(fsize,fsize))
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# ax2.spines['right'].set_visible(False)

# ax.set_xlabel('Days')
# ax2.set_xlabel('Days')
# ll2=ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Baseline Spatial-SIR')[0]
# ll1=ax.plot(np.arange(days),logsb0.mean(axis=1)[:days],'--', color='darkorange',label ='No movement')[0]
# ll3=ax.plot(np.arange(days),logsb20.mean(axis=1)[:days],'--', color='darkorange',label ='20\% Speed')[0]

# days = 213


lastday = 505
b20[320:] = b20[320]
bavg.prtinf = np.pad(bavg.prtinf,(0,lastday))
bavg.prinf = np.pad(bavg.prinf,(0,lastday))

qq1= '$\\frac{R_\\infty}{N}= $'+str(np.round(np.max(b[0].prtinf),2))
pp2= '$\\frac{R_\\infty}{N}= $'+ str(np.round(np.max(dfn[0,:]),2))
pp3= '$\\frac{R_\\infty}{N}= $'+str(np.round(np.max(dfn[1,:]),2))
ax.set_ylim(-0.2,0.4)
ax2.set_ylim(0,0.17)
l1=ax2.plot(np.arange(0,lastday),bavg.prinf[0:lastday],color='olive',label ='Baseline, '+qq1)[0]
l2=ax2.plot(np.arange(0,lastday),b20[4:lastday+4],'--', color='darkorange',label ='20\% speed, '+pp3)[0]
l3=ax2.plot(np.arange(0,lastday),b0[4:lastday+4],':', color='saddlebrown',label ='No movement, ' + pp2)[0]

ax2.set_yticks(np.arange(0,0.17,0.04))

#plt.setp(l.get_title(), multialignment='center')

# ax.set_title('Growth rate')
ax2.set_title('Infected')
ax2.legend(bbox_to_anchor=(.2,.9))
line_labels=['Baseline Spatial-SIR',
    'Speed = 20\%','No Movement']

#fig.legend([ll2,ll1,ll3],line_labels, bbox_to_anchor=(1,.9), fontsize= 12, borderaxespad=0.5,framealpha=1)



# fig,(ax1) = plt.subplots(1,1,figsize=(2*fsize,2*fsize))
# ax2 = ax1.twinx()
# ax2.spines['top'].set_visible(False)
# ax1.spines['top'].set_visible(False)

# ax1.set_xlabel('Days')
# ax1.set_ylabel('Fraction of population')
# ax1.set_ylim(0,1)
# ax2.set_yticks(np.arange(0,0.3,.05))
# ax2.set_ylim(0,0.15)

# dfn[1,320:] = dfn[1,320]
# lastday = 505

# bavg.prtinf = np.pad(bavg.prtinf,(0,lastday))
# bavg.prtinf[bavg.minlastday:]= bavg.prtinf[bavg.minlastday]
# bavg.prinf = np.pad(bavg.prinf,(0,lastday))


# dfn[1,272:] = dfn[1,272]

# ax1.plot(np.arange(0,lastday),bavg.prtinf[0:lastday],c='grey',label='Baseline', linewidth=1)
# ax1.plot(np.arange(0,lastday),dfn[1,4:lastday+4],':', linewidth=2, c='grey',label='Speed=20% baseline')
# ax1.plot(np.arange(0,lastday),dfn[0,4:lastday+4],'--', c='grey',label='No movement')
# ax2.plot(np.arange(0,lastday),bavg.prinf[0:lastday],'',linewidth=1.3, c='olive', label='Baseline')
# ax2.plot(np.arange(0,lastday),df2n[1,4:lastday+4],':',linewidth=1.7, c='saddlebrown', label='Speed=20% baseline')
# ax2.plot(np.arange(0,lastday),df2n[0,4:lastday+4],'--',linewidth=1.3, c='darkorange', label='No movement')

# ax1.legend(bbox_to_anchor=(1, .61), loc='center right',title='Infected + Recovered')
# ax2.legend(bbox_to_anchor=(1, .3),loc='center right',title='Infected \hspace{.02em} (right scale)')




fig.tight_layout()
plt.savefig(imagedir+'short-nomovement-rateslarge.pdf')
plt.close()



#%% Behavioral model, reduction of contacts

# load behavioral spatial
file = gzip.open(outputdir+'basebehLones.pickle.gz','rb')
lones = pickle.load(file)
file.close()
bLon = averageStats(lones)

# behavioral SIR
popSize = benchkwargs['q_popsize']
cluster = 30
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]
bb = SIRmodel(contRate*13.5, recRate, q_popsize=popSize, behModel={'type': 'Lones','phi': 0.01})

fix,ax = plt.subplots(figsize=(fsize*4/3,fsize))
ax.set_ylim(0,0.4)
ax.set_yticks(np.arange(0,0.41,0.1))
days = np.min((bLon.minlastday, bb.fracNotScared.size))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Days ')
ax.set_ylabel('Fraction of population')

ax.plot(np.arange(days), 1-bLon.fracNotScared[:days], color='darkorange', label = 'Behavioral Spatial SIR')
ax.plot(np.arange(days), 1-bb.fracNotScared[:days], '--', color='darkorange', label = 'Behavioral SIR')
ax.legend()
fig.tight_layout()
plt.savefig(imagedir+'SIR_beh_responses.pdf')
plt.close()

#%% Behavioral responses, comparison with benchmark

# load baseline
file = gzip.open(outputdir+'basesim.pickle.gz','rb')
benchplus = pickle.load(file)
file.close()
dfb = pd.DataFrame()
logsb = pd.DataFrame()
for idx,m in enumerate(benchplus) :
    dfb[idx] = np.sum(m.nstates[:,1:3],axis=1)
    dfb[str(idx)+'L1'] = dfb[idx].shift()
    logsb[str(idx)] = np.log(dfb[idx]) - np.log(dfb[str(idx)+'L1'])

# load behavioral
file = gzip.open(outputdir+'basebehLones.pickle.gz','rb')
lones = pickle.load(file)
file.close()
dfLon = pd.DataFrame()
logLon = pd.DataFrame()
for idx,m in enumerate(lones) :
    dfLon[idx] = np.sum(m.nstates[:,1:3],axis=1)
    dfLon[str(idx)+'L1'] = dfLon[idx].shift()
    logLon[str(idx)] = np.log(dfLon[idx]) - np.log(dfLon[str(idx)+'L1'])

# generate standard SIR data
popSize = benchkwargs['q_popsize']
cluster = 30
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]
bb = SIRmodel(contRate*13.5, recRate, q_popsize=popSize, behModel={'type': 'Lones','phi': 0.01})

aa = SIRmodel(contRate*13.5, recRate, q_init=5, q_popsize=popSize)
dfSIR = pd.DataFrame(aa.SIR[:,:3],index=np.arange(aa.day+1),columns=['S','I','R'])
dfSIR['ILag'] = dfSIR['I'].shift()
dfSIR['ld'] = np.log(dfSIR['I']) - np.log(dfSIR['ILag'])
dfSIRld = np.array(dfSIR['ld'])

bb = SIRmodel(contRate*13.5, recRate, q_init=5, q_popsize=popSize,behModel={'type': 'Lones','phi': 0.01},)
dfSIRLon = pd.DataFrame(bb.SIR[:,:3],index=np.arange(bb.day+1),columns=['S','I','R'])
dfSIRLon['ILag'] = dfSIRLon['I'].shift()
dfSIRLon['ld'] = np.log(dfSIRLon['I']) - np.log(dfSIRLon['ILag'])
dfSIRLonld = np.array(dfSIRLon['ld'])

bp = averageStats(benchplus)
bLon = averageStats(lones)
days = 86

fig,(ax,ax2) = plt.subplots(1,2,figsize=(2*fsize,fsize))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Days ')

li1=ax.plot(np.arange(days),dfSIRld[:days],'--',color='olive',label ='SIR')[0]
li2=ax.plot(np.arange(days),dfSIRLonld[:days],'-.',color='darkorange',label ='Behavioral SIR')[0]
li3=ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Spatial SIR')[0]
li4=ax.plot(np.arange(days),logLon.mean(axis=1)[:days],':',linewidth=1.9, color='darkorange',label ='Behavioral spatial SIR')[0]

maSIR = str(np.round(1-aa.SIR[days-1,0]/25600,2))
mbSIR = str(np.round(np.max(ravg.prtinf),2))
mb = str(np.round(np.max(bp.prtinf),2))
mLon = str(np.round(np.max(bLon.prtinf),2))

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Days ')
ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,'--',color='olive',label = maSIR)[0]
ax2.plot(np.arange(days),bb.SIR[:days,1]/25600,'-.',color='darkorange',label = mbSIR)[0]
ax2.plot(np.arange(days),bp.prinf[:days],color='olive',label =mb)[0]
ax2.plot(np.arange(days),bLon.prinf[:days],':',linewidth=1.9, color='darkorange',label = mLon)[0]

ax2.legend(title='$\\frac{R_\\infty}{N}$', bbox_to_anchor=(.52,.38), title_fontsize=15)
ax.set_title('Growth rate')
ax2.set_title('Infected ')

line_labels=['SIR','Behavioral SIR','Baseline Spatial-SIR',
    'Behavioral Spatial-SIR']

fig.legend([li1,li2,li3,li4],line_labels, bbox_to_anchor=(.65,.85), fontsize= 12, framealpha=1)


fig.tight_layout()
plt.savefig(imagedir+'SIR_beh.pdf')
plt.close()

#%% Behavioral model, local, reduction of contacts

# load behavioral spatial
file = gzip.open(outputdir+'basebehLones_local.pickle.gz','rb')
lones = pickle.load(file)
file.close()
bLonL = averageStats(lones)

# load behavioral spatial
file = gzip.open(outputdir+'basebehLones.pickle.gz','rb')
lonesL = pickle.load(file)
file.close()
bLon = averageStats(lonesL)

# behavioral SIR
popSize = benchkwargs['q_popsize']
cluster = 30
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]
bb = SIRmodel(contRate*13.5, recRate, q_popsize=popSize, behModel={'type': 'Lones','phi': 0.01})

fig,(ax,ax2) = plt.subplots(1,2,figsize=(2*fsize,fsize))
#ax.set_ylim(0,0.4)
ax.set_yticks(np.arange(0,0.41,0.05))
days = np.min((bLon.minlastday, bb.fracNotScared.size))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Days ')
#ax.set_ylabel('Fraction of population')

lb1 = ax.plot(np.arange(days), 1-bLonL.fracNotScared[:days],'', color='olive', label = 'Behavioral Spatial SIR (local)')[0]
lb2 = ax.plot(np.arange(days), 1-bLon.fracNotScared[:days],':', linewidth=1.9, color='darkorange', label = 'Behavioral Spatial SIR')[0]

#plt.savefig(imagedir+'SIR_beh_responses_local.pdf')

## %% Behavioral responses, comparison with benchmark

days = 112

mb1 = np.round(np.max(bLonL.prtinf),2)
mb2 = np.round(np.max(bLon.prtinf),2)



ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Days ')
ax2.plot(np.arange(days),bLonL.prinf[:days],color='olive',label =mb1)
ax2.plot(np.arange(days),bLon.prinf[:days],':',linewidth=1.9, color='darkorange',label =mb2)
ax2.set_ylim(0,0.10)

bLon.prtinf = np.append(bLon.prtinf,[0]*300)
bLon.prtinf[bLon.minlastday:] = bLon.prtinf[bLon.minlastday]
bLonL.prtinf[bLonL.minlastday:] = bLonL.prtinf[bLonL.minlastday]

ax.set_title('Reduction of contacts')
ax2.set_title('Infected')

ax2.legend(title='$\\frac{R_\\infty}{N}$', bbox_to_anchor=(.62,.55), title_fontsize=15)

line_labels=['Behavioral Spatial-SIR (local)',
    'Behavioral Spatial-SIR']

fig.legend([lb1,lb2],line_labels, bbox_to_anchor=(.7,.85), fontsize= 12, framealpha=1)

fig.tight_layout()
plt.savefig(imagedir+'SIR_beh_local.pdf')
plt.close()


#%% Figure comparing estimated betas in baseline and behavioral with real beta

#   Before doing this, run sim_estimate_dens.do in stata to generate densbetas.csv
#   and dens-behbetas.csv
betaframe_nobeh = pd.read_csv(outputdir+'densbetas.csv')
density_nobeh = betaframe_nobeh['density']
betalamd_nobeh = betaframe_nobeh['beta']

betaframe = pd.read_csv(outputdir+'dens-beh_pbetas.csv')
density = betaframe['density']
beh_betalamd = betaframe['beta']
betaframe['realbeta']= 0.054*13.5*betaframe['density']

fig,(ax2) = plt.subplots(1,1,figsize=(fsize,fsize))
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Density $d$ relative to benchmark')
ax2.set_ylabel('Coefficient')
#ax2.set_ylim(0.33,1.8)
ax2.plot(density,betaframe['realbeta'],':',color='saddlebrown',label='$\\beta$', marker='s')
ax2.plot(density_nobeh,betalamd_nobeh,color='olive',label='$\\hat{\\beta} $ Baseline Spatial-SIR data', marker='^')
ax2.plot(density,beh_betalamd,color='darkorange',label='$\\hat{\\beta}$ Behavioral Spatial-SIR data', marker='o')
ax2.legend()

fig.tight_layout()
plt.savefig(imagedir+'est_densitybetas.pdf')
plt.close()

#%% Apply policy in SIR using the estimated beta, baseline
if __name__ == "__main__":
    from class_SIRmodel import SIRmodel
    from class_averageStats import averageStats

    thisdensity = 0.5
    shutday = 20

    # importing from estimates in sim_estimate_dens.do
    betaframe_nobeh = pd.read_csv(outputdir+'densbetas.csv')
    betaframe_nobeh['truebeta']= 0.054*13.5*betaframe_nobeh['density']

    baseline = 25600
    cluster = 10
    beta = np.array(betaframe_nobeh[betaframe_nobeh['density']==thisdensity]['beta'])[0]
    truebeta = np.array(betaframe_nobeh[betaframe_nobeh['density']==thisdensity]['truebeta'])[0]
    delta = benchkwargs['p_probr'][0]

    # simulate effect of policy in SIR
    base = SIRmodel(beta, delta ,q_popsize=baseline, q_init=cluster, )
    policy = SIRmodel(beta, delta, q_popsize=baseline, q_init=cluster, p_shutr=[shutday,.25], p_openr=[999,0])

    # compare with effect of policy in Spatial-SIR
    file = gzip.open(outputdir+'dens-20-80-25pc.gz','rb')
    allRandmodels = pickle.load(file)
    file.close()

    modbase = list()
    modpoli = list()
    for model in allRandmodels:
        if abs(model.q_citysize - 1/np.sqrt(thisdensity)) < .000001:
            if model.p_shutr[0] == 999:
                modbase.append(model)
            else:
                modpoli.append(model)
    avbase = averageStats(modbase)
    avpoli = averageStats(modpoli)

    ##############
    thisdensity = 1
    shutday = 20
    beta2 = np.array(betaframe_nobeh[betaframe_nobeh['density']==thisdensity]['beta'])[0]
    delta = benchkwargs['p_probr'][0]
    truebeta2 = np.array(betaframe_nobeh[betaframe_nobeh['density']==thisdensity]['truebeta'])[0]

    # simulate effect of policy in SIR
    base2 = SIRmodel(beta2, delta ,q_popsize=baseline, q_init=cluster, )
    policy2 = SIRmodel(beta2, delta, q_popsize=baseline, q_init=cluster, p_shutr=[shutday,.25], p_openr=[999,0])
    #base.plot(base.day)
    #policy.plot(policy.day)

    # compare with effect of policy in Spatial-SIR
    modbase2 = list()
    modpoli2 = list()
    for model in allRandmodels:
        if abs(model.q_citysize - 1/np.sqrt(thisdensity)) < .000001:
            if model.p_shutr[0] == 999:
                modbase2.append(model)
            else:
                modpoli2.append(model)
    avbase2 = averageStats(modbase2)
    avpoli2 = averageStats(modpoli2)

    ##############
    thisdensity = 1.5
    shutday = 20
    beta3 = np.array(betaframe_nobeh[betaframe_nobeh['density']==thisdensity]['beta'])[0]
    delta = benchkwargs['p_probr'][0]
    truebeta3 = np.array(betaframe_nobeh[betaframe_nobeh['density']==thisdensity]['truebeta'])[0]

    # simulate effect of policy in SIR
    base3 = SIRmodel(beta3, delta ,q_popsize=baseline, q_init=cluster, )
    policy3 = SIRmodel(beta3, delta, q_popsize=baseline, q_init=cluster, p_shutr=[shutday,.25], p_openr=[999,0])
    #base.plot(base.day)
    #policy.plot(policy.day)

    # compare with effect of policy in Spatial-SIR
    modbase3 = list()
    modpoli3 = list()
    for model in allRandmodels:
        if abs(model.q_citysize - 1/np.sqrt(thisdensity)) < .000001:
            if model.p_shutr[0] == 999:
                modbase3.append(model)
            else:
                modpoli3.append(model)
    avbase3 = averageStats(modbase3)
    avpoli3 = averageStats(modpoli3)

    fsize=3.5
    days= 80
    fig = plt.figure(figsize=(3*fsize,fsize*1))
    spec = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)

    ax = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])
    ax3 = fig.add_subplot(spec[0,2])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    ax.set_ylim(-.03,.16)
    ax.set_yticks(np.arange(0,0.16,.05))
    ax.plot(np.arange(20,days),avbase.prinf[20:days]-avpoli.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax.plot(np.arange(20,days),base.SIR[20:days,1]/25600-policy.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.set_ylim(-.03,.16)
    ax2.set_yticks(np.arange(0,0.16,.05))
    ax2.plot(np.arange(20,days),avbase2.prinf[20:days]-avpoli2.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax2.plot(np.arange(20,days),base2.SIR[20:days,1]/25600-policy2.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.set_xlabel('Days ')
    ax3.set_ylim(-.03,.16)
    ax3.set_yticks(np.arange(0,0.16,.05))
    ax3.plot(np.arange(20,days),base3.SIR[20:days,1]/25600-policy3.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')
    ax3.plot(np.arange(20,days),avbase3.prinf[20:days]-avpoli3.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax3.hlines(0,20,80,color='grey',lw=0.3)
    ax2.hlines(0,20,80,color='grey',lw=0.3)
    ax.hlines(0,20,80,color='grey',lw=0.3)

    #first_legend =ax2.legend(bbox_to_anchor=(.48,.35), fontsize=10)
    l2 = ax2.legend(handles=[],frameon=False, title_fontsize=12,
                               title='Density 1 \n $\\hat{\\beta} ='+str(round(beta2,2))
                                +'$',#' \n$  \\beta ='+str(round(truebeta2,2))+'$'
                               loc='upper left')

    la =plt.legend(bbox_to_anchor=(-1.3,.8), fontsize=10, framealpha=1)

    #ax2.add_artist(a3)

    l1 = ax.legend(handles=[], frameon=False, title_fontsize=12,
                         title='Density 0.5 \n $\\hat{\\beta} ='+str(round(beta,2))
                         +'$', #'\n$ \\beta ='+str(round(truebeta,2))+'$',
                         loc='upper left')

    l3 = plt.legend(handles=[], frameon=False, title_fontsize=12,
                         title='Density 1.5 \n $\\hat{\\beta} ='+str(round(beta3,2))
                         +'$', #'\n$ \\beta ='+str(round(truebeta,2))+'$',
                         loc='upper left')
    ax3.add_artist(la)
    ax3.add_artist(l3)
    fig.tight_layout()
    plt.savefig(imagedir+'estimatedbeta_policies.pdf')
    plt.close()

#%% Apply policy in SIR using the estimated beta, behavioral
if __name__ == "__main__":
    from class_SIRmodel import SIRmodel
    from class_averageStats import averageStats

    thisdensity = 0.5
    shutday = 20

    # importing from estimates in sim_estimate_dens.do
    betaframe_beh = pd.read_csv(outputdir+'dens-beh_pbetas.csv')
    betaframe_beh['truebeta']= 0.054*13.5*betaframe_beh['density']

    baseline = 25600
    cluster = 10
    beta = np.array(betaframe_beh[betaframe_beh['density']==thisdensity]['beta'])[0]
    truebeta = np.array(betaframe_beh[betaframe_beh['density']==thisdensity]['truebeta'])[0]
    delta = benchkwargs['p_probr'][0]

    # simulate effect of policy in SIR
    base = SIRmodel(beta, delta ,q_popsize=baseline, q_init=cluster, )
    policy = SIRmodel(beta, delta, q_popsize=baseline, q_init=cluster, p_shutr=[shutday,.25], p_openr=[999,0])

    # compare with effect of policy in Spatial-SIR
    file = gzip.open(outputdir+'dens-beh_p-20-80-25pc.gz','rb')
    allRandmodels = pickle.load(file)
    file.close()

    modbase = list()
    modpoli = list()
    for model in allRandmodels:
        if abs(model.q_citysize - 1/np.sqrt(thisdensity)) < .000001:
            if model.p_shutr[0] == 999:
                modbase.append(model)
            else:
                modpoli.append(model)
    avbase = averageStats(modbase)
    avpoli = averageStats(modpoli)

    ##############
    thisdensity = 1
    shutday = 20
    beta2 = np.array(betaframe_beh[betaframe_beh['density']==thisdensity]['beta'])[0]
    delta = benchkwargs['p_probr'][0]
    truebeta2 = np.array(betaframe_beh[betaframe_beh['density']==thisdensity]['truebeta'])[0]

    # simulate effect of policy in SIR
    base2 = SIRmodel(beta2, delta ,q_popsize=baseline, q_init=cluster, )
    policy2 = SIRmodel(beta2, delta, q_popsize=baseline, q_init=cluster, p_shutr=[shutday,.25], p_openr=[999,0])
    #base.plot(base.day)
    #policy.plot(policy.day)

    # compare with effect of policy in Spatial-SIR
    modbase2 = list()
    modpoli2 = list()
    for model in allRandmodels:
        if abs(model.q_citysize - 1/np.sqrt(thisdensity)) < .000001:
            if model.p_shutr[0] == 999:
                modbase2.append(model)
            else:
                modpoli2.append(model)
    avbase2 = averageStats(modbase2)
    avpoli2 = averageStats(modpoli2)

    ##############
    thisdensity = 1.5
    shutday = 20
    beta3 = np.array(betaframe_beh[betaframe_beh['density']==thisdensity]['beta'])[0]
    delta = benchkwargs['p_probr'][0]
    truebeta3 = np.array(betaframe_beh[betaframe_beh['density']==thisdensity]['truebeta'])[0]

    # simulate effect of policy in SIR
    base3 = SIRmodel(beta3, delta ,q_popsize=baseline, q_init=cluster, )
    policy3 = SIRmodel(beta3, delta, q_popsize=baseline, q_init=cluster, p_shutr=[shutday,.25], p_openr=[999,0])
    #base.plot(base.day)
    #policy.plot(policy.day)

    # compare with effect of policy in Spatial-SIR
    modbase3 = list()
    modpoli3 = list()
    for model in allRandmodels:
        if abs(model.q_citysize - 1/np.sqrt(thisdensity)) < .000001:
            if model.p_shutr[0] == 999:
                modbase3.append(model)
            else:
                modpoli3.append(model)
    avbase3 = averageStats(modbase3)
    avpoli3 = averageStats(modpoli3)


    fsize=3.5
    days= 80
    fig = plt.figure(figsize=(3*fsize,fsize*1))
    spec = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)

    ax = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])
    ax3 = fig.add_subplot(spec[0,2])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    ax.set_ylim(-.03,.16)
    ax.set_yticks(np.arange(0,0.16,.05))
    ax.plot(np.arange(20,days),avbase.prinf[20:days]-avpoli.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax.plot(np.arange(20,days),base.SIR[20:days,1]/25600-policy.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.set_ylim(-.03,.16)
    ax2.set_yticks(np.arange(0,0.16,.05))
    ax2.plot(np.arange(20,days),avbase2.prinf[20:days]-avpoli2.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax2.plot(np.arange(20,days),base2.SIR[20:days,1]/25600-policy2.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.set_xlabel('Days ')
    ax3.set_ylim(-.03,.16)
    ax3.set_yticks(np.arange(0,0.16,.05))
    ax3.plot(np.arange(20,days),base3.SIR[20:days,1]/25600-policy3.SIR[20:days,1]/25600,'--',color='darkorange',label ='Behavioral SIR')
    ax3.plot(np.arange(20,days),avbase3.prinf[20:days]-avpoli3.prinf[20:days],color='olive',label ='Behavioral Spatial-SIR')
    ax3.hlines(0,20,80,color='grey',lw=0.3)
    ax2.hlines(0,20,80,color='grey',lw=0.3)
    ax.hlines(0,20,80,color='grey',lw=0.3)

    l2 = ax2.legend(handles=[],frameon=False, title_fontsize=12,
                               title='Density 1 \n $\\hat{\\beta} ='+str(round(beta2,2))
                                +'$',#' \n$  \\beta ='+str(round(truebeta2,2))+'$'
                               loc='upper left')

    la =plt.legend(bbox_to_anchor=(-1.3,.8), fontsize=10, framealpha=1)

    #ax2.add_artist(a3)

    l1 = ax.legend(handles=[], frameon=False, title_fontsize=12,
                         title='Density 0.5 \n $\\hat{\\beta} ='+str(round(beta,2))
                         +'$', #'\n$ \\beta ='+str(round(truebeta,2))+'$',
                         loc='upper left')
    l3 = plt.legend(handles=[], frameon=False, title_fontsize=12,
                         title='Density 1.5 \n $\\hat{\\beta} ='+str(round(beta3,2))
                         +'$', #'\n$ \\beta ='+str(round(truebeta,2))+'$',
                         loc='upper left')
    ax3.add_artist(la)
    ax3.add_artist(l3)
    fig.tight_layout()
    plt.savefig(imagedir+'estimatedbeta_beh-policies.pdf')
    plt.close()

#%% prediction plots
    import pandas as pd
    fsize=3.5
    popsize = 25600
    nobeh = pd.read_csv(outputdir+'predictions.csv')
    beh = pd.read_csv(outputdir+'predictions_beh.csv')

    fig = plt.figure(figsize=(3*fsize,fsize*1))
    spec = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)

    ax = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])
    ax3 = fig.add_subplot(spec[0,2])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.plot(nobeh['t'],nobeh['active']/popsize, color='olive', )
    ax.plot(nobeh['t'],nobeh['myactive20']/popsize, ':', color='olive', )
    ax.plot(beh['t'],beh['active']/popsize, color='darkorange')
    ax.plot(beh['t'],beh['myactive20']/popsize, ':', color='darkorange')
    ax.set_xlabel('Days ')
    ax.set_yticks(np.arange(0,0.11,0.02))
    ax.set_ylim(0.01,0.08)
    ax.legend(handles=[], title='Infected', title_fontsize=12, frameon=False, loc='lower left')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.plot(nobeh['t'],nobeh['outside'], color='olive', )
    ax2.plot(nobeh['t'],nobeh['myoutside20'], ':', color='olive', )
    ax2.plot(beh['t'],beh['outside'], color='darkorange')
    ax2.plot(beh['t'],beh['myoutside20'], ':', color='darkorange')
    ax2.set_xlabel('Days ')
    ax2.set_yticks(np.arange(0.6,0.751,0.05))
    ax2.set_ylim(0.55,0.76)
    ax2.legend(handles=[], title='Contacts', title_fontsize=12, frameon=False, loc='lower left')

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.plot(nobeh['t'],nobeh['growthi'], color='olive', label='Baseline Spatial-SIR')
    ax3.plot(nobeh['t'],nobeh['mygrowthi20'], ':', color='olive',label='Baseline, prediction' )
    ax3.plot(beh['t'],beh['growthi'], color='darkorange', label='Behavioral Spatial-SIR')
    ax3.plot(beh['t'],beh['mygrowthi20'], ':', color='darkorange', label='Behavioral, prediction')
    ax3.set_xlabel('Days ')
    ax3.set_ylim(-.12,0.025)
    la = ax3.legend(bbox_to_anchor=(0.2,0.48), framealpha=1)
    l3 = plt.legend(handles=[], title='Growth rate', title_fontsize=12, frameon=False,
                    loc='lower center')
    ax3.add_artist(l3)
    ax3.add_artist(la)
    fig.tight_layout()
    plt.savefig(imagedir+'prediction_nobias.pdf')

    fig = plt.figure(figsize=(3*fsize,fsize*1))
    spec = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)
    ax = fig.add_subplot(spec[0,0])
    ax2 = fig.add_subplot(spec[0,1])
    ax3 = fig.add_subplot(spec[0,2])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.plot(nobeh['t'],nobeh['active']/popsize, color='olive')
    ax.plot(nobeh['t'],nobeh['myactive20_dbias']/popsize, ':', color='olive')
    ax.plot(beh['t'],beh['active_dbias']/popsize, color='darkorange')
    ax.plot(beh['t'],beh['myactive20_dbias']/popsize, ':', color='darkorange')
    ax.set_xlabel('Days ')
    ax.set_yticks(np.arange(0,0.11,0.02))
    ax.set_ylim(0,0.113)
    ax.legend(handles=[], title='Infected', title_fontsize=12, frameon=False, loc='lower left')


    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.plot(nobeh['t'],nobeh['outside_dbias'], color='olive', label='Simulation')
    ax2.plot(nobeh['t'],nobeh['myoutside20_dbias'], ':', color='olive', label='Prediction')
    ax2.plot(beh['t'],beh['outside_dbias'], color='darkorange')
    ax2.plot(beh['t'],beh['myoutside20_dbias'], ':', color='darkorange')
    ax2.set_xlabel('Days ')
    ax2.set_yticks(np.arange(0.6,0.751,0.05))
    ax2.set_ylim(0.55,0.76)
    ax2.legend(handles=[], title='Contacts', title_fontsize=12, frameon=False, loc='lower right')

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.plot(nobeh['t'],nobeh['growthi_dbias'], color='olive', label='Baseline Spatial-SIR')
    ax3.plot(nobeh['t'],nobeh['mygrowthi20_dbias'], ':', color='olive',label='Baseline, prediction' )
    ax3.plot(beh['t'],beh['growthi_dbias'], color='darkorange', label='Behavioral Spatial-SIR')
    ax3.plot(beh['t'],beh['mygrowthi20_dbias'], ':', color='darkorange', label='Behavioral, prediction')
    ax3.set_xlabel('Days ')
    ax3.set_ylim(-.12,0.025)
    la = ax3.legend(bbox_to_anchor=(-0.85,0.9), framealpha=1)
    l3 = plt.legend(handles=[], title='Growth rate', title_fontsize=12, frameon=False,
                    loc='lower left')
    ax3.add_artist(l3)
    ax3.add_artist(la)
    fig.tight_layout()
    plt.savefig(imagedir+'prediction_withbias.pdf')
    plt.close()
