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

import labellines
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
    "q_popsize"   : 25600,
    "q_days"      : 300,
    "q_printOption" : 0.6,
    'g_maxinf'    : 1,
}  

imagedir = 'output/images/'+prefix

popSize = benchkwargs['q_popsize']
cluster = 30
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]


# prefix= 'nc6-' 
# outputdir = 'output/'+prefix
# (probc,avgstep)=[0.039607841, 0.072992994]
# benchkwargs = {
#     "q_seed" : 2443,
#     "p_proby"     : [0.0],
#     "p_probr"     : [1/7],
#     "p_probc"     : [[0, probc, 0, 0, 0]], #prob of contagion by type and state
#     "p_probd"     : [0],
#     "p_infradius" : 0.013016,
#     "p_avgstep"   : avgstep,
#     "p_stdstep"   : 0,
#     "p_cluster"   : 'cluster',
#     "q_popsize"   : 25600,
#     "q_days"      : 300,
#     "q_printOption" : 0.6,
#     'g_maxinf'    : 1,
# }  

imagedir = 'output/images/'+prefix

popSize = benchkwargs['q_popsize']
cluster = 30
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]

#%% Spatial progression of infections, benchmarkmodel dates 5, 40, 110

from class_spatialModels import spatialSAYDR    

printdays = [3, 10, 20, 30]
kwargs = deepcopy(benchkwargs)
kwargs['q_days'] = np.max(printdays) + 1

m = spatialSAYDR(**kwargs)
  
for day in range(1,m.q_days):

    print(m.computeStats(day,1))
    if day in printdays:
        m.drawPositions(day,savefig='baseline',S=True,folder=imagedir,ext='.pdf')
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
        m.drawPositions(day,savefig='randomcluster',S=True,folder=imagedir,ext='.pdf')
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
        m.drawPositions(day,savefig='nomove',S=True,folder=imagedir,ext='.pdf')
        pass
    m.aDayInTheLife(day) 
      
    if np.sum(m.state==m.I)+np.sum(m.state==m.Y) <=15:
        break

    #%% Initial location in space, heterogeneous density model

from class_spatialModels import spSAYDR_hetDensity
    
mod1 = spSAYDR_hetDensity(q_lambda=1,**benchkwargs)
mod1.drawPositions(0,{'savefig':'hetdens','folder':imagedir,'S': True, 
                      'colors':['grey','orange','orange','red','green'], 'ext':'.pdf'}, )
    

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
dfSIR = pd.DataFrame(aa.SIR,index=np.arange(aa.day+1),columns=['S','I','R'])
dfSIR['ILag'] = dfSIR['I'].shift()
dfSIR['ld'] = np.log(dfSIR['I']) - np.log(dfSIR['ILag'])
dfSIRld = np.array(dfSIR['ld'])

days=aa.SIR[:,0].size


fig,(ax,ax2) = plt.subplots(1,2,figsize=(7,3.5))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Days ')
ax.set_ylim(-0.2,0.5)
ax.set_xlim(-1,92)
ax.plot(np.arange(len(dfSIRld)),dfSIRld,':', linewidth=2, color='saddlebrown',label ='SIR')
ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Spatial SIR')
ax.plot(np.arange(days),logs.mean(axis=1)[:days],'--', color='darkorange',label ='Sp. SIR with random pos.')
 
ax.legend(title='Infection growth rates')
days=aa.SIR[:,0].size

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Days ')
ax2.set_xlim(-1,92)
ax2.set_ylim(0,0.6)
ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,':', linewidth=2, color='saddlebrown',label ='SIR')
ax2.plot(np.arange(days),bavg.prinf[:days],color='olive',label ='Spatial SIR')
ax2.plot(np.arange(days),ravg.prinf[:days],'--', color='darkorange',label ='Sp. SIR w/ random pos.')   
ax2.set_ylabel('Fraction of population')
ax.set_ylabel('Growth rate')
ax2.legend(title='Infected ')
fig.tight_layout()
 

plt.savefig(imagedir+'density_contagion2.pdf')
plt.show()

#%% Random locations

# Benchmark model
file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
bavg = averageStats(b)
file1.close()

file = gzip.open(outputdir+'random.pickle.gz','rb')
r = pickle.load(file)
ravg = averageStats(r)
    
fig,(ax) = plt.subplots(1,1,figsize=(3.75,3.75))
ax2 = ax.twinx()

ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

minday = 0
g_plotdays = 90
g_maxinf = .6

ravg.prtinf[ravg.minlastday:] = ravg.prtinf[ravg.minlastday]

ax2.set_xticks(np.arange(0,g_plotdays,20))
ax2.set_yticks(np.arange(0,g_maxinf+.1,.2))
ax2.set_ylim(0,g_maxinf)
ax.set_xlabel('Days ')
ax.set_ylabel('Fraction of population')

ax.set_ylim(0,1)
ax.spines['top'].set_visible(False)
ax.plot(np.arange(g_plotdays), bavg.prtinf[0:g_plotdays],'',c='grey', label= 'Baseline        ')
ax.plot(np.arange(g_plotdays), ravg.prtinf[0:g_plotdays],'--' ,c='grey', label= 'Active')

ax2.plot(np.arange(g_plotdays), bavg.prinf[0:g_plotdays], c='olive',label="Baseline        ")
ax2.plot(np.arange(g_plotdays), ravg.prinf[0:g_plotdays], '--', c='darkorange',label="Random")
l = ax.legend(bbox_to_anchor=(1, .74), title='Infected $+$ Recovered', loc='center right')
#plt.setp(l.get_title(), multialignment='center')

#ax2.legend(loc='center right',ncol=1, title= '        Fraction of population \n(black lines: Infected )')
l2 = ax2.legend(title='Infected \\hspace{.25ex} (right scale)', loc='center right')
#plt.setp(l2.get_title(), multialignment='center')

fig.tight_layout()

plt.savefig(imagedir+'short-random-rates.pdf')
plt.show()  
    
print('Peak active, baseline',np.round(max(bavg.prinf),2)
      ,'day ',np.where(max(bavg.prinf)==bavg.prinf)[0])
print('Peak active, random',np.round(max(ravg.prinf),2)
      ,'day ',np.where(max(ravg.prinf)==ravg.prinf)[0])
avgmaxinf = np.average(list(map(lambda x: np.max(x.prtinf),b)))
print('Steady state infected baseline', np.round(avgmaxinf,4))
avgmaxinf2 = np.average(list(map(lambda x: np.max(x.prtinf),r)))
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

fig,(ax1,ax) = plt.subplots(1,2,figsize=(7,3.5))
ax12 = ax1.twinx()
ax2 = ax.twinx()
g_plotdays=150
g_maxinf=0.6

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax12.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax1.set_xticks(np.arange(0,g_plotdays+1,25))
ax1.set_ylim(0,1)
ax12.set_ylim(0,g_maxinf)
ax.set_ylim(0,1)
ax2.set_ylim(0,g_maxinf)
ax2.set_yticks(np.arange(0,g_maxinf+.1,.2))
ax12.set_yticks(np.arange(0,g_maxinf+.1,.2))


ax1.set_xlabel('Days ')
ax1.set_ylabel('Fraction of population')

s_qavg.prtinf[s_qavg.minlastday:] = s_qavg.prtinf[s_qavg.minlastday]
s_bavg.prtinf[s_bavg.minlastday:] = s_bavg.prtinf[s_bavg.minlastday]
   
#spatial SIR
ax12.plot(np.arange(g_plotdays),s_b2avg.prinf[0:g_plotdays],':',c='saddlebrown',linewidth=1.9, label="$4* N$")
ax12.plot(np.arange(g_plotdays),s_bavg.prinf[0:g_plotdays],c='olive',label="Baseline")
ax12.plot(np.arange(g_plotdays),s_qavg.prinf[0:g_plotdays],'--',c='darkorange',label="$1/4*N$")

ax1.plot(np.arange(g_plotdays),s_b2avg.prtinf[0:g_plotdays],':',c='grey',linewidth=1.7, label="$4*N$")
ax1.plot(np.arange(g_plotdays),s_bavg.prtinf[0:g_plotdays],'',c='grey',linewidth=1.3, label="Baseline $N$")
ax1.plot(np.arange(g_plotdays),s_qavg.prtinf[0:g_plotdays],'--',c='grey',linewidth=1.3, label='$1/4*N$')

ax1.legend(bbox_to_anchor=(1, .76), loc='center right',title='Infected + Recovered')
ax12.legend(bbox_to_anchor=(1, .42),loc='center right',ncol=1,title='Infected\hspace{.35em}  (right scale) ')

g_plotdays = 70

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks(np.arange(0,g_plotdays,20))
ax.set_yticks(np.arange(0,1.1,.2))
ax.set_xticks(np.arange(0,201,25))
ax.set_xlabel('Days ')
   
ax2.set_xticks(np.arange(0,201,20))

ax2.plot(np.arange(g_plotdays),b2avg.prinf[0:g_plotdays],':',c='saddlebrown',linewidth=1.9,label="$4*N$")
ax2.plot(np.arange(g_plotdays),bavg.prinf[0:g_plotdays],c='olive',label="Baseline")
ax2.plot(np.arange(g_plotdays),sqavg.prinf[0:g_plotdays],'--',c='darkorange',label="1/$4*N$")

ax.plot(np.arange(g_plotdays),b2avg.prtinf[0:g_plotdays],':',c='grey',linewidth=1.7,label="$4*N$")
ax.plot(np.arange(g_plotdays),bavg.prtinf[0:g_plotdays],'',c='grey',linewidth=1.3,label="Baseline $N$")
ax.plot(np.arange(g_plotdays),sqavg.prtinf[0:g_plotdays],'--',c='grey',linewidth=1.3,label="$1/4*N$")

ax.legend(bbox_to_anchor=(1, .76), loc='center right',title='Infected + Recovered')
ax2.legend(bbox_to_anchor=(1, .42), loc='center right',title='Infected\hspace{0.35em}  (right scale)')

fig.tight_layout()
plt.savefig(imagedir+'SIR-citysize-rates.pdf')
plt.show()           


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

fig,(ax,ax3) = plt.subplots(1,2,figsize=(7,3.5))
ax2 = ax3.twinx()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax3.spines['top'].set_visible(False)

ax.set_xlabel('Days')
ax3.set_xlabel('Days')
ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Baseline')
ax.plot(np.arange(days),logs6.mean(axis=1)[:days],'--', color='darkorange',label ='1/6*density, 6*contagious')

days = 213

ax2.set_ylim(0,1)

ax2.set_xlabel('Days ')
ax.set_ylim(-0.2,0.4)
ax2.set_ylim(0,0.2)
ax3.set_ylim(0,1)

ax2.plot(np.arange(days),bavg.prinf[:days],color='olive',label ='Baseline')
ax2.plot(np.arange(days),b6.prinf[:days],'--', color='darkorange',label ='\\parbox{5em}{1/6*density,\\newline 6*contagious}')

ax2.set_yticks(np.arange(0,0.21,0.05))
ax3.plot(np.arange(days), bavg.prtinf[0:days],'',c='grey', label = 'Baseline')
ax3.plot(np.arange(days), b6.prtinf[0:days],'--',c='grey', label = '\\parbox{6.5em}{1/6*density,\\newline 6*contagious}')
l= ax.legend(loc='upper right', title='Infection growth rate')

ax3.legend(bbox_to_anchor=(1, .75), title='Infected + Recovered', loc='center right')
ax2.legend(bbox_to_anchor=(1, .43), title='Infected \hspace{0.05em} (right scale)', loc='center right')

plt.setp(l.get_title(), multialignment='center')

ax.set_ylabel('Growth rate')
ax3.set_ylabel('Fraction of population')

fig.tight_layout()
plt.savefig(imagedir+'short-density_contagion1.pdf')


#%% Different city density (different size same population) 

file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
bavg = averageStats(b)
file1.close()

file1 = gzip.open(outputdir+'density.pickle.gz','rb')
allmods=  pickle.load(file1)
file1.close()

modh = list()
mod2 = list()
for mod in allmods:
    if np.round(mod.q_citysize,1) == 0.7 :
        modh.append(mod)
    if np.round(mod.q_citysize,1) == 1.4 :
        mod2.append(mod)

b2avg = averageStats(mod2)
havg = averageStats(modh)

      
fig,(ax) = plt.subplots(1,1,figsize=(3.75,3.75))
ax2 = ax.twinx()

g_plotdays = 213

### draw legend with empty data
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax.set_yticks(np.arange(0,1.1,.2))
ax.set_xlabel('Days ')
ax.set_ylabel('Fraction of population')

ax2.set_ylim(0,0.85)
ax.set_ylim(0,1.01)
ax2.set_yticks(np.arange(0,0.81,0.2))
ax.set_xticks(np.arange(0,g_plotdays,50))

bavg.prtinf[bavg.minlastday:] = bavg.prtinf[bavg.minlastday]
havg.prtinf[havg.minlastday:] = havg.prtinf[havg.minlastday]
b2avg.prtinf[b2avg.minlastday:] = b2avg.prtinf[b2avg.minlastday]

ax2.plot(np.arange(g_plotdays),havg.prinf[0:g_plotdays],':',linewidth=2, c='saddlebrown',label="$2*$ density")
ax2.plot(np.arange(g_plotdays),bavg.prinf[0:g_plotdays],c='olive',label="Baseline")
ax2.plot(np.arange(g_plotdays),b2avg.prinf[0:g_plotdays],'--',c='darkorange',label="$1/2*$ density")

ax.plot(np.arange(g_plotdays),havg.prtinf[0:g_plotdays],':',linewidth=1.7 ,c='grey', label="$2*$ density")
ax.plot(np.arange(g_plotdays),bavg.prtinf[0:g_plotdays],'',linewidth=1.3, c='grey', label="Baseline")
ax.plot(np.arange(g_plotdays),b2avg.prtinf[0:g_plotdays],'--',linewidth=1.3, c='grey', label="$1/2*$ density")

ax.legend(bbox_to_anchor=(1, .55), loc='center right',title='Infected + Recovered')
ax2.legend(bbox_to_anchor=(1, .25),loc='center right',title='Infected  \hspace{.1em}  (right scale)')

fig.tight_layout()
plt.savefig(imagedir+'short-3densities.pdf')
plt.show()           

print('Sir DAYS to peak',np.where(max(bavg.prinf)==bavg.prinf)[0]
      ,np.where(max(b2avg.prinf)==b2avg.prinf)[0],
      np.where(max(havg.prinf)==havg.prinf)[0])
print('Peak active, 1*',np.round(max(bavg.prinf),2)
      ,'day ',np.where(max(bavg.prinf)==bavg.prinf)[0])
print('Peak active, 2*',np.round(max(b2avg.prinf),2)
      ,'day ',np.where(max(b2avg.prinf)==b2avg.prinf)[0])
print('Peak active, 1/2*',np.round(max(havg.prinf),2)
      ,'day ',np.where(max(havg.prinf)==havg.prinf)[0])
avgmaxinf = np.average(list(map(lambda x: np.max(x.prtinf),b)))
avgmaxinf = np.max(bavg.prtinf)
print('Steady state infected 1*', np.round(avgmaxinf,4))
avgmaxinf2 = np.average(list(map(lambda x: np.max(x.prtinf),mod2)))
avgmaxinf2 = np.max(b2avg.prtinf)
print('Steady state infected 2*', np.round(avgmaxinf3,4))
avgmaxinf3 = np.max(havg.prtinf)
print('Steady state infected 1/2*', np.round(avgmaxinf2,4))
avgmaxinf3 = np.average(list(map(lambda x: np.max(x.prtinf),modh)))


#%% SIR density

popSize = benchkwargs['q_popsize']
cluster = 30
contRate = benchkwargs['p_probc'][0][1]
recRate = benchkwargs['p_probr'][0]
bavg = SIRmodel(contRate*13.4, recRate, q_popsize=popSize, q_init=cluster)
b2avg = SIRmodel(contRate*13.4*2, recRate, q_popsize=popSize, q_init=cluster)
bhavg = SIRmodel(contRate*13.4/2, recRate, q_popsize=popSize, q_init=cluster)
  
maxday=bavg.day
day2= b2avg.day    
dayh= bhavg.day

fig,(ax) = plt.subplots(1,1,figsize=(3.75,3.75))
ax2 = ax.twinx()
g_plotdays=150
g_maxinf=0.6

ax.set_ylim(0,1.01)
ax2.set_ylim(0,0.85)

g_plotdays =  83

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.set_xticks(np.arange(0,g_plotdays,20))
ax.set_yticks(np.arange(0,1.1,.2))
ax.set_xticks(np.arange(0,201,25))
ax2.set_yticks(np.arange(0,.81,.2))

ax.set_xlabel('Days ')

ax2.plot(np.arange(day2),b2avg.prinf[0:day2],':',c='saddlebrown',linewidth=1.9,label="$2*$ density")
ax2.plot(np.arange(g_plotdays),bavg.prinf[0:g_plotdays],c='olive',label="Baseline")
ax2.plot(np.arange(g_plotdays),bhavg.prinf[0:g_plotdays],'--',c='darkorange',label="$1/2*$ density")

ax.plot(np.arange(day2),b2avg.prtinf[0:day2],':',c='grey',linewidth=1.7,label="$2*$ density")
ax.plot(np.arange(g_plotdays),bavg.prtinf[0:g_plotdays],'',c='grey',linewidth=1.3,label="Baseline")
ax.plot(np.arange(g_plotdays),bhavg.prtinf[0:g_plotdays],'--',c='grey',linewidth=1.3,label="$1/2$ density")

ax.legend(bbox_to_anchor=(1, .7), loc='center right',title='Infected + Recovered')
ax2.legend(bbox_to_anchor=(1, .4), loc='center right',title='Infected \hspace{.02em} (right scale)')

fig.tight_layout()
plt.savefig(imagedir+'SIR-3densities.pdf')
plt.show()           


print('Sir DAYS to peak',np.where(max(bavg.prinf)==bavg.prinf)[0]
      ,np.where(max(b2avg.prinf)==b2avg.prinf)[0],
      np.where(max(bhavg.prinf)==bhavg.prinf)[0])

print('Peak active, 1*  ',np.round(max(bavg.prinf),2)
      ,'day ',np.where(max(bavg.prinf)==bavg.prinf)[0])
print('Peak active, 2*  ',np.round(max(b2avg.prinf),2)
      ,'day ',np.where(max(b2avg.prinf)==b2avg.prinf)[0])
print('Peak active, 1/2*',np.round(max(bhavg.prinf),2)
      ,'day ',np.where(max(bhavg.prinf)==bhavg.prinf)[0])

avgmaxinf = np.max(bavg.prtinf)
print('Steady state infected 1*  ', np.round(avgmaxinf,4))
avgmaxinf2 = np.max(b2avg.prtinf)
print('Steady state infected 2*  ', np.round(avgmaxinf2,4))
avgmaxinf3 = np.max(bhavg.prtinf)
print('Steady state infected 1/2*', np.round(avgmaxinf3,4))




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

days = 600
prtinf = pd.DataFrame(list(map(lambda x: x.prtinf[0:days],models)))
prinf = pd.DataFrame(list(map(lambda x: x.prinf[0:days],models)))

df1 = pd.concat([df,prtinf],axis=1)
df2 = pd.concat([df,prinf],axis=1)

dfav = df1.groupby('speed').mean()
df2av = df2.groupby('speed').mean()

dfn = dfav.to_numpy()[:,5:]
df2n = df2av.to_numpy()[:,5:]

fig,(ax1) = plt.subplots(1,1,figsize=(3.75,3.75))
ax2 = ax1.twinx()
ax2.spines['top'].set_visible(False)
ax1.spines['top'].set_visible(False)

ax1.set_xlabel('Days')
ax1.set_ylabel('Fraction of population')
ax1.set_ylim(0,1)
ax2.set_yticks(np.arange(0,0.3,.05))
ax2.set_ylim(0,0.15)

dfn[1,320:] = dfn[1,320]
lastday = 505

bavg.prtinf = np.pad(bavg.prtinf,(0,lastday))
bavg.prtinf[bavg.minlastday:]= bavg.prtinf[bavg.minlastday]
bavg.prinf = np.pad(bavg.prinf,(0,lastday))


dfn[1,272:] = dfn[1,272]

ax1.plot(np.arange(0,lastday),bavg.prtinf[0:lastday],c='grey',label='Baseline', linewidth=1)
ax1.plot(np.arange(0,lastday),dfn[1,4:lastday+4],':', linewidth=2, c='grey',label='Speed=20% baseline')
ax1.plot(np.arange(0,lastday),dfn[0,4:lastday+4],'--', c='grey',label='No movement')
ax2.plot(np.arange(0,lastday),bavg.prinf[0:lastday],'',linewidth=1.3, c='olive', label='Baseline')
ax2.plot(np.arange(0,lastday),df2n[1,4:lastday+4],':',linewidth=1.7, c='saddlebrown', label='Speed=20% baseline')
ax2.plot(np.arange(0,lastday),df2n[0,4:lastday+4],'--',linewidth=1.3, c='darkorange', label='No movement')

ax1.legend(bbox_to_anchor=(1, .61), loc='center right',title='Infected + Recovered')
ax2.legend(bbox_to_anchor=(1, .3),loc='center right',title='Infected \hspace{.02em} (right scale)')

fig.tight_layout()
plt.savefig(imagedir+'short-nomovement-rateslarge.pdf')
plt.show() 


#%% Heterogeneous density

file1 = gzip.open(outputdir+'basesim.pickle.gz','rb')
b = pickle.load(file1)
file1.close()
bavg = averageStats(b)

fileq = 'output/hetdens-lambda-1.pickle.gz'
file = gzip.open(fileq,'rb')
mod1 = pickle.load(file)
file.close()
b1avg = averageStats(mod1)

fig,ax = plt.subplots(figsize=(3.75,3.75))
ax2 = ax.twinx()
ax2.set_ylim(0,0.45)
ax.set_ylim(0,1)
ax2.set_yticks(np.arange(0,0.5,0.1))
g_plotdays = b1avg.minlastday
bavg.prtinf[bavg.minlastday:] = bavg.prtinf[bavg.minlastday]

### draw legend with empty data
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.set_xticks(np.arange(0,g_plotdays+10,20))
ax.set_xlabel('Days ')
ax.set_ylabel('Fraction of population')

ax.plot(np.arange(g_plotdays),bavg.prtinf[0:g_plotdays],c='grey', label="Baseline", linewidth=1.2)
ax.plot(np.arange(g_plotdays),b1avg.prtinf[0:g_plotdays],'--',c='grey',linewidth=1.2, label="Heterogeneous density")
ax2.plot(np.arange(g_plotdays),bavg.prinf[0:g_plotdays],'',c='olive',linewidth=1.3, label="Baseline")
ax2.plot(np.arange(g_plotdays),b1avg.prinf[0:g_plotdays],'--',c='darkorange',linewidth=1.7, label="Heterogeneous density")
ax.legend(bbox_to_anchor=(.95,.5), loc='center right',ncol=1,title='        Fraction of population \n(black lines: Infected )')

ax.legend(bbox_to_anchor=(1, .68), loc='center right',title='Infected + Recovered')
ax2.legend(bbox_to_anchor=(1, .44),loc='center right',title='Infected  (right scale)')


fig.tight_layout()
plt.savefig(imagedir+'hetdens1.pdf')
plt.show()    
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

fix,ax = plt.subplots(figsize=(4,3.5))
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
dfSIR = pd.DataFrame(aa.SIR,index=np.arange(aa.day+1),columns=['S','I','R'])
dfSIR['ILag'] = dfSIR['I'].shift()
dfSIR['ld'] = np.log(dfSIR['I']) - np.log(dfSIR['ILag'])
dfSIRld = np.array(dfSIR['ld'])

bb = SIRmodel(contRate*13.5, recRate, q_init=5, q_popsize=popSize,behModel={'type': 'Lones','phi': 0.01},)
dfSIRLon = pd.DataFrame(bb.SIR,index=np.arange(bb.day+1),columns=['S','I','R'])
dfSIRLon['ILag'] = dfSIRLon['I'].shift()
dfSIRLon['ld'] = np.log(dfSIRLon['I']) - np.log(dfSIRLon['ILag'])
dfSIRLonld = np.array(dfSIRLon['ld'])

bp = averageStats(benchplus)
bLon = averageStats(lones)
days = 86

fig,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(10.5,3.5))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Days ')

ax.plot(np.arange(days),dfSIRld[:days],'--',color='olive',label ='SIR')
ax.plot(np.arange(days),dfSIRLonld[:days],'-.',color='darkorange',label ='Behavioral SIR')
ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Spatial SIR')
ax.plot(np.arange(days),logLon.mean(axis=1)[:days],':',linewidth=1.9, color='darkorange',label ='Behavioral spatial SIR')

ax.legend(title='Infected growth rate')

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Days ')
ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,'--',color='olive',label ='SIR')
ax2.plot(np.arange(days),bb.SIR[:days,1]/25600,'-.',color='darkorange',label ='Behavioral SIR')
ax2.plot(np.arange(days),bp.prinf[:days],color='olive',label ='Spatial SIR')
ax2.plot(np.arange(days),bLon.prinf[:days],':',linewidth=1.9, color='darkorange',label ='Behavioral spatial SIR')
 
days=86
bp.prtinf = np.append(bp.prtinf,[0]*300)
bp.prtinf[bp.minlastday:] = bp.prtinf[bp.minlastday]
bLon.prtinf[bLon.minlastday:] = bLon.prtinf[bLon.minlastday]
aa.SIR = np.append(aa.SIR,[aa.SIR[aa.day,:]]*300,axis=0)
bb.SIR = np.append(bb.SIR,[bb.SIR[bb.day,:]]*300,axis=0)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.set_xlabel('Days ')
ax3.plot(np.arange(days),1-aa.SIR[:days,0]/25600,'--',color='olive',label ='SIR')
ax3.plot(np.arange(days),1-bb.SIR[:days,0]/25600,'-.',color='darkorange',label ='Behavioral SIR')
ax3.plot(np.arange(days),bp.prtinf[:days],color='olive',label ='Spatial SIR')
ax3.plot(np.arange(days),bLon.prtinf[:days],':',linewidth=1.9, color='darkorange',label ='Behavioral spatial SIR')
ax3.legend(title='Infected + Recovered')
ax.set_ylabel('Growth rate')
ax2.set_ylabel('Fraction of population')

ax2.legend(title='Infected ')

fig.tight_layout()    
plt.savefig(imagedir+'SIR_beh.pdf')

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

fix,ax = plt.subplots(figsize=(3.5,3.5))
#ax.set_ylim(0,0.4)
ax.set_yticks(np.arange(0,0.41,0.1))
days = np.min((bLon.minlastday, bb.fracNotScared.size))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Days ')
ax.set_ylabel('Fraction of population')

ax.plot(np.arange(days), 1-bLonL.fracNotScared[:days],'', color='olive', label = 'Behavioral Spatial SIR (local)')
ax.plot(np.arange(days), 1-bLon.fracNotScared[:days],':', linewidth=1.9, color='darkorange', label = 'Behavioral Spatial SIR')
#ax.plot(np.arange(days), 1-bb.fracNotScared[:days], '--', color='darkorange', label = 'Behavioral SIR')
ax.legend(bbox_to_anchor=(.5, .64), loc='center')
 
plt.savefig(imagedir+'SIR_beh_responses_local.pdf')

#%% Behavioral responses, comparison with benchmark

days = 112

fig,(ax2,ax3) = plt.subplots(1,2,figsize=(7,3.5))

ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xlabel('Days ')
ax2.plot(np.arange(days),bLonL.prinf[:days],color='olive',label ='Behavioral spatial SIR (local)')
ax2.plot(np.arange(days),bLon.prinf[:days],':',linewidth=1.9, color='darkorange',label ='Behavioral spatial SIR')
#ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,'--',color='olive',label ='SIR')
#ax2.plot(np.arange(days),bb.SIR[:days,1]/25600,'-.',color='darkorange',label ='Behavioral SIR')
 
days=112
bLon.prtinf = np.append(bLon.prtinf,[0]*300)
bLon.prtinf[bLon.minlastday:] = bLon.prtinf[bLon.minlastday]
bLonL.prtinf[bLonL.minlastday:] = bLonL.prtinf[bLonL.minlastday]
# aa.SIR = np.append(aa.SIR,[aa.SIR[aa.day,:]]*300,axis=0)
# bb.SIR = np.append(bb.SIR,[bb.SIR[bb.day,:]]*300,axis=0)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.set_xlabel('Days ')
ax3.plot(np.arange(days),bLonL.prtinf[:days],color='olive',label ='Behavioral spatial SIR (local)')
ax3.plot(np.arange(days),bLon.prtinf[:days],':',linewidth=1.9, color='darkorange',label ='Behavioral spatial SIR')
#ax3.plot(np.arange(days),1-aa.SIR[:days,0]/25600,'--',color='olive',label ='SIR')
#ax3.plot(np.arange(days),1-bb.SIR[:days,0]/25600,'-.',color='darkorange',label ='Behavioral SIR')
ax3.legend(title='Infected + Recovered')
#ax.set_ylabel('Growth rate')
ax2.set_ylabel('Fraction of population')
ax3.set_ylabel('Fraction of population')

ax2.legend(loc='center', title='Infected ')

fig.tight_layout()    
plt.savefig(imagedir+'SIR_beh_local.pdf')

 
