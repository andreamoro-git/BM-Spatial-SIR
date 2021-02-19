#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 18:28:27 2020

@author: moroa
"""

from class_averageStats import averageStats
import numpy as np
import pickle,gzip
from class_simul import contagionModel
import matplotlib.pyplot as plt

class simul_policy(contagionModel) :
    ''' class simulating a lockdown
        shutr: shut down shutr[1]% of people i if fraction in state I >shutr[0]
    '''
    
    def __init__(self, p_shutr, p_openr,
                 shutOnce = False,
                 *args,**kwargs,) :
        contagionModel.__init__(self,*args,**kwargs)
        self.p_shutr = p_shutr #this is to avoid hitting exactly the threshold
        self.p_openr = p_openr
        self.shutOnce = shutOnce
        self.reopened = False
        
        # draws shutdown-specific ids to lock down the same workers / firms
        self.shutid = np.arange(self.q_popsize)
        np.random.shuffle(self.shutid)
        
        # define vecors of people shut down by govt
        self.shut = np.zeros(self.q_popsize,dtype=bool)
        self.shutdownMode = False
  
    def shutdown(self,day) :
        # detect people shut down by government
        
        self.shut = np.zeros(self.q_popsize,dtype=bool)
        
        #### WARNING -2 vs -1
        fY = self.prstates[day,self.I]+self.prstates[day,self.Y]
        fYprev = self.prstates[day-2,self.I]+self.prstates[day-2,self.Y]
        
        # shut either when allowed multiple times, or when never reopened
        if (self.shutOnce == False or self.reopened == False) :
        
            # if going over the shutdown threshold from below: shutdown
            if (fY>=self.p_shutr[0]) & (fYprev<self.p_shutr[0]) & (self.shutdownMode==False):
                self.shutdownMode=True
                
        # if going over the reopen threshold from above: reopen
        if (fY<self.p_openr[0]) and (fYprev>=self.p_openr[0]):
            if self.shutdownMode==True:
                self.shutdownMode = False
                self.reopened = True
        
        ## now define masks of who is staying home
        # shut down (p_shutr[1]%) of city if asymptomatics > p_shutr[0]% 
        if self.shutdownMode :
            if self.p_shutr[1]==1 :
                self.shut[:] = True
            else:
                lockedDown = self.p_shutr[1]*self.q_popsize
                self.shut = (self.shutid<lockedDown)
                    
        if self.q_printOption>=2:
            print('\t\t Shutdownmodes: ',self.shutdownMode)
            
        return 
 
    def aDayInTheLife(self,day):
    
        # figure out who is scared and save their positions
        self.shutdown(day)
        self.savepos = np.copy(self.pos)
        
        # move scared out of the box
        self.pos[self.shut,1] =  self.q_citysize*(2*(self.id[self.shut]+1))
        
        # generate infections, deaths and recoveries
        self.infections(day)
        self.symptoms(day)
        self.deaths(day)
        self.recoveries(day)      ,
                    
        # return people to their position before moving them arund 
        self.pos[self.shut,1] =  self.savepos[self.shut,1]

        # move people around
        self.pos = self.movepeople(self.pos)
        self.lastday = day

          
       
def simulatePool(kwargs) :
    ''' function that simulates the model once, to be used in conjunction
        with multiprocessing.Pool to get multiple replications quickly 
    '''
    m = simul_policy(**kwargs)
    m.simulateOnce()
    m.summaryStats(m.lastday,nondef=kwargs,savefile='output/results.txt',clear='a')

    return m

#%%  simulate 10-5 lockdown policy from base model
if __name__ == "__main__":

    
    from multiprocessing import Pool
    
    benchkwargs = {"q_seed"      : 2443,
        "p_probc"     : [[0,0.038130388,0.038130388,0,0]], #prob of contagion by type and state
        "q_popsize"   : 25600,
        "p_infradius" : 0.013016,
        "p_avgstep"   : 0.033993284,
        "p_stdstep"   : 0,
        "q_printOption" : 0.6,
        'q_days'      : 500,
        'p_cluster'   : 'cluster',
        'p_shutr'     : [0.10,0.30],
        'p_openr'     : [0.05,1],
    }
    
    # a = simul_policy(**benchkwargs)
    # a.simulateOnce()
    
    kwargs = benchkwargs.copy()
    nboots = 4
    nprocs = 4

    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    list(map(lambda x: delattr(x,'pos'),allRandmodels))
    file = gzip.open('output/basePolicy-10-5-30pc.pickle.gz','wb')
    pickle.dump(allRandmodels,file)

#%% Figures 10-5 base model 30pc
    import pandas as pd
    from class_SIRmodel import SIRmodel
    
    
    # benchmark   
    file = open('output/zbasesim.pickle','rb')
    bench = pickle.load(file)
    file.close()
    dfb = pd.DataFrame()
    logsb = pd.DataFrame()    
    for idx,m in enumerate(bench) :
        dfb[idx] = np.sum(m.nstates[:,1:3],axis=1)
        dfb[str(idx)+'L1'] = dfb[idx].shift()
        logsb[str(idx)] = np.log(dfb[idx]) - np.log(dfb[str(idx)+'L1'])        

    file = gzip.open('output/basePolicy-10-5-30pc.pickle.gz','rb')
    pol = pickle.load(file)
    file.close()
    dfpol = pd.DataFrame()
    logpol = pd.DataFrame()    
    for idx,m in enumerate(pol) :
        dfpol[idx] = np.sum(m.nstates[:,1:3],axis=1)
        dfpol[str(idx)+'L1'] = dfpol[idx].shift()
        logpol[str(idx)] = np.log(dfpol[idx]) - np.log(dfpol[str(idx)+'L1'])        
    
    aa = SIRmodel(.03813*13.5,0.05,q_popsize=25600,q_init=30)
    dfSIR = pd.DataFrame(aa.SIR,index=np.arange(aa.day+1),columns=['S','I','R'])
    dfSIR['ILag'] = dfSIR['I'].shift()
    dfSIR['ld'] = np.log(dfSIR['I']) - np.log(dfSIR['ILag'])
    dfSIRld = np.array(dfSIR['ld'])
    
    bb = SIRmodel(.03813*13.5,0.05,q_popsize=25600,
                    q_init=30,
                    behModel = {'type':'Lockdown','p_shutr':[0.1,.3],'p_openr':[0.05,1],'p_shutOnce': False}
                    )
    dfSIRp = pd.DataFrame(bb.SIR,index=np.arange(bb.day+1),columns=['S','I','R'])
    dfSIRp['ILag'] = dfSIRp['I'].shift()
    dfSIRp['ld'] = np.log(dfSIRp['I']) - np.log(dfSIRp['ILag'])
    dfSIRpld = np.array(dfSIRp['ld'])
 
        
    bp = averageStats(bench)
    bpol = averageStats(pol)
    days = 180
    
    #first plot: compare SIR with benchmark 
    fig,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')

    ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Spatial-SIR')
    ax.plot(np.arange(days),logpol.mean(axis=1)[:days],':',linewidth=1.9, color='darkorange',label ='Spatial-SIR with lockdown')
    ax.plot(np.arange(days),dfSIRld[:days],'--',color='olive',label ='SIR')
    ax.plot(np.arange(days),dfSIRpld[:days],'-.',color='darkOrange',label ='SIR with lockdown')

    ax.legend(title='Infection growth rates')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(days),bp.prinf[:days],color='olive',label ='Spatial-SIR')
    ax2.plot(np.arange(days),bpol.prinf[:days],':',linewidth=1.9, color='darkorange',label ='Spatial-SIR with lockdown')
    ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,'--',color='olive',label ='SIR')
    ax2.plot(np.arange(days),bb.SIR[:days,1]/25600,'-.',color='darkOrange',label ='SIR with lockdown')
    # ax2.axhline(pol[0].p_shutr[0],color='olive',ls='dotted')
    # ax2.axhline(pol[0].p_openr[0],color='olive',ls='-.')

    days=230
    bp.prtinf = np.append(bp.prtinf,[0]*300)
    bp.prtinf[140:] = bp.prtinf[140]
    aa.SIR = np.append(aa.SIR,[aa.SIR[aa.day,:]]*300,axis=0)

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.set_xlabel('Days ')
    ax3.plot(np.arange(days),bp.prtinf[:days],color='olive',label ='Spatial-SIR')
    ax3.plot(np.arange(days),bpol.prtinf[:days],':',linewidth=1.9, color='darkorange',label ='Spatial-SIR with lockdown')
    ax3.plot(np.arange(days),1-aa.SIR[:days,0]/25600,'--',color='olive',label ='SIR')
    ax3.plot(np.arange(days),1-bb.SIR[:days,0]/25600,'-.',color='darkOrange',label ='SIR with lockdown')
    ax3.legend(title='Total cases',loc='center right')

    ax2.legend(title='Active cases')

    fig.tight_layout()    
    plt.savefig('../Covid BisinMoro/write/images/base_policy-10-5-30pc-comp.png')

#%% Compute 10-5-50pc

    benchkwargs = {"q_seed"      : 2443,
        "p_probc"     : [[0,0.038130388,0.038130388,0,0]], #prob of contagion by type and state
        "q_popsize"   : 25600,
        "p_infradius" : 0.013016,
        "p_avgstep"   : 0.033993284,
        "p_stdstep"   : 0,
        "q_printOption" : 0.6,
        'q_days'      : 500,
        'p_cluster'   : 'cluster',
        'p_shutr'     : [0.10,0.50],
        'p_openr'     : [0.05,1],
    }
    
    # a = simul_policy(**benchkwargs)
    # a.simulateOnce()
    
    kwargs = benchkwargs.copy()
    nboots = 4
    nprocs = 4

    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    list(map(lambda x: delattr(x,'pos'),allRandmodels))
    file = gzip.open('output/basePolicy-10-5-50pc.pickle.gz','wb')
    pickle.dump(allRandmodels,file)

#%% Figures 10-5 base model 30pc
    import pandas as pd
    from class_SIRmodel import SIRmodel
    
    
    # benchmark   
    file = open('output/zbasesim.pickle','rb')
    bench = pickle.load(file)
    file.close()
    dfb = pd.DataFrame()
    logsb = pd.DataFrame()    
    for idx,m in enumerate(bench) :
        dfb[idx] = np.sum(m.nstates[:,1:3],axis=1)
        dfb[str(idx)+'L1'] = dfb[idx].shift()
        logsb[str(idx)] = np.log(dfb[idx]) - np.log(dfb[str(idx)+'L1'])        

    file = gzip.open('output/basePolicy-10-5-50pc.pickle.gz','rb')
    pol = pickle.load(file)
    file.close()
    dfpol = pd.DataFrame()
    logpol = pd.DataFrame()    
    for idx,m in enumerate(pol) :
        dfpol[idx] = np.sum(m.nstates[:,1:3],axis=1)
        dfpol[str(idx)+'L1'] = dfpol[idx].shift()
        logpol[str(idx)] = np.log(dfpol[idx]) - np.log(dfpol[str(idx)+'L1'])        
    
    aa = SIRmodel(.03813*13.5,0.05,q_popsize=25600,q_init=30)
    dfSIR = pd.DataFrame(aa.SIR,index=np.arange(aa.day+1),columns=['S','I','R'])
    dfSIR['ILag'] = dfSIR['I'].shift()
    dfSIR['ld'] = np.log(dfSIR['I']) - np.log(dfSIR['ILag'])
    dfSIRld = np.array(dfSIR['ld'])
    
    bb = SIRmodel(.03813*13.5,0.05,q_popsize=25600,
                    q_init=30,
                    behModel = {'type':'Lockdown','p_shutr':[0.1,.5],'p_openr':[0.05,1],'p_shutOnce': False}
                    )
    dfSIRp = pd.DataFrame(bb.SIR,index=np.arange(bb.day+1),columns=['S','I','R'])
    dfSIRp['ILag'] = dfSIRp['I'].shift()
    dfSIRp['ld'] = np.log(dfSIRp['I']) - np.log(dfSIRp['ILag'])
    dfSIRpld = np.array(dfSIRp['ld'])
 
        
    bp = averageStats(bench)
    bpol = averageStats(pol)
    days = 180
    
    #first plot: compare SIR with benchmark 
    fig,(ax,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')

    ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Spatial-SIR')
    ax.plot(np.arange(days),logpol.mean(axis=1)[:days],':',linewidth=1.9, color='darkorange',label ='Spatial-SIR with lockdown')
    ax.plot(np.arange(days),dfSIRld[:days],'--',color='olive',label ='SIR')
    ax.plot(np.arange(days),dfSIRpld[:days],'-.',color='darkOrange',label ='SIR with lockdown')

    ax.legend(title='Infection growth rates')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(days),bp.prinf[:days],color='olive',label ='Spatial-SIR')
    ax2.plot(np.arange(days),bpol.prinf[:days],':',linewidth=1.9, color='darkorange',label ='Spatial-SIR with lockdown')
    ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,'--',color='olive',label ='SIR')
    ax2.plot(np.arange(days),bb.SIR[:days,1]/25600,'-.',color='darkOrange',label ='SIR with lockdown')
    # ax2.axhline(pol[0].p_shutr[0],color='black',ls='dotted',linewidth=1)
    # ax2.axhline(pol[0].p_openr[0],color='black',ls='-.',linewidth=.8)

    days=230
    bp.prtinf = np.append(bp.prtinf,[0]*300)
    bp.prtinf[140:] = bp.prtinf[140]
    aa.SIR = np.append(aa.SIR,[aa.SIR[aa.day,:]]*300,axis=0)

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.set_xlabel('Days ')
    ax3.plot(np.arange(days),bp.prtinf[:days],color='olive',label ='Spatial-SIR')
    ax3.plot(np.arange(days),bpol.prtinf[:days],':',linewidth=1.9, color='darkorange',label ='Spatial-SIR with lockdown')
    ax3.plot(np.arange(days),1-aa.SIR[:days,0]/25600,'--',color='olive',label ='SIR')
    ax3.plot(np.arange(days),1-bb.SIR[:days,0]/25600,'-.',color='darkOrange',label ='SIR with lockdown')
    ax3.legend(title='Total cases', loc='center right')

    ax2.legend(title='Active cases')

    fig.tight_layout()    
    plt.savefig('../Covid BisinMoro/write/images/base_policy-10-5-50pc-comp.png')

#%% another figure option
    # benchmark   
    file = open('output/zbasesim.pickle','rb')
    bench = pickle.load(file)
    file.close()

    file = gzip.open('output/basePolicy-10-5-50pc.pickle.gz','rb')
    pol = pickle.load(file)
    file.close()

    file = gzip.open('output/basepolicy-10-5-30pc.pickle.gz','rb')
    pol30 = pickle.load(file)
    file.close()
    
    bp = averageStats(bench)
    bpol = averageStats(pol)
    bpol30 = averageStats(pol30)

    
    aa = SIRmodel(.03813*13.5,0.05,q_popsize=25600,q_init=30)
    
    bb = SIRmodel(.03813*13.5,0.05,q_popsize=25600,
                    q_init=30,
                    behModel = {'type':'Lockdown','p_shutr':[0.1,.3],'p_openr':[0.05,1],'p_shutOnce': False}
                    )

    cc = SIRmodel(.03813*13.5,0.05,q_popsize=25600,
                    q_init=30,
                    behModel = {'type':'Lockdown','p_shutr':[0.1,.5],'p_openr':[0.05,1],'p_shutOnce': False}
                    )
    
    #first plot: compare SIR with benchmark 
    days=250
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(10,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    ax.set_ylabel('Percent infected')
    ax.plot(np.arange(days),bp.prinf[:days],color='olive',label ='Benchmark')
    ax.plot(np.arange(days),bpol30.prinf[:days],':',linewidth=1.9, color='saddlebrown',label ='30% lockdown')
    ax.plot(np.arange(days),bpol.prinf[:days],'--', color='darkorange',label ='50% lockdown')

    ax.legend(title='Spatial-SIR, Active cases')

    days=130
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,'',color='olive',label ='SIR')
    ax2.plot(np.arange(days),bb.SIR[:days,1]/25600,':',linewidth=1.9,color='saddlebrown',label ='SIR, 30% lockdown')
    ax2.plot(np.arange(days),cc.SIR[:days,1]/25600,'--',color='darkorange',label ='SIR, 50% lockdown')
    # ax2.axhline(pol[0].p_shutr[0],color='black',ls='dotted',linewidth=1)
    # ax2.axhline(pol[0].p_openr[0],color='black',ls='-.',linewidth=.8)
    ax2.legend(title='SIR model, Active cases')
    
    fig.tight_layout()
    plt.savefig('../Covid BisinMoro/write/images/base_policy-10-5-active-.png')

