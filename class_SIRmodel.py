#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 23:28:30 2020

"""
import numpy as np
import matplotlib.pyplot as plt
import pickle


class SIRmodel():
    ''' 
        simulates a (behavioral) standard SIR model
    
    '''
    def __init__(self, beta, gamma, 
                 q_init=1, 
                 q_popsize=25600, 
                 behModel=False, # use behavioral adjustments?
                 p_shutr=[0,0],  # day of lockdown, fraction locked down
                 p_openr=[0,0],  # day of reopen, fraction put back (check??)
                 printOption=0
                 ) :
        self.beta = beta
        self.gamma = gamma
        self.behModel = behModel
        self.q_popsize = q_popsize
        self.printOption = printOption
        self.q_init = q_init
        self.p_shutr = p_shutr
        self.p_openr = p_openr
        
        # array of susceptibles, infected, recovered, lockeddown
        self.SIR = np.array([[self.q_popsize-self.q_init, self.q_init, 0, 0]])
        self.fracNotScared = np.array([1])
        self.shutdownMode = False
        self.simulate()
        self.computeStats()
        
    def simulate(self):
        
        day = 0
        
        # loop while susceptibles are positive and infected are less than popsize (why?)
        while np.round(self.SIR[day,1])>0 and np.round(self.SIR[day,1])<self.q_popsize:
            
            day += 1
            if self.behModel != False:
                beta = self.beta*self.scaredycats(day)
            else :
                beta = self.beta
            
            dI = beta*self.SIR[day-1, 0]*self.SIR[day-1, 1]/self.q_popsize
            dR = self.gamma*self.SIR[day-1, 1]
            S = self.SIR[day-1, 0] - dI
            I = self.SIR[day-1, 1] + dI - dR
            R = self.SIR[day-1, 2] + dR
            L = self.SIR[day-1, 3] #locked-down
            
            #lockedout or reopen
            if day==self.p_shutr[0]: 
                
                # correct for the relative size of S
                L = S*self.p_shutr[1] * S/self.q_popsize
                S = S-L
            if day==self.p_openr[0]:
                S = S + L*self.p_openr[1]
                L = L*(1-self.p_openr[1])
            
            #correction just in case
            if S<0:
                I = I-S
                S = 0
                
            self.SIR = np.row_stack((self.SIR, [S, I, R, L]))
            
            if I<0:
                break
            
            if self.printOption>=1:
                print(day,np.round(self.SIR[day,:],0))
                     
        self.day = day

    def plot(self,maxday,log=False):
        fig,ax = plt.subplots()

        if log : 
            plotData = np.log(self.SIR)
            #ax.set_ylim(-2,12)
        else :
            plotData = self.SIR
            #ax.set_ylim(0,self.q_popsize)
        ax.plot(np.arange(maxday),plotData[:maxday,0],c='saddlebrown',label='S')
        ax.plot(np.arange(maxday),plotData[:maxday,1],c='darkorange',label='I')
        ax.plot(np.arange(maxday),plotData[:maxday,2],c='olive',label='R')
        ax.plot(np.arange(maxday),plotData[:maxday,3],':',c='black',label='L')
        ax.legend()
        plt.show()
        
    def scaredycats(self,day) :
        '''
        Detects people scared to go outside so that infections() can 
        put them in a land far and far away

        '''
        par = self.behModel
        
        # figure out how many are infected and compute how many are scared      
        if par['type'] == 'Lones':
            
            fY = self.SIR[day-1,1] / self.q_popsize
            if fY <= par['phi'] :
                reducedContagion = 1.
            else:
                reducedContagion = (par['phi']/fY)**0.12        
                
        elif par['type'] == 'Jesus':
            
            b0 = par['beta0']
            bstar = par['betastar']
            lambd = par['lambda']
            reducedContagion = b0 * np.exp(-par['lambda']*day) + bstar * (1-np.exp(-lambd*day))
        
        elif par['type'] == 'Lockdown':
            p_shutr = par['p_shutr']
            p_openr = par['p_openr']
            p_shutOnce = par['p_shutOnce']
            reducedContagion = 1
        
            #### WARNING -2 vs -1
            fY = self.SIR[day-1,1] / self.q_popsize
            fYprev = self.SIR[day-2,1] / self.q_popsize
        
            # shut either when allowed multiple times, or when never reopened
            if (p_shutOnce == False or self.reopened == False) :
        
                # if going over the shutdown threshold from below: shutdown
                if (fY>=p_shutr[0]) & (fYprev<p_shutr[0]) & (self.shutdownMode==False):
                    self.shutdownMode=True
                
            # if going over the reopen threshold from above: reopen
            if (fY<p_openr[0]) and (fYprev>=p_openr[0]):
                if self.shutdownMode==True:
                    self.shutdownMode = False
                    self.reopened = True
        
            ## now define masks of who is staying home
            # shut down (p_shutr[1]%) of city if asymptomatics > p_shutr[0]% 
            if self.shutdownMode :
                if p_shutr[1]==1 :
                    reducedContagion = 0
                else:
                    reducedContagion = (1-p_shutr[1])
             
        self.fracNotScared = np.append(self.fracNotScared,reducedContagion)

        return reducedContagion

    def computeStats(self):
        
        self.prtinf = 1-self.SIR[:,0]/self.q_popsize
        self.prinf = self.SIR[:,1]/self.q_popsize
        self.R0 = self.beta/self.gamma
        self.Ipeak = -1*np.log(self.R0*(self.q_popsize-self.q_init)/self.q_popsize)/self.R0  - 1/self.R0 + 1


#%% testing the class
if __name__ == '__main__':
    
    baseline = 25600
    cluster = 10
    beta = .15 #.054*13.5,
    delta = .02 #1/6.5
    bavg = SIRmodel(beta, delta ,q_popsize=baseline,
                    q_init=cluster,
                    )
    bavg.plot(bavg.day)
    
    b2 = SIRmodel(beta, delta, q_popsize=baseline,
                    q_init=cluster, p_shutr=[20,.3], p_openr=[999,1]
                    )
    b2.plot(b2.day)


    
#%% double cont rate and 2* density in benchmark plus
if __name__ == '__main__':
    
    from class_spatialModels import simulatePool
    kwargs=benchkwargs.copy()
    kwargs['q_citysize'] = 1*np.sqrt(2)
    kwargs['p_probc'] = 2* np.array([[0,0.038130388,0.038130388,0,0]])

    
    nboots = 10
    nprocs = 5

    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandBigger = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    list(map(lambda x: delattr(x,'pos'),allRandBigger))
    file = open('output/benchmark_noY_bigger_morecont.pickle','wb')
    pickle.dump(allRandBigger,file)
    
    kwargs=benchkwargs.copy()
    kwargs['q_citysize'] = 1
    kwargs['p_probc'] = np.array([[0,0.038130388,0.038130388,0,0]])
    
    nboots = 10
    nprocs = 10

    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandBigger = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    list(map(lambda x: delattr(x,'pos'),allRandBigger))
    file = open('output/benchmark_noY.pickle','wb')
    pickle.dump(allRandBigger,file)


    kwargs=benchkwargs.copy()
    kwargs['q_citysize'] = 1*np.sqrt(6)
    kwargs['p_probc'] = 6* np.array([[0,0.038130388,0.038130388,0,0]])

    
    nboots = 10
    nprocs = 5

    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandBigger = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    list(map(lambda x: delattr(x,'pos'),allRandBigger))
    file = open('output/benchmark_noY_huge.pickle','wb')
    pickle.dump(allRandBigger,file)
    


#%%
if __name__ == '__main__':
    import pandas as pd
    

    file = open('../base25.6K/output/zbasesim.pickle','rb')
    bench = pickle.load(file)
    file.close()

    file = open('output/allrandmodels.pickle','rb')
    rand = pickle.load(file)
    file.close()

    file = open('output/allrandmodelsbig.pickle','rb')
    randb = pickle.load(file)
    file.close()

    file = open('output/benchmark_noY.pickle','rb')
    benchplus = pickle.load(file)
    file.close()
    
    file = open('output/benchmark_noY_bigger_morecont.pickle','rb')
    benchplusbig = pickle.load(file)
    file.close()
    
    file = open('output/benchmark_noY_huge.pickle','rb')
    benchplushuge  = pickle.load(file)
    file.close()
    

    # true benchmark    
    dfbench = pd.DataFrame()
    logsbench = pd.DataFrame()    
    for idx,m in enumerate(bench):
        dfbench[idx] = np.sum(m.nstates[:,1:3],axis=1)
        dfbench[str(idx)+'L1'] = dfbench[idx].shift()
        logsbench[str(idx)] = np.log(dfbench[idx]) - np.log(dfbench[str(idx)+'L1'])

    # random matching in geo
    df = pd.DataFrame()
    logs = pd.DataFrame()
    for idx,m in enumerate(rand) :
        df[idx] = np.sum(m.nstates[:,1:3],axis=1)
        df[str(idx)+'L1'] = df[idx].shift()
        logs[str(idx)] = np.log(df[idx]) - np.log(df[str(idx)+'L1'])

    # random matching in geo bigger (6*)
    dfr6 = pd.DataFrame()
    logsr6 = pd.DataFrame()
    for idx,m in enumerate(randb) :
        dfr6[idx] = np.sum(m.nstates[:,1:3],axis=1)
        dfr6[str(idx)+'L1'] = dfr6[idx].shift()
        logsr6[str(idx)] = np.log(dfr6[idx]) - np.log(dfr6[str(idx)+'L1'])

    # benchmark plus (Y contagious)   
    dfb = pd.DataFrame()
    logsb = pd.DataFrame()    
    for idx,m in enumerate(benchplus) :
        dfb[idx] = np.sum(m.nstates[:,1:3],axis=1)
        dfb[str(idx)+'L1'] = dfb[idx].shift()
        logsb[str(idx)] = np.log(dfb[idx]) - np.log(dfb[str(idx)+'L1'])        
        
    # benchmark plus bigger 
    df2 = pd.DataFrame()
    logs2 = pd.DataFrame()    
    for idx,m in enumerate(benchplusbig) :
        df2[idx] = np.sum(m.nstates[:,1:3],axis=1)
        df2[str(idx)+'L1'] = df2[idx].shift()
        logs2[str(idx)] = np.log(df2[idx]) - np.log(df2[str(idx)+'L1'])

    # benchmark plus huge 
    df6 = pd.DataFrame()
    logs6 = pd.DataFrame()    
    for idx,m in enumerate(benchplushuge) :
        df6[idx] = np.sum(m.nstates[:,1:3],axis=1)
        df6[str(idx)+'L1'] = df6[idx].shift()
        logs6[str(idx)] = np.log(df6[idx]) - np.log(df6[str(idx)+'L1'])



    aa = SIRmodel(.03813*13.5,0.05,q_popsize=25600)
    dfSIR = pd.DataFrame(aa.SIR,index=np.arange(aa.day+1),columns=['S','I','R'])
    dfSIR['ILag'] = dfSIR['I'].shift()
    dfSIR['ld'] = np.log(dfSIR['I']) - np.log(dfSIR['ILag'])
    dfSIRld = np.array(dfSIR['ld'])

    days = 150

    from class_averageStats import averageStats

    ravg = averageStats(rand)
    ravgb = averageStats(randb)
    bp = averageStats(benchplus)
    b6 = averageStats(benchplushuge)

    #first plot: compare SIR with benchmark plus
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(10,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    #ax.plot(np.arange(days),logsbench.mean(axis=1)[:days],color='saddlebrown',label ='True benchmark')
    #ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='darkorange',label ='Benchmark plus')
    #ax.plot(np.arange(days),logs2.mean(axis=1)[:days],color='olive',label ='B. plus 1/2 dense 2*contagious')
    #ax.plot(np.arange(days),logs6.mean(axis=1)[:days],color='pink',label ='B. plus 1/6 dense 6*contagious')
    ax.plot(np.arange(days),dfSIRld[:days],color='saddlebrown',label ='(i) SIR')
    ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='darkorange',label ='(ii) Spatial SIR model')

    ax.legend(title='Infection growth rates')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,color='saddlebrown',label ='SIR')
    ax2.plot(np.arange(days),bp.prinf[:days],color='darkorange',label ='(ii) Spatial SIR model')
 
    ax2.legend(title='Active cases')
    fig.tight_layout()
    
    plt.savefig('../Covid BisinMoro/write/images/SIR_randomGEO_comp1.png')

    # second plot: add bench plus with less density
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(10,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    ax.plot(np.arange(days),dfSIRld[:days],color='saddlebrown',label ='(i) SIR')
    ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='darkorange',label ='(ii) Spatial SIR model')
    ax.plot(np.arange(days),logs6.mean(axis=1)[:days],color='olive',label ='(iii) 1/6 dense 6*contagious')


    ax.legend(title='Infection growth rates')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,color='saddlebrown',label ='SIR')
    ax2.plot(np.arange(days),bp.prinf[:days],color='darkorange',label ='(ii) Spatial SIR model')
    ax2.plot(np.arange(days),b6.prinf[:days],color='olive',label ='(iii) 1/6 dense 6 *contagious')
    ax2.legend(title='Total and active cases')
    fig.tight_layout()
    
    plt.savefig('../Covid BisinMoro/write/images/SIR_randomGEO_comp2.png')

    #third plot: add the effect of non random matching
    days = 150
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(10,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    ax.set_ylim(-0.15,0.4)
    ax.plot(np.arange(days),dfSIRld[:days],color='saddlebrown',label ='(i) SIR')
    ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='darkorange',label ='(ii) Geospatial model')
    ax.plot(np.arange(days),logs.mean(axis=1)[:days],color='lightsteelblue',label ='(iii) Geo with random positions')
 
    ax.legend(title='Infection growth rates')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(days),aa.SIR[:days,1]/25600,color='saddlebrown',label ='(i) SIR')
    ax2.plot(np.arange(days),bp.prinf[:days],color='darkorange',label ='(ii) Geospatial model')
    ax2.plot(np.arange(days),ravg.prinf[:days],color='lightsteelblue',label ='(iii) Geo with random positions')   

    ax2.legend(title='Active cases')
    fig.tight_layout()
    plt.savefig('../Covid BisinMoro/write/images/SIR_random_noMatching_comp1.png')


    #fourth plot: compare benchmark with one with less density
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(10,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    ax.set_ylim(-0.15,0.4)
    ax.plot(np.arange(days),dfSIRld[:days],color='saddlebrown',label ='(i) SIR')
    ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='darkorange',label ='(ii) Spatial SIR model')
    ax.plot(np.arange(days),logs6.mean(axis=1)[:days],color='olive',label ='(iii) 1/6 dense 6*contagious')
    ax.plot(np.arange(days),logs.mean(axis=1)[:days],color='lightsteelblue',label ='(iv) Sp. SIR with random positions')
 
    ravg = averageStats(rand)
    ravgb = averageStats(randb)
    bp = averageStats(benchplus)
    b6 = averageStats(benchplushuge)
    ax.legend(title='Infection growth rates')
    days=250
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(aa.day),aa.SIR[:aa.day,1]/25600,color='saddlebrown',label ='(i) SIR')
    ax2.plot(np.arange(days),bp.prinf[:days],color='darkorange',label ='(ii) Spatial SIR model')
    ax2.plot(np.arange(days),b6.prinf[:days],color='olive',label ='(iii) 1/6 dense 6 *contagious')
    ax2.plot(np.arange(days),ravg.prinf[:days],color='lightsteelblue',label ='(iv) Sp. SIR with random positions')   
    
    ax2.legend(title='Active cases')
    fig.tight_layout()
 

    plt.savefig('../Covid BisinMoro/write/images/SIR_random_noMatching_comp2.png')
    plt.show()
#%% fifth plot: compare benchmark with one with less density
    days = 150
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(10,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    ax.plot(np.arange(days),dfSIRld[:days],color='saddlebrown',label ='(i) SIR')
    ax.plot(np.arange(days),logsr6.mean(axis=1)[:days],color='lightsteelblue',label ='(iii) Sp. SIR with random positions')
    ax.plot(np.arange(days),logs.mean(axis=1)[:days],color='pink',label ='(iv) Sp. SIR, random pos, 1/6 dense')
 

    ravg = averageStats(rand)
    ravgb = averageStats(randb)
    bp = averageStats(benchplus)
    b6 = averageStats(benchplushuge)
    ax.legend(title='Infection growth rates')
    days=250
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(aa.day),aa.SIR[:aa.day,1]/25600,color='saddlebrown',label ='(i) SIR')
    ax2.plot(np.arange(days),ravg.prinf[:days],color='lightsteelblue',label ='(iii) Sp. SIR, random pos, 1/6 dense')   
    ax2.plot(np.arange(days),ravgb.prinf[:days],color='pink',label ='(iv) Sp. SIR with random positions')   
    
    ax2.legend(title='Active cases')
    fig.tight_layout()
 

    plt.savefig('../Covid BisinMoro/write/images/SIR_random_matching_comp2.png')
    plt.show()

#%% for paper

    # first plot: add bench plus with less density
    days = 150
    fig,(ax,ax2) = plt.subplots(1,2,figsize=(10,5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlabel('Days ')
    ax.plot(np.arange(days),logsb.mean(axis=1)[:days],color='olive',label ='Benchmark')
    ax.plot(np.arange(days),logs6.mean(axis=1)[:days],'--', color='darkorange',label ='(ii) 1/6*density, 6*contagious')


    ax.legend(title='Infection growth rates')
    days = 250

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(days),bp.prinf[:days],color='olive',label ='Benchmark')
    ax2.plot(np.arange(days),b6.prinf[:days],'--', color='darkorange',label ='1/6*density, 6*contagious')
    ax2.legend(title='Active cases')
    fig.tight_layout()
    
    plt.savefig('../Covid BisinMoro/write/images/density_contagion1.png')

