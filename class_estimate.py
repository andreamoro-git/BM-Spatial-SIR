#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 08:37:42 2020

@author: moroa
"""
from class_spatialModels import spatialSAYDR
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool
from datetime import datetime
from copy import deepcopy

# these two functions are needed for my constrained optimization technique
def squeeze(x,bottom,top) :
    if (x<=700 and x>=-700) :    # stupid way to avoid overflow errors
        squeeze=(top-bottom)*1.2**x/(1.2**x+1)+bottom
    elif (x>700) :
        squeeze=top
    else :
        squeeze=bottom
    return squeeze

def blowup(x,bottom,top) :
    if (x<bottom+1.E-10) : # again to avoid overflow errors
         blowup=-700.
    elif (x>top-1.E-10) :
        blowup=700
    else :
        blowup=np.log((x-bottom)/(top-x))/np.log(1.2)
    return blowup


class estimatePars(spatialSAYDR) :
    def __init__(self,
                 benchkwargs,
                 est_days = 34,
                 nboots=10,
                 print_screen=1,
                 parallelize=3,     # number of processes
                 shortcut = True,
                 ):
        
#        spatialSAYDR.__init__(self,*args,**kwargs)
        self.nboots = nboots
        self.print_screen = print_screen
        self.parallelize = parallelize
        self.shortcut = shortcut
        self.minimum = float('inf')

        self.df = pd.DataFrame()
        self.benchkwargs = benchkwargs
        
        self.est_days = est_days
        self.mldr = np.zeros(self.est_days)
        milan = open('input/drlombardia.txt','r')
        for idx,line in enumerate(milan) :
            if idx < self.est_days: 
                self.mldr[idx] = line
        
    def computeBoot(self,kwargs):

        m = spatialSAYDR(**kwargs)
        m.computeStats(0,m.q_printOption)
        for day in range(1,self.est_days):

            cs = m.computeStats(day,m.q_printOption)
            if m.q_printOption >=3 :
                m.drawRates(day)
            if m.q_printOption >= 1: 
                print(cs)
            m.aDayInTheLife(day)    
            
            if np.sum(m.state==m.I)+np.sum(m.state==m.Y) <=15:
                break
            
            # hack to speed up some wrong guesses
            # not entirely working if different boots have wildly different
            # predictions (TODO) but there's certainly something wrong with
            # the guess anyways
            if self.shortcut :
                if day==8 :
                    df = pd.DataFrame(m.igrowth[2:])        
                    ma = df.rolling(3,win_type='triang').mean()
                    ma = ma[2:]
                    loss = ((ma[0:day-2]-self.mldr[0:day-2][0])**2).sum()
                    if loss[0] > 5 * self.minimum :
                        return m.igrowth[2:]
                if day==15 :
                    df = pd.DataFrame(m.igrowth[2:])        
                    ma = df.rolling(3,win_type='triang').mean()
                    ma = ma[2:]
                    loss = ((ma[0:day-2]-self.mldr[0:day-2][0])**2).sum()
                    if loss[0] > 2 * self.minimum :
                        return m.igrowth[2:]

        # m.drawRates(day)
        return m.igrowth[2:self.est_days+2]
        
    def computeMoments(self,pars):

        timestart=datetime.now()
        df = pd.DataFrame()
        
        # add seed number 
#        pars = np.insert(pars, 0, 0)
        kw=deepcopy(self.benchkwargs)
        kw['p_probc'][0][1] = pars[0]
        kw['p_avgstep'] = pars[1]
        
        if self.parallelize <= 1 :
        ## serialize bootstraps
            for boot in range(self.nboots) :
                newkw = deepcopy(kw)
                newkw['q_seed'] += boot
                growth = self.computeBoot(kwargs=newkw)
                df[boot] = growth
        else :
            ## parallelize bootstraps
            iternumbers = []
            for bootn in range(self.parallelize):
                newkw = deepcopy(kw)
                newkw['q_seed'] = kw['q_seed'] + bootn
                iternumbers.append(newkw) 

            
            # iternumbers = np.zeros((self.nboots,len(pars)))
            # iternumbers[:,0] = np.arange(self.nboots) # add seed numbers
            # iternumbers[:,1:] = pars[1:]
            pool = Pool(processes=self.parallelize)
            growth = pool.map(self.computeBoot,(iternumbers))
            pool.close()
            pool.join()
            df = pd.DataFrame(np.transpose(growth[0:self.est_days]))        
            ## end or parallelization
        
        self.ma = df.rolling(5).mean()
        self.ma = self.ma[4:]
        self.ma['lom'] = self.mldr[0:len(self.ma)]
        meanBoots = self.ma[self.ma.columns[0:self.nboots]].mean(axis=1)
        loss = ((meanBoots-self.ma['lom'])**2).sum()

        #draw a nice figure
        if self.print_screen >= 2 :
            fix, ax = plt.subplots(figsize=(5,5))
            ax.set_ylim(0,0.5)
            ax.plot(self.ma.index,meanBoots,c='chocolate',label='Infection growth, baseline')
            ax.plot(self.ma.index,self.ma['lom'],c='olive',label='Fatalities growth, lombardy')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            if self.print_screen>=2:
                ax.legend(title=str(np.round(pars,4))+' Loss:'+str(np.round(loss,5)))
            else:
                ax.legend()
            plt.show()

        if loss<self.minimum:
            self.minimum = loss

        #and print some stuff
        if self.print_screen >= 1 :
            print('[%.7f %.7f] Loss %.7f, Min %8.7f' % (pars[0],pars[1],loss,self.minimum)
                  ,str(datetime.now()-timestart)[:-7])
            
        return loss
    
    def computeMomentsWithSd(self,pars):
        
        global bottom, top
        timestart=datetime.now()
        normpars = np.copy(pars)
#        normpars[2] = squeeze(pars[2],bottom,top)  # make sure stedv is positiv
        loss = self.computeMoments(normpars)
        
        if self.print_screen == 2 :
            print('[%.7f %.7f %.7f] Loss %8.7f, Min %8.7f' % (pars[0],pars[1],pars[2],loss,self.minimum)
              ,str(datetime.now()-timestart)[:-7])

        return loss
#%% 
if __name__ == "__main__":

    benchkwargs = {
        "q_seed" : 2443,
        "p_proby"     : [0.0],
        "p_probr"     : [1/7],
        "p_probc"     : [[0, 0.048, 0, 0, 0]], #prob of contagion by type and state
        "p_probd"     : [0],
        "p_infradius" : 0.013016,
        "p_avgstep"   : 0.045,
        "p_stdstep"   : 0,
        "p_cluster"   : 'cluster',
        "q_popsize"   : 25600,
        "q_days"      : 42,
        "q_printOption" : 0.6,
        'g_maxinf'    : 1,
    }  
    # kwargs = deepcopy(benchkwargs)
    # est = estimatePars(benchkwargs = kwargs, parallelize= 6,
    #                     nboots=6, print_screen=2, shortcut=False,)
    # start = [0.06, 0.034]
    # est.computeMoments(start)
    
    # timestart = datetime.now()
    # est = estimatePars(benchkwargs,nboots=2,print_screen=3,parallelize=0)
    # a = est.computeMoments([0.0384902, 0.0344617])
    # print(a)
    # print(datetime.now()-timestart)
    
    # m = spatialSAYDR(**benchkwargs)
    # m.simulateOnce()
    #%%
if __name__ == "__main__":
    
    from scipy.optimize import minimize,basinhopping
    global bottom 
    global top 
    bottom = 0
    top = 2
    print('')

    kwargs = deepcopy(benchkwargs)
    est = estimatePars(benchkwargs = kwargs, est_days=42, parallelize= 6,
                       nboots=6, print_screen=2, shortcut=False,)
    start = np.array([0.040262457, 0.058719231])
    est.computeMoments(start)

    # define initial simplex
    nonzdelt = 0.15
    zdelt = 0.01
    sim = np.zeros((3, 2), dtype=start[0].dtype)
    sim[0] = start
    for k in range(2):
        y = np.array(start, copy=True)
        if y[k] != 0:
            y[k] = (1 + nonzdelt)*y[k]
        else:
            y[k] = zdelt
        sim[k + 1] = y
    ret = minimize(est.computeMoments,start,method='nelder-mead',options={'disp':True,'adaptive': True,'initial_simplex':sim})

    print("global minimum: x = [%.9f, %.9f], f(x0) = %.5f" % (ret.x[0],ret.x[1],ret.fun))
    
#%%
    kwargs = deepcopy(benchkwargs)
    est = estimatePars(benchkwargs = kwargs, parallelize= 6,
                       nboots=20, print_screen=2, shortcut=False,)
    start = np.array([0.045, 0.04])
    est.computeMoments(start)

    # start = np.array([0.039659132, 0.054135346])
    # Optimization terminated successfully.
    #      Current function value: 0.187858
    #      Iterations: 37
    #      Function evaluations: 95
    # global minimum: x = [0.039659132, 0.054135346], f(x0) = 0.18786
    
    ## start = [0.04, 0.03] #probc, avgstep
    ## global minimum: x = [0.03808820, 0.03402238], f(x0) = 0.12846
    # start = [0.03808820, 0.03402238]
    # global minimum: x = [0.03808823, 0.03402237], f(x0) = 0.12781
    # ret = minimize(est.computeMoments,start,method='nelder-mead',options={'disp':True,'adaptive': True})
    # start = [0.03808823, 0.03402237]
    # bounds = ((0.01,0.045),(0.01,0.5))
    # est = estimatePars(nboots=3,print_screen=1,shortcut=False)
    # ret = minimize(est.computeMoments,start,method='L-BFGS-B',bounds=bounds,options={'disp':True})
    # print("global minimum: x = [%.8f, %.8f], f(x0) = %.5f" % (ret.x[0],ret.x[1],ret.fun))
    # # minimizer_kwargs = {"method":"L-BFGS-B"}
    # ret = basinhopping(est.computeMoments,start,minimizer_kwargs=minimizer_kwargs,disp=True)
    #  ret = minimize(est.computeMoments, start, method='SLSQP',bounds = bounds, options={'disp':True})    
   

    # bounds = ((0.03, 0.05), (0.02, 0.05))
    # minimizer_kwargs = {"method":"SLSQP", 'bounds': bounds}
    # ret = basinhopping(est.computeMoments,start,minimizer_kwargs=minimizer_kwargs,disp=True)
    # basinhopp qing step 100: f 0.354081 trial_f 0.354081 accepted 1  lowest_f 0.161787
    # global minimum: x = [0.038088230, 0.034022370], f(x0) = 0.16179

    # ## compute with sd
#     est = estimatePars(nboots=8,print_screen=2,shortcut=False)
# #     startPars = [0.04, 0.03,blowup(0.03,bottom,top)]
# #     startPars = [0.037907, 0.030653, 0.03008]
# #     bounds = ((0.025, 0.045), (0.02, 0.4), (0, .3))

#     startPars = [0.0379070,0.0306536, 0.0300808] # 0.1147765
#     bounds = ((0.04, 0.042), (0.025, 0.35), (0, .1))
#     ret = minimize(est.computeMomentsWithSd, startPars, method='SLSQP',bounds = bounds, options={'disp':True})    
#     print("global minimum: x = [%.9f, %.9f %.9f], f(x0) = %.5f" % (ret.x[0],ret.x[1],ret[2],ret.fun))

    # minimizer_kwargs = {"method":"L-BFGS-B"}
    # ret = basinhopping(est.computeMomentsWithSd,startPars,minimizer_kwargs=minimizer_kwargs,disp=True)

    # start = [0.03656626, 0.02797655]
    # minimizer_kwargs = {"method":"L-BFGS-B"}
    

   
# print(ret)
#  final_simplex: (array([[0.03656626, 0.02797655],
#        [0.03656625, 0.02797655],
#        [0.03656626, 0.02797655]]), array([0.07917224, 0.07917815, 0.07924045]))
#            fun: 0.07917223591941705
#        message: 'Optimization terminated successfully.'
#           nfev: 90
#            nit: 37
#         status: 0
#        success: True
#              x: array([0.03656626, 0.02797655])
