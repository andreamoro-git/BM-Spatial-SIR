#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 11:41:03 2020

"""

import numpy as np
import gzip, pickle
from copy import deepcopy
from class_spatialModels import simulatePool
from multiprocessing import Pool

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

#%% get replications with different n. of clusters
if __name__ == "__main__":

    kwargs = deepcopy(benchkwargs)
    
    nboots = 4
    nprocs = 6
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    # define speeds as fraction of baseline speed
    basespeed = kwargs['p_avgstep']
     
    allclusters = [2,3,5,6,10]
    for clusters in allclusters :
        for boot in range(nboots) :
            kwargs = deepcopy(benchkwargs)
            kwargs['q_seed'] += boot
            if clusters == 2 :
                kwargs['q_firstinf'] = [[0.25,0.25],[0.75,0.75]]
                kwargs['p_ininf'] = 15
            elif clusters == 3 :
                kwargs['q_firstinf'] = [[0.25,0.25],[0.75,0.75],[0.5,0.5]]
                kwargs['p_ininf'] = 10
            elif clusters == 5 :
                kwargs['q_firstinf'] = [[0.25,0.25],[0.25,0.75],[0.75,0.25],[0.5,0.5],[0.75,0.75]]
                kwargs['p_ininf'] = 6
            elif clusters == 6 :
               kwargs['q_firstinf'] = [[0.20,0.20],[0.2,0.8],[0.8,0.2],[0.33,0.5],[0.67,0.5],[0.8,0.8]]
               kwargs['p_ininf'] = 5
            elif clusters == 10 :
               kwargs['q_firstinf'] = [[0.2,0.2],[0.2,0.8],[0.8,0.2],[0.8,0.8],[0.2,0.5],[0.5,0.2],[0.8,0.5],[0.5,0.8],[0.33,0.3],[0.67,0.67]]
               kwargs['p_ininf'] = 3
            elif clusters == 15 :
               kwargs['q_firstinf'] = [[0.20,0.20],[0.2,0.8],[0.8,0.2],[0.8,0.8],[0.2,0.5],[0.5,0.2],[0.8,0.5],[0.5,0.8],[0.33,0.3],[0.67,0.67],[0.33,06.],[0.67,0.33]]
               kwargs['p_ininf'] = 2
              
            iterlist.append(kwargs)
        
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allmodels))
        
    # save all objects    
    file = gzip.open(outputdir+'allclusters.pickle.gz','wb')
    pickle.dump(allmodels,file)
    file.close()    
        
#%% replication for many city sizes
if __name__ == "__main__":

    nboots = 4
    nprocs = 6
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    kwargs = benchkwargs
    basespeed = kwargs['p_avgstep']
    
    
    sizes = np.append(np.arange(10000,40000,5000),np.arange(40000,91000,10000))
    for size in sizes :
        for bootn in range(nboots):
            newkw = kwargs.copy()
            newkw['q_seed'] = kwargs['q_seed'] + bootn
            newkw['q_popsize'] = size
            newkw['q_citysize'] = np.sqrt(size)/160
            newkw['q_days'] = 400
            iterlist.append((newkw)) 
    
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    list(map(lambda x: delattr(x,'pos'),allmodels))
    file = gzip.open(outputdir+'allsizes.pickle.gz','wb')
    pickle.dump(allmodels,file)
    file.close()

#%% Replication of more densities

if __name__ == "__main__":
    
    kwargs = deepcopy(benchkwargs)
    nboots = 4
    nprocs = 6
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    # change seeds
    for size in np.arange(0.5,2.1,0.1) :
        if size in [0.7,1.4] :
            continue
        for bootn in range(nboots):
            newkw = kwargs.copy()
            newkw['q_seed'] = kwargs['q_seed'] + bootn
            newkw['q_citysize'] = size
            iterlist.append((newkw)) 
       
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allmodels))
    
    # save all objects    
    with gzip.open(outputdir+'app-density.pickle.gz', 'wb') as outputfile:
        pickle.dump(allmodels,outputfile)
    outputfile.close()
  
