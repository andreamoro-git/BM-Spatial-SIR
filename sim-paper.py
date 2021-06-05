#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Apr 11 18:29:52 2020

Use this file to replicate simulations. 
Results are saved to in output/*.pickle.gz for further analysis
"""

import numpy as np
from class_spatialModels import simulatePool
import pickle,gzip
from multiprocessing import Pool
import time
from copy import deepcopy

prefix= 'nc5-' 
start_time = time.time()
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


#%% Simulate the basic model twenty times
if __name__ == "__main__":
 
    kwargs = deepcopy(benchkwargs)
    
    nboots = 20
    nprocs = 7
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allmodels))
    
    # save all objects
    with gzip.open(outputdir+'basesim.pickle.gz', 'wb') as outputfile:
        pickle.dump(allmodels,outputfile)
    outputfile.close()
      
    print('Time spent ', time.time()-start_time)

#%% Random locations every day

from class_spatialModels import simulateRandPool
if __name__ == '__main__':

    kwargs = deepcopy(benchkwargs)   
    nboots = 10
    nprocs = 7

    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandmodels = pool.map(simulateRandPool,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allRandmodels))
    
    # save all objects    
    file = gzip.open(outputdir+'allrandmodels.pickle.gz','wb')
    pickle.dump(allRandmodels,file)
    file.close()

    print('Time spent ', time.time()-start_time)
    
#%% Random initial location

    kwargs = deepcopy(benchkwargs)
    kwargs['p_cluster']   = 'random'
    
    nboots = 20
    nprocs = 7
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allmodels))
    
    # save all objects
    with gzip.open(outputdir+'random.pickle.gz', 'wb') as outputfile:
        pickle.dump(allmodels,outputfile)

    print('Time spent ', time.time()-start_time)


#%% 1/4 the city size

if __name__ == "__main__":
    
    kwargs = deepcopy(benchkwargs)
    kwargs["q_popsize"]   = np.int(0.25*160**2)
    kwargs["q_citysize"]  = 1/2
    
    nboots = 20
    nprocs = 7
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()

    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allmodels))
    
    # save all objects
    with gzip.open(outputdir+'quartersim.pickle.gz', 'wb') as outputfile:
        pickle.dump(allmodels,outputfile)
    outputfile.close()

    print('Time spent ', time.time()-start_time)

#%% 6x contagion, 1/6 density

if __name__ == "__main__":
    
    nboots = 20
    nprocs = 7
    
    kwargs=deepcopy(benchkwargs)
    kwargs['q_citysize'] = 1 * np.sqrt(6)
    kwargs['p_probc'][0][1] = 6 * benchkwargs['p_probc'][0][1]
    
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

    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allRandBigger))
    
    # save all objects    
    file = gzip.open(outputdir+'benchmark_6x_cont.pickle.gz','wb')
    pickle.dump(allRandBigger,file)
    file.close()
    
    print('Time spent ', time.time()-start_time)

# #%% No movement

# if __name__ == "__main__":
    
#     from class_spatialModels import simulatePool

#     kwargs = deepcopy(benchkwargs)
#     kwargs['p_avgstep']   = 0
    
#     nboots = 20
#     nprocs = 7
    
#     # constructt list to pass to processes
#     iterboots = np.arange(nboots)
#     iterlist = []
    
#     # change seeds
#     for bootn in range(nboots):
#         newkw = kwargs.copy()
#         newkw['q_seed'] = kwargs['q_seed'] + bootn
#         iterlist.append((newkw)) 
    
#     # spawn nboots processes, nprocs at once    
#     pool = Pool(processes=nprocs)
#     allmodels = pool.map(simulatePool,(iterlist))
#     pool.close()
#     pool.join()
        
#     file = gzip.open(outputdir+'nomovement.pickle.gz','wb')
#     pickle.dump(allmodels,file)
#     file.close()

#     print('Time spent ', time.time()-start_time)

#%% Different city density (different size same population) 

if __name__ == "__main__":
    
    kwargs = deepcopy(benchkwargs)
    nboots = 10
    nprocs = 7
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    # change seeds
    for size in [0.7,1.4] :
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
    with gzip.open(outputdir+'density.pickle.gz', 'wb') as outputfile:
        pickle.dump(allmodels,outputfile)
    outputfile.close()
    
    print('Time spent ', time.time()-start_time)
             
    
#%% Different movement speed 

    nboots = 10
    nprocs = 7
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    kwargs = deepcopy(benchkwargs)

    # define speeds as fraction of baseline speed
    basespeed = kwargs['p_avgstep']
    
    #speeds = [0,0.04,0.08,0.12,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,2]
    speeds = [0,0.2]
    for speed in speeds :
        for bootn in range(nboots):
            newkw = kwargs.copy()
            newkw['q_seed'] = kwargs['q_seed'] + bootn
            newkw['p_avgstep'] = basespeed * speed
            if speed<=0.4:
                newkw['q_days'] = 600
            iterlist.append((newkw)) 
        
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allmodels))
    
    # save all objects    
    file = gzip.open(outputdir+'allmovements.pickle.gz','wb')
    pickle.dump(allmodels,file)
    file.close()

    print('Time spent ', time.time()-start_time)


#%% Heterogeneous density

if __name__ == "__main__":

    from class_spatialModels import simulateHetPool
    
    kwargs = deepcopy(benchkwargs)
    kwargs['p_cluster'] = 'cluster'
    kwargs['q_lambda'] = 1

    #simulatePool(kwargs)
    nboots = 20
    nprocs = 4
  
    # construct list of parameters to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulateHetPool,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allmodels))
        
    # save all objects    
    file = gzip.open(outputdir+'hetdens-lambda-1.pickle.gz','wb')
    pickle.dump(allmodels,file)
    file.close()
    
    print('Time spent ', time.time()-start_time)



#%% Spatial behavioral SIR

from class_spatialModels import simulateBehPool
if __name__ == '__main__':

    kwargs = deepcopy(benchkwargs)
    nboots = 20
    nprocs = 7
    kwargs['behModel'] = {'type': 'Lones','phi': 0.01}

    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allBehmodels = pool.map(simulateBehPool,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allBehmodels))
    
    # save all objects    
    file = gzip.open(outputdir+'basebehLones.pickle.gz','wb')
    pickle.dump(allBehmodels,file)
    file.close()
    
    print('Time spent ', time.time()-start_time)

#%% Spatial SIR with twice the city side and population
 
    kwargs = deepcopy(benchkwargs) 
    kwargs["q_popsize"]   = 25600*4
    kwargs["q_citysize"]  = 2
    kwargs['p_cluster']   = 'cluster'
    
    nboots = 10
    nprocs = 7
    
    # construct list to pass to processes
    iterboots = np.arange(nboots)
    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
        
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allmodels = pool.map(simulatePool,(iterlist))
    pool.close()
    pool.join()
        
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allmodels))
    
    # save all objects
    with gzip.open(outputdir+'basesize2.pickle.gz', 'wb') as outputfile:
        pickle.dump(allmodels,outputfile)
    outputfile.close()

    print('Time spent ', time.time()-start_time)

#%% Spatial behavioral SIR. Local

from class_spatialModels import simulateBehPool_local
if __name__ == '__main__':

    kwargs = deepcopy(benchkwargs)
    nboots = 14
    nprocs = 7
    kwargs['behModel'] = {'type': 'Lones','phi': 0.01}

    iterlist = []
    
    # change seeds
    for bootn in range(nboots):
        newkw = kwargs.copy()
        newkw['q_seed'] = kwargs['q_seed'] + bootn
        iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allBehmodels = pool.map(simulateBehPool_local,(iterlist))
    pool.close()
    pool.join()
    
    # we don't need to save positions (saves some space)
    list(map(lambda x: delattr(x,'pos'),allBehmodels))
    
    # save all objects    
    file = gzip.open(outputdir+'basebehLones_local.pickle.gz','wb')
    pickle.dump(allBehmodels,file)
    file.close()
    
    print('Time spent ', time.time()-start_time)


    #%%
    
    prefix= 'nc5-' 
    start_time = time.time()
    outputdir = 'output/'+prefix
    imagedir = 'output/images/'+prefix
    
    
    
    estbenchkwargs = {
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
        "q_days"      : 81,
        "q_printOption" : 2.6,
        'g_maxinf'    : 1,
        'p_shutr'     : [20,0.30],
        'p_openr'     : [999,1],
    }  


    #%% Spatial-SIR, generate data from several densities, no benhavior

    from class_simul_policy import policyTimePool
    estbenchkwargs = {
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
        "q_days"      : 80,
        "q_printOption" : 0.6,
        'g_maxinf'    : 1,
        'p_shutr'     : [20,0.30],
        'p_openr'     : [999,1],
        'behModel'    : {'type': 'None'}
    }  
    
    nboots = 5
    nprocs = 10
    
    iterlist = []
    
    # change seeds
    for citysize in (1/np.sqrt(np.append(np.arange(.5,.6,.02),np.arange(.6,1.6,.1)))):
        for shutdown in ( [20, 0.25], [999, 0] ) :
            for bootn in range(nboots):
                newkw = estbenchkwargs.copy()
                newkw['q_citysize'] = citysize
                newkw['p_shutr'] = shutdown
                newkw['q_seed'] = estbenchkwargs['q_seed'] + bootn
                iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandmodels = pool.map(policyTimePool,(iterlist))
    pool.close()
    pool.join()
          
    list(map(lambda x: delattr(x,'pos'),allRandmodels))
    
    file = gzip.open(outputdir+'dens-20-80-25pc.gz','wb')
    pickle.dump(allRandmodels,file)
    file.close()

#%% Spatial-SIR, densities only, no behavior, data save

    import pandas as pd
    
    file = gzip.open(outputdir+'dens-20-80-25pc.gz','rb')
    allRandmodels = pickle.load(file)
    file.close()

    moddata = pd.DataFrame()
    for idx, model in enumerate(allRandmodels):
        modeldf = pd.DataFrame([idx]*model.q_days,columns=['naics'])
        modeldf['active']= model.nstates[:,1]
        modeldf['susceptible'] = model.nstates[:,0]
        modeldf['citysize'] = model.q_citysize
        modeldf['total'] = model.nstates[-1,1] + model.nstates[-1,4]
        modeldf['density'] = model.q_popsize/(model.q_citysize**2)*1/benchkwargs['q_popsize']
        modeldf['q_popsize'] = model.q_popsize
        modeldf['npi_date']= model.p_shutr[0]
        modeldf['npi_size']= model.p_shutr[1]
        modeldf['outside']= model.outside
        modeldf['bootn']= model.q_seed
        moddata = moddata.append(modeldf)
        
    moddata.to_csv(outputdir+'dens-20-80-25pc.csv',index_label='t')
 
#%% Spatial-SIR, different densities with behavior

    from class_simul_policy import policyTimePool
    
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
        "q_days"      : 80,
        "q_printOption" : 0.6,
        'g_maxinf'    : 1,
        'p_shutr'     : [20,0.25],
        'p_openr'     : [999,0],
        'behModel'    : {'type': 'Lones','phi': 0.01}
    }  
    
    kwargs = benchkwargs.copy()
    nboots = 5
    nprocs = 10
    
    iterlist = []
    
    # change seeds
    for citysize in (1/np.sqrt(np.append(np.arange(.5,.6,.02),np.arange(.6,1.6,.1)))):
        for shutdown in ( [20, 0.25], [999, 0] ) :
            for bootn in range(nboots):
                newkw = benchkwargs.copy()
                newkw['q_citysize'] = citysize
                newkw['p_shutr'] = shutdown
                newkw['q_seed'] = kwargs['q_seed'] + bootn
                iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandmodels = pool.map(policyTimePool,(iterlist))
    pool.close()
    pool.join()
          
    list(map(lambda x: delattr(x,'pos'),allRandmodels))
    
    file = gzip.open(outputdir+'dens-beh_p-20-80-25pc.gz','wb')
    pickle.dump(allRandmodels,file)
    file.close()
    
#   create csv data
    file = gzip.open(outputdir+'dens-beh_p-20-80-25pc.gz','rb')
    allRandmodels = pickle.load(file)
    file.close()
    
    moddata = pd.DataFrame()
    for idx, model in enumerate(allRandmodels):
        modeldf = pd.DataFrame([idx]*model.q_days,columns=['naics'])
        modeldf['active']= model.nstates[:,1]
        modeldf['susceptible'] = model.nstates[:,0]
        modeldf['total'] = model.nstates[-1,1] + model.nstates[-1,4]
        modeldf['citysize'] = model.q_citysize
        modeldf['density'] = model.q_popsize/(model.q_citysize**2)*1/benchkwargs['q_popsize']
        modeldf['q_popsize'] = model.q_popsize
        modeldf['bootn']= model.q_seed
        modeldf['outside']= model.outside
        modeldf['npi_date'] = model.p_shutr[0]
        modeldf['npi_size']= model.p_shutr[1]

        moddata = moddata.append(modeldf)
      
    moddata.to_csv(outputdir+'dens-beh_p-20-80-25pc.csv',index_label='t')

#%% different densities without behavior, policies at different t

    from class_simul_policy import policyTimePool
    
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
        "q_days"      : 80,
        "q_printOption" : 0.6,
        'g_maxinf'    : 1,
        'p_shutr'     : [20,0.25],
        'p_openr'     : [999,1],
        'behModel'    : {'type': 'None'}
    }  
    
    kwargs = benchkwargs.copy()
    nboots = 5
    nprocs = 10
    
    iterlist = []
    
    # change seeds
    for citysize in (1/np.sqrt(np.arange(.5, 1.6, .5))):
        for shutdate in (15,40,999) :
            for shutdown in ( [shutdate, 0.25], [999, 0] ) :
                newkw = benchkwargs.copy()
                newkw['q_citysize'] = citysize
                newkw['p_shutr'] = shutdown
                #newkw['q_seed'] = kwargs['q_seed'] + bootn
                iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandmodels2 = pool.map(policyTimePool,(iterlist))
    pool.close()
    pool.join()
          
    list(map(lambda x: delattr(x,'pos'),allRandmodels2))
    
    file = gzip.open(outputdir+'dens-times-25pc.gz','wb')
    pickle.dump(allRandmodels2,file)
    file.close()
    
#%% create csv data
    file = gzip.open(outputdir+'dens-times-25pc.gz','rb')
    allRandmodels = pickle.load(file)
    file.close()
    
    moddata = pd.DataFrame()
    for idx, model in enumerate(allRandmodels):
        modeldf = pd.DataFrame([idx]*model.q_days,columns=['naics'])
        modeldf['active']= model.nstates[:,1]
        modeldf['susceptible'] = model.nstates[:,0]
        modeldf['citysize'] = model.q_citysize
        modeldf['total'] = model.nstates[-1,1] + model.nstates[-1,4]
        modeldf['density'] = model.q_popsize/(model.q_citysize**2)*1/benchkwargs['q_popsize']
        modeldf['q_popsize'] = model.q_popsize
        modeldf['bootn']= model.q_seed
        modeldf['outside']= model.outside
        modeldf['npi_date'] = model.p_shutr[0]
        modeldf['npi_size']= model.p_shutr[1]
        moddata = moddata.append(modeldf)
      
    moddata.to_csv(outputdir+'dens-times-25pc.csv',index_label='t')



#%% different densities with behavior, policies at different t

    from class_simul_policy import policyTimePool
    
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
        "q_days"      : 80,
        "q_printOption" : 0.6,
        'g_maxinf'    : 1,
        'p_shutr'     : [20,0.25],
        'p_openr'     : [999,1],
        'behModel'    : {'type': 'Lones','phi': 0.01}
    }  
    
    kwargs = benchkwargs.copy()
    nboots = 5
    nprocs = 7
    
    iterlist = []
    
    # change seeds
    for citysize in (1/np.sqrt(np.arange(.5, 1.6, .1))):
        for shutdate in np.arange(15,41,5) :
            for shutdown in ( [shutdate, 0.25], [999, 0] ) :
                newkw = benchkwargs.copy()
                newkw['q_citysize'] = citysize
                newkw['p_shutr'] = shutdown
                #newkw['q_seed'] = kwargs['q_seed'] + bootn
                iterlist.append((newkw)) 
    
    # spawn nboots processes, nprocs at once    
    pool = Pool(processes=nprocs)
    allRandmodels2 = pool.map(policyTimePool,(iterlist))
    pool.close()
    pool.join()
          
    list(map(lambda x: delattr(x,'pos'),allRandmodels2))
    
    file = gzip.open(outputdir+'dens-beh_p-times-25pc.gz','wb')
    pickle.dump(allRandmodels2,file)
    file.close()
    
# create csv data
    file = gzip.open(outputdir+'dens-beh_p-times-25pc.gz','rb')
    allRandmodels = pickle.load(file)
    file.close()
    
    moddata = pd.DataFrame()
    for idx, model in enumerate(allRandmodels):
        modeldf = pd.DataFrame([idx]*model.q_days,columns=['naics'])
        modeldf['active']= model.nstates[:,1]
        modeldf['susceptible'] = model.nstates[:,0]
        modeldf['citysize'] = model.q_citysize
        modeldf['total'] = model.nstates[-1,1] + model.nstates[-1,4]
        modeldf['density'] = model.q_popsize/(model.q_citysize**2)*1/benchkwargs['q_popsize']
        modeldf['q_popsize'] = model.q_popsize
        modeldf['bootn']= model.q_seed
        modeldf['outside']= model.outside
        modeldf['npi_date'] = model.p_shutr[0]
        modeldf['npi_size']= model.p_shutr[1]
        moddata = moddata.append(modeldf)
      
    moddata.to_csv(outputdir+'dens-beh_p-times-25pc.csv',index_label='t')

