# -*- coding: utf-8 -*-
"""
Created on Mon, Mar 1, 2021

This file simulates data for comparison with estimation
"""

import numpy as np
from class_simul_policy import simul_policyByTime, policyTimePool
import pickle,gzip
from multiprocessing import Pool
import time
from copy import deepcopy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

prefix= 'nc5-' 
start_time = time.time()
outputdir = 'output/'+prefix
imagedir = 'output/images/'+prefix

fsize = 3.3

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
    "q_days"      : 81,
    "q_printOption" : 2.6,
    'g_maxinf'    : 1,
    'p_shutr'     : [20,0.30],
    'p_openr'     : [999,1],
}  

import matplotlib
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"]})
    
# #%% Generate data from SIR model

# if __name__ == "__main__":
#     from class_SIRmodel import SIRmodel
#     baseline = 25600
#     cluster = 10
#     beta = .15  #.054*13.5,
#     delta = .02 #1/6.5
    
#     models = list()
    
#     base = SIRmodel(beta, delta ,q_popsize=baseline,
#                 q_init=cluster,
#                 )
#     base.plot(base.day)
#     models.append(base)

#     pol = SIRmodel(beta, delta, q_popsize=baseline,
#                 q_init=cluster, p_shutr=[20,.3], p_openr=[999,1]
#                 )
#     pol.plot(pol.day)
#     models.append(pol)
    
#     dsize = base.day-1
#     dpsize = pol.day-1
#     dp = base.SIR[-dsize:, 1] - base.SIR[:dsize, 1]
#     dpol = pol.SIR[-dpsize:, 1] - pol.SIR[:dpsize, 1]
    
#     #
#     days = 100
#     fig,(ax2) = plt.subplots(1,1,figsize=(fsize,fsize))
#     ax2.spines['right'].set_visible(False)
#     ax2.spines['top'].set_visible(False)
#     ax2.set_xlabel('Days ')

#     # ax2.plot(np.arange(days),bp.prinf[:days],color='olive',label ='Untreated')
#     # ax2.plot(np.arange(days),bpol.prinf[:days],':',linewidth=1.9, color='darkorange',label ='Treated')
#     ax2.plot(np.arange(days),dp[:days],color='olive',label ='Untreated')
#     ax2.plot(np.arange(days),dpol[:days],':',linewidth=1.9, color='darkorange',label ='Treated')


#     ax2.legend(title='Cases daily change')
#     fig.tight_layout() 
    
    
#     moddata = pd.DataFrame()
#     days = 80
#     for idx, model in enumerate(models):
#         modeldf = pd.DataFrame([idx]*days,columns=['naics'])
#         modeldf['Susceptible']= model.SIR[:days,0]
#         modeldf['Active']= model.SIR[:days,1]
#         modeldf['Total'] = model.q_popsize-model.SIR[:days,0]-model.SIR[:days,3]
#         modeldf['Locked'] = model.SIR[:days,3]
#         modeldf['npi_date']= model.p_shutr[0]
#         modeldf['npi_size']= model.p_shutr[1]
#         moddata = moddata.append(modeldf)
        
#     moddata.to_csv(outputdir+'SIR-policy-20-30pc.csv',index_label='t')

#%% Spatial-SIR, generate data from several densities, no benhavior

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
        'p_shutr'     : [20,0.30],
        'p_openr'     : [999,1],
        'behModel'    : {'type': 'None'}
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
    
    file = gzip.open(outputdir+'dens-20-80-25pc.gz','wb')
    pickle.dump(allRandmodels,file)
    file.close()

#%% Spatial-SIR, densities only, no behavior, data save

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
 
#%% dSpatial-SIR, different densities with behavior

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


#%% Figure comparing estimated betas in baseline and behavioral with real beta

#   Before doing this, run sim_estimate_dens.do in stata to generate densbetas.csv
#   and dens-behbetas.csv
    betaframe_nobeh = pd.read_csv(outputdir+'densbetas.csv')
    density_nobeh = betaframe_nobeh['density']
    betalamd_nobeh = betaframe_nobeh['beta']

    betaframe = pd.read_csv(outputdir+'dens-behbetas.csv')
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
    #ax2.axhline(y=1, color='grey',linewidth=.5 )    
    
    fig.tight_layout()
    plt.savefig(imagedir+'est_densitybetas.pdf')

    # only data before peak
    # betaframe = pd.read_csv(outputdir+'dens-behbetas-beforepeak.csv')
    # density = betaframe['density']
    # beh_betalamd = betaframe['beta']
    
    # fig,(ax2) = plt.subplots(1,1,figsize=(fsize,fsize))
    # ax2.spines['right'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # ax2.set_xlabel('Density')
    # ax2.set_ylabel('Coefficient')
    # ax2.set_ylim(0.33,1.8)
    # ax2.plot(density,beh_betalamd/beh_betalamd[9],color='darkorange',label='behavioral', marker='o')
    # ax2.plot(density,betalamdbp/betalamdbp[9],color='saddlebrown',label='benchmark', marker='^')
    # #ax2.plot(density,density,':',color='olive',label='density')
    # ax2.legend(loc='lower right')
    # ax2.axhline(y=1, color='grey',linewidth=.5 )    

    # plt.savefig(imagedir+'est_densitybetas_beforepeak.pdf')

    
    
# #%% Figure comparing high/low density with and without policy

#     file = gzip.open(outputdir+'dens-20-80-25pc.gz','rb')
#     allRandmodels = pickle.load(file)
#     file.close()
#     import matplotlib.pyplot as plt
#     from class_averageStats import averageStats
    
#     bpol = averageStats(allRandmodels[0:4])
#     bp = averageStats(allRandmodels[5:9])
#     bpol1 = averageStats(allRandmodels[140:144])
#     bp1 = averageStats(allRandmodels[145:149])


#     # first plot: compare SIR with benchmark 
#     dsize = bp.prinf.size-1
#     dpol = bpol.prinf[-dsize:] - bpol.prinf[:dsize]
#     dp = bp.prinf[-dsize:] - bp.prinf[:dsize]
#     # bpol=allRandmodels[0]
#     # bp = allRandmodels[5]
#     days=80
#     fig,(ax2) = plt.subplots(1,1,figsize=(fsize*1.7,fsize))
#     ax2.spines['right'].set_visible(False)
#     ax2.spines['top'].set_visible(False)
#     ax2.set_xlabel('Days ')

#     avp = (bpol1.prinf*3+bpol.prinf*3)/10
#     av = (bp1.prinf*3+bp.prinf*7)/10
#     # avp[21:] = (bpol1.prinf[21:]*3+bpol.prinf[21:]*7)/10
#     # av[21:] = (bp1.prinf[21:]*3+bp.prinf[21:]*7)/10
#     bvp = (bpol1.prinf*7+bpol.prinf*3)/10
#     bv = (bp1.prinf*7+bp.prinf*3)/10


#     ax2.spines['right'].set_visible(False)
#     ax2.spines['top'].set_visible(False)
#     ax2.set_xlabel('Days ')
#     ax2.plot(np.arange(days),bp1.prinf[:days],color='olive',lw=.5,label ='high density')
#     ax2.plot(np.arange(days),bpol1.prinf[:days],':',linewidth=1, color='olive',)
#     ax2.plot(np.arange(days),bp.prinf[:days],color='darkorange', lw=.5, label ='low density')
#     ax2.plot(np.arange(days),bpol.prinf[:days],':',linewidth=1, color='darkorange',)
#     ax2.plot(np.arange(days),bv[:days],color='black', label='true 70-30% avg')
#     ax2.plot(np.arange(days),bvp[:days],':',linewidth=1.9, color='black',label= 'true 70-30% treatment')
#     ax2.plot(np.arange(days),avp[:days],color='saddlebrown',label ='treated (30% low density)')
    
#     # ax2.plot(np.arange(days-1),dp[:days-1],color='olive',label ='Untreated')
#     # ax2.plot(np.arange(days-1),dpol[:days-1],':',linewidth=1.9, color='darkorange',label ='Treated')

#     ax2.axvline(x=20, color='grey',linewidth=.5 )    
#     ax2.legend(title='Cases daily % change')
#     fig.tight_layout()    
    
#     plt.savefig(imagedir+'ddcomp.pdf')

    
#%% different population only
    pass
#     benchkwargs = {
#         "q_seed" : 2443,
#         "p_proby"     : [0.0],
#         "p_probr"     : [1/6.5],
#         "p_probc"     : [[0, 0.054, 0, 0, 0]], #prob of contagion by type and state
#         "p_probd"     : [0],
#         "p_infradius" : 0.013016,
#         "p_avgstep"   : 0.034,
#         "p_stdstep"   : 0,
#         "p_cluster"   : 'cluster',
#         "q_popsize"   : 25600,
#         "q_days"      : 81,
#         "q_printOption" : 0.6,
#         'g_maxinf'    : 1,
#         'p_shutr'     : [20,0.30],
#         'p_openr'     : [80,1],
#     }  
    
#     kwargs = benchkwargs.copy()
#     nboots = 5
#     nprocs = 10
    
#     iterlist = []
    
#     # change seeds
#     for popscale in (np.append(np.arange(.5,.6,.02),np.arange(.6,1.6,.1))):
#         for shutdown in ( [20, 0.25], [999, 0] ) :
#             for bootn in range(nboots):
#                 newkw = benchkwargs.copy()
#                 newkw['q_popsize'] = np.int(np.round(popscale*benchkwargs['q_popsize']))
#                 newkw['p_shutr'] = shutdown
#                 newkw['q_seed'] = kwargs['q_seed'] + bootn
#                 iterlist.append((newkw)) 
    
#     # spawn nboots processes, nprocs at once    
#     pool = Pool(processes=nprocs)
#     allRandmodels = pool.map(policyTimePool,(iterlist))
#     pool.close()
#     pool.join()
          
#     list(map(lambda x: delattr(x,'pos'),allRandmodels))
    
#     file = gzip.open(outputdir+'pop-20-80-25pc.gz','wb')
#     pickle.dump(allRandmodels,file)
#     file.close()

        
#     moddata = pd.DataFrame()
#     for idx, model in enumerate(allRandmodels):
#         modeldf = pd.DataFrame([idx]*model.q_days,columns=['naics'])
#         modeldf['active']= model.nstates[:,1]
#         modeldf['susceptible'] = model.nstates[:,0]
#         modeldf['citysize'] = model.q_citysize
#         modeldf['density'] = model.q_popsize/(model.q_citysize**2)*1/benchkwargs['q_popsize']
#         modeldf['q_popsize'] = model.q_popsize
#         modeldf['npi_date']= model.p_shutr[0]
#         modeldf['npi_size']= model.p_shutr[1]
#         modeldf['bootn']= model.q_seed
#         moddata = moddata.append(modeldf)
        
#     moddata.to_csv(outputdir+'pop-20-80-25pc.csv',index_label='t')
  
# #%% 
#     betaframe = pd.read_csv(outputdir+'popbetas.csv')
#     density = betaframe['density']
#     betalam = betaframe['beta']
    
#     fig,(ax2) = plt.subplots(1,1,figsize=(fsize,fsize))
#     ax2.spines['right'].set_visible(False)
#     ax2.spines['top'].set_visible(False)
#     ax2.set_xlabel('Density')
#     ax2.set_ylabel('Coefficient')
#     ax2.plot(density,betalam/betalam[9],color='darkorange',label='coefficient', marker='o')
#     ax2.plot(density,density,color='olive',label='density')
#     ax2.legend()
    

# #%% different population+initial cluster only

#     benchkwargs = {
#         "q_seed" : 2443,
#         "p_proby"     : [0.0],
#         "p_probr"     : [1/6.5],
#         "p_probc"     : [[0, 0.054, 0, 0, 0]], #prob of contagion by type and state
#         "p_probd"     : [0],
#         "p_infradius" : 0.013016,
#         "p_avgstep"   : 0.034,
#         "p_stdstep"   : 0,
#         "p_cluster"   : 'cluster',
#         "q_popsize"   : 25600,
#         "q_days"      : 81,
#         "q_printOption" : 0.6,
#         'g_maxinf'    : 1,
#         'p_shutr'     : [20,0.30],
#         'p_openr'     : [80,1],
#         'p_ininf'     : 30,
#     }  
    
#     kwargs = benchkwargs.copy()
#     nboots = 5
#     nprocs = 10
    
#     iterlist = []
    
#     # change seeds
#     for popscale in (np.append(np.arange(.5,.6,.05),np.arange(.6,1.6,.1))):
#         for shutdown in ( [20, 0.25], [999, 0] ) :
#             for bootn in range(nboots):
#                 newkw = benchkwargs.copy()
#                 newkw['q_popsize'] = np.int(np.round(popscale*benchkwargs['q_popsize']))
#                 newkw['p_ininf'] = np.int(np.round(popscale*benchkwargs['p_ininf']))
#                 newkw['p_shutr'] = shutdown
#                 newkw['q_seed'] = kwargs['q_seed'] + bootn
#                 iterlist.append((newkw)) 
    
#     # spawn nboots processes, nprocs at once    
#     pool = Pool(processes=nprocs)
#     allRandmodels = pool.map(policyTimePool,(iterlist))
#     pool.close()
#     pool.join()
          
#     list(map(lambda x: delattr(x,'pos'),allRandmodels))
    
#     file = gzip.open(outputdir+'pop_cl-20-80-25pc.gz','wb')
#     pickle.dump(allRandmodels,file)
#     file.close()

        
#     moddata = pd.DataFrame()
#     for idx, model in enumerate(allRandmodels):
#         modeldf = pd.DataFrame([idx]*model.q_days,columns=['naics'])
#         modeldf['active']= model.nstates[:,1]
#         modeldf['susceptible'] = model.nstates[:,0]
#         modeldf['citysize'] = model.q_citysize
#         modeldf['density'] = model.q_popsize/(model.q_citysize**2)*1/benchkwargs['q_popsize']
#         modeldf['q_popsize'] = model.q_popsize
#         modeldf['npi_date']= model.p_shutr[0]
#         modeldf['npi_size']= model.p_shutr[1]
#         modeldf['bootn']= model.q_seed
#         moddata = moddata.append(modeldf)
        
#     moddata.to_csv(outputdir+'pop_cl-20-80-25pc.csv',index_label='t')
  
# #%% Figure
#     import matplotlib.pyplot as plt
#     from class_averageStats import averageStats
    
#     bpol = averageStats(allRandmodels[95:99])
#     bp = averageStats(allRandmodels[115:119])


#     # first plot: compare SIR with benchmark 
#     dsize = bp.prinf.size-1
#     dpol = bpol.prinf[-dsize:] - bpol.prinf[:dsize]
#     dp = bp.prinf[-dsize:] - bp.prinf[:dsize]
#     # bpol=allRandmodels[0]
#     # bp = allRandmodels[5]
#     days=81
#     fig,(ax2) = plt.subplots(1,1,figsize=(fsize,fsize))
#     ax2.spines['right'].set_visible(False)
#     ax2.spines['top'].set_visible(False)
#     ax2.set_xlabel('Days ')

#     ax2.spines['right'].set_visible(False)
#     ax2.spines['top'].set_visible(False)
#     ax2.set_xlabel('Days ')
#     ax2.plot(np.arange(days),bp.prinf[:days],color='olive',label ='Untreated')
#     ax2.plot(np.arange(days),bpol.prinf[:days],':',linewidth=1.9, color='darkorange',label ='Treated')
#     # ax2.plot(np.arange(days-1),dp[:days-1],color='olive',label ='Untreated')
#     # ax2.plot(np.arange(days-1),dpol[:days-1],':',linewidth=1.9, color='darkorange',label ='Treated')


#     ax2.legend(title='Cases daily change')
#     fig.tight_layout() 
# #%% 
#     betaframe = pd.read_csv(outputdir+'pop_clbetas.csv')
#     density = betaframe['density']
#     betalam = betaframe['beta']
    
#     fig,(ax2) = plt.subplots(1,1,figsize=(fsize,fsize))
#     ax2.spines['right'].set_visible(False)
#     ax2.spines['top'].set_visible(False)
#     ax2.set_xlabel('Density')
#     ax2.set_ylabel('Coefficient')
#     ax2.plot(density,betalam/betalam[6],color='darkorange',label='coefficient', marker='o')
#     ax2.plot(density,density,color='olive',label='density')
#     ax2.legend()
    
    
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

#%% Apply policy in SIR using the estimated beta, behavioral
if __name__ == "__main__":
    from class_SIRmodel import SIRmodel
    from class_averageStats import averageStats
    
    thisdensity = 0.5
    shutday = 20
    
    # importing from estimates in sim_estimate_dens.do
    betaframe_beh = pd.read_csv(outputdir+'dens-behbetas.csv')
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
    plt.savefig(imagedir+'estimatedbeta_beh-policies.pdf')

#%% prediction plots
    import pandas as pd
    fsize=3.5
    popsize = 25600
    nobeh = pd.read_csv('output/predictions.csv')
    beh = pd.read_csv('output/predictions_beh.csv')

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

#%% SIR, many betas

if __name__ == "__main__":
    from class_SIRmodel import SIRmodel
    baseline = 25600
    cluster = 10
    beta = .15  #.054*13.5,
    delta = .02 #1/6.5
    
    models = list()
    
    for scale in np.arange(.5,1.6,0.1):
        base = SIRmodel(beta*scale, delta ,q_popsize=baseline,
                q_init=cluster,
                )
        base.plot(base.day)
        models.append(base)    
    
        moddata = pd.DataFrame()

    days = 80
    for idx, model in enumerate(models):
        modeldf = pd.DataFrame([idx]*days,columns=['naics'])
        modeldf['Susceptible']= model.SIR[:days,0]
        modeldf['Density'] = model.beta/beta
        modeldf['Active']= model.SIR[:days,1]
        modeldf['Total'] = model.q_popsize-model.SIR[:days,0]-model.SIR[:days,3]
        modeldf['Locked'] = model.SIR[:days,3]
        modeldf['npi_date']= model.p_shutr[0]
        modeldf['npi_size']= model.p_shutr[1]
        moddata = moddata.append(modeldf)
        
    moddata.to_csv(outputdir+'SIR-dens.csv',index_label='t')


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
    nprocs = 10
    
    iterlist = []
    
    # change seeds
    for citysize in (1/np.sqrt(np.arange(.5, 1.6, .1))):
        for shutdate in np.arange(15,40,5) :
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
    
#%% create csv data
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

#%% Figure

    file = gzip.open(outputdir+'dens-beh_p-times-25pc.gz','rb')
    allRandmodels = pickle.load(file)
    file.close()
    import matplotlib.pyplot as plt
    from class_averageStats import averageStats
    
    bpol = averageStats(allRandmodels[10:11])
    bp = averageStats(allRandmodels[11:12])
    bpol1 = averageStats(allRandmodels[120:121])
    bp1 = averageStats(allRandmodels[121:122])


    # first plot: compare SIR with benchmark 
    dsize = bp.prinf.size-1
    dpol = bpol.prinf[-dsize:] - bpol.prinf[:dsize]
    dp = bp.prinf[-dsize:] - bp.prinf[:dsize]
    # bpol=allRandmodels[0]
    # bp = allRandmodels[5]
    days=80
    fig,(ax2) = plt.subplots(1,1,figsize=(fsize*1.7,fsize))
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')

    avp = (bpol1.prinf*3+bpol.prinf*3)/10
    av = (bp1.prinf*3+bp.prinf*7)/10
    # avp[21:] = (bpol1.prinf[21:]*3+bpol.prinf[21:]*7)/10
    # av[21:] = (bp1.prinf[21:]*3+bp.prinf[21:]*7)/10
    bvp = (bpol1.prinf*7+bpol.prinf*3)/10
    bv = (bp1.prinf*7+bp.prinf*3)/10


    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.plot(np.arange(days),bp1.prinf[:days],color='olive',lw=.5,label ='high density')
    ax2.plot(np.arange(days),bpol1.prinf[:days],':',linewidth=1, color='olive',)
    ax2.plot(np.arange(days),bp.prinf[:days],color='darkorange', lw=.5, label ='low density')
    ax2.plot(np.arange(days),bpol.prinf[:days],':',linewidth=1, color='darkorange',)
    ax2.plot(np.arange(days),bv[:days],color='black', label='true 70-30% avg')
    ax2.plot(np.arange(days),bvp[:days],':',linewidth=1.9, color='black',label= 'true 70-30% treatment')
    ax2.plot(np.arange(days),avp[:days],color='saddlebrown',label ='treated (30% low density)')
    
    # ax2.plot(np.arange(days-1),dp[:days-1],color='olive',label ='Untreated')
    # ax2.plot(np.arange(days-1),dpol[:days-1],':',linewidth=1.9, color='darkorange',label ='Treated')

    ax2.axvline(x=15, color='grey',linewidth=.5 )    
    ax2.axvline(x=40, color='grey',linewidth=.5 )    
    ax2.legend(title='Cases daily % change')
    fig.tight_layout()    
    
    plt.savefig(imagedir+'ddcomp.pdf')
    
    
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


#%% TEST- Apply policy in SIR using the estimated beta, levels 
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
    ax.set_ylim(0,.35)
    ax.set_yticks(np.arange(0,0.16,.05))
    ax.plot(np.arange(20,days),avbase.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax.plot(np.arange(20,days),base.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')
    ax.plot(np.arange(20,days),avpoli.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax.plot(np.arange(20,days),policy.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')
    
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('Days ')
    ax2.set_ylim(0,.35)
    ax2.set_yticks(np.arange(0,0.16,.05))
    ax2.plot(np.arange(20,days),avbase2.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax2.plot(np.arange(20,days),base2.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')
    ax2.plot(np.arange(20,days),avpoli2.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax2.plot(np.arange(20,days),policy2.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.set_xlabel('Days ')
    ax3.set_ylim(0,.35)
    ax3.set_yticks(np.arange(0,0.16,.05))
    ax3.plot(np.arange(20,days),avbase3.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax3.plot(np.arange(20,days),base3.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')
    ax3.plot(np.arange(20,days),avpoli3.prinf[20:days],color='olive',label ='Spatial-SIR')
    ax3.plot(np.arange(20,days),policy3.SIR[20:days,1]/25600,'--',color='darkorange',label ='SIR')
    ax3.hlines(0,20,80,color='grey',lw=0.3)
    ax2.hlines(0,20,80,color='grey',lw=0.3)
    ax.hlines(0,20,80,color='grey',lw=0.3)
    
    #first_legend =ax2.legend(bbox_to_anchor=(.48,.35), fontsize=10)
    l2 = ax2.legend(handles=[],frameon=False, title_fontsize=12, 
                               title='Density 1 \n $\\hat{\\beta} ='+str(round(beta2,2))
                                +'$',#' \n$  \\beta ='+str(round(truebeta2,2))+'$'
                               loc='upper left')
    
    la =plt.legend(bbox_to_anchor=(-.12,1), fontsize=10, framealpha=1)

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
