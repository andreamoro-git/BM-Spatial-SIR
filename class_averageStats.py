#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 10:00:34 2020

This file contains the averageStats class used to compute average statistics 
for several replications of each simulation. 

The class accepts as input a list of objects of the spatialSAYDR class
"""
import numpy as np
from datetime import datetime

class averageStats() :
    def __init__(self,models,caption='') :
        self.classname = type(models[0]).__name__
        self.caption = caption
        self.models = models
        
        m = self.models
        
        self.popsize = m[0].q_popsize
        
        # these are legacy codes
        try :
            self.lastday = list(map(lambda x: x.lastday,m))
        except : # legacy
            for model in m :
                model.lastday = np.max(np.where(model.prtinf==max(model.prtinf)))
            self.lastday = list(map(lambda x: x.lastday,m))

        if (max(self.lastday)-min(self.lastday)>40) :
            print ('warning: check last days very different')

        self.minlastday = min(self.lastday)
        self.R0 = np.average(list(map(lambda x: np.max(x.R0),m)))
        self.R0std = np.std(list(map(lambda x: np.max(x.R0),m)))
        self.igrowth = np.average(list(map(lambda x: x.igrowth,m)),axis=0)    
        self.prtinf = np.average(list(map(lambda x: x.prtinf,m)),axis=0)    
        self.prinf = np.average(list(map(lambda x: x.prinf,m)),axis=0)   
        activeByType = np.array(list(map( lambda x: np.sum(x.nstatesByType[:,:,1:3],axis=2),m)))
        self.frActiveByType = np.average(activeByType,axis=0)/m[0].q_typeprobs/m[0].q_popsize
        self.maxCases = np.average(list(map(lambda x: max(x.prinf),m)),axis=0)
        self.maxCasesByType = np.max(self.frActiveByType,axis=0)
        try: 
            self.prasymp = np.average(list(map(lambda x: x.prstates[:,x.I],m)),axis=0) # fraction of asymptomatics each day
            self.prsymp = np.average(list(map(lambda x: x.prstates[:,x.Y],m)),axis=0)  # fraction of symptomatics each day
        except: #legacy
            try : #legacy stats with old name
                self.frasymp = np.average(list(map(lambda x: x.frasymp,m)),axis=0)    
            except : 
                print('frasymp not found')
                self.frasymp = np.average(list(map(lambda x: x.prasymp,m)),axis=0)    
            self.prasymp = self.frasymp * self.prinf # fraction of asymptomatics each day
            self.prsymp = self.prinf - self.prasymp  # fraction of symptomatics each day

        self.ninstate = np.average(list(map(lambda x: x.nstates[x.lastday]/x.q_popsize, m)),axis=0)
        self.ninstateByType = np.average(list(map(
            lambda x: 
                np.transpose(np.transpose(x.nstatesByType[x.lastday])/np.sum(x.nstatesByType[x.lastday]
                                                                ,axis=1)), 
                    m)),axis=0)
        try: 
            IFR = list(map(lambda x: np.max(x.prstates[:,x.D]),m))
            np.average(list(map(lambda i: IFR[i]/m[i].prtinf[m[i].lastday],np.arange(len(m)))))
        except: 
            #self.IFR = list(map(lambda x: np.max(x.prdead),m))
            pass
        self.maxasy = np.average(list(map(lambda x: x.maxinf,m)))/self.popsize
        self.maxdays = np.average(list(map(lambda x: np.max((x.tdea,x.trec)),m)))
            
        if self.classname != 'spatialSAYDR' and self.classname !='spSAYDR_randLoc' and self.classname!='spSAYDR_hetDensity' and self.classname != 'class_averageStats' and self.classname != 'simul_policy':
            
            self.fracNotScared = np.average(list(map(lambda x: x.fracNotScared,m)),axis=0)    
            try: 
                self.whereInfected = np.average(list(map(lambda x: x.whereInfected,m)),axis=0)
                self.whereInfByTypeLast = np.average(list(map(
                    lambda x: x.whereInfByType[x.lastday],m )),axis = 0)  
                self.whereInfectedLast =  np.average(list(map(
                    lambda x: x.whereInfected[x.lastday],m )),axis = 0)     
                
                # report only conditional on being infected
                self.whereInfByTypeLast = np.transpose(np.transpose(self.whereInfByTypeLast[:,1:])/np.sum(self.whereInfByTypeLast[:,1:],axis=1))
                self.whereInfectedLast =  self.whereInfectedLast[1:] /np.sum(self.whereInfectedLast[1:])

            except:
                print('no where')
                pass
            try:
                self.p_oldhomesize = m[0].p_oldhomesize
            except:
                self.p_oldhomesize = 0

        
    def printout(self):
        print('ha')
        returndic = {
                 'igrowth': self.igrowth,
                 'ninstate': self.ninstate,
                 'prinf': self.prinf,
                 'prtinf': self.prtinf,
                 'prasymp': self.prasymp,
                 'maxCases': self.maxCases,
                 'IFR': self.IFR,
                 'popsize': self.popsize,
                 'R0': self.R0,
                 'R0std': self.R0std,
                 'maxAsy': self.maxasy,
                 'maxdays': self.maxdays,
                 'p_oldhomesize' : self.p_oldhomesize,
                 }
        if self.classname != 'spatialSAYDR':
            print('self.classname')
            returndic.update({'fracNotScaredFirms': self.fracNotScaredFirms,
                 'fracNotScared': self.fracNotScared,
                 'whereInfected': self.whereInfected,
                 })
        return (returndic)

    def printLaTeX(self):
        # prepare table
        texstring = '%% \\usepackage{array,booktabs,threeparttable}'
        texstring += '\n%%%\n\\begin{table}'
        texstring += '\n\\caption{Counterfactual: '+self.caption+' '
        texstring += str(datetime.now()) + '} \label{tab:'+str(datetime.now())+'}'
        texstring += '\n \\centering \n \\begin{threeparttable} \n\t \\begin{tabular}'
        texstring += '{m{0.20\linewidth}>{\\raggedright}m{0.08\linewidth}>{\centering}m{0.08\linewidth}>{\centering}m{0.08\linewidth}>{\centering}m{0.08\linewidth}>{\centering}m{0.08\linewidth}>{\centering \\arraybackslash}m{0.08\linewidth}}' 
        texstring += '\n\t\\toprule \n'
        texstring += '& \multicolumn{3}{c}{Infection location} & \multicolumn{3}{c}{Steady-state outcomes}\\tabularnewline \n'
        texstring += '& City & W/S & Home & D & R & Peak A+Y \\tabularnewline \n'
        texstring += '\t\\midrule \n'
        
        # aggregate outcome
        statespr = self.ninstate
        locinf = self.whereInfectedLast
        texstring += '\t\tAll \t& %4.3f \t& %4.3f \t& %4.3f \t& %4.3f \t& %4.3f \t& %4.3f \\tabularnewline \n\t\\midrule' % (
            locinf[0],locinf[1],locinf[2],statespr[3],statespr[4],self.maxCases
            )
        
        # type outcomes
        for ttype in (['Young',0],['Not employed',1],['Old',2]) :
            statespr = self.ninstateByType[ttype[1],:]
            locinf = self.whereInfByTypeLast[ttype[1]]
            texstring += '\n \t\t '+ ttype[0]+' \t\t& %4.3f \t& %4.3f \t& %4.3f \t& %4.3f \t& %4.3f \t& %4.3f \\tabularnewline ' % (
                locinf[0],locinf[1],locinf[2],statespr[3],statespr[4],self.maxCasesByType[ttype[1]]
                )
        #     states = self.nstatesByType[day,1,:]/np.sum(self.type==1)
        # locinf = self.whereInfByType[day,1,1:]/(1-self.whereInfByType[day,1,0])
        # texstring += 'Not employed & %4.3f & %4.3f & %4.3f & %4.3f & %4.3f \\tabularnewline \n' % (locinf[0],locinf[1],locinf[2],statespr[3],statespr[4])
        # states = self.nstatesByType[day,1,:]/np.sum(self.type==1)
        # locinf = self.whereInfByType[day,2,1:]/(1-self.whereInfByType[day,2,0])
        # texstring += 'Old & %4.3f & %4.3f & %4.3f & %4.3f & %4.3f \\tabularnewline \n' % (locinf[0],locinf[1],locinf[2],statespr[3],statespr[4])
        
        # close down table
        texstring += '\n\t \\bottomrule'
        texstring += '\n\t \\end{tabular} \n\\end{threeparttable} \n\\end{table} \n%%%\n'
     
        return texstring
   
#%%#--------------------------------------------------------------------------#
if __name__ == "__main__":

    import pickle
    import matplotlib.pyplot as plt
    file = open('output/simple-periods.pickle', 'rb')
    b = pickle.load(file)
    
    bavg = averageStats(b)
    
    print(bavg.printout())
    fig,ax = plt.subplots(figsize=(5,5))

    g_plotdays = 205
    g_maxinf = 1.3
    ### draw legend with empty data
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks(np.arange(0,g_plotdays,20))
    ax.set_yticks(np.arange(0,1,.2))
    ax.set_ylim(0,g_maxinf)
    ax.set_xlabel('Days ')
    ax.set_ylabel('Percent infected')
    
    bavg.prtinf[120:] = bavg.prtinf[120]
    print('\n*********\nWarning: hand-changed prtinf value aaa\n\n*********')
    
    # draw the actual chart
    # draw the actual chart. Rescale position so that it fits the square properly
    ax.plot(np.arange(g_plotdays),bavg.prtinf[0:g_plotdays],c='olive',label="Baseline")
    ax.plot(np.arange(g_plotdays),b[0].prtinf[0:g_plotdays],c='chocolate',label="Baseline")
    
    ax.plot(np.arange(g_plotdays),bavg.prinf[0:g_plotdays],'--',c='olive')
    ax.plot(np.arange(g_plotdays),b[0].prinf[0:g_plotdays],'--',c='chocolate')
    
    #ax.plot(np.arange(g_plotdays),bavg['igrowth'][0:g_plotdays],'-.',c='chocolate')
    # ax.plot(np.arange(g_plotdays),b2avg['igrowth'][0:g_plotdays],'-.',c='olive')
    #ax.plot(np.arange(g_plotdays),havg['igrowth'][0:g_plotdays],'-.',c='DarkKhaki')
    
    #ax.plot(np.arange(40),m2.mldr[0:40],'-.',c='gray')
    ax.plot([],[],'--',c='black',label='Active')
    #ax.plot([],[],'-.',c='black',label='Growth rate')
    
    ax.legend(loc='upper left',ncol=1,title='Percent Infected')
    

def avgDist(radius,pos,maxc='') :
    
    q_popsize = pos.shape[0]
    if type(maxc) == str :
        maxc = q_popsize
    
    np.random.shuffle(pos)
    sampsize = min(q_popsize,maxc)
    ncontacts = np.zeros(sampsize)
    distance = np.zeros(sampsize)
    for i in range(sampsize) :
        if np.mod(i,1000)==0 :
#            print(i)
            pass
        disti = np.linalg.norm(pos-pos[i],axis=1)
        distance[i] = np.average(disti)
        ncontacts[i] = np.sum(disti <= radius) - 1
    return (np.sum(distance)/sampsize,np.sum(ncontacts)/sampsize)

vecavgDist = np.vectorize(avgDist,excluded=['pos','maxc'])
