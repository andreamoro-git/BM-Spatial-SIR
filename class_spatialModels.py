#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 10:42:38 2020

Defines the following classes:

1) spatialSAYDR 
    Contains all methods necessary to simulate a 5-state (behavioral) spatial-SIR model, 
    and generates object attributes containing contagion statistics. 
    Simulation data is saved in object attributes.
2) spSAYDR_randLoc (child of spatialSAYDR)
    Simulates spatial-SIR with agents moved in randomly drawn locations every day
3) spSAYDR_hetDensity (child of spatialSAYDR)
    Simulates spatial-SIR with initial density of agents decreasing from center of city
4) spSAYDR_behav (child of spatialSAYDR)
    Simulates spatial-SIR with behavioral responses
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from multiprocessing import freeze_support

class spatialSAYDR() :
    
    """ 
        Class for simulating a behavioral spatialSAYRD model 
        with multiple types
        
    """
    
    def __init__(self,
                 
        # basic parameters 
        q_printOption = 1,              # global print option
        q_debug     = 0,
        q_days      = 300,              # max number of days to converge
        q_popsize   = 160**2,           # size of population
        q_citysize  = 1,                # side of city square
        
        q_firstinf  = [[0.25,0.25]],    # location first cluster of infection
        q_typeprobs = [1],              # fraction of people in each type (must sum to 1)
        q_seed      = 2443,             # randomization seed
        
        # graphic parameters
        g_plotdays  = 200,              # of days in plot
        g_maxinf    = 1.,               # scale of y axis in plot
        
        # other model parameters
        p_avgstep   = 3.5*.00805,       # mean travel distance
        p_stdstep   = .0,               # sd travel distance 
        p_infradius = 0.00805,          # min distance for contagion
        p_timerec = 100,                # days until recover for sure
        p_timey = 1,                    # min days to become symptomatic
        p_probc = [[0, 0.14, 0, 0, 0]], # probability of contagion by type (and status)
        p_proby = [0.09],               # prob to become symptomatic by type
        p_probr = [0.05],               # recover probab by type
        p_timedea = 3,                  # min n. days to die after symptoms
        p_probd = [0.001/7.5],          # death probab by type
        p_ininf = 30,                   # number of initial infected
        p_cluster = 'cluster',          # position of initial infected (cluster, special or random)
        p_clustermask = [False],        # position of initial infected if p_cluster=special
        p_agilityS = np.array([1,1,1,0,1]),  # limit ability to move by state
        p_agilityT = np.array([1]),          # limit ability to move by type
    ):
        
        ## import parameters
        self.q_printOption = q_printOption
        self.q_debug = q_debug
        self.q_days = q_days
        self.q_popsize = q_popsize
        self.q_citysize = q_citysize
        self.p_avgstep = p_avgstep
        self.p_stdstep = p_stdstep
        self.q_firstinf = q_firstinf
        self.q_typeprobs = np.asarray(q_typeprobs)
        self.q_seed = q_seed
        if g_plotdays > q_days:
            self.g_plotdays = q_days
        else :
            self.g_plotdays = g_plotdays
        self.g_maxinf = g_maxinf
        self.p_infradius = p_infradius
        self.p_timerec = p_timerec
        self.p_timedea = p_timedea
        self.p_probc = np.array(p_probc)
        self.p_proby = np.array(p_proby)
        self.p_probd = np.array(p_probd)
        self.p_probr = np.array(p_probr)
        self.p_ininf = p_ininf
        self.p_timey = p_timey
        self.p_cluster = p_cluster
        self.p_agilityS = np.array(p_agilityS)
        self.p_agilityT = np.array(p_agilityT)
        
        ##define variables       
        np.random.seed(self.q_seed)
        
        # state codes
        self.S = 0      #susceptibles
        self.I = 1      #infected
        self.Y = 2      #symptomatic
        self.D = 3      #dead
        self.R = 4      #recovered
        
        # id
        self.id = np.arange(self.q_popsize)
        np.random.shuffle(self.id)
        
        # assign types
        self.ntypes = len(self.q_typeprobs)
        self.type = np.random.choice(self.ntypes,q_popsize,p=self.q_typeprobs)
        
        # this is in case one forgets to set the agility by type vector
        if self.ntypes > 1 :
            if len(self.p_agilityT)==1 :
                self.p_agilityT = np.array([self.p_agilityT[0]]*self.ntypes)
            if len(self.p_probc)==1 :
                self.p_probc = np.array([self.p_probc[0]]*self.ntypes)
            if len(self.p_probd)==1 :
                self.p_probd = np.array([self.p_probd[0]]*self.ntypes)
            if len(self.p_probr)==1 :
                self.p_probr = np.array([self.p_probr[0]]*self.ntypes)
            if len(self.p_proby)==1 :
                self.p_proby = np.array([self.p_proby[0]]*self.ntypes)
        
        # initial setup of people's position
        self.pos = np.random.rand(q_popsize,2)*q_citysize

        # threshold for when people  are too close to the boundary
        self.q_bthr = 0.5*self.p_avgstep*1/np.sqrt(self.q_popsize)
    
        # state value of each individual
        self.state = np.zeros(q_popsize,dtype=int)

        # time of infection, symptoms, death and recovery
        self.tinf = np.zeros(self.q_popsize,dtype=int)-9
        self.tsym = np.zeros(self.q_popsize,dtype=int)-9
        self.tdea = np.zeros(q_popsize,dtype=int)-9
        self.trec = np.zeros(q_popsize,dtype=int)-9
        
        # statistics
        self.prinf = np.zeros(q_days)                       #proportion of infected (out of living)
        self.prtinf = np.zeros(q_days)                      #proportion of cases (out of total pop)
        self.dgrowth = np.zeros(q_days)                     #death growth
        self.igrowth = np.zeros(q_days)                     #infected growth
        self.nstates = np.zeros((q_days,5),dtype=int)       #n. in each state
        self.prstates = np.zeros((q_days,5))                #fraction in each state
        self.nstatesByType = np.zeros((q_days,self.ntypes,5),dtype=int)
        self.maxinf = 0                                     #max n. of infected
        self.ninf = np.zeros(q_popsize,dtype=int)           # number of ppl infected by each person
        self.R0 = 0                                         #estimate of R0
        
        ## initial setup of infected 
        # set their location depending on model choice 
        if self.p_cluster == 'cluster' :
            for firstpositions in self.q_firstinf :
                distance = np.linalg.norm(self.pos-firstpositions,axis=1)
                argdist = np.argsort(distance)
                self.state[argdist[0:self.p_ininf]] = self.I
        elif self.p_cluster == 'special' :
            self.state[self.clustermask] = self.I
        else :  # random position
            self.state[np.random.choice(self.q_popsize,self.p_ininf)] = self.I
        
        # set their initial time of infection and state
        self.tinf[self.state==self.I] = 0
        self.nstates[0,0:2] = [np.sum(self.state==self.S),np.sum(self.state==self.I)]
        
        ## read lombardy's data
        self.mldr = np.zeros(self.q_days)
        milan = open('input/drlombardia.txt','r')
        for idx,line in enumerate(milan) :
            if idx < self.q_days: 
                self.mldr[idx] = line
        return
    
    def movepeople(self,positions):
        ''' move people around from their initial positions 
            arguments 
                positions: initial position
            returns 
                array of new positions
        '''
        
        q_bthr = self.q_bthr
        pos = positions
        newpos = np.zeros((self.q_popsize,2))
        
        # compute distance and correct by mobility
        traveldist = np.random.normal(self.p_avgstep,self.p_stdstep,self.q_popsize)
        traveldist = traveldist*self.p_agilityS[self.state]
        traveldist = traveldist*self.p_agilityT[self.type]
        
        traveldir = np.random.rand(self.q_popsize)*2*np.pi #draw angle of direction   
        
        # change direction if too close to boundary
        ubound = self.q_citysize
        mask = (pos[:,0]<q_bthr)
        traveldir[mask] = np.random.rand(sum(mask))*np.pi-np.pi/2
        mask = (pos[:,0]>ubound-q_bthr)
        traveldir[mask] = np.random.rand(sum(mask))*np.pi+np.pi/2
        mask = (pos[:,1]<q_bthr)
        traveldir[mask] = np.random.rand(sum(mask))*np.pi
        mask = (pos[:,1]>ubound-q_bthr)
        traveldir[mask] = np.random.rand(sum(mask))*np.pi+2*np.pi/2
                  
        newpos[:,0] = pos[:,0] + traveldist[:]*np.cos(traveldir[:])
        newpos[:,1] = pos[:,1] + traveldist[:]*np.sin(traveldir[:])
       
        # correct if out of bounds (this may be handled better but it's good enough)
        mask = (newpos[:,0]<0)  
        newpos[:,0][mask] = np.minimum(ubound,-newpos[:,0][mask])
        mask = (newpos[:,0]>ubound)
        newpos[:,0][mask] = np.maximum(0,2*ubound-newpos[:,0][mask])
        mask = (newpos[:,1]<0)
        newpos[:,1][mask] = np.minimum(ubound,-newpos[:,1][mask])
        mask = (newpos[:,1]>ubound)
        newpos[:,1][mask] = np.maximum(0,2*ubound-newpos[:,1][mask])
            
        return newpos     

    def infections(self,today,maskdata=['']) :

        ''' who gets infected today? changes state and tinf 
            arguments
                today: what day is today
                maskdata: array of ids of who can get infected (default: everyone)
            returns
                array of size q_popsize =1 if person is newly infected, 0 otherwise
        '''
        
        pos = self.pos
        state = self.state
        ttype = self.type
        
        newinf = np.zeros(self.q_popsize,dtype=int)
        
        # get the list of id to loop through. Either everyone or only the list passed on
        if maskdata[0] == '' :
            infectives = np.arange(self.q_popsize)
        else :
            infectives = np.ma.nonzero(maskdata)
            
        #each individual i that is contagious infects his neighbors
        for i in infectives :
            p_probc = self.p_probc[ttype[i],state[i]]
            if p_probc > 0 :

                #compute distance from all individuals
                dist_i = np.linalg.norm(pos-pos[i],axis=1)
                
                #draw infection only from the closest, and if they are not infected already
                mask = ((dist_i<self.p_infradius) & (state==self.S) & (newinf==False))
                draw = np.random.choice(2,np.sum(mask),p=[1-p_probc,p_probc])
                newinf[mask]= self.I * draw
                self.ninf[i] += np.sum(draw)
                        
        if (sum((newinf)&(self.state>0))) >0 :
            print ('WARNING: check if some newly infected is not susceptible')
        
        # set state and time of infection
        self.tinf[newinf==self.I] = today 
        self.state[newinf==self.I] = newinf[newinf==self.I]
        return newinf

    def symptoms(self,today) :
        
        ''' who becomes symptomatic today? changes state and tsym '''
        
        for ttype in range(self.ntypes):
 
            # select symptomatic of each type such that min days for death since inf
            typemask = (self.type == ttype)
       
            p_proby = self.p_proby[ttype] 
        
            mask = (typemask & (self.state==self.I) & (today-self.tinf>self.p_timey))
            draw = np.random.choice(2,np.sum(mask),p=[1-p_proby,p_proby])
        
            self.state[mask] = self.Y*draw + self.state[mask]*(1-draw)
            self.tsym[mask] = today*draw + (-9)*(1-draw)
    
        return 

    def deaths(self,today) :
        
        ''' who dies today? changes state and tdea '''
      
        for ttype in range(self.ntypes):
            # select symptomatic of each type such that min days for death since inf
            typemask = (self.type == ttype)
       
            probd = self.p_probd[ttype]
            mask = (typemask & (self.state==self.Y) & (today-self.tsym>self.p_timedea))
            draw = np.random.choice(2,np.sum(mask),p=[1-probd,probd])
            
            self.state[mask] = self.D*draw + self.state[mask]*(1-draw)
            self.tdea[mask] = today*draw + (-9)*(1-draw)
        
        return 

    def recoveries(self,today):
        
        ''' who recovers today? changes state and trec '''
        
        for status in (self.I,self.Y): #only infected or symptomatics
            # recover for sure after p_timerec days
            mask = ((self.state==status) & (today-self.tinf>self.p_timerec))
            self.state[mask] = self.R
            self.trec[mask] = today
        
            # ... or recover with prob p_probr before 
            for ttype in range(self.ntypes):
                # select symptomatic of each type such that min days for death since inf
                typemask = (self.type == ttype)
           
                p_probr = self.p_probr[ttype]
            
                mask2 = (typemask & (self.state==status))
                draw = np.random.choice(2,np.sum(mask2),p=[1-p_probr,p_probr])
    
                self.state[mask2] = self.R*draw + self.state[mask2]*(1-draw)
                self.trec[mask2] = today*draw + (-9)*(1-draw)
            
        return

    def drawPositions(self,day,maskdata='',
                      savefig='',folder='output/images/',
                      S=False,nojunk=False,
                      colors=['lightgrey','orange','orange','red','green'],
                      focus=[],
                      ext='.png') :

        ''' draw the city with everyone's position and color-coded state '''
        pos = self.pos
        state = self.state
        focus = np.array(focus)
        
        # allow for plotting subsets of data
        if type(maskdata)==str :
            maskdata = np.array([True]*self.q_popsize)
        
        if day>self.g_plotdays:
            self.g_plotdays = day
          
        fig,ax = plt.subplots(figsize=(3.3,3.3))
        
        ### draw legend with empty data 
        if nojunk==False :
            ax.set_ylim(0,self.q_citysize)
            #ax.set_xlabel('Day '+str(day))
            ax.set_xlim(0,self.q_citysize)
            ax.set_xticks([])
            ax.set_yticks([])
        if focus.size>1:
            ax.set_xlim(focus[0,0],focus[0,1])
            ax.set_ylim(focus[1,0],focus[1,1])

        # draw the actual chart. Rescale position so that it fits the square properly
        dpos = np.copy(pos[maskdata])
        dstate = np.copy(state[maskdata])

        if S==True :       
            ax.scatter([],[],c=colors[0],s=5,label='Susceptible ') #+ str(np.round(self.prstates[day,0],3)))
            ax.scatter(dpos[:,0][dstate==self.S],dpos[:,1][dstate==self.S],c=colors[0],s=.1)
        ax.scatter(dpos[:,0][dstate==self.I],dpos[:,1][dstate==self.I],c=colors[1],s=1,) # + str(np.round(self.prstates[day,1],3)))
        ax.scatter(dpos[:,0][dstate==self.R],dpos[:,1][dstate==self.R],c=colors[4],s=1,) # + str(np.round(self.prstates[day,4],3)))
        ax.scatter(dpos[:,0][dstate==self.Y],dpos[:,1][dstate==self.Y],c=colors[2],s=1)
        ax.scatter(dpos[:,0][dstate==self.D],dpos[:,1][dstate==self.D],c=colors[3],s=1)
        ax.scatter([],[],c=colors[1],s=5,label='Infected')
        ax.scatter([],[],c=colors[4],s=5,label='Recovered')
        ax.legend(loc='upper left',title='Day '+str(day),title_fontsize=14,fontsize=12,framealpha=0.85,labelspacing=0.14)
        fig.tight_layout()

        if savefig != '' :
            plt.savefig(folder+savefig+'_pos-day'+str(day)+ext)
        plt.show()    
  
    def drawRates(self,day,lines='',ylim='',aaamaskdata=[''],savefig='',folder='output/images/') :

        ''' chart with growth rates '''
        
        if day>self.g_plotdays:
            self.g_plotdays = day
            
        fig,ax = plt.subplots(figsize=(5,5))
        
        ### draw legend with empty data
        ax.plot([],[],label='Lombardy',c='olive')
        ax.plot([],[],label='Growth rate',c='DarkOrange')
        ax.legend(loc='upper right', title='Day '+str(day))
        fig.tight_layout()
        if day<120 :
            tickrange = 10
        elif day<200 :
            tickrange = 20
        else :
            tickrange = 50
        ax.set_xticks(np.arange(0,max(self.g_plotdays,day),tickrange))
        ax.set_yticks(np.arange(0,self.g_maxinf,.1))
        ax.set_ylim(0,self.g_maxinf)
        ax.set_xlabel('Day')
        ax.set_ylabel('Growth rate or fraction infected')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)


        df = pd.DataFrame(self.igrowth)
        ma = df.rolling(5).mean()

        ax.plot(np.arange(30),self.mldr[0:30],c='olive')
        ax.plot(np.arange(day-2),self.igrowth[2:day],c='DarkOrange')
        if day>=5 :
             ax.plot(np.arange(day-4),np.array(ma[4:day][0]))      
             ax.plot(np.arange(day),self.prinf[0:day])
        ax.plot(np.arange(day),self.prtinf[0:day])
        fig.tight_layout()
        if savefig != '' :
            plt.savefig(folder+savefig+'_rates-day'+str(day)+'.png')
        plt.show()      

    def drawChart(self,day,maskdata=[''],savefig='',folder='output/images/') :
        
        ''' figure combining rates and positions '''
        
        pos = self.pos
        state = self.state
        
        # allow for plotting subsets of data
        if maskdata[0]=='' :
            maskdata = np.array([True]*self.q_popsize)
        
        if day>self.g_plotdays:
            self.g_plotdays = day
            
        fig,ax = plt.subplots(figsize=(5,5))
        
        ### draw legend with empty data
        ax2 = ax.twinx()
        ax2.scatter([],[],c='orange',s=3,label='Infected')
        ax2.scatter([],[],c='red',s=3,label='Dead')
        ax2.scatter([],[],c='green',s=3,label='Recovered')
        prtinf = np.copy(self.prtinf)
        prtinf[day:] = np.nan
        ax.plot(np.arange(self.g_plotdays),prtinf[0:self.g_plotdays],c='black')
        
        if day<120 :
            tickrange = 10
        elif day<200 :
            tickrange = 20
        else:
            tickrange = 50
        ax.set_xticks(np.arange(0,max(self.g_plotdays,day),tickrange))
        ax.set_ylim(0,self.g_maxinf)
        ax.set_xlim(0,self.g_plotdays)
        ax.set_xlabel('Days')

        ax2.set_ylabel('% Cases')
        ax2.set_ylim(0,self.g_maxinf)

        # draw the actual chart. Rescale position so that it fits the square properlyns
        dpos = np.copy(pos[maskdata])
        dpos[:,0] = dpos[:,0]*max(day,self.g_plotdays)/self.q_citysize
        dpos[:,1] = dpos[:,1]*self.g_maxinf/self.q_citysize
        dstate = np.copy(state[maskdata])
        
#        ax.scatter(dpos[:,0][dstate==self.S],dpos[:,1][dstate==self.S],c='black',s=.1)
        ax2.scatter(dpos[:,0][dstate==self.I],dpos[:,1][dstate==self.I],c='orange',s=1)
        ax2.scatter(dpos[:,0][dstate==self.Y],dpos[:,1][dstate==self.Y],c='orange',s=1)
        ax2.scatter(dpos[:,0][dstate==self.D],dpos[:,1][dstate==self.D],c='red',s=3)
        ax2.scatter(dpos[:,0][dstate==self.R],dpos[:,1][dstate==self.R],c='green',s=.2)
  
        ax2.plot([],[],c='black',label='Percent total cases')

        ax.set_yticks(np.arange(0,0.1,0.7))
        ax2.legend(loc='upper left',title='Day '+str(day),framealpha=0.8)
        if savefig != '' :
            plt.savefig(folder+savefig+'_pos-day'+str(day)+'.png')
        plt.show()
    
    def computeStats(self, day, printOption=0) :
        
        ''' compute and save a bunch of stats, print them if requested '''
        state = self.state

        nstates = [np.sum(state==0),np.sum(state==1),np.sum(state==2),np.sum(state==3),np.sum(state==4)]
        self.nstates[day] = nstates
        self.prstates[day] = self.nstates[day]/self.q_popsize
        self.maxinf = max(self.maxinf,np.sum(state==self.I))
        self.prinf[day] = sum(self.prstates[day,1:3])
        
        for tt in range(self.ntypes) :
            statet = self.state[self.type==tt]
            nstatest = [np.sum(statet==0),np.sum(statet==1),np.sum(statet==2),np.sum(statet==3),np.sum(statet==4)]
            self.nstatesByType[day,tt,:] = nstatest    
        
        if day>=1 :
            self.prtinf[day] = 1-self.prstates[day,0]
            prdead = self.prstates[:,self.D]/(1-self.prtinf[day])
            if day>= 2: 
                self.igrowth[day] = (self.prtinf[day]-self.prtinf[day-1])/self.prtinf[day-1]
            if day>=(self.p_timedea+1) :
                if prdead[day-1] > 0 :             
                    self.dgrowth[day] = (prdead[day]-prdead[day-1])/prdead[day-1]   
                else : 
                    self.dgrowth[day] = np.nan
            try : 
                self.R0 = np.average(self.ninf[(self.tinf==0) & (self.tinf<=5)])
            except :
                self.R0 = 0
        prtstring = ''
        if printOption >=1 :
            prtstring += '\nDay %3d: R0: %3.2f, Active %4.3f, Total %4.3f, CDR %4.3f' % (day,self.R0,self.prinf[day],self.prtinf[day],self.prstates[day][self.D])
            prtstring += ("\n         [%5d %5d %5d %5d %5d] [S I Y D R]" 
                  % (nstates[self.S],nstates[self.I],nstates[self.Y],nstates[self.D],nstates[self.R]))
            prtstring += '\n         ' + str(np.round(self.nstates[day]/self.q_popsize,3)) + ' [S I Y D R]'
            print(prtstring)
        return

    def summaryStats(self,day,nondef='',savefile='',clear='a') :
        
        ''' print summary stats '''
        
        string = ''
        if nondef != '' :
            string += "\n_______________________________________________________\n"
            string += datetime.today().strftime('%d/%m/%Y %H:%M:%S')
            string += "\nModel "+type(self).__name__
            string += ("\nConvergence complete with non-default arguments: %s" % nondef)
        state = self.state
        string += ("\n\nModel "+type(self).__name__+" summary statistics day %d: \nMax # of infected: %d, Infected: %d, Deaths: %d \nR0: %3.2f, CDR: %5.4f"
               % (day,self.maxinf,self.q_popsize-np.sum(state==self.S),np.sum(state==self.D),
                 self.R0,self.prstates[day][self.D]))
        string += ('          State rates: ')+str(np.round(self.prstates[day],3))
    
        print(string)
        if savefile :
            text_file = open(savefile, clear) 
            text_file.write(string)
            text_file.close()
            
        return string
    
    def aDayInTheLife(self,day):
        
        ''' what happens in a given day '''
        
        # get infected, become symptomatic, deaths and recoveries
        self.infections(day)
        self.symptoms(day)
        self.deaths(day)
        self.recoveries(day)      ,
                    
        # move people around
        self.pos = self.movepeople(self.pos)
        self.lastday = day

    def simulateOnce(self) :
        
        ''' simulates the model for q_days, to be used in conjunction
            with multiprocessing.Pool to get multiple replications quickly 
        '''
        
        stats = self.computeStats(0)
        if self.q_printOption >= 1 :
            print(stats)
        
        for day in range(1,self.q_days):

            stats = self.computeStats(day,self.q_printOption)
            self.aDayInTheLife(day)

            if self.q_printOption >= 1: 
                print(stats)
            if self.q_printOption >= 2:
                self.drawRates(day,lines=['Y','I','active'],ylim=0.4)
            if self.q_printOption >= 3:
                self.drawPositions(day)

            if np.sum(self.state==self.I)+np.sum(self.state==self.Y) <= 15:
                break
        
        if self.q_printOption>=0.5:
            print('Ended '+type(self).__name__+' lastday: '+str(self.lastday))        
        
        summary = (self.summaryStats(day))
        return summary

class spSAYDR_randLoc(spatialSAYDR):
    
    ''' this class simulates the contagion model but puts people
        in random locations every day
    '''

    def movepeople(self,positions):
        newpos = np.random.rand(self.q_popsize,2)*self.q_citysize

        return newpos

class spSAYDR_hetDensity(spatialSAYDR):
    
    ''' this class simulates the contagion model where people are 
        placed initially in a city where density decreases from center
        to periphery
        
        parameter: q_lambda, standard deviation of distance from center
    '''
    
    def __init__(self,q_lambda,*args,**kwargs) :
        
        self.q_lambda = q_lambda
        kwargs['q_firstinf']=[[-0,0]]
        
        # since __init__ is redefined we need to initialize the parent class
        # to load all attributes (typically defined as kwargs) from the parent class
        spatialSAYDR.__init__(self,
                                *args,**kwargs)
        
        # move to position with distance from center normally distributed
        # iterate to make sure everyone is inside the unit square
        # for comparison with baseline model
        outsideSquare = [True]*self.q_popsize
        while sum(outsideSquare) > 0:

            # draw distance and direction
            traveldist = np.random.normal(loc=0,scale=self.q_lambda,size=self.q_popsize)
            traveldir = np.random.rand(self.q_popsize)*2*np.pi #draw angle of direction   
       
            # position people (center is [0.5,0.5] for conmparison with baseline)
            self.pos[outsideSquare,0] = 0.5 + traveldist[outsideSquare]*np.cos(traveldir[outsideSquare])
            self.pos[outsideSquare,1] = 0.5 + traveldist[outsideSquare]*np.sin(traveldir[outsideSquare])
            
            # determine who is outside the unit aquare
            outsideSquare = (self.pos[:,0]>1) | (self.pos[:,0]<0) | (self.pos[:,1]>1) | (self.pos[:,1]<0)
                   
        # set up initially infected near the center
        self.state[:] = self.S
        distance = np.linalg.norm(self.pos-[0.5,0.5],axis=1)
        argdist = np.argsort(distance)
        self.state[argdist[0:self.p_ininf]] = self.I
        self.tinf[self.state==self.I] = 0
 
    def drawPositions(self,day,kwargs={'focus':[[0,1],[0,1]], 'S': True}):
        super().drawPositions(day=day,**kwargs)


class spSAYDR_behav(spatialSAYDR):
    
    ''' extends spatialSAYDR to allow for behavioral effect 
        arguments:
            behModel: dictionary with types and params
            behModel['type'] = 'Lones','Jesus'
    '''

    def __init__(self,
                 behModel='',
                 *args,
                 **kwargs
                 ) :
        
        # since __init__ is redefined we need to initialize the parent class
        # to load all attributes (typically defined as kwargs) from the parent class
        spatialSAYDR.__init__(self, 
                                *args,**kwargs
                                )
            
        self.behModel = behModel
        
        # extra vars and stats to keep track of
        self.fracNotScared = np.zeros(self.q_days)+1        #nunmber of scared
        self.scared = np.zeros(self.q_popsize,dtype=bool)   #who is scared
        
    def scaredycats(self,day) :
        
        '''
        Detects people scared to go outside so that infections() can 
        put them in a land far and far away

        '''
        par = self.behModel
        
        # figure out how many are infected and compute how many are scared 
        
        if par['type'] == 'Lones':
            
            fY = self.prstates[day,self.I] + self.prstates[day,self.Y]
            if fY <= par['phi'] :
                reducedContagion = 1.
            else:
                reducedContagion = (par['phi']/fY)**0.12        
                
        elif par['type'] == 'Jesus':
            
            b0 = par['beta0']
            bstar = par['betastar']
            lambd = par['lambda']
            reducedContagion = b0 * np.exp(-par['lambda']*day) + bstar * (1-np.exp(-lambd*day))
        
        self.fracNotScared[day] = reducedContagion
        scaredcats = np.round((1-self.fracNotScared[day])*self.q_popsize)
        self.scared = (self.id<scaredcats) #strict otherwise picks id=0 when scared are 0

        return scaredcats
      
    def aDayInTheLife(self,day):
          
        # figure out who is scared and save their positions
        self.scaredycats(day)

        # move scared out of the box and far from each other
        self.savepos = np.copy(self.pos)
        self.pos[self.scared,1] = self.q_citysize*(2*(self.id[self.scared]+1))
        
        # generate infections, deaths and recoveries
        self.infections(day)
        self.symptoms(day)
        self.deaths(day)
        self.recoveries(day)      ,
                    
        # return people to their position before moving them arund 
        self.pos[self.scared,1] =  self.savepos[self.scared,1]

        # move people around
        self.pos = self.movepeople(self.pos)
        self.lastday = day
 
class spSAYDR_behav_local(spSAYDR_behav):

    def scaredycats(self,day) :
        
        '''
        Detects people scared to go outside so that infections() can 
        put them in a land far and far away

        '''
        par = self.behModel
        pos = self.pos
        fY = np.zeros(self.q_popsize)
        reducedContagion = np.ones(self.q_popsize)
        
        self.scared = np.array([False]*self.q_popsize)
        for i in np.arange(self.q_popsize) :
            if self.state[i] == self.S :

                #compute distance from all individuals
                dist_i = np.linalg.norm(pos-pos[i],axis=1)
                
                # figure out who is infected in the neighborhood
                mask = (dist_i<self.p_infradius) 
                neighbors = np.sum(mask)
                # infectedn = np.sum(mask[self.state==self.I])
                    
                fY[i] = np.sum(self.state[mask]==self.I)/np.sum(mask)
                if fY[i]>par['phi']:
                    reducedContagion = (par['phi']/fY[i])**0.12    
                else:
                    reducedContagion = 1
                
                # check who has id among the correct percentile defined
                # by the reducedContagion among the neighbors
                self.scared[i] = np.where(np.sort(self.id[mask])==self.id[i])/neighbors > reducedContagion
                pass
        
        self.fracNotScared[day] = 1 - np.average(self.scared)
        # code for no risk aversion
        # for i in np.arange(self.q_popsize) :
        #     if self.state[i] == self.S :

        #         #compute distance from all individuals
        #         dist_i = np.linalg.norm(pos-pos[i],axis=1)
                
        #         # figure out who is infected in the neighborhood
        #         mask = (dist_i<self.p_infradius) 
        #         # neighbors = np.sum(mask)
        #         # infectedn = np.sum(mask[self.state==self.I])
                    
        #         fY[i] = np.sum(self.state[mask]==self.I)/np.sum(mask)

        
        # # figure out how many are infected and compute how many are scared         
        # if par['type'] == 'Lones':
            
        #     reducedContagion[fY<=par['phi']] = 1.
        #     if np.sum(fY>par['phi'])>0:
        #         reducedContagion[fY>par['phi']] = (par['phi']/fY[fY>par['phi']])**0.12       
                
        # elif par['type'] == 'Jesus':
        #     exit('Jesus not implemented') 
            
                    
        # self.fracNotScared[day] = np.average(reducedContagion)
        # prob = (1-reducedContagion)
        
        # self.scared = (np.random.sample(self.q_popsize) < prob)
        return np.sum(self.scared)

    
# the following functions are needed only if using the multiprocessing package
# to produce multiple replications of the simulations at once
def simulatePool(kwargs) :
    ''' function that simulates the model once, to be used in conjunction
        with multiprocessing.Pool to get multiple replications quickly 
    '''
    m = spatialSAYDR(**kwargs)
    m.simulateOnce()
    m.summaryStats(m.lastday,nondef=kwargs,savefile='output/results.txt',clear='a')

    return m

def simulateRandPool(kwargs) :
    ''' function that simulates the model once, to be used in conjunction
        with multiprocessing.Pool to get multiple replications quickly 
    '''
    m = spSAYDR_randLoc(**kwargs)
    m.simulateOnce()
    m.summaryStats(m.lastday,nondef=kwargs,savefile='output/results.txt',clear='a')

    return m

def simulateHetPool(kwargs) :
    ''' function that simulates the model once, to be used in conjunction
        with multiprocessing.Pool to get multiple replications quickly 
    '''
    m = spSAYDR_hetDensity(**kwargs)
    m.simulateOnce()
    m.summaryStats(m.lastday,nondef=kwargs,savefile='output/results.txt',clear='a')

    return m

def simulateBehPool(kwargs) :
    ''' function that simulates the model once, to be used in conjunction
        with multiprocessing.Pool to get multiple replications quickly 
    '''
    m = spSAYDR_behav(**kwargs)
    m.simulateOnce()
    m.summaryStats(m.lastday,nondef=kwargs,savefile='output/results.txt',clear='a')

    return m

def simulateBehPool_local(kwargs) :
    ''' function that simulates the model once, to be used in conjunction
        with multiprocessing.Pool to get multiple replications quickly 
    '''
    m = spSAYDR_behav_local(**kwargs)
    m.simulateOnce()
    m.summaryStats(m.lastday,nondef=kwargs,savefile='output/results.txt',clear='a')

    return m

#%%#--------------------------------------------------------------------------#
if __name__ == "__main__":
    
    from multiprocessing import freeze_support
    freeze_support()

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
        "q_days"      : 134,
        "q_printOption" : 1,
        'g_maxinf'    : 1,
    }   

    kwargs = benchkwargs.copy()
    kwargs['behModel'] = {'type': 'Lones','phi': 0.01}

#    m = spatialSAYDR(**kwargs)
    m = spSAYDR_behav_local(**kwargs)
  
    m.computeStats(0,m.q_printOption)
    #m.drawPositions(0,savefig='baseline')
#   m.drawChart(0)
    for day in range(1,m.q_days):
    
        m.computeStats(day,1)
        if np.mod(day,5)==0 :
#            m.drawPositions(day,savefig='baseline',S=True,) #,folder='../Covid BisinMoro/write/images/')
            #m.drawChart(day)#,savefig='baseline',folder='../Covid BisinMoro/write/images/')
            pass
        m.drawRates(day)
        m.aDayInTheLife(day) 
          
        if np.sum(m.state==m.I)+np.sum(m.state==m.Y) <=15:
            break

    print(m.summaryStats(day,nondef=kwargs,savefile='output/results.txt',clear='a'))

       