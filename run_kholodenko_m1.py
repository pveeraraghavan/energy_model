#==============================================================================
# RUN THE MODEL AND PLOT
#==============================================================================

from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
import kholodenko_m1

#define run time
tspan = np.linspace(0, 10, 1000)
#create solver
sim = ScipyOdeSimulator(kholodenko_m1.model, tspan)

#==============================================================================
# SIMPLE SIMULATION WITH 4 SETUP FOR F AND G
#==============================================================================
#==============================================================================
# 
# #define f,g paramter set to simulate
# ff=[1.0, 0.01, 1.0, 0.01]
# gg=[1.0, 1.0, 100, 100];
# 
# #define I inhibitor concentration dilution
# I_dil = np.logspace(-2, 4, 20)
# 
# #run simulations with different f,g paramters over an inihibitor I dilution to
# #generate steady state plots of B, BB, BI, BBI and IBBI concentrations
# for i in range(len(ff)):
#     #build steady state response vector
#     ss_v = np.empty((len(I_dil), len(kholodenko_m1.model.observables)))
#     for j in range(len(I_dil)):
#         res = sim.run(param_values={'conc_I_0':I_dil[j], 'f':ff[i], 'g':gg[i]})
#         #res.dataframe.loc[:, 'B_obs':'IBBI_obs'].plot()
#         #plt.show()
#         ss_v[j] = list(res.observables[-1])
#         ss_v
#     plt.figure()    
#     plt.plot(I_dil, ss_v[:,2:8], marker='o')
#     plt.xscale('log')
#     str_leg=[d[0] for d in res.observables.dtype.descr]
#     plt.legend(str_leg[2:8])
#==============================================================================

#==============================================================================
# COMPLEX SIMULATION WITH SCAN FOR F AND G OVER MULTIPLE VALUES
#==============================================================================

#set param scan
ff=np.logspace(-3, 3, 7, base=10.0)
gg=np.logspace(-3, 3, 7, base=10.0)

#define I inhibitor concentration dilution
I_dil = np.logspace(-2, 4, 20)

#set dimension multiple subplots
nrow=len(ff)
ncol=len(gg)
#create figure with subplot
fig, axarray = plt.subplots(nrow,ncol)
#font size
fs=6
#ind_obs=range(2,9)
#ind_obs=range(7,9)
#perform simulations
for i in range(len(ff)):
    for j in range(len(gg)):
        #build steady state response vector
        ss_v = np.empty((len(I_dil), len(kholodenko_m1.model.observables)))
        for z in range(len(I_dil)):
            res = sim.run(param_values={'conc_I_0':I_dil[z], 'f':ff[i], 'g':gg[j]})
            ss_v[z] = list(res.observables[-1])
        ax=axarray[i,j]
        ax.plot(I_dil, ss_v[:,2:9], marker='o', markersize=3)
        ax.set_xscale('log')
        ax.tick_params(labelsize=fs)
        ax.set_title("f:%.1e,g:%.1e"%(ff[i],gg[j]), fontsize=fs)
        str_leg=[d[0] for d in res.observables.dtype.descr]
        if (i==1) and (j==1):
           ax.legend(str_leg[2:9], fontsize=fs)
    print "Simulation %d of %d done.." % (i+1,len(ff))        