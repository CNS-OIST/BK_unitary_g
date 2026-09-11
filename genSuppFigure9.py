import steps.interface

from steps.saving import *

from matplotlib import pyplot as plt
import numpy as np

import scipy
from math import pi

SEED='222' # is lin space
caconces = np.arange(0, 101e-6, 1e-6)

caconces = np.flip(caconces)

lw=3
BK_maxes=[]
CaP_maxes=[]

BK_gs = [50, 150, 250]

"""
indexes, i, in Currents.data[0,:,i]
0 = Kv3
1 = NaP
2 = Rsg
3 = CaP
4 = CaT
5 = BK
6 = SK
7 = Ih
8 = L
label = Currents.labels[5].split('.')[-2], the name of the current
"""

"""
BKstates index and state
0 memb.BKchan[BKCa, BKCa, BKCa, BKCa, BKclose].Count
1 memb.BKchan[BKCa, BKCa, BKCa, BKCa, BKopen].Count
2 memb.BKchan[BK, BKCa, BKCa, BKCa, BKclose].Count
3 memb.BKchan[BK, BKCa, BKCa, BKCa, BKopen].Count
4 memb.BKchan[BK, BK, BKCa, BKCa, BKclose].Count
5 memb.BKchan[BK, BK, BKCa, BKCa, BKopen].Count
6 memb.BKchan[BK, BK, BK, BKCa, BKclose].Count
7 memb.BKchan[BK, BK, BK, BKCa, BKopen].Count
8 memb.BKchan[BK, BK, BK, BK, BKclose].Count
9 memb.BKchan[BK, BK, BK, BK, BKopen].Count

BKstates index and bound Ca
0 Ca, Ca, Ca, Ca
1 Ca, Ca, Ca, Ca
2 Ca, Ca, Ca
3 Ca, Ca, Ca
4 Ca, Ca
5 Ca, Ca
6 Ca
7 Ca
8
9
"""

counter=1

for BK_g in BK_gs:
    plt.figure(1)
    plt.subplot(1,len(BK_gs)+1,counter)
    counter+=1
    acts=[]
    plotcnt=10
    for caconc in caconces:
        dataset = f'BKonlymodel_'+SEED+f'_10.0ms_Cabind_Cylinder3dia1umL10um_TEMP34.0_capacfac50_BK{BK_g}p_Caconc{caconc}_pot0.0'
        with HDF5Handler('STEPS/data/'+dataset) as hdf:
            Currents, BKstates, CaConcs, Pot =  hdf['BKonlymodel'].results
        
            BKdat = BKstates.data
            countTot = np.sum(BKdat[0,:,:], axis=1)
            actTot = BKdat[0,:,1] + BKdat[0,:,3] + BKdat[0,:,5] + BKdat[0,:,7] + BKdat[0,:,9]
            inactTot = BKdat[0,:,0] + BKdat[0,:,2] + BKdat[0,:,4] + BKdat[0,:,6] + BKdat[0,:,8]
            if plotcnt==10:
                if counter ==4:
                    plt.plot(1e3 * BKstates.time[0], 100 * actTot/countTot, linewidth=lw, label = f'{int(round(caconc*1e6))}µM')
                else:
                    plt.plot(1e3 * BKstates.time[0], 100 * actTot/countTot, linewidth=lw)
                plotcnt=1
            else: plotcnt+=1

            acts.append(100*np.mean(actTot[-200:-1])/countTot[-1]) # averaging the last 2ms

    plt.xlabel('Time (ms)')
    plt.ylabel(f'Percentage of activated BK channels, {BK_g}pS')
    plt.ylim(0, 100)
    if counter ==4: plt.legend()

    plt.figure(2)
    plt.plot(caconces*1e6, acts, label = f'{BK_g}pS', linewidth=lw )


plt.figure(1)
plt.subplot(1,len(BK_gs)+1,counter)
acts=[]

plotcnt=10
for caconc in caconces:
    # seed 333 has no calcium binding
    dataset = f'BKonlymodel_333_10.0ms_Canobind_Cylinder3dia1umL10um_TEMP34.0_capacfac50_BK250p_Caconc{caconc}_pot0.0'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKonlymodel'].results
    
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        actTot = BKdat[0,:,1] + BKdat[0,:,3] + BKdat[0,:,5] + BKdat[0,:,7] + BKdat[0,:,9]
        inactTot = BKdat[0,:,0] + BKdat[0,:,2] + BKdat[0,:,4] + BKdat[0,:,6] + BKdat[0,:,8]

        if plotcnt==10:
            plt.plot(1e3 * BKstates.time[0], 100 * actTot/countTot, linewidth=lw)
            plotcnt=1
        else: plotcnt+=1

        acts.append(100*np.mean(actTot[-200:-1])/countTot[-1]) # averaging the last 2ms

plt.xlabel('Time (ms)')
plt.ylabel(f'Percentage of activated BK channels, no calcium binding')
plt.ylim(0, 100)

plt.figure(2)
plt.plot(caconces*1e6, acts, label = 'No binding', linewidth=lw )




plt.figure(1)
fig = plt.gcf()
fig.set_size_inches(9*1.5, 4.5*1.5)
plt.subplots_adjust(wspace=0.4, hspace=0.3)
plt.savefig('Figures_Sup/FigureS9a.pdf', dpi=300)


plt.figure(2)
plt.xlabel('Calcium concentration (uM)')
plt.ylabel('Percentage of activated BK channels')
fig = plt.gcf()
fig.set_size_inches(9*1.5, 4.5*1.5)
plt.subplots_adjust(wspace=0.4, hspace=0.3)
plt.legend()
plt.savefig('Figures_Sup/FigureS9b.pdf', dpi=300)

plt.close()




