import steps.interface

from steps.saving import *

from matplotlib import pyplot as plt
#plt.style.use('./figures.naturestyle')
import numpy as np

SEED='123'

BK_gs = [50,  250]

lw=3
BK_maxes=[]
CaP_maxes=[]

mesh = 'Cylinder3dia1umL10um'

mode = 'perc_all' # if not this it'll be percentage of active

BKfac=4

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
plt.subplot(251)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        inactTot = BKdat[0,:,0] + BKdat[0,:,2] + BKdat[0,:,4] + BKdat[0,:,6] + BKdat[0,:,8]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,8]/countTot, linewidth=lw  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,8]/inactTot, linewidth=lw  )
plt.xlabel('Time (ms)')
if mode == 'perc_all':plt.ylabel('Percentage of all BK channels in state')
else: plt.ylabel('Percentage of inactive BK channels in state')
plt.xlim(1,4)
plt.title('C0')
plt.ylim(0,100)

plt.subplot(252)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        inactTot = BKdat[0,:,0] + BKdat[0,:,2] + BKdat[0,:,4] + BKdat[0,:,6] + BKdat[0,:,8]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,6]/countTot, linewidth=lw  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,6]/inactTot, linewidth=lw  )
plt.xlabel('Time (ms)')
plt.xlim(1,4)
plt.title('C1')
plt.ylim(0,100)

plt.subplot(253)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        inactTot = BKdat[0,:,0] + BKdat[0,:,2] + BKdat[0,:,4] + BKdat[0,:,6] + BKdat[0,:,8]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,4]/countTot, linewidth=lw , label=f"{BK_g}pS"  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,4]/inactTot, linewidth=lw , label=f"{BK_g}pS"  )
plt.xlabel('Time (ms)')
plt.xlim(1,4)
plt.title('C2')
plt.ylim(0,100)
plt.legend()

plt.subplot(254)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        inactTot = BKdat[0,:,0] + BKdat[0,:,2] + BKdat[0,:,4] + BKdat[0,:,6] + BKdat[0,:,8]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,2]/countTot, linewidth=lw  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,2]/inactTot, linewidth=lw  )
plt.xlabel('Time (ms)')
plt.xlim(1,4)
plt.title('C3')
plt.ylim(0,100)

plt.subplot(255)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        inactTot = BKdat[0,:,0] + BKdat[0,:,2] + BKdat[0,:,4] + BKdat[0,:,6] + BKdat[0,:,8]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,0]/countTot, linewidth=lw )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,0]/inactTot, linewidth=lw  )
plt.xlabel('Time (ms)')
plt.xlim(1,4)
plt.title('C4')
plt.ylim(0,100)

plt.subplot(256)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        actTot = BKdat[0,:,1] + BKdat[0,:,3] + BKdat[0,:,5] + BKdat[0,:,7] + BKdat[0,:,9]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,9]/countTot, linewidth=lw  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,9]/actTot, linewidth=lw  )
plt.xlabel('Time (ms)')
if mode == 'perc_all':plt.ylabel('Percentage of all BK channels in state')
else: plt.ylabel('Percentage of active BK channels in state')
plt.xlim(1,4)
plt.title('O0')
plt.ylim(0,4)


plt.subplot(257)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        actTot = BKdat[0,:,1] + BKdat[0,:,3] + BKdat[0,:,5] + BKdat[0,:,7] + BKdat[0,:,9]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,7]/countTot, linewidth=lw  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,7]/actTot, linewidth=lw  )
plt.xlabel('Time (ms)')
plt.xlim(1,4)
plt.title('O1')
plt.ylim(0,4)

plt.subplot(258)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        actTot = BKdat[0,:,1] + BKdat[0,:,3] + BKdat[0,:,5] + BKdat[0,:,7] + BKdat[0,:,9]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,5]/countTot, linewidth=lw  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,5]/actTot, linewidth=lw  )
plt.xlabel('Time (ms)')
plt.xlim(1,4)
plt.title('O2')
plt.ylim(0,4)

plt.subplot(259)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        actTot = BKdat[0,:,1] + BKdat[0,:,3] + BKdat[0,:,5] + BKdat[0,:,7] + BKdat[0,:,9]
        if mode == 'perc_all': plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,3]/countTot, linewidth=lw  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,3]/actTot, linewidth=lw  )
plt.xlabel('Time (ms)')
plt.xlim(1,4)
plt.title('O3')
plt.ylim(0,4)

plt.subplot(2,5,10)
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_10.0ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac{BKfac}_CaPfac1'
    with HDF5Handler('STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        BKdat = BKstates.data
        countTot = np.sum(BKdat[0,:,:], axis=1)
        actTot = BKdat[0,:,1] + BKdat[0,:,3] + BKdat[0,:,5] + BKdat[0,:,7] + BKdat[0,:,9]
        if mode == 'perc_all':plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,1]/countTot, linewidth=lw  )
        else: plt.plot(1e3 * BKstates.time[0], 100*BKdat[0,:,1]/actTot, linewidth=lw  )
plt.xlabel('Time (ms)')
plt.xlim(1,4)
plt.title('O4')
plt.ylim(0,4)


fig = plt.gcf()
fig.set_size_inches(1.3*6.8, 1.3*5.5)
plt.subplots_adjust(wspace=0.42, hspace=0.35)
plt.savefig(f'Figures_Sup/FigureS4.pdf', dpi=300)
plt.close()
