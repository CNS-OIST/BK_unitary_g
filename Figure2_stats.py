import steps.interface
from steps.saving import *

import matplotlib.pyplot as plt
plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np


rat50=[]
rat250=[]

for SEED in range (1,101):
    SEED = str(SEED)

    BK_facs = CaP_facs = [1, 4]
    BK_facs=np.array(BK_facs)

    mesh = 'Cylinder3dia1umL10um'
    #mesh = 'Cylinder3dia5umL10um'
    capacfac = 50
    endtime = 3.0

    CaP_G = 0.419 # CaP approx single-channel conductance

    fig, (axis1, axis2) = plt.subplots(2, 3, figsize=[7.08, 5.2])
    ax1, ax2, ax3=axis1[0], axis1[1], axis1[2]
    ax4, ax5, ax6=axis2[0], axis2[1], axis2[2]

    BK250_maxes=[]
    BK50_maxes=[]
    CaP_maxes=[]


    BK_G = '250' # note, this choice basically doesn't matter for the CaP data. 50 and 250 look very similar, as you would expect.

    for BK_fac in BK_facs:
        dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_fac}_CaPfac1'
        with HDF5Handler('./STEPS/data/'+dataset) as hdf:
            Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
            actTot = BKstates.data[0,:,1] + BKstates.data[0,:,3] + BKstates.data[0,:,5] + BKstates.data[0,:,7] + BKstates.data[0,:,9]
            ax2.plot(1e3 * BKstates.time[0], actTot*250e-3, label = BK_fac )
            BK250_maxes.append(np.max(actTot*250))
            ax5.plot(1e3 * BKstates.time[0], 100*actTot/np.sum(BKstates.data[0,:,:], axis=1), label = BK_fac )

    BK_G = '50'
    for BK_fac in BK_facs:
        dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_fac}_CaPfac1'
        with HDF5Handler('./STEPS/data/'+dataset) as hdf:
            Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
            actTot = BKstates.data[0,:,1] + BKstates.data[0,:,3] + BKstates.data[0,:,5] + BKstates.data[0,:,7] + BKstates.data[0,:,9]
            ax3.plot(1e3 * BKstates.time[0], actTot*50e-3, label = BK_fac )
            BK50_maxes.append(np.max(actTot*50))
            ax6.plot(1e3 * BKstates.time[0], 100*actTot/np.sum(BKstates.data[0,:,:], axis=1), label = BK_fac )
            
            
    print("Seed:", SEED, "250 rat:", BK250_maxes[1]/BK250_maxes[0], "50 rat:", BK50_maxes[1]/BK50_maxes[0])
    
    rat50.append(BK50_maxes[1]/BK50_maxes[0])
    rat250.append(BK250_maxes[1]/BK250_maxes[0])

print ("Mean 50:", np.mean(rat50))
print ("Std 50:", np.std(rat50))
print ("Range 50:", np.min(rat50), np.max(rat50))

print ("Mean 250:", np.mean(rat250))
print ("Std 250:", np.std(rat250))
print ("Range 250:", np.min(rat250), np.max(rat250))

