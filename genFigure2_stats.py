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

"""

axins1 = inset_axes(ax1, width=0.7, height=0.7, bbox_to_anchor=(0.65, .1, .3, .5), bbox_transform=ax1.transAxes)
axins1.plot(range(1,6), range(1,6), 'k--')

axins2 = inset_axes(ax2, width=0.7, height=0.7, bbox_to_anchor=(0.67, .5, .3, .5), bbox_transform=ax2.transAxes)
axins2.plot(range(1,6), range(1,6), 'k--')
for i in range(len(CaP_facs)): axins2.plot(BK_facs[i]/BK_facs[0], BK250_maxes[i]/BK250_maxes[0], 'o')

axins3 = inset_axes(ax3, width=0.7, height=0.7, bbox_to_anchor=(0.65, .5, .3, .5), bbox_transform=ax3.transAxes)
axins3.plot(range(1,6), range(1,6), 'k--')
for i in range(len(CaP_facs)): axins3.plot(BK_facs[i]/BK_facs[0], BK50_maxes[i]/BK50_maxes[0], 'o')


ax1.legend(CaP_facs, loc='upper right',facecolor='white', framealpha=1)
ax1.set_ylabel('CaP conductance (nS)')
ax1.set_xlim(1.1,2.5)
ax1.set_xticks((1.5, 2))

for ax in (ax2, ax3):
    ax.legend(BK_facs, loc='lower right',facecolor='white', framealpha=1)
    ax.set_xlim(1.2,4.2)
    ax.set_ylim(0,350)
    #ax.set_ylim(0,1200)
ax2.set_ylabel(f'BK250 conductance (nS)')
ax3.set_ylabel(f'BK50 conductance (nS)')


ax4.legend(loc='best')
ax4.set_xlabel('Time (ms)')
ax4.set_ylabel('% of activated CaP channels')
ax4.set_xlim(1.1,2.3)
#ax4.set_ylim(0,105)
ax4.set_xticks((1.5, 2))

for ax in (ax5, ax6):
    ax.legend(loc='best')
    ax.set_xlabel('Time (ms)')
    ax.set_xlim(1.2,4.2)
    ax.set_ylim(0,9)

ax5.set_ylabel(f'% of activated BK250 channels')
ax6.set_ylabel(f'% of activated BK50 channels')

for axins in (axins1, axins2, axins3):
    axins.set_ylabel('Peak cond')
    axins.set_xlabel('Increase factor')
    axins.set_yticks(())
    axins.set_xticks(())
    axins.set_ylim(0.8, 4.3)
    axins.set_xlim(0.7, 4.3)
axins1.set_ylim(0.7, 4.3)

ax1.text(-0.15,1.05,'a', fontweight = 'bold', transform=ax1.transAxes)
ax2.text(-0.15,1.05,'b', fontweight = 'bold', transform=ax2.transAxes)
ax3.text(-0.15,1.05,'c', fontweight = 'bold', transform=ax3.transAxes)
ax4.text(-0.15,1.05,'d', fontweight = 'bold', transform=ax4.transAxes)
ax5.text(-0.15,1.05,'e', fontweight = 'bold', transform=ax5.transAxes)
ax6.text(-0.15,1.05,'f', fontweight = 'bold', transform=ax6.transAxes)

plt.subplots_adjust(wspace=0.4, hspace=0.3)
#plt.savefig('Figures/Figure2.pdf', dpi=300)
plt.show()
plt.close()
"""
