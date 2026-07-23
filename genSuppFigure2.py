import steps.interface
from steps.saving import *

import matplotlib.pyplot as plt
plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np


SEED='123'

BK_facs = CaP_facs = [1, 2, 3, 4]
BK_facs=np.array(BK_facs)

mesh = 'Cylinder3dia1umL10um'

capacfac = 50
endtime = 10.0

CaP_G = 0.419 # CaP approx single-channel conductance

fig, (axis1, axis2) = plt.subplots(2, 2, figsize=[5.2, 5.2])
ax1, ax2 = axis1[0], axis1[1]
ax3, ax4 =axis2[0], axis2[1]

CaP_BK250_maxes=[]
CaP_BK50_maxes=[]


BK_G = '250' # note, this choice basically doesn't matter for the CaP data. 50 and 250 look very similar, as you would expect.
for CaP_fac in CaP_facs:
    dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_facs[0]}_CaPfac{CaP_fac}'
    with HDF5Handler('./STEPS/data_CaPstates/'+dataset) as hdf:
        Currents, CaPstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        ax1.plot(1e3 * CaPstates.time[0], CaPstates.data[0,:,1]*CaP_G*1e-3, label = CaP_fac )
        CaP_BK250_maxes.append(np.max(CaPstates.data[0,:,1]*CaP_G))
        ax3.plot(1e3 * CaPstates.time[0], 100*CaPstates.data[0,:,1]/np.sum(CaPstates.data[0,:,:], axis=1), label = CaP_fac )

BK_G = '50'

for CaP_fac in CaP_facs:
    dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_facs[0]}_CaPfac{CaP_fac}'
    with HDF5Handler('./STEPS/data_CaPstates/'+dataset) as hdf:
        Currents, CaPstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        ax2.plot(1e3 * CaPstates.time[0], CaPstates.data[0,:,1]*CaP_G*1e-3, label = CaP_fac )
        CaP_BK50_maxes.append(np.max(CaPstates.data[0,:,1]*CaP_G))
        ax4.plot(1e3 * CaPstates.time[0], 100*CaPstates.data[0,:,1]/np.sum(CaPstates.data[0,:,:], axis=1), label = CaP_fac )


axins1 = inset_axes(ax1, width=0.7, height=0.7, bbox_to_anchor=(0.65, .1, .3, .5), bbox_transform=ax1.transAxes)
axins1.plot(range(1,6), range(1,6), 'k--')
for i in range(4): axins1.plot(CaP_facs[i], CaP_BK250_maxes[i]/CaP_BK250_maxes[0], 'o')

axins2 = inset_axes(ax2, width=0.7, height=0.7, bbox_to_anchor=(0.65, .1, .3, .5), bbox_transform=ax2.transAxes)
axins2.plot(range(1,6), range(1,6), 'k--')
for i in range(4): axins2.plot(CaP_facs[i], CaP_BK50_maxes[i]/CaP_BK50_maxes[0], 'o')


for ax in (ax1, ax2):
    ax.legend(CaP_facs, loc='upper right',facecolor='white', framealpha=1)
    ax.set_ylabel('CaP conductance (nS)')
    ax.set_xlim(1.1,2.5)
    ax.set_xticks((1.5, 2))


for ax in (ax3, ax4):
    ax.legend(loc='best')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('% of activated CaP channels')
    ax.set_xlim(1.1,2.3)
    ax.set_xticks((1.5, 2))


for axins in (axins1, axins2):
    axins.set_ylabel('Peak cond')
    axins.set_xlabel('Increase factor')
    axins.set_yticks(())
    axins.set_xticks(())
    axins.set_ylim(0.8, 4.3)
    axins.set_xlim(0.7, 4.3)
axins1.set_ylim(0.7, 4.3)
axins1.set_ylabel('Peak perm')

ax1.text(-0.15,1.05,'a', fontweight = 'bold', transform=ax1.transAxes)
ax1.set_title('250pS BK', fontsize=12)
ax2.text(-0.15,1.05,'b', fontweight = 'bold', transform=ax2.transAxes)
ax2.set_title('50pS BK', fontsize=12)
ax3.text(-0.15,1.05,'c', fontweight = 'bold', transform=ax3.transAxes)
ax4.text(-0.15,1.05,'d', fontweight = 'bold', transform=ax4.transAxes)

plt.subplots_adjust(wspace=0.4, hspace=0.3)
plt.savefig(f'Figures_Sup/FigureS2.pdf', dpi=300)
plt.close()
