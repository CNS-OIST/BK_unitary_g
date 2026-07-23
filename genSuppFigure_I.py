import steps.interface
from steps.saving import *

import matplotlib.pyplot as plt
plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np


SEED='123'


mesh = 'Cylinder3dia1umL10um'
#mesh = 'Cylinder3dia5umL10um'
capacfac = 50
endtime = 10.0

CaP_P = 2.5e-2 # CaP single-channel permeability, µm3/s

fig, (axis1, axis2) = plt.subplots(2, 2, figsize=[7.08, 5.2])
ax1, ax2 = axis1[0], axis1[1]
ax3, ax4 = axis2[0], axis2[1]



BK_G = '250'
BK_fac = 1

dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_fac}_CaPfac1'
with HDF5Handler('./STEPS/data/'+dataset) as hdf:
    Currents, BKstates, CaConcs, Pot = hdf['BKmodel_axononlySim'].results

    ax1.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,5], label = 'BK')
    ax1.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,0], label = 'Kv3')
    ax1.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,1], label = 'NaP')
    ax1.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,2], label = 'NaR')
    ax1.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,3], label = 'CaP')
    ax1.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,4], label = 'CaT')
    ax1.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,7], label = 'Ih')
    ax1.set_title(f"{BK_G}pS, factor={BK_fac}", fontsize=10)

BK_G = '50'
BK_fac = 1

dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_fac}_CaPfac1'
with HDF5Handler('./STEPS/data/'+dataset) as hdf:
    Currents, BKstates, CaConcs, Pot = hdf['BKmodel_axononlySim'].results

    ax2.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,5], label = 'BK')
    ax2.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,0], label = 'Kv3')
    ax2.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,1], label = 'NaP')
    ax2.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,2], label = 'NaR')
    ax2.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,3], label = 'CaP')
    ax2.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,4], label = 'CaT')
    ax2.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,7], label = 'Ih')
    ax2.set_title(f"{BK_G}pS, factor={BK_fac}", fontsize=10)

BK_G = '250'
BK_fac = 4

dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_fac}_CaPfac1'
with HDF5Handler('./STEPS/data/'+dataset) as hdf:
    Currents, BKstates, CaConcs, Pot = hdf['BKmodel_axononlySim'].results

    ax3.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,5], label = 'BK')
    ax3.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,0], label = 'Kv3')
    ax3.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,1], label = 'NaP')
    ax3.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,2], label = 'NaR')
    ax3.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,3], label = 'CaP')
    ax3.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,4], label = 'CaT')
    ax3.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,7], label = 'Ih')
    ax3.set_title(f"{BK_G}pS, factor={BK_fac}", fontsize=10)

BK_G = '50'
BK_fac = 4

dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_fac}_CaPfac1'
with HDF5Handler('./STEPS/data/'+dataset) as hdf:
    Currents, BKstates, CaConcs, Pot = hdf['BKmodel_axononlySim'].results

    ax4.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,5], label = 'BK')
    ax4.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,0], label = 'Kv3')
    ax4.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,1], label = 'NaP')
    ax4.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,2], label = 'NaR')
    ax4.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,3], label = 'CaP')
    ax4.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,4], label = 'CaT')
    ax4.plot(1e3 * Currents.time[0], 1e9 * Currents.data[0,:,7], label = 'Ih')
    ax4.set_title(f"{BK_G}pS, factor={BK_fac}", fontsize=10)


ax1.legend(loc='center right',facecolor='white', framealpha=1)
for ax in (ax1, ax3):
    ax.set_ylabel('Current (nA)')

for ax in (ax3, ax4):
    ax.set_xlabel('Time (ms)')

for ax in (ax1, ax2, ax3, ax4):
    ax.set_xlim(0.8, 2.6)
    ax.set_ylim(-6.5, 6.5)

#ax1.set_xticks((1.5, 2))
#ax1.set_ylim(-3,0.1)



ax1.text(-0.15,1.05,'a', fontweight = 'bold', transform=ax1.transAxes)
ax2.text(-0.15,1.05,'b', fontweight = 'bold', transform=ax2.transAxes)
ax3.text(-0.15,1.05,'c', fontweight = 'bold', transform=ax3.transAxes)
ax4.text(-0.15,1.05,'d', fontweight = 'bold', transform=ax4.transAxes)

plt.subplots_adjust(wspace=0.4, hspace=0.3)
#plt.savefig('Figures/Figure2.pdf', dpi=300)
plt.savefig(f'Figures_Sup/FigureS_I.pdf', dpi=300)
plt.show()
plt.close()
