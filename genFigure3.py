import steps.interface
from steps.saving import *

import matplotlib.pyplot as plt
plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import scipy

SEED='1234'

BK_gs = [50, 100, 150, 200, 250]

mesh = 'Cylinder3dia1umL10um'

endtime = 100.0

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

fig, (axis1, axis2) = plt.subplots(2, 3, figsize=[7.08, 5.2])
ax, ax2, ax3=axis1[0], axis1[1], axis1[2]
axins = inset_axes(ax, width=0.5, height=0.7, bbox_to_anchor=(0.2, .4, .3, .5), bbox_transform=ax.transAxes)

ax4, ax5, ax6=axis2[0], axis2[1], axis2[2]
ax7 = ax6.twinx()

BK_maxes=[]
CaP_maxes=[]
peaks_v = []
peaks_ca = []
widths_ca = []

for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac1_CaPfac1'
    with HDF5Handler('./STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        actTot = BKstates.data[0,:,1] + BKstates.data[0,:,3] + BKstates.data[0,:,5] + BKstates.data[0,:,7] + BKstates.data[0,:,9]
        ax2.plot(1e3 * (BKstates.time[0]-np.argmax(actTot)*1e-5), actTot*BK_g*1e-3, label = f'{BK_g}pS')
        BK_maxes.append(np.max(actTot*BK_g*1e-3))
        ax.plot(1e3 * Pot.time[0], 1e3 * Pot.data[0])
        axins.plot(1e3 * Pot.time[0], 1e3 * Pot.data[0] )
        ps = scipy.signal.find_peaks(1e3 * Pot.data[0,:,0], prominence=15)[0]
        peaks_v.append(np.take(1e3 * Pot.time[0], ps))

        BKdat = BKstates.data
        CaTot = 4*BKdat[0,:,0] + 4*BKdat[0,:,1] + 3*BKdat[0,:,2] + 3*BKdat[0,:,3] + 2*BKdat[0,:,4] + 2*BKdat[0,:,5] + BKdat[0,:,6] + BKdat[0,:,7]
        ax4.plot(1e3 * (BKstates.time[0]-np.argmax(CaTot)*1e-5), CaTot,label = f'{BK_g}pS')
        ax5.plot(1e3 * (CaConcs.time[0]-np.argmax(CaConcs.data[0,:,1])*1e-5), 1e6 * CaConcs.data[0,:,1])
        
        ps = scipy.signal.find_peaks(1e6 * CaConcs.data[0,:,1], prominence=10)[0]
        ws = 1e-2*np.array(scipy.signal.peak_widths(1e6 * CaConcs.data[0,:,1], ps, rel_height=0.9))[0]
        peaks_ca.append(np.mean(np.take(1e6 * CaConcs.data[0,:,1], ps)))
        widths_ca.append(np.mean(ws))
        

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Membrane potential (mV)')
axins.set_xlim(1,2)
axins.set_xticks((1,2))
axins.set_yticks(())

ax2.legend(loc='best')
ax2.set_xlabel('Time rel to peak (ms)')
ax2.set_ylabel('BK conductance (nS)')
ax2.set_xlim(-0.5,0.8)

ax3.plot(BK_gs, BK_maxes, 'o-')
ax3.set_xlabel('Unitary BK conductance (pS)')
ax3.set_ylabel('Max BK conductance (nS)')

ax4.set_xlabel('Time rel to peak (ms)')
ax4.set_ylabel('Ca ions bound to BK channels')
ax4.set_xlim(1.1-1.6,2.9-1.6)
ax4.legend(loc='best')
ax4.set_ylim(0)

ax5.set_xlabel('Time rel to peak (ms)')
ax5.set_ylabel('Submemb Ca concentration (µM)')
ax5.set_xlim(-0.2, 0.7)
ax5.set_ylim(0)

ax6.plot(BK_gs, peaks_ca, 'bo-')
ax6.set_ylabel('Submemb Ca peak heights (µM)')
ax6.set_xlabel('Unitary BK conductance (pS)')
ax6.yaxis.label.set_color('blue')
ax6.tick_params(axis='y', colors='blue')

ax7.plot(BK_gs, widths_ca, 'ro-')
ax7.set_ylabel('Submemb Ca peak width (ms)')
ax7.yaxis.label.set_color('red')
ax7.tick_params(axis='y', colors='red')

if mesh == 'Cylinder3dia1umL10um':
    ax.set_xlim(85,95)
    ax.set_xticks((86,88,90,92,94))
    ax6.set_ylim(44, 61)
    ax6.set_yticks((45,50,55,60))
    ax7.set_ylim(0.6, 0.8)
    ax7.set_yticks((0.6, 0.7, 0.8))
else:
    ax.set_xlim(79,87)
    ax6.set_ylim(40, 60)
    ax6.set_yticks((40, 50, 60))
    ax7.set_yticks((0.5, 0.6, 0.7))

text_x = -0.15
text_x = -0.4
text_y = 1.06
text_y = 1.1
ax.text(text_x,text_y,'a', fontweight = 'bold', transform=ax.transAxes)
ax2.text(text_x,text_y,'b', fontweight = 'bold', transform=ax2.transAxes)
ax3.text(text_x,text_y,'c', fontweight = 'bold', transform=ax3.transAxes)
ax4.text(text_x,text_y,'d', fontweight = 'bold', transform=ax4.transAxes)
ax5.text(text_x,text_y,'e', fontweight = 'bold', transform=ax5.transAxes)
ax6.text(text_x,text_y,'f', fontweight = 'bold', transform=ax6.transAxes)


plt.subplots_adjust(wspace=0.6, hspace=0.4)
plt.savefig('Figures/Figure3A-F.pdf', dpi=300)
plt.close()

