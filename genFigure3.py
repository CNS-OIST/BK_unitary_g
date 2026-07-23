import steps.interface
from steps.saving import *

import matplotlib.pyplot as plt
plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import scipy

SEED='1234'
NA = 6.02214076e23 # Avogadro's number

BK_gs = [50, 100, 150, 200, 250]
BK_gs_red = [50, 150, 250]
colours=['blue', 'orange', 'green']

mesh = 'Cylinder3dia1umL10um'
submemb_tets_vol = 1.8853612004593083e-15 #L
cyto_tets_vol = 7.573408483497916e-15 #L

endtime = 5.0

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

fig, (axis1, axis2) = plt.subplots(2, 4, figsize=[9.08, 5.2])
ax, ax2, ax3, ax4=axis1[0], axis1[1], axis1[2], axis1[3]
axins = inset_axes(ax, width=0.4, height=0.5, bbox_to_anchor=(0.7, .45, .3, .5), bbox_transform=ax.transAxes)

ax7, ax6, ax5, ax8 =axis2[0], axis2[1], axis2[2],  axis2[3]
axins4 = inset_axes(ax4, width=0.3, height=0.4, bbox_to_anchor=(0.4, -0.13, .3, .5), bbox_transform=ax4.transAxes)
ax9 = ax8.twinx()

BK_maxes=[]
CaP_maxes=[]
peaks_v = []
peaks_ca = []
widths_ca = []

relT=False


idx=-1
for BK_g in BK_gs:
    dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_SKfac0.0_TEMP34.0_capacfac50_BK{BK_g}p_BKfac1_CaPfac1'
    with HDF5Handler('./STEPS/data/'+dataset) as hdf:
        Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
        actTot = BKstates.data[0,:,1] + BKstates.data[0,:,3] + BKstates.data[0,:,5] + BKstates.data[0,:,7] + BKstates.data[0,:,9]
        
        if BK_g in BK_gs_red:
            if relT: ax2.plot(1e3 * (BKstates.time[0]-np.argmax(actTot)*1e-5), actTot*BK_g*1e-3, label = f'{BK_g}pS')
            else: ax2.plot(1e3 * BKstates.time[0], actTot*BK_g*1e-3, label = f'{BK_g}pS')
        
            ax.plot(1e3 * Pot.time[0], 1e3 * Pot.data[0])
            axins.plot(1e3 * (Pot.time[0]-np.argmax(Pot.data[0])*1e-5), 1e3 * Pot.data[0] )

        BK_maxes.append(np.max(actTot*BK_g*1e-3))
        ps = scipy.signal.find_peaks(1e3 * Pot.data[0,:,0], prominence=15)[0]
        peaks_v.append(np.take(1e3 * Pot.time[0], ps))

        BKdat = BKstates.data
        CaTot = 4*BKdat[0,:,0] + 4*BKdat[0,:,1] + 3*BKdat[0,:,2] + 3*BKdat[0,:,3] + 2*BKdat[0,:,4] + 2*BKdat[0,:,5] + BKdat[0,:,6] + BKdat[0,:,7]
        AllTot = BKdat[0,:,0] + BKdat[0,:,1] + BKdat[0,:,2] + BKdat[0,:,3] + BKdat[0,:,4] + BKdat[0,:,5] + BKdat[0,:,6] + BKdat[0,:,7] + BKdat[0,:,8] + BKdat[0,:,9]
        AllTot = np.sum(BKdat[0,:,:], axis=1)
        
        if BK_g in BK_gs_red:
            idx+=1
            if relT: ax4.plot(1e3 * (BKstates.time[0]-np.argmax(CaTot)*1e-5), CaTot/AllTot,label = f'{BK_g}pS', color=colours[idx])
            else: ax4.plot(1e3 * BKstates.time[0], CaTot/AllTot,label = f'{BK_g}pS', color=colours[idx])
            axins4.plot(1e3 * (BKstates.time[0]-np.argmax(CaTot)*1e-5), CaTot/AllTot,label = f'{BK_g}pS')


            """  # concentration 
            if relT: ax5.plot(1e3 * (CaConcs.time[0]-np.argmax(CaConcs.data[0,:,1])*1e-5), 1e6 * CaConcs.data[0,:,1], label = f'{BK_g}pS subm', color=colours[idx])
            else: ax5.plot(1e3 * CaConcs.time[0], 1e6 * CaConcs.data[0,:,1], label = f'{BK_g}pS subm', color=colours[idx])
            ax5.plot(1e3 * BKstates.time[0], 1e6*CaTot/(submemb_tets_vol*NA), '--', label = f'{BK_g}pS chan', color=colours[idx]) # 'concentration' bound to channel
            ax5.plot(1e3 * CaConcs.time[0], 1e6 * CaConcs.data[0,:,0], ':', label = f'{BK_g}pS comp', color=colours[idx])
            """
            if BK_g == 50: axx = ax5
            elif BK_g == 150: axx = ax6
            else: axx = ax7
            axx.plot(1e3 * BKstates.time[0], 0.001*CaTot, '-', label = f'chan', color=colours[idx]) # 'concentration' bound to channel
            axx.plot(1e3 * CaConcs.time[0], 0.001*submemb_tets_vol*NA*CaConcs.data[0,:,1], '--', label = f'subm', color=colours[idx])
            axx.plot(1e3 * CaConcs.time[0], 0.001*((cyto_tets_vol*NA * CaConcs.data[0,:,0])-(submemb_tets_vol*NA*CaConcs.data[0,:,1])), ':', label = f'core', color=colours[idx])
            axx.set_title(f"{BK_g}pS", fontsize=10)
        
        ps = scipy.signal.find_peaks(1e6 * CaConcs.data[0,:,1], prominence=10)[0]
        ws = 1e-2*np.array(scipy.signal.peak_widths(1e6 * CaConcs.data[0,:,1], ps, rel_height=0.9))[0]
        peaks_ca.append(np.mean(np.take(1e6 * CaConcs.data[0,:,1], ps)))
        widths_ca.append(np.mean(ws))
        

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Membrane potential (mV)')
ax.set_xlim(1,2)

axins.set_xlim(-0.05,0.05)
axins.set_ylim(16,20)
axins.set_xticks(())

ax2.legend(loc='best')
if relT:
    ax2.set_xlabel('Time rel to peak (ms)')
    ax2.set_xlim(-0.5,0.8)
else:
    ax2.set_xlabel('Time (ms)')
    ax2.set_xlim(1.2,3.0)
ax2.set_ylabel('BK conductance (nS)')

ax3.plot(BK_gs, BK_maxes, 'o-')
ax3.set_xlabel('Unitary BK conductance (pS)')
ax3.set_ylabel('Max BK conductance (nS)')
ax3.set_xticks((50,150,250))

if relT:
    ax4.set_xlabel('Time rel to peak (ms)')
    ax4.set_xlim(-0.2, 0.3)
else:
    ax4.set_xlabel('Time (ms)')
    ax4.set_xlim(1.2, 2.3)
ax4.set_ylabel('Ca ions bound per BK channel')

axins4.set_xlim(-0.1,0.15)
axins4.set_ylim(3.05,3.4)
axins4.set_xticks(())

for axx in [ax5, ax6, ax7]:
    if relT:
        axx.set_xlabel('Time rel to peak (ms)')
        axx.set_xlim(-0.2, 0.7)
        axx.set_ylim(0, 130)
    else:
        axx.set_xlabel('Time (ms)')
        axx.set_xlim(1.2,2.5)
        axx.set_ylim(0, 130)

ax7.set_ylabel('Ca ion number (/1000)')
ax6.legend(loc='best')

ax8.plot(BK_gs, peaks_ca, 'bo-')
ax8.set_ylabel('Submemb Ca peak heights (µM)')
ax8.set_xlabel('Unitary BK conductance (pS)')
ax8.yaxis.label.set_color('blue')
ax8.tick_params(axis='y', colors='blue')
ax8.set_xticks((50,150,250))

ax9.plot(BK_gs, widths_ca, 'ro-')
ax9.set_ylabel('Submemb Ca peak width (ms)')
ax9.yaxis.label.set_color('red')
ax9.tick_params(axis='y', colors='red')

ax8.set_ylim(42, 60)
ax8.set_yticks((45,50,55,60))
ax9.set_ylim(0.58, 0.78)
ax9.set_yticks((0.6, 0.7, 0.8))

text_x = -0.15
text_x = -0.4
text_y = 1.06
text_y = 1.1
ax.text(text_x,text_y,'a', fontweight = 'bold', transform=ax.transAxes)
ax2.text(text_x,text_y,'b', fontweight = 'bold', transform=ax2.transAxes)
ax3.text(text_x,text_y,'c', fontweight = 'bold', transform=ax3.transAxes)
ax4.text(text_x,text_y,'d', fontweight = 'bold', transform=ax4.transAxes)

ax7.text(text_x,text_y,'e', fontweight = 'bold', transform=ax7.transAxes)

ax8.text(text_x,text_y,'f', fontweight = 'bold', transform=ax8.transAxes)


plt.subplots_adjust(wspace=0.6, hspace=0.5)

plt.savefig(f'Figures/Figure3A-F.pdf', dpi=300)
plt.close()

