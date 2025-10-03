import steps.interface
from steps.saving import *

import matplotlib.pyplot as plt
plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import scipy

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Create a figure
fig = plt.figure(figsize=(7.08, 10))

# Define grid specifications for the layout
gs = gridspec.GridSpec(5, 15, figure=fig)

# Create subplots based on the grid specifications
# Row 1 (5 columns)
ax1 = fig.add_subplot(gs[0, :3])
ax2 = fig.add_subplot(gs[0, 3:6])
ax3 = fig.add_subplot(gs[0, 6:9])
ax4 = fig.add_subplot(gs[0, 9:12])
ax5 = fig.add_subplot(gs[0, 12:15])

# Row 2 (5 columns)
ax6 = fig.add_subplot(gs[1, :3])
ax7 = fig.add_subplot(gs[1, 3:6])
ax8 = fig.add_subplot(gs[1, 6:9])
ax9 = fig.add_subplot(gs[1, 9:12])
ax10 = fig.add_subplot(gs[1, 12:15])

# Row 3 (5 columns)
ax11 = fig.add_subplot(gs[2, :3])
ax12 = fig.add_subplot(gs[2, 3:6])
ax13 = fig.add_subplot(gs[2, 6:9])
ax14 = fig.add_subplot(gs[2, 9:12])
ax15 = fig.add_subplot(gs[2, 12:15])

# Row 4 (3 columns)
ax16 = fig.add_subplot(gs[3, :4])
ax17 = fig.add_subplot(gs[3, 5:9])
ax18 = fig.add_subplot(gs[3, 10:14])

# Row 5 (3 columns)
ax19 = fig.add_subplot(gs[4, :4])
ax20 = fig.add_subplot(gs[4, 5:9])
ax21 = fig.add_subplot(gs[4, 10:14])

BK_axs = (ax1, ax2, ax3, ax4, ax5, ax16, ax19, ax6, ax7, ax8, ax9, ax10, ax17, ax20)
CaP_axs = (ax11, ax12, ax13, ax14, ax15, ax18, ax21)


BK_facs = [1,2,3,4,5]
CaP_facs = [1,2,3,4,5]

mesh = 'Cylinder3dia1umL10um'
capacfac = 50
endtime = 10.0

BK_Gs = [50, 100, 150, 200, 250]

SEEDS={}
SEEDS['Cabind']  = ['123', '234', '345', '456', '567', '678', '789', '890']
SEEDS['Canobind']= ['123', '234', '345', '456', '567', '678', '789', '890']

axidx = 0

for bindmode in ['Cabind', 'Canobind']:
    for i in range(len(BK_facs)):
        ax = BK_axs[axidx+i]
        
        BK_maxes=np.zeros([len(BK_Gs), len(SEEDS[bindmode])])
        Time_to_peak=[]
        bk_idx=0
        for BK_G in BK_Gs:
            cond_data = []
            seed_idx=0
            time_data=None
            for SEED in SEEDS[bindmode]:
                dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_{bindmode}_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK{BK_G}p_BKfac{BK_facs[i]}_CaPfac1'
                with HDF5Handler('./STEPS/data/'+dataset) as hdf:
                    Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
                    time_data=1e3 * Currents.time[0]
                    BKcond = BK_G*1e-3*(BKstates.data[0,:,1] + BKstates.data[0,:,3] + BKstates.data[0,:,5] + BKstates.data[0,:,7] + BKstates.data[0,:,9])
                    cond_data.append(BKcond)
                    BK_maxes[bk_idx][seed_idx] = np.max(BKcond)
                seed_idx+=1
            cond_data = np.mean(np.array(cond_data), axis=0)
            ax.plot(time_data, cond_data, label = BK_G)
            bk_idx+=1

            t_idx = 0
            for c in cond_data:
                if c < 0.1: t_idx+=1
                else: break
            t_start = time_data[t_idx]
            t_start=0.0 #override
            Time_to_peak.append(time_data[np.argmax(cond_data)]-t_start)

        BK_maxes = np.mean(BK_maxes, axis=1)

        ax.set_xlabel('Time (ms)')
        ax.set_xlim(1,3.2)
        if i==0:
            ax.legend(loc='best')
            ax.set_xlim(1,3.99)
            if bindmode == 'Cabind': ax.set_ylabel('BK g (Ca buffered) (nS)')
            else: ax.set_ylabel('BK g (Ca not buffered) (nS)')
        else: ax.set_yticklabels(())
        if bindmode == 'Cabind' or 1: ax.set_title(f'G_BK increase {BK_facs[i]}X', fontweight="bold")
        
        ax.set_ylim(0,500)
        
        if bindmode == 'Cabind': ax.set_ylim(0,400)
        else: ax.set_ylim(0,475)
        
        ax = BK_axs[axidx+5]
        ax.plot(BK_Gs, BK_maxes/BK_maxes[-1], '-', label=f'{BK_facs[i]}X', color='midnightblue', alpha=(0.2+0.2*i))

        ax = BK_axs[axidx+6]
        ax.plot(BK_Gs, Time_to_peak, '-', label=f'{BK_facs[i]}X', color='darkmagenta', alpha=(0.2+0.2*i))

    ax = BK_axs[axidx+5]
    ax.set_ylabel('Peak BK g (normed)')
    ax.set_xlabel('Unitary conductance (pS)')
    ax.set_yticks((0.5,1))
    ax.set_ylim(0.25,1.05)
    if bindmode == 'Cabind':
        ax.set_title(f'BK (Ca buffered)', fontweight="bold")
    else:
        ax.set_title(f'BK (Ca not buffered)', fontweight="bold")
        ax.legend(loc='best')

    ax = BK_axs[axidx+6]

    ax.set_ylabel('Time to peak BK g (ms)')
    ax.set_xlabel('Unitary conductance (pS)')
    ax.set_ylim(0.36+1.2,0.66+1.2)
    if bindmode == 'Cabind':  ax.set_title(f'BK (Ca buffered)', fontweight="bold")
    else:
        ax.set_title(f'BK (Ca not buffered)', fontweight="bold")
        ax.legend()
    
    axidx+=7


seed_to_p = { ('1','11','111','1111','11111','111111','1111111','11111111'):2.5e-2, ('2','22','222','2222','22222','222222','2222222','22222222'):5.0e-2, ('3','33','333','3333','33333','333333','3333333','33333333'):7.5e-2, ('4','44','444','4444','44444','444444','4444444','44444444'):10.0e-2, ('5','55','555','5555','55555','555555','5555555','55555555'):12.5e-2}

CaP_ps = seed_to_p.values()

axidx = 0

for i in range(len(CaP_facs)):
    ax = CaP_axs[axidx+i]
    CaP_maxes=[]
    Time_to_peak=[]
    for seeds in seed_to_p.keys():
        CaP_P = seed_to_p[seeds]
        Cap_max = []
        tpeak=[]
        plotted=False # Can only plot one, because the mean gives misleading results because of small timing differences
        for seed in seeds:
            dataset = f'BKmodel_axononly_'+seed+f'_{endtime}ms_Canobind_{mesh}_SKfac0.0_TEMP34.0_capacfac{capacfac}_BK250p_BKfac1_CaPfac{CaP_facs[i]}'
            with HDF5Handler('./STEPS/data_CaPstates/'+dataset) as hdf:
                Currents, CaPstates, CaConcs, Pot =  hdf['BKmodel_axononlySim'].results
                Cap_max.append(np.max(CaPstates.data[0,:,1]*CaP_P))
                if not plotted:
                    ax.plot(1e3 * Currents.time[0], CaPstates.data[0,:,1]*CaP_P, label = seed_to_p[seeds])
                    plotted=True
                tpeak.append(1e3 * Currents.time[0][np.argmax(CaPstates.data[0,:,1])]-0.0) 
        CaP_maxes.append(np.mean(Cap_max))
        Time_to_peak.append(np.mean(tpeak))

        if i==0:
            ax.legend(loc='best')
            ax.set_ylabel('CaP p ($µm^3s^{-1}$)')
        else:
            ax.set_yticks(())
        ax.set_xlabel('Time (ms)')
        ax.set_xlim(1,2)
        ax.set_ylim(0,3500)
        ax.set_title(f'P_CaP increase {CaP_facs[i]}X', fontweight="bold")
    
    ax18.plot(CaP_ps, CaP_maxes/np.min(CaP_maxes), '-', label=f'{CaP_facs[i]}X', color='midnightblue', alpha=(0.2+0.2*i))
    ax18.set_ylabel('Peak CaP p (normed)')
    ax18.set_xlabel('Unitary permeability ($µm^3s^{-1}$)')
    ax18.set_ylim(0.25,1.05)
    ax18.set_yticks((0.5,1))
    ax18.set_title(f'CaP', fontweight="bold")

    ax21.plot(CaP_ps, Time_to_peak, '-', label=f'{CaP_facs[i]}X', color='darkmagenta', alpha=(0.2+0.2*i))
    ax21.set_ylabel('Time to peak CaP p (ms)')
    ax21.set_xlabel('Unitary permeability ($µm^3s^{-1}$)')
    ax21.set_ylim(0.05+1.2,0.35+1.2)
    ax21.set_title(f'CaP', fontweight="bold")


ax1.text(-0.25,1.2,'a', fontweight = 'bold', transform=ax1.transAxes)
ax6.text(-0.25,1.2,'b', fontweight = 'bold', transform=ax6.transAxes)
ax11.text(-0.25,1.2,'c', fontweight = 'bold', transform=ax11.transAxes)

ax16.text(-0.15,1.13,'d', fontweight = 'bold', transform=ax16.transAxes)
ax19.text(-0.15,1.13,'e', fontweight = 'bold', transform=ax19.transAxes)

# Show the plot
plt.subplots_adjust(wspace=2, hspace=0.7)
fig = ax.get_figure()
plt.savefig('Figures/Figure4.pdf', dpi=300)
plt.close()
