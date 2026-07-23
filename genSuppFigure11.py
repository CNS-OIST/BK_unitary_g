import matplotlib.pyplot as plt
#plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from math import pi
plt.style.use('./figures.naturestyle')

dx = '200nm'


BK_facs = CaP_facs = [1, 2, 3, 4]
BK_facs=np.array(BK_facs)


fig, (axis1, axis2) = plt.subplots(2, 3, figsize=[7.08*1.3, 5.2*1.3])
ax1, ax2, ax3=axis1[0], axis1[1], axis1[2]
ax4, ax5, ax6=axis2[0], axis2[1], axis2[2]


BK250_maxes=[]
BK50_maxes=[]
CaP_maxes=[]


BK_G = '250'
for CaP_fac in CaP_facs:
        neuronfile = open(f'NEURON/simple_model_{BK_G}pS/data_{dx}/BK1_CaP{CaP_fac}/CaP_I.dat', 'r')
        neuronfile_lines = neuronfile.readlines()
        neuronv=[]
        neuront=[]

        for line in neuronfile_lines[2:]:
            line=line.split()
            neuront.append(float(line[0]))
            neuronv.append(float(line[1])*(pi/10))
        
        ax1.plot(neuront, neuronv, label = CaP_fac )
        CaP_maxes.append(np.max(neuronv))
        
        
        neuronfile = open(f'NEURON/simple_model_{BK_G}pS/data_{dx}/BK1_CaP{CaP_fac}/CaP_m.dat', 'r')
        neuronfile_lines = neuronfile.readlines()
        neuronv=[]
        neuront=[]

        for line in neuronfile_lines[2:]:
            line=line.split()
            neuront.append(float(line[0]))
            neuronv.append(float(line[1]))
            
        ax4.plot(neuront, 100*np.array(neuronv), label = CaP_fac )



for BK_fac in BK_facs:
        neuronfile = open(f'NEURON/simple_model_{BK_G}pS/data_{dx}/BK{BK_fac}_CaP1/BK_I.dat', 'r')
        neuronfile_lines = neuronfile.readlines()
        neuronv=[]
        neuront=[]

        for line in neuronfile_lines[2:]:
            line=line.split()
            neuront.append(float(line[0]))
            neuronv.append(float(line[1])*(pi/10))

        ax2.plot(neuront, neuronv, label = BK_fac )
        BK250_maxes.append(np.max(neuronv))
        
        
        neuronfile = open(f'NEURON/simple_model_{BK_G}pS/data_{dx}/BK{BK_fac}_CaP1/BK_g.dat', 'r')
        neuronfile_lines = neuronfile.readlines()
        neuronv=[]
        neuront=[]

        for line in neuronfile_lines[2:]:
            line=line.split()
            neuront.append(float(line[0]))
            neuronv.append(float(line[1])/(6*BK_fac)) # 6 is base g
            
        ax5.plot(neuront, 100*np.array(neuronv), label = BK_fac )

BK_G = '50'

for BK_fac in BK_facs:
        neuronfile = open(f'NEURON/simple_model_{BK_G}pS/data_{dx}/BK{BK_fac}_CaP1/BK_I.dat', 'r')
        neuronfile_lines = neuronfile.readlines()
        neuronv=[]
        neuront=[]

        for line in neuronfile_lines[2:]:
            line=line.split()
            neuront.append(float(line[0]))
            neuronv.append(float(line[1])*(pi/10))

        ax3.plot(neuront, neuronv, label = BK_fac )
        BK50_maxes.append(np.max(neuronv))
        
        
        neuronfile = open(f'NEURON/simple_model_{BK_G}pS/data_{dx}/BK{BK_fac}_CaP1/BK_g.dat', 'r')
        neuronfile_lines = neuronfile.readlines()
        neuronv=[]
        neuront=[]

        for line in neuronfile_lines[2:]:
            line=line.split()
            neuront.append(float(line[0]))
            neuronv.append(float(line[1])/(6*BK_fac)) # 6 is base g
            
        ax6.plot(neuront, 100*np.array(neuronv), label = BK_fac )

axins1 = inset_axes(ax1, width=0.7, height=0.7, bbox_to_anchor=(0.65, .1, .3, .5), bbox_transform=ax1.transAxes)
axins1.plot(range(1,6), range(1,6), 'k--')
for i in range(4): axins1.plot(CaP_facs[i], CaP_maxes[i]/CaP_maxes[0], 'o')

axins2 = inset_axes(ax2, width=0.7, height=0.7, bbox_to_anchor=(0.67, .5, .3, .5), bbox_transform=ax2.transAxes)
axins2.plot(range(1,6), range(1,6), 'k--')
for i in range(4): axins2.plot(BK_facs[i]/BK_facs[0], BK250_maxes[i]/BK250_maxes[0], 'o')

axins3 = inset_axes(ax3, width=0.7, height=0.7, bbox_to_anchor=(0.65, .5, .3, .5), bbox_transform=ax3.transAxes)
axins3.plot(range(1,6), range(1,6), 'k--')
for i in range(4): axins3.plot(BK_facs[i]/BK_facs[0], BK50_maxes[i]/BK50_maxes[0], 'o')


ax1.legend(CaP_facs, loc='upper right',facecolor='white', framealpha=1)
ax1.set_ylabel('CaP current (nA)')
ax1.set_xlim(1.4,2.4)
ax1.set_xticks((1.5, 2))

for ax in (ax2, ax3):
    ax.legend(BK_facs, loc='lower right',facecolor='white', framealpha=1)
    ax.set_xlim(1.5,3.0)
    ax.set_xticks((1.5, 2, 2.5, 3))
    ax.set_ylim(0,6)
ax2.set_ylabel(f'BK250 current (nA)')
ax3.set_ylabel(f'BK50 current (nA)')


ax4.legend(loc='best')
ax4.set_xlabel('Time (ms)')
ax4.set_ylabel('% of activated CaP channels')
ax4.set_xlim(1.4,2.4)
#ax4.set_ylim(0,105)
ax4.set_xticks((1.5, 2))

for ax in (ax5, ax6):
    ax.legend(loc='best')
    ax.set_xlabel('Time (ms)')
    ax.set_xlim(1.4,4.0)
    ax.set_ylim(0,9)

ax5.set_ylabel(f'% of activated BK250 channels')
ax6.set_ylabel(f'% of activated BK50 channels')

for axins in (axins1, axins2, axins3):
    axins.set_ylabel('Peak I')
    axins.set_xlabel('Increase factor')
    axins.set_yticks(())
    axins.set_xticks(())
    axins.set_ylim(0.6, 2.3)
    axins.set_xlim(0.7, 4.3)
axins1.set_ylim(0.7, 4.3)
axins1.set_ylabel('Peak I')

ax1.text(-0.15,1.05,'a', fontweight = 'bold', transform=ax1.transAxes)
ax2.text(-0.15,1.05,'b', fontweight = 'bold', transform=ax2.transAxes)
ax3.text(-0.15,1.05,'c', fontweight = 'bold', transform=ax3.transAxes)
ax4.text(-0.15,1.05,'d', fontweight = 'bold', transform=ax4.transAxes)
ax5.text(-0.15,1.05,'e', fontweight = 'bold', transform=ax5.transAxes)
ax6.text(-0.15,1.05,'f', fontweight = 'bold', transform=ax6.transAxes)

plt.subplots_adjust(wspace=0.4, hspace=0.3)
plt.savefig(f'Figures_Sup/FigureS11.pdf', dpi=300)
plt.close()

