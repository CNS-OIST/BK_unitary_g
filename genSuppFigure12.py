import matplotlib.pyplot as plt
import numpy as np
from numpy import array
import scipy

gs = [50, 100, 150, 200, 250]

fig, (axis1, axis2) = plt.subplots(2, 2, figsize=[7.08, 5.2])
ax1, ax2 = axis1[0], axis1[1]
ax3, ax4 = axis2[0], axis2[1]
ax5 = ax4.twinx()

BK_gs=gs

data_dir = './NEURON/purkinje_pub_BK_10nm/data/'

comp='ais'
BK_maxes=[]
peaks_ca = []
widths_ca = []

for g in gs:
    neuronO0=[]
    neuronO1=[]
    neuronO2=[]
    neuronO3=[]
    neuronO4=[]

    neuront=[]

    neuronfileI = open(data_dir+f'{g}pS/ais_O0.dat', 'r')
    neuronfileI_lines = neuronfileI.readlines()
    for line in neuronfileI_lines[2:]:
        line=line.split()
        neuront.append(float(line[0]))
        neuronO0.append(float(line[1]))
    neuronfileI = open(data_dir+f'{g}pS/ais_O1.dat', 'r')
    neuronfileI_lines = neuronfileI.readlines()
    for line in neuronfileI_lines[2:]:
        line=line.split()
        neuronO1.append(float(line[1]))

    neuronfileI = open(data_dir+f'{g}pS/ais_O2.dat', 'r')
    neuronfileI_lines = neuronfileI.readlines()
    for line in neuronfileI_lines[2:]:
        line=line.split()
        neuronO2.append(float(line[1]))

    neuronfileI = open(data_dir+f'{g}pS/ais_O3.dat', 'r')
    neuronfileI_lines = neuronfileI.readlines()
    for line in neuronfileI_lines[2:]:
        line=line.split()
        neuronO3.append(float(line[1]))

    neuronfileI = open(data_dir+f'{g}pS/ais_O4.dat', 'r')
    neuronfileI_lines = neuronfileI.readlines()
    for line in neuronfileI_lines[2:]:
        line=line.split()
        neuronO4.append(float(line[1]))
    
    neuronG = 6*(array(neuronO0)+array(neuronO1)+array(neuronO2)+array(neuronO3)+array(neuronO4)) # S/cm2
    ax1.plot(neuront, neuronG, label=f'{g}pS', linewidth=3)

    neuronfileca = open(data_dir+f'{g}pS/{comp}_ca.dat', 'r')
    neuronfileca_lines = neuronfileca.readlines()
    neuronca=[]
    neuront=[]
    for line in neuronfileca_lines[2:]:
        line=line.split()
        neuronca.append(float(line[1])*1e3)
        neuront.append(float(line[0]))

    ax3.plot(neuront, neuronca, label=f'{g}pS', linewidth=3)
    
    BK_maxes.append(np.max(neuronG))

    ps = scipy.signal.find_peaks(neuronca, prominence=0.5)[0]
    ws = 0.02*np.array(scipy.signal.peak_widths(neuronca, ps, rel_height=0.9))[0]
    peaks_ca.append(np.mean(np.take(neuronca, ps)))
    widths_ca.append(np.mean(ws))
    
ax4.plot(BK_gs, peaks_ca, 'bo-')
ax4.set_ylabel('Submemb Ca peak heights (µM)')
ax4.set_xlabel('Unitary BK conductance (pS)')
ax4.yaxis.label.set_color('blue')
ax4.tick_params(axis='y', colors='blue')

ax5.plot(BK_gs, widths_ca, 'ro-')
ax5.set_ylabel('Submemb Ca peak width (ms)')
ax5.yaxis.label.set_color('red')
ax5.tick_params(axis='y', colors='red')

ax2.plot(BK_gs, BK_maxes, 'o-', label = comp)
ax2.set_xlabel('Unitary BK conductance (pS)')
ax2.set_ylabel('Max BK conductance (S/$cm^2$)')

ax1.set_xlim(20.3, 25.0)
ax1.set_ylabel('BK conductance (S/$cm^2$)')
ax1.set_xlabel('Time (ms)')

ax3.legend()

ax3.set_xlim(20.5, 25.5)
ax3.set_ylabel('Submemb Ca (µM)')
ax3.set_xlabel('Time (ms)')


plt.subplots_adjust(wspace=0.6, hspace=0.4)
plt.savefig(f'Figures_Sup/FigureS12.pdf', dpi=300)
plt.close()
