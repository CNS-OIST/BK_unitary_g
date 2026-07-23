import steps.interface
from steps.saving import *
from steps.geom import *

import matplotlib.pyplot as plt
plt.style.use('./figures.naturestyle')
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import scipy

SEED='123'
NA = 6.02214076e23 # Avogadro's number

BK_gs=[250, ]


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

"""

mesh_file =   ['./STEPS/meshes/Cylinder_dia1um_L1um_size0.2_602tets.msh', 'CylinderL1um']
mesh_file =   ['./STEPS/meshes/Cylinder_dia1um_L5um_size0.2_2485tets.msh', 'CylinderL5um']
mesh_file =   ['./STEPS/meshes/Cylinder_dia1um_L10um_size0.2_4922tets.msh', 'CylinderL10um']
mesh_file =   ['./STEPS/meshes/Cylinder_dia1um_L15um_size0.2_7306tets.msh', 'CylinderL15um']
mesh_file =   ['./STEPS/meshes/Cylinder_dia1um_L18.72um_size0.2_9035tets.msh', 'CylinderL18.72um']

mesh = TetMesh.LoadGmsh(mesh_file[0], 1e-6)

with mesh:
    cyto = Compartment.Create(mesh.tets)
    
    if mesh_file[1][:3] == 'Cyl':
        ends = [cyto.bbox.min.z, cyto.bbox.max.z]
        memb_tris = TriList(tri for tri in cyto.surface if tri.center.z not in ends)
    else:
        memb_tris = cyto.surface    
    
    memb = Patch.Create(memb_tris, cyto)
    print ("Area: ", memb.Area)
exit()
"""


# For cylinders, size 0.2
L_to_area = {1: 3.1264643140246244e-12, 5:1.5629165712911503e-11, 10:3.125759585902539e-11, 15:4.6885405505897675e-11, 18.72:5.85129770004717e-11}


idx=-1

shift=0.0
CV_at_peak={}
CV_at_peak[50]=[]
CV_at_peak[250]=[]

STD_at_peak={}
STD_at_peak[50]=[]
STD_at_peak[250]=[]

Q = {}
Q[50]=[]
Q[250]=[]


lengths=[1,5,10,15,18.72]


indices = lengths
    
for i in indices:
    mesh = f'CylinderL{i}um'
        
    for BK_g in BK_gs:
        dataset = f'BKmodel_axononly_'+SEED+f'_{endtime}ms_Cabind_{mesh}_TEMP34.0_capacfac1_BK{BK_g}p'
        with HDF5Handler('./STEPS/data_BKonly/'+dataset) as hdf:
            Currents, BKstates, CaConcs, Pot =  hdf['BKmodel_BKonlySim'].results
            
            area= L_to_area[i]
            lab = 'Length'

            # 0.1 converts amp/m2 to mA/cm2
            plt.errorbar(1e3 * Currents.time[0]+shift, 0.1 * np.mean(Currents.data[:,:,0], axis=0)/area, 0.1 * np.std(Currents.data[:,:,0], axis=0)/area, label = f'{lab}={i}µm, {BK_g}pS BK current mean and std')
            plt.ylabel('Current (mA/cm2)')

            #CV
            #plt.plot(1e3 * Currents.time[0]+shift, np.std(Currents.data[:,:,0], axis=0)/np.mean(Currents.data[:,:,0], axis=0), label = f'Rad={R}µm, {BK_g}pS')
            
            # STD
            #plt.plot(1e3 * Currents.time[0]+shift, 0.1 * np.std(Currents.data[:,:,0], axis=0)/R_to_area[R], label = f'Rad={R}µm, {BK_g}pS')
            #plt.ylabel('Current std (mA/cm2)')

            peak_idx =  np.argmax(np.mean(Currents.data[:,:,0], axis=0))

            CV_at_peak[BK_g].append(np.std(Currents.data[:,:,0], axis=0)[peak_idx]/np.mean(Currents.data[:,:,0], axis=0)[peak_idx])
            
            STD_at_peak[BK_g].append(0.1 * np.max(np.std(Currents.data[:,:,0], axis=0))/area)
            
            Q[BK_g].append(np.std(np.sum(Currents.data[:,:,0], axis=1)/area))
            
        shift+=0.0001
    shift+=0.0004

#plt.savefig(f'Figures/Figure3A-F_{mesh}.pdf', dpi=300)
#plt.ylabel('CV')
plt.xlim(1.2, 2.0)
plt.xlabel('Time (ms)')
plt.legend()
plt.show()
plt.close()

for BK_g in BK_gs:
    plt.plot(indices, CV_at_peak[BK_g], label = f'{BK_g}pS BK current')
plt.ylabel('CV at current peak')
plt.xlabel('Cylinder length (µm)')
plt.legend()
plt.show()

for BK_g in BK_gs:
    plt.plot(indices, STD_at_peak[BK_g], label = f'{BK_g}pS BK')
plt.ylabel('STD at current peak (mA/cm2)')
plt.xlabel('Cylinder length (µm)')
plt.legend()
plt.show()

for BK_g in BK_gs:
    plt.plot(indices, Q[BK_g], label = f'{BK_g}pS BK')

plt.ylabel('STD of total integrated current (arb)')
plt.xlabel('Cylinder length (µm)')
plt.legend()
plt.show()


