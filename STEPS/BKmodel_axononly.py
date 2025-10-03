# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import steps.interface

from steps.geom import *
from steps.model import *
from steps.rng import *
from steps.saving import *
from steps.sim import *
from steps.utils import *
SetVerbosity(10)

import math
import sys

from myconstants import *

import time
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###########################################################
# Simulation Parameters
###########################################################

SEED = int(sys.argv[1])
print (sys.argv)
BK_fac = int(sys.argv[2]) if len(sys.argv)>2 else 1
CaP_fac = int(sys.argv[3]) if len(sys.argv)>3 else 1

Camode = 'Cabind'
#Camode = 'Canobind'

CaPstates = False

NBRUNS = 1

EF_DT = 2.0e-6

DT =  1.0e-5

ENDT = 10.0e-3
ENDT = 100.0e-3


###########################################################
# Mesh Parameters
###########################################################

mesh_file =  ['./meshes/Cylinder3_dia1um_L10um_noouter_0.3shell_0.3size_13617tets_adaptive.inp', 'Cylinder3dia1umL10um'] # mesh file and label

###########################################################
# Biochemical model
###########################################################

mdl = Model()
with mdl:
    # Species
    Pump, CaPump, InfluxPump, CaInfluxPump, PV, PVMg, PVCa, Mg = Species.Create()
    Ca = Species.Create(valence=2)

    # Calbindin
    CBs, CBf, CBsCa, CBfCa, CBmob, CBimmob = SubUnitState.Create()
    CBsSU, CBfSU, CBmobSU = SubUnit.Create(
        [CBs, CBsCa], [CBf, CBfCa], [CBmob, CBimmob]
    )
    CB = Complex.Create([CBsSU, CBfSU, CBmobSU], statesAsSpecies=True)

    # Channels
    Kv3c, Kv3o = SubUnitState.Create()
    Kv3_SU = SubUnit.Create([Kv3c, Kv3o])
    Kv3chan = Channel.Create([Kv3_SU, Kv3_SU, Kv3_SU, Kv3_SU])

    NaPc, NaPo = SubUnitState.Create()
    NaP_SU = SubUnit.Create([NaPc, NaPo])
    NaPchan = Channel.Create([NaP_SU, NaP_SU, NaP_SU])

    Rsg_o, Rsg_c, Rsg_C, Rsg_I, Rsg_on, Rsg_off, Rsg_b = SubUnitState.Create()
    RsgocSU, RsgCISU, RsgonoffSU = SubUnit.Create(
        [Rsg_o, Rsg_c],[Rsg_C, Rsg_I],[Rsg_on, Rsg_off, Rsg_b]
    )
    Rsgchan = Channel.Create([RsgocSU, RsgocSU, RsgocSU, RsgocSU, RsgCISU, RsgonoffSU])
    Rsgchan_C5 = Rsgchan[Rsg_o, Rsg_o, Rsg_o, Rsg_o, Rsg_C, Rsg_off]
    Rsgchan_O  = Rsgchan[Rsg_o, Rsg_o, Rsg_o, Rsg_o, Rsg_C, Rsg_on]
    Rsgchan_B  = Rsgchan[Rsg_o, Rsg_o, Rsg_o, Rsg_o, Rsg_C, Rsg_b]
    Rsgchan_I5 = Rsgchan[Rsg_o, Rsg_o, Rsg_o, Rsg_o, Rsg_I, Rsg_off]
    Rsgchan_I6 = Rsgchan[Rsg_o, Rsg_o, Rsg_o, Rsg_o, Rsg_I, Rsg_on]
    RsgStates = Rsgchan[..., Rsg_off] | Rsgchan_O | Rsgchan_B | Rsgchan_I6
    
    CaPc, CaPo = SubUnitState.Create()
    CaPchan = Channel.Create([CaPc, CaPo])
    
    CaTmc, CaTmo, CaThc, CaTho = SubUnitState.Create()
    CaTm_SU = SubUnit.Create([CaTmc, CaTmo])
    CaTh_SU = SubUnit.Create([CaThc, CaTho])
    CaTchan = Channel.Create([CaTm_SU, CaTm_SU, CaTm_SU, CaTh_SU])

    BKslw, BKslwCa, BKslwopen, BKslwclose = SubUnitState.Create()
    BKslwCaSU = SubUnit.Create([BKslw, BKslwCa])
    BKslwocSU = SubUnit.Create([BKslwopen, BKslwclose])
    BKslwchan = Channel.Create([BKslwCaSU, BKslwCaSU, BKslwCaSU, BKslwCaSU, BKslwocSU])

    BK, BKCa, BKopen, BKclose = SubUnitState.Create()
    BKCaSU = SubUnit.Create([BK, BKCa])
    BKocSU = SubUnit.Create([BKopen, BKclose])
    BKchan = Channel.Create([BKCaSU, BKCaSU, BKCaSU, BKCaSU, BKocSU])

    SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2 = SubUnitState.Create()
    SKchan = Channel.Create([SK_C1, SK_C2, SK_C3, SK_C4, SK_O1, SK_O2])

    Ihc, Iho = SubUnitState.Create()
    Ihchan = Channel.Create([Ihc, Iho])

    Leak = SubUnitState.Create()
    L = Channel.Create([Leak])
    
    r = ReactionManager()
    
    vsys = VolumeSystem.Create()
    with vsys:
        # PVCa
        PV + Ca <r[1]> PVCa
        r[1].K = PVca_f, PVca_b

        # PVMg
        PV + Mg <r[1]> PVMg
        r[1].K = PVmg_f_kcst, PVmg_b_kcst

        # Ca Influx Pump
        InfluxPump + Ca <r[1]> CaInfluxPump >r[2]> InfluxPump
        r[1].K = Inf_f, Inf_b
        r[2].K = Inf_k

        with CB[...]:
            # Fast binding
            CBf + Ca <r[1]> CBfCa
            r[1].K = CBf_f_kcst, CBf_b_kcst

            # Slow binding
            CBs + Ca <r[1]> CBsCa
            r[1].K = CBs_f_kcst, CBs_b_kcst

        diff_Ca   = Diffusion.Create(Ca, DCST)
        diff_CB   = Diffusion.Create(CB[:, :, CBmob], DCB)
        diff_InfluxPump   = Diffusion.Create(InfluxPump, DCB)
        diff_PV   = Diffusion.Create(PV, DPV)
        diff_PVCa = Diffusion.Create(PVCa, DPV)
        diff_PVMg = Diffusion.Create(PVMg, DPV)
    
    ssys = SurfaceSystem.Create()
    with ssys:
        
        # Ca Pump
        Pump.s + Ca.i <r[1]> CaPump.s >r[2]> Pump.s
        r[1].K = P_f, P_b
        r[2].K = P_k

        # Kv3 channel
        with Kv3chan[...]:
            Kv3c.s <r[1]> Kv3o.s
            r[1].K = a_n, b_n
        OC_Kv3 = OhmicCurr.Create(Kv3chan[Kv3o, Kv3o, Kv3o, Kv3o], Kv3_G, Kv3_rev)

        # NaP channel
        with NaPchan[...]:
            NaPc.s <r[1]> NaPo.s
            r[1].K = a_m, b_m
        OC_NaP = OhmicCurr.Create(NaPchan[NaPo, NaPo, NaPo], NaP_G, NaP_rev)

        # Rsg channel
        with Rsgchan[..., Rsg_C, Rsg_off]:
            Rsg_c.s <r[1]> Rsg_o.s
            r[1].K = rsgC_f, rsgC_b

        with Rsgchan[..., Rsg_I, Rsg_off]:
            Rsg_c.s <r[1]> Rsg_o.s
            r[1].K = rsgI_f, rsgI_b

        with Rsgchan[..., Rsg_off]:
            Rsg_C.s <r[1]> Rsg_I.s

            rsg_f = CompDepRate.Create(lambda state: rsgCI_f * ((alfac_rsg) ** state.Count(Rsg_o)), [Rsgchan])
            rsg_b = CompDepRate.Create(lambda state: rsgCI_b * ((btfac_rsg) ** state.Count(Rsg_o)), [Rsgchan])
            r[1].K = rsg_f, rsg_b
        
        Rsgchan_C5.s <r[1]> Rsgchan_O.s #C5 <-> O
        Rsgchan_I5.s <r[2]> Rsgchan_I6.s #I5 <-> I6
        Rsgchan_O.s <r[3]> Rsgchan_I6.s #O <-> I6
        Rsgchan_O.s <r[4]> Rsgchan_B.s #O <-> B
        r[1].K = rsgCO_f, rsgCO_b
        r[2].K = rsgIn_f, rsgIn_b
        r[3].K = rsgCIn_f, rsgCIn_b
        r[4].K = rsgOB_f, rsgOB_b
            
        OC_Rsg = OhmicCurr.Create(Rsgchan[..., Rsg_C, Rsg_on], Rsg_G, Rsg_rev)
    
        # CaP channel
        with CaPchan[...]:
            CaPc.s <r[1]> CaPo.s
            r[1].K = alpha_cap, beta_cap
        OC_CaP = GHKCurr.Create(
            CaPchan[CaPo], Ca, CaP_P,
            computeflux=True,
            virtual_oconc=Ca_oconc,
        )
        
        # CaT channel
        with CaTchan[...]:
            CaTmc.s <r[1]> CaTmo.s
            r[1].K = alpham_cat, betam_cat
            
            CaThc.s <r[1]> CaTho.s
            r[1].K = alphah_cat, betah_cat
        OC_CaT = GHKCurr.Create(
            CaTchan[CaTmo, CaTmo, CaTmo, CaTho], Ca, CaT_P,
            computeflux=False,
            virtual_oconc=Ca_oconc,
        )

        # BKslw channel # ignore Ca
        with BKslwchan[..., BKslwclose]:
            BKslw.s + Ca.i >r[1]> BKslwCa.s + Ca.i
            BKslwCa.s >r[2]> BKslw.s
            r[1].K = BKslw_f
            r[2].K = BKslwc_b

        with BKslwchan[..., BKslwopen]:
            BKslw.s + Ca.i >r[1]> BKslwCa.s + Ca.i
            BKslwCa.s >r[2]> BKslw.s
            r[1].K = BKslw_f
            r[2].K = BKslwo_b

        with BKslwchan[...]:
            BKslwclose.s <r[1]> BKslwopen.s

            BKslw_f = CompDepRate.Create(lambda s: BKslw_oc_f[s.Count(BKslwCa)], [BKslwchan])
            BKslw_b = CompDepRate.Create(lambda s: BKslw_oc_b[s.Count(BKslwCa)], [BKslwchan])
            r[1].K = BKslw_f, BKslw_b
        OC_BKslw = OhmicCurr.Create(BKslwchan[..., BKslwopen], BKslw_G, BKslw_rev)
        
        # BK channel # ignore Ca
        
        if Camode == 'Cabind':
            with BKchan[..., BKclose]:
                BK.s + Ca.i <r[1]> BKCa.s
                r[1].K = BK_f, BKc_b

            with BKchan[..., BKopen]:
                BK.s + Ca.i <r[1]> BKCa.s
                r[1].K = BK_f, BKo_b

        elif Camode == 'Canobind':
            with BKchan[..., BKclose]:
                BK.s + Ca.i >r[1]> BKCa.s + Ca.i
                BKCa.s >r[2]> BK.s
                r[1].K = BK_f
                r[2].K = BKc_b

            with BKchan[..., BKopen]:
                BK.s + Ca.i >r[1]> BKCa.s + Ca.i
                BKCa.s >r[2]> BK.s
                r[1].K = BK_f
                r[2].K = BKo_b
        
        else: raise Exception("Unkown Ca bind mode!")
        
        with BKchan[...]:
            BKclose.s <r[1]> BKopen.s

            BK_f = CompDepRate.Create(lambda s: BK_oc_f[s.Count(BKCa)], [BKchan])
            BK_b = CompDepRate.Create(lambda s: BK_oc_b[s.Count(BKCa)], [BKchan])
            r[1].K = BK_f, BK_b
        OC_BK = OhmicCurr.Create(BKchan[..., BKopen], BK_G, BK_rev)

        # SK channel ignore Ca
        with SKchan[...]:
            
            if Camode == 'Cabind':

                ((SK_C1.s + Ca.i <r[1]> SK_C2.s)\
                        + Ca.i <r[2]> SK_C3.s)\
                        + Ca.i <r[3]> SK_C4.s
                r[1].K = dirc2_t, invc1_t
                r[2].K = dirc3_t, invc2_t
                r[3].K = dirc4_t, invc3_t
        
            elif Camode == 'Canobind':

                SK_C1.s + Ca.i >r[1]> SK_C2.s + Ca.i
                SK_C2.s + Ca.i >r[2]> SK_C3.s + Ca.i
                SK_C3.s + Ca.i >r[3]> SK_C4.s + Ca.i
                SK_C4.s >r[6]> SK_C3.s
                SK_C3.s >r[5]> SK_C2.s
                SK_C2.s >r[4]> SK_C1.s
                r[1].K = dirc2_t
                r[2].K = dirc3_t
                r[3].K = dirc4_t
                r[4].K = invc1_t
                r[5].K = invc2_t
                r[6].K = invc3_t
            
            else: raise Exception("Unkown Ca bind mode!")

            SK_C3.s <r[1]> SK_O1.s
            SK_C4.s <r[2]> SK_O2.s
            r[1].K = diro1_t, invo1_t
            r[2].K = diro2_t, invo2_t
        OC_SK = OhmicCurr.Create(SKchan[SK_O1|SK_O2], SK_G, SK_rev)

        # Ih current channel
        with Ihchan[...]:
            Ihc.s <r[1]> Iho.s
            r[1].K = alphan_ih, betan_ih
        OC_Ih = OhmicCurr.Create(Ihchan[Iho], Ih_G, Ih_rev)

        # Leak current channel
        OC_L = OhmicCurr.Create(L[Leak], L_G, L_rev)

###########################################################
# Mesh and compartmentalization
###########################################################

mesh = TetMesh.LoadAbaqus(mesh_file[0], 1e-6)

with mesh:

    cyto = Compartment.Create(mesh.tets, vsys)

    ends = [cyto.bbox.min.z, cyto.bbox.max.z]
    memb_tris = TriList(tri for tri in cyto.surface if tri.center.z not in ends)
    memb = Patch.Create(memb_tris, cyto, None, ssys)
    
    submemb_tets = TetList()
    for tri in memb.tris:
            submemb_tets |= tri.tetNeighbs
    
    membrane = Membrane.Create([memb])
    
###########################################################
# Simulation
###########################################################

rng = RNG('mt19937', 512, SEED)

part = LinearMeshPartition(mesh, 1, 1, MPI.nhosts)

sim = Simulation('TetOpSplit', mdl, mesh, rng, part, calcMembPot=MPI.EF_DV_BDSYS)

rs = ResultSelector(sim)

Currents = rs.SUM(rs.TRIS(memb.tris).OC_Kv3.I) << \
           rs.SUM(rs.TRIS(memb.tris).OC_NaP.I) << \
           rs.SUM(rs.TRIS(memb.tris).OC_Rsg.I) << \
           rs.SUM(rs.TRIS(memb.tris).OC_CaP.I) << \
           rs.SUM(rs.TRIS(memb.tris).OC_CaT.I) << \
           rs.SUM(rs.TRIS(memb.tris).OC_BK.I) << \
           rs.SUM(rs.TRIS(memb.tris).OC_SK.I) << \
           rs.SUM(rs.TRIS(memb.tris).OC_Ih.I) << \
           rs.SUM(rs.TRIS(memb.tris).OC_L.I)


if CaPstates:
    Chanstates = rs.memb.LIST(*CaPchan[...]).Count
else:
    Chanstates = rs.memb.LIST(*BKchan[...]).Count


CaConcs = rs.cyto.Ca.Conc << \
         (rs.SUM(rs.TETS(submemb_tets).Ca.Count) / (AVOGADRO * submemb_tets.Vol * 1e3))

Pot = rs.TET(0).V

sim.toSave(Currents, Chanstates, CaConcs, Pot, dt=DT)

#Capump off, ignore Ca on BK and SK

if CaPstates:
    datadir = 'data_CaPstates'
else:
    datadir = 'data'


filename = f'{datadir}/BKmodel_axononly_{SEED}_{ENDT*1e3}ms_{Camode}_{mesh_file[1]}_SKfac{SK_fac}_TEMP{TEMPERATURE-273.15}_capacfac{capac_fac}_BK{BK_cond}p_BKfac{BK_fac}_CaPfac{CaP_fac}'

with HDF5Handler(filename) as hdf:
    sim.toDB(hdf, 'BKmodel_axononlySim')
    for i in range(NBRUNS):

                sim.newRun()
        
                # Setting initial conditions
                area = Parameter(memb.Area, 'm^2')
                
                for s in Kv3chan[...]:
                    sim.memb.LIST(s).Count = round(Kv3_ro*area*Kv3_p[s.Count(Kv3o)])
                
                for s in NaPchan[...]:
                    sim.memb.LIST(s).Count = round(NaP_ro*area*NaP_p[s.Count(NaPo)])
                
                for s in Rsgchan[..., Rsg_off]:
                    ocCnt, CICnt = s.Count(Rsg_o), s.Count(Rsg_I)
                    sim.memb.LIST(s).Count = round(Rsg_ro*area*Rsg_p[CICnt][ocCnt])
                    
                sim.memb.LIST(Rsgchan_O).Count = round(Rsg_ro*area*Rsg_p[0][5])
                sim.memb.LIST(Rsgchan_B).Count = round(Rsg_ro*area*Rsg_p[0][6])
                sim.memb.LIST(Rsgchan_I6).Count = round(Rsg_ro*area*Rsg_p[1][5])

                for s in CaPchan[...]:
                    sim.memb.LIST(s).Count = round(CaP_fac*CaP_ro*area*CaP_p[s.Count(CaPo)])
                
                for s in CaTchan[...]:
                    hCnt, mCnt = s.Count(CaTho), s.Count(CaTmo)
                    sim.memb.LIST(s).Count = round(CaT_ro*area*CaT_p[hCnt][mCnt])

                for s in BKchan[...]:
                    isOpen, nbCa = s.Count(BKopen), s.Count(BKCa)
                    sim.memb.LIST(s).Count = round(BK_fac*BK_ro*area*BK_p[isOpen][nbCa])

                for s in Ihchan[...]:
                    sim.memb.LIST(s).Count = round(Ih_ro*area*Ih_p[s.Count(Iho)])
                
                sim.memb.L[Leak].Count = round(L_ro * area)
                
                sim.cyto.Ca.Conc = Ca_iconc
                sim.cyto.Mg.Conc = Mg_conc
        
                
                sim.cyto.CB[CBs,   CBf,   CBimmob].Conc = iCBsf_conc
                sim.cyto.CB[CBsCa, CBf,   CBimmob].Conc = iCBCaf_conc
                sim.cyto.CB[CBs,   CBfCa, CBimmob].Conc = iCBsCa_conc
                sim.cyto.CB[CBsCa, CBfCa, CBimmob].Conc = iCBCaCa_conc
        
                sim.cyto.CB[CBs,   CBf,   CBmob].Conc = CBsf_conc
                sim.cyto.CB[CBsCa, CBf,   CBmob].Conc = CBCaf_conc
                sim.cyto.CB[CBs,   CBfCa, CBmob].Conc = CBsCa_conc
                sim.cyto.CB[CBsCa, CBfCa, CBmob].Conc = CBCaCa_conc
        
                sim.cyto.PV.Conc = PV_conc
                sim.cyto.PVCa.Conc = PVCa_conc
                sim.cyto.PVMg.Conc = PVMg_conc

                sim.cyto.InfluxPump.Conc = InfluxPump_compconc_soma
                
                sim.EfieldDT = EF_DT
        
                sim.ALL(Membrane).Potential = init_pot
                sim.membrane.VolRes = Ra
                sim.membrane.Capac = memb_capac
        
                # Set temperature for ghk reactions
                sim.Temp = TEMPERATURE
            
                btime=time.time()
                
                for j in range(1001):
                    t = ENDT * j / 1000
                    if MPI.rank == 0:
                        print(f'run {i}: CaP {CaP_fac}: BK {BK_fac}: {t} / {ENDT}s', time.time()-btime, flush=True)
                    sim.run(round(t,6))
