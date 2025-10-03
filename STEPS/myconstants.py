# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import steps.interface

from steps.model import *
from steps.utils import *

import math

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###########################################################
# Model Parameters
###########################################################

TEMPERATURE = 34.0 + 273.15
Q10 = 3

FARADAY = 96485.3365     # C/mol
R = 8.3144621            # J/mol K
AVOGADRO = 6.02214129e23 # /mol

Qt = math.pow(Q10, (TEMPERATURE - (23 + 273.15)) / 10)
Qt_kv3 = math.pow(2.7, (TEMPERATURE - (22 + 273.15)) / 10)
Qt_nap = math.pow(2.7, (TEMPERATURE - (22 + 273.15)) / 10)
Qt_rsg = math.pow(2.7, (TEMPERATURE - (22 + 273.15)) / 10)
Qt_CaP = math.pow(Q10, (TEMPERATURE - (22 + 273.15)) / 10)
Qt_CaT = math.pow(1.0913, (TEMPERATURE - (32 + 273.15)) / 10)
Qt_mslo = math.pow(2.0, (TEMPERATURE - (25 + 273.15)) / 10)
Qt_sk = math.pow(2.7, (TEMPERATURE - (23 + 273.15)) / 10)
Qt_ih = math.pow(Q10, (TEMPERATURE - (22 + 273.15)) / 10)

#######################################
# Membrane Parameters
#######################################
capac_fac = 50

init_pot = Parameter(-60, 'mV', Description='Initial membrane potential')
Ra = Parameter(120.0*1.0e-2, 'ohm m', Description='Bulk resistivity')
memb_capac = Parameter(capac_fac*0.64e-2, 'F m^-2', Description='Membrane capacitance')

#######################################
# Kv3 channels parameters
#######################################

Kv3_G = Parameter(16, 'pS', Description='Kv3 single channel conductance')
gbar_Kv3 = Parameter((0.3*6*2*2*2*4*2/20)*1e4, 'pS um^-2', Description='Kv3 total conductance')
#gbar_Kv3 = Parameter((0.3*6*2*2*2*4*2)*1e4, 'pS um^-2', Description='Kv3 total conductance')
Kv3_ro = Parameter(gbar_Kv3/Kv3_G, 'um^-2', Description='Kv3 channels density')
Kv3_rev = Parameter(-77, 'mV', Description='Kv3 channel reversal potential')

# Reaction rates

def alphaFktkv3(Vm): #(1/ms)->(1/s)
    ca = 0.22e3 #(1/ms)->(1/s)
    cva = 16 #(mV)
    vshift_kv3 = 4.0
    cka = -26.5 #(mV)
    return ca*math.exp(-(Vm+cva+vshift_kv3)/cka)

def betaFktkv3(Vm): #(1/ms)->(1/s)
    cb = 0.22e3 #(1/ms)->(1/s)
    cvb = 16 #(mV)
    vshift_kv3 = 4.0
    ckb = 26.5 #(mV)
    return cb*math.exp(-(Vm+cvb+vshift_kv3)/ckb)

a_n = VDepRate(lambda V: Qt_kv3 * alphaFktkv3(V*1e3))
b_n = VDepRate(lambda V: Qt_kv3 * betaFktkv3(V*1e3))

# Initial conditions
Kv3_p = [0.8262961138601734, 0.16147973214386666, 0.011834000905909248, 0.00038544521437743485, 4.70787567326776e-6]

#######################################
# NaP channels parameters
#######################################

NaP_G = Parameter(0.37, 'pS', Description='NaP single channel conductance')
gbar_NaP = Parameter(0.0023 * 1e4, 'pS um^-2', Description='NaP total conductance')
#NaP_ro = Parameter(0.6116216216216216216e9, 'um^-2', Description='NaP channels density')
NaP_ro = Parameter(gbar_NaP/NaP_G, 'um^-2', Description='NaP channels density')
NaP_rev = Parameter(63, 'mV', Description='NaP channel reversal potential')

# Reaction rates

def alphamnap(Vm): #(1/ms)->(1/s), Vm (mV)
    Aamp = 17.235e3 #(1/ms)->(1/s) : A for alpha m persis
    Bamp = 27.58 #(mV)
    Camp = -11.47 #(mV)
    return Aamp/(1.0 + math.exp((Vm + Bamp)/Camp))

def betamnap(Vm): #(1/ms)->(1/s), Vm (mV)
    Abmp = 17.235e3 #(1/ms)->(1/s) : A for beta m persis
    Bbmp = 86.2 #(mV)
    Cbmp = 19.8 #(mV)
    return Abmp/(1.0 + math.exp((Vm + Bbmp)/Cbmp))

a_m = VDepRate(lambda V: 1 / (1 + math.exp(-(V*1e3 + 66) / 5)) / (1/(alphamnap(V*1e3) + betamnap(V*1e3))) * Qt_nap)
b_m = VDepRate(lambda V: (1 - 1 / (1 + math.exp(-(V*1e3 + 66) / 5))) / (1/(alphamnap(V*1e3) + betamnap(V*1e3))) * Qt_nap)

# Initial conditions
NaP_p = [0.01240262169113613, 0.1235344624891214, 0.41014885945129553, 0.4539140563684469]

#######################################
# Rsg sodium channels parameters
#######################################

Rsg_G = Parameter(0.764, 'pS', Description='Rsg single channel conductance')
gbar_Rsg = Parameter(0.8 * 0.7 * 1e4, 'pS um^-2', Description='Rsg total conductance')
#Rsg_ro = Parameter(1.96335e9, 'um^-2', Description='Rsg channels density')
Rsg_ro = Parameter(gbar_Rsg/Rsg_G, 'um^-2', Description='Rsg channels density')
Rsg_rev = Parameter(63, 'mV', Description='Rsg channel reversal potential')

# Reaction rates

#Units (M)
alpha_rsg = 150.0 #(*/ms*)(*commented in Neuron file: activation*)
beta_rsg = 3. #(*/ms*)(*commented in Neuron file: deactivation*)
gamma_rsg = 150.0 #(*/ms*)(*commented in Neuron file: opening*)
delta_rsg = 40. #(*/ms*)(*commented in Neuron file: closing,greater than BEAN/KUO=0.2*)
epsilon_rsg = 1.75 #(*/ms*)(*open->Iplus for tau=0.3 ms at+30 with x5*)
zeta_rsg = 0.03 #(*/ms*)(*commented in Neuron file: Iplus->open for tau=25 ms at-30 with x6*)
x1_rsg = 20.0 #(*mV*)(*commented in Neuron file: Vdep of activation (alpha)*)
x2_rsg = -20. #(*mV*) (*commented in Neuron file: Vdep of deactivation (beta)*)
x3_rsg = 1e12 #(*mV*)(*commented in Neuron file: Vdep of activation (gamma)*)
x4_rsg = -1e12 #(*mV*) (*commented in Neuron file: Vdep of deactivation (delta)*)
x5_rsg = 1e12 #(*mV*)(*commented in Neuron file: Vdep of activation (epsilon)*)
x6_rsg = -25 #(*mV*)(*commented in Neuron file: Vdep out of Ipos (zeta)*)

vshifta_rsg_axon = 15. #(*mV*) (*commented in Neuron file: steady state activation curve*)
vshiftk_rsg_axon = 5.
vshifti_rsg_axon = -5. #(*mV*)(*commented in Neuron file: steady state inactivation curve*)

Oon_rsg = 0.75 #(*/ms*)(*commented in Neuron file: open->Ineg transition*)
Con_rsg = 0.005 #(*/ms*)(*commented in Neuron file: closed->inactivated transitions*)
alfac_rsg = math.pow((Oon_rsg/Con_rsg),(1.0/4.0))

Ooff_rsg = 0.005 #(*/ms*)(*commented in Neuron file: Ineg->open transition,set to 0 to eliminate the persistent*)
Coff_rsg = 0.5 #(*/ms*)(*commented in Neuron file: inactivated->closed transitions*)
btfac_rsg = math.pow((Ooff_rsg/Coff_rsg),(1.0/4.0))

#(*forward along C*)
rsgC_f = VDepRate.Create(lambda V: alpha_rsg*math.exp(((V*1e3)/(x1_rsg)))*Qt_rsg*1e3)
rsgCO_f = gamma_rsg*Qt_rsg*1e3
rsgOB_f = epsilon_rsg*Qt_rsg*1e3 #(*commented in Neuron file: *exp(v/x5)*)

#(*backward along C*)
rsgC_b = VDepRate.Create(lambda V: beta_rsg*math.exp((((V*1e3) + vshifta_rsg_axon)/(x2_rsg + vshiftk_rsg_axon)))*Qt_rsg*1e3)
rsgCO_b = delta_rsg*Qt_rsg*1e3
rsgOB_b = VDepRate.Create(lambda V: zeta_rsg*math.exp(((V*1e3)/x6_rsg))*Qt_rsg*1e3)

#(*forward along I*)
rsgI_f = VDepRate.Create(lambda V: alpha_rsg*alfac_rsg*math.exp((((V*1e3) + vshifti_rsg_axon)/x1_rsg))*Qt_rsg*1e3)
rsgIn_f = gamma_rsg*Qt_rsg*1e3 #(*f1n*)(*commented in Neuron file: *exp(v/x3) dito*)

#(*backward along I*)
rsgI_b = VDepRate.Create(lambda V: beta_rsg*btfac_rsg*math.exp((((V*1e3) + vshifti_rsg_axon)/x2_rsg))*Qt_rsg*1e3)
rsgIn_b = delta_rsg*Qt_rsg*1e3 #(*1n*)(*commented in Neuron file: *exp(v/x4) dito*)

#(*forward along CI*)
rsgCI_f = Con_rsg*Qt_rsg*1e3
rsgCIn_f = Oon_rsg*Qt_rsg*1e3

#(*backward along CI*)
rsgCI_b = Coff_rsg*Qt_rsg*1e3
rsgCIn_b = Ooff_rsg*Qt_rsg*1e3

# Initial conditions
Rsg_p = [
    [0.5484840504151031,0.2718795626236479,0.05052711033631685,0.004169841737733529,0.0001277736943237481,0.000472392875100075,0.002499846695358826],
    [0.009298378119390719,0.031019754518843094,0.038750650374045514,0.02150371122375334,0.0044758250332146135,0.016791102353168746]
]
'''
Rsg_p = [
    [C1,C2,C3,C4,C5,O,B],
    [I1,I2,I3,I4,I5,I6]
]
'''

#######################################
# CaP channels parameters
#######################################

CaP_P = Parameter(2.5e-2, 'um^3 s^-1', Description='CaP single channel permeability')
pbar_CaP = Parameter(0.95e-4 * 24 * 1e4, 'um s^-1', Description='CaP total permeability')
#CaP_ro = Parameter(38, 'um^-2', Description='CaP channels density')
CaP_ro = Parameter(pbar_CaP/CaP_P, 'um^-2', Description='CaP channels density')

# Reaction rates

def minf_cap(Vm):
    cv = 30.5 #(mV)
    ck = 4.113 #(mV)
    vshift_cap = -5.0
    return 1/(1 + math.exp(-(Vm+cv+vshift_cap)/ck))

def tau_cap(Vm):
    vshift_cap = -5.0
    return (0.2 + 0.7031 * math.exp(-(((Vm + 30.0 + vshift_cap)/14) ** 2)))/Qt_CaP

alpha_cap = VDepRate.Create(lambda V: (minf_cap(V * 1e3) / tau_cap(V * 1e3)) *1e3)
beta_cap = VDepRate.Create(lambda V: (1 - minf_cap(V * 1e3)) / tau_cap(V * 1e3) *1e3)

# Initial conditions
#CaP_m0_p = 1-minf_cap(init_pot)
#CaP_m1_p = minf_cap(init_pot)
CaP_p = [0.9997724784220603, 0.0002275215779397552]

#######################################
# CaT channels parameters
#######################################

CaT_P = Parameter(1.65e-5, 'um^3 s^-1', Description='CaT single channel permeability')
pbar_CaT = Parameter((6.4e-5) * 2 * 1e4, 'um s^-1', Description='CaT total permeability')
#CaT_ro = Parameter(3.7576, 'um^-2', Description='CaT channels density')
CaT_ro = Parameter(pbar_CaT/CaT_P, 'um^-2', Description='CaT channels density')

# Reaction rates

def minf_cat(Vm):
    vhalfm = -42.206
    vshift_cat = 0
    cvm = -4.7056
    return math.pow(1 / (1 + math.exp((Vm - vhalfm) / cvm)), 1 / 3)

def taum_cat(Vm):
    C_tau_m = 1.2757
    A_tau_m = -2.3199
    v0_tau_m1 = -48.048
    vshift_cat = 0.0
    k_tau_m1 = 30.655
    B_tau_m = 2.5712
    v0_tau_m2 = -28.386
    k_tau_m2 = 9.6306
    return 1/(C_tau_m+A_tau_m/(1+math.exp((v0_tau_m1-Vm-vshift_cat)/k_tau_m1))+B_tau_m/(1+math.exp((v0_tau_m2-Vm-vshift_cat)/k_tau_m2)))

def hinf_cat(Vm):
    v0_h_inf = -75.118
    vshift_cat = 0.0
    k_h_inf = 6.4635
    return 1/(1+math.exp((Vm-v0_h_inf-vshift_cat)/k_h_inf))

def tauh_cat(Vm):
    C_tau_h = 0.0076
    A_tau_h = 0.17746
    v0_tau_h1 = -58.535
    vshift_cat = 0.0
    k_tau_h1 = 6.2692
    B_tau_h = 0.13402
    v0_tau_h2= -101.436
    k_tau_h2 = -5.5845
    return 1/(C_tau_h+A_tau_h/(1+math.exp((v0_tau_h1-Vm-vshift_cat)/k_tau_h1))+B_tau_h/(1+math.exp((v0_tau_h2-Vm-vshift_cat)/k_tau_h2)))

alpham_cat = VDepRate.Create(lambda V: minf_cat(V * 1e3) / taum_cat(V * 1e3) * Qt_CaT * 1e3)
betam_cat = VDepRate.Create(lambda V: (1 - minf_cat(V * 1e3)) / taum_cat(V * 1e3) * Qt_CaT * 1e3)

alphah_cat = VDepRate.Create(lambda V: hinf_cat(V * 1e3) / tauh_cat(V * 1e3) * Qt_CaT * 1e3)
betah_cat = VDepRate.Create(lambda V: (1 - hinf_cat(V * 1e3)) / tauh_cat(V * 1e3) * Qt_CaT * 1e3)

# Initial conditions
CaT_p = [
    [0.33844731242243303,   0.397593135955208,   0.1556916896234537, 0.020322200448766186],   # h0
    [0.03263508711711888, 0.03833827645489287, 0.015012711485003546, 0.001959586493123762], # h1
]

#######################################
# BK(slow) channels parameters
#######################################

BKslw_G = Parameter(210, 'pS', Description='BKslw single channel conductance')
gbar_BKslw = Parameter(1.5 * 0.7 * 1e4, 'pS um^-2', Description='BKslw total conductance')
#BKslw_ro = Parameter(0.4761904761904761905, 'um^-2', Description='BKslw channels density')
BKslw_ro = Parameter(gbar_BKslw/BKslw_G, 'um^-2', Description='BKslw channels density')
BKslw_rev = Parameter(-77, 'mV', Description='BKslw channel reversal potential')

# Reaction rates

#Units (1)
Qo_slw = 1.0264 # fitted
Qc_slw = -0.0016661 # fitted

#Units (/s)
pf0_slw = 0.00095197*1e3 # fitted
pf1_slw = 0.012711*1e3 # fitted
pf2_slw = 0.15735*1e3 # fitted
pf3_slw = 1.2613*1e3 # fitted
pf4_slw = 0.077161*1e3 # fitted

pb0_slw = 11.837*1e3 # fitted
pb1_slw = 0.053557*1e3 # fitted
pb2_slw = 0.0050029*1e3 # fitted
pb3_slw = 0.22836*1e3 # fitted
pb4_slw = 0.00012541*1e3 # fitted

#Units(/M)
k1_slw = 1.0e6 #ok

#Units(/s)
onoffrate_slw = 1.0e3 #ok

L0_slw = 1806 #ok

#Units (M)
Kc_slw = 0.018748*1e-3 # fitted
Ko_slw = 0.0001028*1e-3 # fitted

BKslw_f = k1_slw*onoffrate_slw*Qt_mslo
BKslwo_b = Ko_slw*k1_slw*onoffrate_slw*Qt_mslo
BKslwc_b = Kc_slw*k1_slw*onoffrate_slw*Qt_mslo

BKslw_f0 = VDepRate.Create(
    lambda V: pf0_slw * Qt_mslo * (math.exp((Qo_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_f1 = VDepRate.Create(
    lambda V: pf1_slw * Qt_mslo * (math.exp((Qo_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_f2 = VDepRate.Create(
    lambda V: pf2_slw * Qt_mslo * (math.exp((Qo_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_f3 = VDepRate.Create(
    lambda V: pf3_slw * Qt_mslo * (math.exp((Qo_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_f4 = VDepRate.Create(
    lambda V: pf4_slw * Qt_mslo * (math.exp((Qo_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_oc_f = [BKslw_f0, BKslw_f1, BKslw_f2, BKslw_f3, BKslw_f4]

BKslw_b0 = VDepRate.Create(
    lambda V: pb0_slw * Qt_mslo * (math.exp((Qc_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_b1 = VDepRate.Create(
    lambda V: pb1_slw * Qt_mslo * (math.exp((Qc_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_b2 = VDepRate.Create(
    lambda V: pb2_slw * Qt_mslo * (math.exp((Qc_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_b3 = VDepRate.Create(
    lambda V: pb3_slw * Qt_mslo * (math.exp((Qc_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_b4 = VDepRate.Create(
    lambda V: pb4_slw * Qt_mslo * (math.exp((Qc_slw * FARADAY * V) / (R * TEMPERATURE)))
)
BKslw_oc_b = [BKslw_b0, BKslw_b1, BKslw_b2, BKslw_b3, BKslw_b4]

# Initial conditions
BKslw_p = [
    [0.990358,      0.0095082,  0.0000342597, 7.64949e-8, 4.68805e-11],
    [6.06811e-6, 0.0000558796,  0.0000318605, 5.34511e-6,  5.84768e-7],
]

#######################################
# BK(mslo) channels parameters
#######################################
BK_cond = 250
BK_G = Parameter(BK_cond, 'pS', Description='BK single channel conductance')
#BK_G = Parameter(240, 'pS', Description='BK single channel conductance')
gbar_BK = Parameter(3 * 2 * 1e4, 'pS um^-2', Description='BK total conductance')
#BK_ro = Parameter(0.4761904761904761905, 'um^-2', Description='BK channels density')
BK_ro = Parameter(gbar_BK/BK_G, 'um^-2', Description='BK channels density')
BK_rev = Parameter(-77, 'mV', Description='BK channel reversal potential')

# Reaction rates

#Units (1)
Qo = 0.73
Qc = -0.58

#Units (/s)
pf0 = 5.5
pf1 = 8
pf2 = 2
pf3 = 884
pf4 = 900

pb0 = 8669
pb1 = 1127
pb2 = 25.2
pb3 = 1013
pb4 = 125.7

#Units(/M.s)
k1 = 1.0e9

#Units(/s)
kc_off = 11917
ko_off = 1065

L0 = 1576

BK_f = k1*Qt_mslo
BKo_b = ko_off*Qt_mslo
BKc_b = kc_off*Qt_mslo



BK_f0 = VDepRate.Create(
    lambda V: pf0 * Qt_mslo * (math.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_f1 = VDepRate.Create(
    lambda V: pf1 * Qt_mslo * (math.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_f2 = VDepRate.Create(
    lambda V: pf2 * Qt_mslo * (math.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_f3 = VDepRate.Create(
    lambda V: pf3 * Qt_mslo * (math.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_f4 = VDepRate.Create(
    lambda V: pf4 * Qt_mslo * (math.exp((Qo * FARADAY * V) / (R * TEMPERATURE)))
)
BK_oc_f = [BK_f0, BK_f1, BK_f2, BK_f3, BK_f4]

BK_b0 = VDepRate.Create(
    lambda V: pb0 * Qt_mslo * (math.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_b1 = VDepRate.Create(
    lambda V: pb1 * Qt_mslo * (math.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_b2 = VDepRate.Create(
    lambda V: pb2 * Qt_mslo * (math.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_b3 = VDepRate.Create(
    lambda V: pb3 * Qt_mslo * (math.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_b4 = VDepRate.Create(
    lambda V: pb4 * Qt_mslo * (math.exp((Qc * FARADAY * V) / (R * TEMPERATURE)))
)
BK_oc_b = [BK_b0, BK_b1, BK_b2, BK_b3, BK_b4]

# Initial conditions
BK_p = [
    [0.979391,     0.0204276, 0.000159776,  5.5542e-7, 7.24043e-10],
    [0.000016302, 4.47363e-6,  4.59886e-7, 2.10196e-8, 3.60282e-10],
]

#######################################
# SK channels parameters
#######################################

SK_fac = 0.0
SK_G = Parameter(10, 'pS', Description='SK single channel conductance')
gbar_SK = Parameter((SK_fac*0.02) / 1.8 * 2.5  * 1e4, 'pS um^-2', Description='SK total conductance')
#gbar_SK = Parameter(0.02 / 1.8 * 2.5 * 1.5 * 1e4, 'pS um^-2', Description='SK total conductance')
#SK_ro = Parameter(0.31, 'um^-2', Description='SK channels density')
SK_ro = Parameter(gbar_SK/SK_G, 'um^-2', Description='SK channels density')
SK_rev = Parameter(-77, 'mV', Description='SK channel reversal potential')

# Reaction rates
diff_sk = 1.0
scal_sk = 1.0

#Units (/s)
invc1 = 80
invc2 = 80
invc3 = 200

invo1 = 1000
invo2 = 100

diro1 = 160
diro2 = 1200

#Units ( /s M)

dirc2 = 200e6/diff_sk
dirc3 = 160e6/diff_sk
dirc4 = 80e6/diff_sk

invc1_t = invc1*Qt_sk*scal_sk
invc2_t = invc2*Qt_sk*scal_sk
invc3_t = invc3*Qt_sk*scal_sk

invo1_t = invo1*Qt_sk*scal_sk
invo2_t = invo2*Qt_sk*scal_sk

diro1_t = diro1*Qt_sk*scal_sk
diro2_t = diro2*Qt_sk*scal_sk

dirc2_t = dirc2*Qt_sk*scal_sk
dirc3_t = dirc3*Qt_sk*scal_sk
dirc4_t = dirc4*Qt_sk*scal_sk

# Intital conditions
SK_C1_p= 0.887615
SK_C2_p= 0.0998567
SK_C3_p= 0.0089871
SK_C4_p= 0.000161768

SK_O1_p= 0.00143794
SK_O2_p= 0.00194121

#######################################
# Ih channels parameters
#######################################

Ih_G = Parameter(9.7e-4, 'pS', Description='Ih single channel conductance')
gbar_Ih = Parameter(3.6e-4 * 0.3 * 1e4, 'pS um^-2', Description='Ih total conductance')
#Ih_ro = Parameter(0.10309278350515463918, 'um^-2', Description='Ih channels density')
Ih_ro = Parameter(gbar_Ih/Ih_G, 'um^-2', Description='Ih channels density')
Ih_rev = Parameter(-30, 'mV', Description='Ih channel reversal potential')

# Reaction rates
def ninf_ih(Vm): #(mV)
    return (1.0/(1.0+math.exp((Vm+90.3+3.0+3.0)/9.67)))

def taun_ih(Vm): #(ms)
    return (1000.0/(0.62*(math.exp((Vm+68.0)/-22.0)+math.exp((Vm+68.0)/7.14)))/Qt_ih/1.3)

alphan_ih = VDepRate.Create(lambda V: ninf_ih(V*1e3)/taun_ih(V*1e3) * 1e3)
betan_ih = VDepRate.Create(lambda V: (1.0-ninf_ih(V*1e3))/taun_ih(V*1e3) * 1e3)

# Intital conditions
Ih_p = [0.9771095274411775, 0.02289047255882251]

#######################################
# Leak channels parameters
#######################################

L_G = Parameter(0.04e-2, 'pS', Description='Leak single channel conductance')
#gbar_L = Parameter(1 / 120236 * 1e4, 'pS um^-2', Description='Leak total conductance')
gbar_L = Parameter(0.001 * 1e4, 'pS um^-2', Description='Leak total conductance')
#L_ro = Parameter(12.5, 'um^-2', Description='Leak channels density')
L_ro = Parameter(gbar_L/L_G, 'um^-2', Description='Leak channels density')
L_rev = Parameter(-65, 'mV', Description='Leak channel reversal potential')

#######################################
# Ca pump channels parameters
#######################################

P_ro = Parameter(6.022141, 'um^-2', Description='Ca2+ pump density')

# Reaction rates
P_f = 3e9
P_b = 1.75e4
P_k = 7.255e4

#######################################
# Ca Influx pump channels parameters
#######################################

InfluxPump_compconc_soma = 0.00559768876562243*math.pow(10.0,-3.0)

Inf_f = 3077.65379846454*math.pow(10.0,6.0) 
Inf_b = 1.0*1.0e-15*math.pow(10.0,3.0)
Inf_k = 17.7595788877003*math.pow(10.0,3.0)
#Inf_k = 17.7595788877003*math.pow(10.0,6.0)

#######################################
# Calcium buffering parameters
#######################################

# Ca concentrations

Ca_oconc = 2e-3
Ca_iconc = 45e-9

# Mg concentrations

Mg_conc = 590e-6

# Diffusion constants

DCST = 0.233e-9 #0.223e-9 # Ca
DCB = 0.028e-9  # Calbindin (CB)
DPV = 0.043e-9  # Parvalbumin (PV)

# Reaction rates

CBf_f_kcst = 4.35e7 #nf1
CBf_b_kcst = 35.8 #nf2

CBs_f_kcst = 0.55e7 #ns1
CBs_b_kcst = 2.6 #ns2

PVca_f = 10.7e7 #m1
PVca_b = 0.95 #m2

PVmg_f_kcst = 0.8e6 #p1
PVmg_b_kcst = 25 #p2

# Buffer concentrations

mobile_ratio = 0.8
imobile_ratio = (1.0-mobile_ratio)

CBnull = 0.08e-3 #(M) #CBnull = 0.08 (mM)
PVnull = 0.04e-3 #(M) #PVnull = 0.04 (mM)

def kdf(cainull,nf1,nf2): #(1)
    return ((cainull*nf1)/nf2)

def kds(cainull,ns1,ns2): #(1)
    return ((cainull*ns1)/ns2)

def ssCB(cainull,nf1,nf2,ns1,ns2,CBnull): #(mM) -> (M)
    return (CBnull/(1.0+kdf(cainull,nf1,nf2)+kds(cainull,ns1,ns2)+(kdf(cainull,nf1,nf2)*kds(cainull,ns1,ns2))))

def ssCBfast(cainull,nf1,nf2,ns1,ns2,CBnull): #(mM) -> (M)
    return ((CBnull*kds(cainull,ns1,ns2))/(1.0+kdf(cainull,nf1,nf2)+kds(cainull,ns1,ns2)+(kdf(cainull,nf1,nf2)*kds(cainull,ns1,ns2))))

def ssCBslow(cainull,nf1,nf2,ns1,ns2,CBnull): #(mM) -> (M)
    return ((CBnull*kdf(cainull,nf1,nf2))/(1.0+kdf(cainull,nf1,nf2)+kds(cainull,ns1,ns2)+(kdf(cainull,nf1,nf2)*kds(cainull,ns1,ns2))))

def ssCBca(cainull,nf1,nf2,ns1,ns2,CBnull): #(mM) -> (M)
    return ((CBnull*kdf(cainull,nf1,nf2)*kds(cainull,ns1,ns2))/(1.0+kdf(cainull,nf1,nf2)+kds(cainull,ns1,ns2)+(kdf(cainull,nf1,nf2)*kds(cainull,ns1,ns2))))

def kdc(cainull,m1,m2): #(1)
    return ((cainull*m1)/m2)

def kdm(mginull,p1,p2): #(1)
    return ((mginull*p1)/p2)

def ssPV(cainull,m1,m2,mginull,p1,p2,PVnull): #(mM) -> (M)
    return (PVnull/(1.0+kdc(cainull,m1,m2)+kdm(mginull,p1,p2)))

def ssPVca(cainull,m1,m2,mginull,p1,p2,PVnull): #(mM) -> (M)
    return ((PVnull*kdc(cainull,m1,m2))/(1.0+kdc(cainull,m1,m2)+kdm(mginull,p1,p2)))

def ssPVmg(cainull,m1,m2,mginull,p1,p2,PVnull): #(mM) -> (M)
    return ((PVnull*kdm(mginull,p1,p2))/(1.0+kdc(cainull,m1,m2)+kdm(mginull,p1,p2)))

CBsf_conc = mobile_ratio*ssCB(Ca_iconc,CBf_f_kcst,CBf_b_kcst,CBs_f_kcst,CBs_b_kcst,CBnull) #(M) #110.82e-6
CBCaf_conc = mobile_ratio*ssCBfast(Ca_iconc,CBf_f_kcst,CBf_b_kcst,CBs_f_kcst,CBs_b_kcst,CBnull) #(M) #10.549e-6
CBsCa_conc = mobile_ratio*ssCBslow(Ca_iconc,CBf_f_kcst,CBf_b_kcst,CBs_f_kcst,CBs_b_kcst,CBnull) #(M) #6.0595e-6
CBCaCa_conc = mobile_ratio*ssCBca(Ca_iconc,CBf_f_kcst,CBf_b_kcst,CBs_f_kcst,CBs_b_kcst,CBnull) #(M) #0.57682e-6

iCBsf_conc = imobile_ratio*ssCB(Ca_iconc,CBf_f_kcst,CBf_b_kcst,CBs_f_kcst,CBs_b_kcst,CBnull) #(M) #27.704e-6
iCBCaf_conc = imobile_ratio*ssCBfast(Ca_iconc,CBf_f_kcst,CBf_b_kcst,CBs_f_kcst,CBs_b_kcst,CBnull) #(M) #2.6372e-6
iCBsCa_conc = imobile_ratio*ssCBslow(Ca_iconc,CBf_f_kcst,CBf_b_kcst,CBs_f_kcst,CBs_b_kcst,CBnull) #(M) #1.5148e-6
iCBCaCa_conc = imobile_ratio*ssCBca(Ca_iconc,CBf_f_kcst,CBf_b_kcst,CBs_f_kcst,CBs_b_kcst,CBnull) #(M) #0.14420e-6

PV_conc = ssPV(Ca_iconc,PVca_f,PVca_b,Mg_conc,PVmg_f_kcst,PVmg_b_kcst,PVnull) #(M) #3.2066e-6
PVCa_conc = ssPVca(Ca_iconc,PVca_f,PVca_b,Mg_conc,PVmg_f_kcst,PVmg_b_kcst,PVnull) #(M) #16.252e-6
PVMg_conc = ssPVmg(Ca_iconc,PVca_f,PVca_b,Mg_conc,PVmg_f_kcst,PVmg_b_kcst,PVnull) #(M) #60.541e-6
