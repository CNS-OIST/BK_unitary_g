/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__cdpAIS
#define _nrn_initial _nrn_initial__cdpAIS
#define nrn_cur _nrn_cur__cdpAIS
#define _nrn_current _nrn_current__cdpAIS
#define nrn_jacob _nrn_jacob__cdpAIS
#define nrn_state _nrn_state__cdpAIS
#define _net_receive _net_receive__cdpAIS 
#define factors factors__cdpAIS 
#define rates rates__cdpAIS 
#define state state__cdpAIS 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg(int);
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define ica_pmp _p[0]
#define ica_pmp_columnindex 0
#define ca1 _p[1]
#define ca1_columnindex 1
#define ca2 _p[2]
#define ca2_columnindex 2
#define ca3 _p[3]
#define ca3_columnindex 3
#define ca4 _p[4]
#define ca4_columnindex 4
#define ca5 _p[5]
#define ca5_columnindex 5
#define ca9 _p[6]
#define ca9_columnindex 6
#define ca (_p + 7)
#define ca_columnindex 7
#define mg (_p + 17)
#define mg_columnindex 17
#define CB (_p + 27)
#define CB_columnindex 27
#define CB_f_ca (_p + 37)
#define CB_f_ca_columnindex 37
#define CB_ca_s (_p + 47)
#define CB_ca_s_columnindex 47
#define CB_ca_ca (_p + 57)
#define CB_ca_ca_columnindex 57
#define iCB (_p + 67)
#define iCB_columnindex 67
#define iCB_f_ca (_p + 77)
#define iCB_f_ca_columnindex 77
#define iCB_ca_s (_p + 87)
#define iCB_ca_s_columnindex 87
#define iCB_ca_ca (_p + 97)
#define iCB_ca_ca_columnindex 97
#define PV (_p + 107)
#define PV_columnindex 107
#define PV_ca (_p + 117)
#define PV_ca_columnindex 117
#define PV_mg (_p + 127)
#define PV_mg_columnindex 127
#define pump _p[137]
#define pump_columnindex 137
#define pumpca _p[138]
#define pumpca_columnindex 138
#define C0 _p[139]
#define C0_columnindex 139
#define C1 _p[140]
#define C1_columnindex 140
#define C2 _p[141]
#define C2_columnindex 141
#define C3 _p[142]
#define C3_columnindex 142
#define C4 _p[143]
#define C4_columnindex 143
#define O0 _p[144]
#define O0_columnindex 144
#define O1 _p[145]
#define O1_columnindex 145
#define O2 _p[146]
#define O2_columnindex 146
#define O3 _p[147]
#define O3_columnindex 147
#define O4 _p[148]
#define O4_columnindex 148
#define BK_ro _p[149]
#define BK_ro_columnindex 149
#define BK_conc _p[150]
#define BK_conc_columnindex 150
#define dr_two _p[151]
#define dr_two_columnindex 151
#define rad_outer _p[152]
#define rad_outer_columnindex 152
#define rad_inner _p[153]
#define rad_inner_columnindex 153
#define vol_shell _p[154]
#define vol_shell_columnindex 154
#define ica _p[155]
#define ica_columnindex 155
#define ibg _p[156]
#define ibg_columnindex 156
#define parea _p[157]
#define parea_columnindex 157
#define cai _p[158]
#define cai_columnindex 158
#define mgi _p[159]
#define mgi_columnindex 159
#define c01 _p[160]
#define c01_columnindex 160
#define c12 _p[161]
#define c12_columnindex 161
#define c23 _p[162]
#define c23_columnindex 162
#define c34 _p[163]
#define c34_columnindex 163
#define o01 _p[164]
#define o01_columnindex 164
#define o12 _p[165]
#define o12_columnindex 165
#define o23 _p[166]
#define o23_columnindex 166
#define o34 _p[167]
#define o34_columnindex 167
#define f0 _p[168]
#define f0_columnindex 168
#define f1 _p[169]
#define f1_columnindex 169
#define f2 _p[170]
#define f2_columnindex 170
#define f3 _p[171]
#define f3_columnindex 171
#define f4 _p[172]
#define f4_columnindex 172
#define c10 _p[173]
#define c10_columnindex 173
#define c21 _p[174]
#define c21_columnindex 174
#define c32 _p[175]
#define c32_columnindex 175
#define c43 _p[176]
#define c43_columnindex 176
#define o10 _p[177]
#define o10_columnindex 177
#define o21 _p[178]
#define o21_columnindex 178
#define o32 _p[179]
#define o32_columnindex 179
#define o43 _p[180]
#define o43_columnindex 180
#define b0 _p[181]
#define b0_columnindex 181
#define b1 _p[182]
#define b1_columnindex 182
#define b2 _p[183]
#define b2_columnindex 183
#define b3 _p[184]
#define b3_columnindex 184
#define b4 _p[185]
#define b4_columnindex 185
#define Dca (_p + 186)
#define Dca_columnindex 186
#define Dmg (_p + 196)
#define Dmg_columnindex 196
#define DCB (_p + 206)
#define DCB_columnindex 206
#define DCB_f_ca (_p + 216)
#define DCB_f_ca_columnindex 216
#define DCB_ca_s (_p + 226)
#define DCB_ca_s_columnindex 226
#define DCB_ca_ca (_p + 236)
#define DCB_ca_ca_columnindex 236
#define DiCB (_p + 246)
#define DiCB_columnindex 246
#define DiCB_f_ca (_p + 256)
#define DiCB_f_ca_columnindex 256
#define DiCB_ca_s (_p + 266)
#define DiCB_ca_s_columnindex 266
#define DiCB_ca_ca (_p + 276)
#define DiCB_ca_ca_columnindex 276
#define DPV (_p + 286)
#define DPV_columnindex 286
#define DPV_ca (_p + 296)
#define DPV_ca_columnindex 296
#define DPV_mg (_p + 306)
#define DPV_mg_columnindex 306
#define Dpump _p[316]
#define Dpump_columnindex 316
#define Dpumpca _p[317]
#define Dpumpca_columnindex 317
#define DC0 _p[318]
#define DC0_columnindex 318
#define DC1 _p[319]
#define DC1_columnindex 319
#define DC2 _p[320]
#define DC2_columnindex 320
#define DC3 _p[321]
#define DC3_columnindex 321
#define DC4 _p[322]
#define DC4_columnindex 322
#define DO0 _p[323]
#define DO0_columnindex 323
#define DO1 _p[324]
#define DO1_columnindex 324
#define DO2 _p[325]
#define DO2_columnindex 325
#define DO3 _p[326]
#define DO3_columnindex 326
#define DO4 _p[327]
#define DO4_columnindex 327
#define DBK_ro _p[328]
#define DBK_ro_columnindex 328
#define DBK_conc _p[329]
#define DBK_conc_columnindex 329
#define Ddr_two _p[330]
#define Ddr_two_columnindex 330
#define Drad_outer _p[331]
#define Drad_outer_columnindex 331
#define Drad_inner _p[332]
#define Drad_inner_columnindex 332
#define Dvol_shell _p[333]
#define Dvol_shell_columnindex 333
#define _g _p[334]
#define _g_columnindex 334
#define _ion_cao	*_ppvar[0]._pval
#define _ion_cai	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _style_ca	*((int*)_ppvar[3]._pvoid)
#define diam	*_ppvar[4]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_factors(void);
 static void _hoc_kdm(void);
 static void _hoc_kdc(void);
 static void _hoc_kds(void);
 static void _hoc_kdf(void);
 static void _hoc_rates(void);
 static void _hoc_ssPVmg(void);
 static void _hoc_ssPVca(void);
 static void _hoc_ssPV(void);
 static void _hoc_ssCBca(void);
 static void _hoc_ssCBslow(void);
 static void _hoc_ssCBfast(void);
 static void _hoc_ssCB(void);
 static void _hoc_u(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_cdpAIS", _hoc_setdata,
 "factors_cdpAIS", _hoc_factors,
 "kdm_cdpAIS", _hoc_kdm,
 "kdc_cdpAIS", _hoc_kdc,
 "kds_cdpAIS", _hoc_kds,
 "kdf_cdpAIS", _hoc_kdf,
 "rates_cdpAIS", _hoc_rates,
 "ssPVmg_cdpAIS", _hoc_ssPVmg,
 "ssPVca_cdpAIS", _hoc_ssPVca,
 "ssPV_cdpAIS", _hoc_ssPV,
 "ssCBca_cdpAIS", _hoc_ssCBca,
 "ssCBslow_cdpAIS", _hoc_ssCBslow,
 "ssCBfast_cdpAIS", _hoc_ssCBfast,
 "ssCB_cdpAIS", _hoc_ssCB,
 "u_cdpAIS", _hoc_u,
 0, 0
};
#define kdm kdm_cdpAIS
#define kdc kdc_cdpAIS
#define kds kds_cdpAIS
#define kdf kdf_cdpAIS
#define ssPVmg ssPVmg_cdpAIS
#define ssPVca ssPVca_cdpAIS
#define ssPV ssPV_cdpAIS
#define ssCBca ssCBca_cdpAIS
#define ssCBslow ssCBslow_cdpAIS
#define ssCBfast ssCBfast_cdpAIS
#define ssCB ssCB_cdpAIS
#define u u_cdpAIS
 extern double kdm( );
 extern double kdc( );
 extern double kds( );
 extern double kdf( );
 extern double ssPVmg( double , double );
 extern double ssPVca( double , double );
 extern double ssPV( double , double );
 extern double ssCBca( double , double );
 extern double ssCBslow( double , double );
 extern double ssCBfast( double , double );
 extern double ssCB( double , double );
 extern double u( double , double );
 /* declare global and static user variables */
#define BK_g BK_g_cdpAIS
 double BK_g = 24;
#define BK_sing_chan_g BK_sing_chan_g_cdpAIS
 double BK_sing_chan_g = 50;
#define CBnull CBnull_cdpAIS
 double CBnull = 0.08;
#define Dpar Dpar_cdpAIS
 double Dpar = 0.043;
#define Dcbd2 Dcbd2_cdpAIS
 double Dcbd2 = 0;
#define Dcbd1 Dcbd1_cdpAIS
 double Dcbd1 = 0.028;
#define Ddmnpe Ddmnpe_cdpAIS
 double Ddmnpe = 0.08;
#define Dbtc Dbtc_cdpAIS
 double Dbtc = 0.007;
#define DCa DCa_cdpAIS
 double DCa = 0.233;
#define Ko Ko_cdpAIS
 double Ko = 0.001065;
#define Kc Kc_cdpAIS
 double Kc = 0.011917;
#define Kp Kp_cdpAIS
 double Kp = 0.0027;
#define L0 L0_cdpAIS
 double L0 = 1576;
#define PVnull PVnull_cdpAIS
 double PVnull = 0.04;
#define Qc Qc_cdpAIS
 double Qc = -0.58;
#define Qo Qo_cdpAIS
 double Qo = 0.73;
#define TotalPump TotalPump_cdpAIS
 double TotalPump = 1e-15;
#define beta beta_cdpAIS
 double beta = 1;
#define cainull cainull_cdpAIS
 double cainull = 4.5e-05;
#define k1 k1_cdpAIS
 double k1 = 1000;
#define kpmp3 kpmp3_cdpAIS
 double kpmp3 = 72.55;
#define kpmp2 kpmp2_cdpAIS
 double kpmp2 = 17.5;
#define kpmp1 kpmp1_cdpAIS
 double kpmp1 = 3000;
#define m2 m2_cdpAIS
 double m2 = 0.00095;
#define m1 m1_cdpAIS
 double m1 = 107;
#define mginull mginull_cdpAIS
 double mginull = 0.59;
#define ns2 ns2_cdpAIS
 double ns2 = 0.0026;
#define ns1 ns1_cdpAIS
 double ns1 = 5.5;
#define nf2 nf2_cdpAIS
 double nf2 = 0.0358;
#define nf1 nf1_cdpAIS
 double nf1 = 43.5;
#define onoffrate onoffrate_cdpAIS
 double onoffrate = 1;
#define pb4 pb4_cdpAIS
 double pb4 = 0.1257;
#define pb3 pb3_cdpAIS
 double pb3 = 1.013;
#define pb2 pb2_cdpAIS
 double pb2 = 0.0252;
#define pb1 pb1_cdpAIS
 double pb1 = 1.127;
#define pb0 pb0_cdpAIS
 double pb0 = 8.669;
#define pf4 pf4_cdpAIS
 double pf4 = 0.9;
#define pf3 pf3_cdpAIS
 double pf3 = 0.884;
#define pf2 pf2_cdpAIS
 double pf2 = 0.002;
#define pf1 pf1_cdpAIS
 double pf1 = 0.008;
#define pf0 pf0_cdpAIS
 double pf0 = 0.0055;
#define p2 p2_cdpAIS
 double p2 = 0.025;
#define p1 p1_cdpAIS
 double p1 = 0.8;
#define vmax vmax_cdpAIS
 double vmax = 0.1;
#define vrat vrat_cdpAIS
 double vrat[10];
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "cainull_cdpAIS", "mM",
 "mginull_cdpAIS", "mM",
 "DCa_cdpAIS", "um2/ms",
 "Dbtc_cdpAIS", "um2/ms",
 "Ddmnpe_cdpAIS", "um2/ms",
 "Dcbd1_cdpAIS", "um2/ms",
 "Dcbd2_cdpAIS", "um2/ms",
 "Dpar_cdpAIS", "um2/ms",
 "CBnull_cdpAIS", "mM",
 "nf1_cdpAIS", "/ms",
 "nf2_cdpAIS", "/ms",
 "ns1_cdpAIS", "/ms",
 "ns2_cdpAIS", "/ms",
 "PVnull_cdpAIS", "mM",
 "m1_cdpAIS", "/ms",
 "m2_cdpAIS", "/ms",
 "p1_cdpAIS", "/ms",
 "p2_cdpAIS", "/ms",
 "kpmp1_cdpAIS", "/mM-ms",
 "kpmp2_cdpAIS", "/ms",
 "kpmp3_cdpAIS", "/ms",
 "beta_cdpAIS", "1",
 "Kp_cdpAIS", "mM",
 "k1_cdpAIS", "/mM",
 "onoffrate_cdpAIS", "/ms",
 "Kc_cdpAIS", "mM",
 "Ko_cdpAIS", "mM",
 "pf0_cdpAIS", "/ms",
 "pf1_cdpAIS", "/ms",
 "pf2_cdpAIS", "/ms",
 "pf3_cdpAIS", "/ms",
 "pf4_cdpAIS", "/ms",
 "pb0_cdpAIS", "/ms",
 "pb1_cdpAIS", "/ms",
 "pb2_cdpAIS", "/ms",
 "pb3_cdpAIS", "/ms",
 "pb4_cdpAIS", "/ms",
 "vrat_cdpAIS", "1",
 "ca_cdpAIS", "mM",
 "mg_cdpAIS", "mM",
 "CB_cdpAIS", "mM",
 "CB_f_ca_cdpAIS", "mM",
 "CB_ca_s_cdpAIS", "mM",
 "CB_ca_ca_cdpAIS", "mM",
 "iCB_cdpAIS", "mM",
 "iCB_f_ca_cdpAIS", "mM",
 "iCB_ca_s_cdpAIS", "mM",
 "iCB_ca_ca_cdpAIS", "mM",
 "PV_cdpAIS", "mM",
 "PV_ca_cdpAIS", "mM",
 "PV_mg_cdpAIS", "mM",
 "pump_cdpAIS", "mol/cm2",
 "pumpca_cdpAIS", "mol/cm2",
 "ica_pmp_cdpAIS", "mA/cm2",
 0,0
};
 static double BK_conc0 = 0;
 static double BK_ro0 = 0;
 static double C40 = 0;
 static double C30 = 0;
 static double C20 = 0;
 static double C10 = 0;
 static double C00 = 0;
 static double CB_ca_ca0 = 0;
 static double CB_ca_s0 = 0;
 static double CB_f_ca0 = 0;
 static double CB0 = 0;
 static double O40 = 0;
 static double O30 = 0;
 static double O20 = 0;
 static double O10 = 0;
 static double O00 = 0;
 static double PV_mg0 = 0;
 static double PV_ca0 = 0;
 static double PV0 = 0;
 static double ca0 = 0;
 static double delta_t = 0.01;
 static double dr_two0 = 0;
 static double iCB_ca_ca0 = 0;
 static double iCB_ca_s0 = 0;
 static double iCB_f_ca0 = 0;
 static double iCB0 = 0;
 static double mg0 = 0;
 static double pumpca0 = 0;
 static double pump0 = 0;
 static double rad_inner0 = 0;
 static double rad_outer0 = 0;
 static double vol_shell0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "cainull_cdpAIS", &cainull_cdpAIS,
 "mginull_cdpAIS", &mginull_cdpAIS,
 "DCa_cdpAIS", &DCa_cdpAIS,
 "Dbtc_cdpAIS", &Dbtc_cdpAIS,
 "Ddmnpe_cdpAIS", &Ddmnpe_cdpAIS,
 "Dcbd1_cdpAIS", &Dcbd1_cdpAIS,
 "Dcbd2_cdpAIS", &Dcbd2_cdpAIS,
 "Dpar_cdpAIS", &Dpar_cdpAIS,
 "CBnull_cdpAIS", &CBnull_cdpAIS,
 "nf1_cdpAIS", &nf1_cdpAIS,
 "nf2_cdpAIS", &nf2_cdpAIS,
 "ns1_cdpAIS", &ns1_cdpAIS,
 "ns2_cdpAIS", &ns2_cdpAIS,
 "PVnull_cdpAIS", &PVnull_cdpAIS,
 "m1_cdpAIS", &m1_cdpAIS,
 "m2_cdpAIS", &m2_cdpAIS,
 "p1_cdpAIS", &p1_cdpAIS,
 "p2_cdpAIS", &p2_cdpAIS,
 "kpmp1_cdpAIS", &kpmp1_cdpAIS,
 "kpmp2_cdpAIS", &kpmp2_cdpAIS,
 "kpmp3_cdpAIS", &kpmp3_cdpAIS,
 "TotalPump_cdpAIS", &TotalPump_cdpAIS,
 "beta_cdpAIS", &beta_cdpAIS,
 "vmax_cdpAIS", &vmax_cdpAIS,
 "Kp_cdpAIS", &Kp_cdpAIS,
 "Qo_cdpAIS", &Qo_cdpAIS,
 "Qc_cdpAIS", &Qc_cdpAIS,
 "k1_cdpAIS", &k1_cdpAIS,
 "onoffrate_cdpAIS", &onoffrate_cdpAIS,
 "L0_cdpAIS", &L0_cdpAIS,
 "Kc_cdpAIS", &Kc_cdpAIS,
 "Ko_cdpAIS", &Ko_cdpAIS,
 "pf0_cdpAIS", &pf0_cdpAIS,
 "pf1_cdpAIS", &pf1_cdpAIS,
 "pf2_cdpAIS", &pf2_cdpAIS,
 "pf3_cdpAIS", &pf3_cdpAIS,
 "pf4_cdpAIS", &pf4_cdpAIS,
 "pb0_cdpAIS", &pb0_cdpAIS,
 "pb1_cdpAIS", &pb1_cdpAIS,
 "pb2_cdpAIS", &pb2_cdpAIS,
 "pb3_cdpAIS", &pb3_cdpAIS,
 "pb4_cdpAIS", &pb4_cdpAIS,
 "BK_sing_chan_g_cdpAIS", &BK_sing_chan_g_cdpAIS,
 "BK_g_cdpAIS", &BK_g_cdpAIS,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 "vrat_cdpAIS", vrat_cdpAIS, 10,
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_synonym(int, double**, Datum**);
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"cdpAIS",
 0,
 "ica_pmp_cdpAIS",
 "ca1_cdpAIS",
 "ca2_cdpAIS",
 "ca3_cdpAIS",
 "ca4_cdpAIS",
 "ca5_cdpAIS",
 "ca9_cdpAIS",
 0,
 "ca_cdpAIS[10]",
 "mg_cdpAIS[10]",
 "CB_cdpAIS[10]",
 "CB_f_ca_cdpAIS[10]",
 "CB_ca_s_cdpAIS[10]",
 "CB_ca_ca_cdpAIS[10]",
 "iCB_cdpAIS[10]",
 "iCB_f_ca_cdpAIS[10]",
 "iCB_ca_s_cdpAIS[10]",
 "iCB_ca_ca_cdpAIS[10]",
 "PV_cdpAIS[10]",
 "PV_ca_cdpAIS[10]",
 "PV_mg_cdpAIS[10]",
 "pump_cdpAIS",
 "pumpca_cdpAIS",
 "C0_cdpAIS",
 "C1_cdpAIS",
 "C2_cdpAIS",
 "C3_cdpAIS",
 "C4_cdpAIS",
 "O0_cdpAIS",
 "O1_cdpAIS",
 "O2_cdpAIS",
 "O3_cdpAIS",
 "O4_cdpAIS",
 "BK_ro_cdpAIS",
 "BK_conc_cdpAIS",
 "dr_two_cdpAIS",
 "rad_outer_cdpAIS",
 "rad_inner_cdpAIS",
 "vol_shell_cdpAIS",
 0,
 0};
 static Symbol* _morphology_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 335, _prop);
 	/*initialize range parameters*/
 	_prop->param = _p;
 	_prop->param_size = 335;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[4]._pval = &prop_ion->param[0]; /* diam */
 prop_ion = need_memb(_ca_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[1]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for ca */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 "mg_cdpAIS", 1e-07,
 "pump_cdpAIS", 1e-15,
 "pumpca_cdpAIS", 1e-15,
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _cdp_AIS_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_morphology_sym = hoc_lookup("morphology");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 335, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "#ca_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 4, "diam");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_synonym(_mechtype, _ode_synonym);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cdpAIS /Users/iain/GIT/BK_unitary_g/NEURON/simple_model_50pS/mod/cdp_AIS.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define FARADAY _nrnunit_FARADAY[_nrnunit_use_legacy_]
static double _nrnunit_FARADAY[2] = {0x1.34c0c8b92a9b7p+3, 9.64853}; /* 9.64853321233100125 */
 
#define PI _nrnunit_PI[_nrnunit_use_legacy_]
static double _nrnunit_PI[2] = {0x1.921fb54442d18p+1, 3.14159}; /* 3.14159265358979312 */
 
#define R _nrnunit_R[_nrnunit_use_legacy_]
static double _nrnunit_R[2] = {0x1.0a1013e8990bep+3, 8.3145}; /* 8.3144626181532395 */
 static double q10 = 2;
 static double N_A = 6.02214076e23;
 static double cao = 2;
 static double _zfactors_done ;
 static double _zradii [ 10 ] ;
 static double _zfrat [ 10 ] ;
 static double _zdsq , _zdsqvol ;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int factors();
static int rates(double, double);
 extern double *_getelm(int, int);
 
#define _MATELM1(_row,_col)	*(_getelm(_row + 1, _col + 1))
 
#define _RHS1(_arg) _coef1[_arg + 1]
 static double *_coef1;
 
#define _linmat1  0
 static void* _sparseobj1;
 static void* _cvsparseobj1;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[140], _dlist1[140]; static double *_temp1;
 static int state();
 
static int  factors (  ) {
   double _lr , _ldr2 , _ldr3 ;
 _lr = 1.0 / 2.0 ;
   _ldr2 = 0.2 / diam ;
   _ldr3 = ( _lr - _ldr2 ) / ( 10.0 - 1.0 ) ;
   _zradii [ 0 ] = _lr ;
   _zradii [ 1 ] = _lr - _ldr2 ;
   {int  _li ;for ( _li = 2 ; _li <= 10 - 1 ; _li ++ ) {
     _zradii [ _li ] = _zradii [ _li - 1 ] - _ldr3 ;
     printf ( "%f\n" , _zradii [ _li ] ) ;
     } }
   vrat [ 0 ] = 0.0 ;
   _zfrat [ 0 ] = 2.0 * _lr ;
   {int  _li ;for ( _li = 0 ; _li <= 10 - 2 ; _li ++ ) {
     vrat [ _li ] = PI * ( ( _zradii [ _li ] * _zradii [ _li ] ) - ( _zradii [ _li + 1 ] * _zradii [ _li + 1 ] ) ) ;
     } }
   vrat [ 10 - 1 ] = PI * _zradii [ 10 - 1 ] * _zradii [ 10 - 1 ] ;
   {int  _li ;for ( _li = 1 ; _li <= 10 - 1 ; _li ++ ) {
     if ( ((double) _li )  == 1.0 ) {
       _zfrat [ _li ] = 2.0 * PI * _zradii [ _li ] / ( _ldr2 + ( _ldr3 / 2.0 ) ) ;
       }
     else if ( ((double) _li ) > 1.0  && ((double) _li ) < ( 10.0 - 1.0 ) ) {
       _zfrat [ _li ] = 2.0 * PI * _zradii [ _li ] / _ldr3 ;
       }
     else if ( ((double) _li )  == ( 10.0 - 1.0 ) ) {
       _zfrat [ _li ] = 2.0 * PI * _zradii [ _li ] / ( ( _ldr3 / 2.0 ) + _zradii [ _li ] ) ;
       }
     } }
    return 0; }
 
static void _hoc_factors(void) {
  double _r;
   _r = 1.;
 factors (  );
 hoc_retpushx(_r);
}
 
static int state ()
 {_reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<140;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 5) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 5, _i + 5) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 15) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 15, _i + 15) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 25) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 25, _i + 25) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 35) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 35, _i + 35) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 50) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 50, _i + 50) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 60) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 60, _i + 60) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 70) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 70, _i + 70) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 80) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 80, _i + 80) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 90) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 90, _i + 90) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 100) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 100, _i + 100) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 110) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 110, _i + 110) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 120) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 120, _i + 120) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 130) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 130, _i + 130) *= ( diam * diam * vrat [ ((int) _i ) ]);  } }
 /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
     ca mg CB CB_f_ca CB_ca_s CB_ca_ca iCB iCB_f_ca iCB_ca_s iCB_ca_ca PV PV_ca PV_mg }
   */
 /* COMPARTMENT ( 1e10 ) * parea {
     pump pumpca }
   */
 /* ~ ca [ 0 ] < < ( - ica * PI * diam / ( 2.0 * FARADAY ) )*/
 f_flux = b_flux = 0.;
 _RHS1( 80 +  0) += (b_flux =   ( - ica * PI * diam / ( 2.0 * FARADAY ) ) );
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 10 - 1 ; _li ++ ) {
     /* ~ ca [ _li ] < < ( - beta * vmax * vrat [ _li ] * ca [ _li ] / ( ca [ _li ] + kpmp2 / kpmp1 ) )*/
 f_flux = b_flux = 0.;
 _RHS1( 80 +  _li) += (b_flux =   ( - beta * vmax * vrat [ _li ] * ca [ _li ] / ( ca [ _li ] + kpmp2 / kpmp1 ) ) );
 /*FLUX*/
  } }
   {int  _li ;for ( _li = 0 ; _li <= 10 - 2 ; _li ++ ) {
     /* ~ ca [ _li ] <-> ca [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 f_flux =  DCa * _zfrat [ _li + 1 ] * ca [ _li] ;
 b_flux =  DCa * _zfrat [ _li + 1 ] * ca [ _li + 1] ;
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li + 1) += (f_flux - b_flux);
 
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li + 1 ,80 +  _li)  -= _term;
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 80 +  _li ,80 +  _li + 1)  -= _term;
 _MATELM1( 80 +  _li + 1 ,80 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ mg [ _li ] <-> mg [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 f_flux =  DCa * _zfrat [ _li + 1 ] * mg [ _li] ;
 b_flux =  DCa * _zfrat [ _li + 1 ] * mg [ _li + 1] ;
 _RHS1( 130 +  _li) -= (f_flux - b_flux);
 _RHS1( 130 +  _li + 1) += (f_flux - b_flux);
 
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 130 +  _li ,130 +  _li)  += _term;
 _MATELM1( 130 +  _li + 1 ,130 +  _li)  -= _term;
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 130 +  _li ,130 +  _li + 1)  -= _term;
 _MATELM1( 130 +  _li + 1 ,130 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ CB [ _li ] <-> CB [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 f_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB [ _li] ;
 b_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB [ _li + 1] ;
 _RHS1( 35 +  _li) -= (f_flux - b_flux);
 _RHS1( 35 +  _li + 1) += (f_flux - b_flux);
 
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 35 +  _li ,35 +  _li)  += _term;
 _MATELM1( 35 +  _li + 1 ,35 +  _li)  -= _term;
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 35 +  _li ,35 +  _li + 1)  -= _term;
 _MATELM1( 35 +  _li + 1 ,35 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ CB_f_ca [ _li ] <-> CB_f_ca [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 f_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_f_ca [ _li] ;
 b_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_f_ca [ _li + 1] ;
 _RHS1( 25 +  _li) -= (f_flux - b_flux);
 _RHS1( 25 +  _li + 1) += (f_flux - b_flux);
 
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 25 +  _li ,25 +  _li)  += _term;
 _MATELM1( 25 +  _li + 1 ,25 +  _li)  -= _term;
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 25 +  _li ,25 +  _li + 1)  -= _term;
 _MATELM1( 25 +  _li + 1 ,25 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ CB_ca_s [ _li ] <-> CB_ca_s [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 f_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_ca_s [ _li] ;
 b_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_ca_s [ _li + 1] ;
 _RHS1( 15 +  _li) -= (f_flux - b_flux);
 _RHS1( 15 +  _li + 1) += (f_flux - b_flux);
 
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 15 +  _li ,15 +  _li)  += _term;
 _MATELM1( 15 +  _li + 1 ,15 +  _li)  -= _term;
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 15 +  _li ,15 +  _li + 1)  -= _term;
 _MATELM1( 15 +  _li + 1 ,15 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ CB_ca_ca [ _li ] <-> CB_ca_ca [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 f_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_ca_ca [ _li] ;
 b_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_ca_ca [ _li + 1] ;
 _RHS1( 5 +  _li) -= (f_flux - b_flux);
 _RHS1( 5 +  _li + 1) += (f_flux - b_flux);
 
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 5 +  _li ,5 +  _li)  += _term;
 _MATELM1( 5 +  _li + 1 ,5 +  _li)  -= _term;
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 5 +  _li ,5 +  _li + 1)  -= _term;
 _MATELM1( 5 +  _li + 1 ,5 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ PV [ _li ] <-> PV [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 f_flux =  Dpar * _zfrat [ _li + 1 ] * PV [ _li] ;
 b_flux =  Dpar * _zfrat [ _li + 1 ] * PV [ _li + 1] ;
 _RHS1( 70 +  _li) -= (f_flux - b_flux);
 _RHS1( 70 +  _li + 1) += (f_flux - b_flux);
 
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 70 +  _li ,70 +  _li)  += _term;
 _MATELM1( 70 +  _li + 1 ,70 +  _li)  -= _term;
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 70 +  _li ,70 +  _li + 1)  -= _term;
 _MATELM1( 70 +  _li + 1 ,70 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ PV_ca [ _li ] <-> PV_ca [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 f_flux =  Dpar * _zfrat [ _li + 1 ] * PV_ca [ _li] ;
 b_flux =  Dpar * _zfrat [ _li + 1 ] * PV_ca [ _li + 1] ;
 _RHS1( 60 +  _li) -= (f_flux - b_flux);
 _RHS1( 60 +  _li + 1) += (f_flux - b_flux);
 
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 60 +  _li ,60 +  _li)  += _term;
 _MATELM1( 60 +  _li + 1 ,60 +  _li)  -= _term;
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 60 +  _li ,60 +  _li + 1)  -= _term;
 _MATELM1( 60 +  _li + 1 ,60 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ PV_mg [ _li ] <-> PV_mg [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 f_flux =  Dpar * _zfrat [ _li + 1 ] * PV_mg [ _li] ;
 b_flux =  Dpar * _zfrat [ _li + 1 ] * PV_mg [ _li + 1] ;
 _RHS1( 50 +  _li) -= (f_flux - b_flux);
 _RHS1( 50 +  _li + 1) += (f_flux - b_flux);
 
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 50 +  _li ,50 +  _li)  += _term;
 _MATELM1( 50 +  _li + 1 ,50 +  _li)  -= _term;
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 50 +  _li ,50 +  _li + 1)  -= _term;
 _MATELM1( 50 +  _li + 1 ,50 +  _li + 1)  += _term;
 /*REACTION*/
  } }
   _zdsq = diam * diam ;
   {int  _li ;for ( _li = 0 ; _li <= 10 - 1 ; _li ++ ) {
     _zdsqvol = _zdsq * vrat [ _li ] ;
     /* ~ ca [ _li ] + CB [ _li ] <-> CB_ca_s [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 f_flux =  nf1 * _zdsqvol * CB [ _li] * ca [ _li] ;
 b_flux =  nf2 * _zdsqvol * CB_ca_s [ _li] ;
 _RHS1( 35 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 15 +  _li) += (f_flux - b_flux);
 
 _term =  nf1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 35 +  _li ,35 +  _li)  += _term;
 _MATELM1( 80 +  _li ,35 +  _li)  += _term;
 _MATELM1( 15 +  _li ,35 +  _li)  -= _term;
 _term =  nf1 * _zdsqvol * CB [ _li] ;
 _MATELM1( 35 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 15 +  _li ,80 +  _li)  -= _term;
 _term =  nf2 * _zdsqvol ;
 _MATELM1( 35 +  _li ,15 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,15 +  _li)  -= _term;
 _MATELM1( 15 +  _li ,15 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + CB [ _li ] <-> CB_f_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 f_flux =  ns1 * _zdsqvol * CB [ _li] * ca [ _li] ;
 b_flux =  ns2 * _zdsqvol * CB_f_ca [ _li] ;
 _RHS1( 35 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 25 +  _li) += (f_flux - b_flux);
 
 _term =  ns1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 35 +  _li ,35 +  _li)  += _term;
 _MATELM1( 80 +  _li ,35 +  _li)  += _term;
 _MATELM1( 25 +  _li ,35 +  _li)  -= _term;
 _term =  ns1 * _zdsqvol * CB [ _li] ;
 _MATELM1( 35 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 25 +  _li ,80 +  _li)  -= _term;
 _term =  ns2 * _zdsqvol ;
 _MATELM1( 35 +  _li ,25 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,25 +  _li)  -= _term;
 _MATELM1( 25 +  _li ,25 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + CB_f_ca [ _li ] <-> CB_ca_ca [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 f_flux =  nf1 * _zdsqvol * CB_f_ca [ _li] * ca [ _li] ;
 b_flux =  nf2 * _zdsqvol * CB_ca_ca [ _li] ;
 _RHS1( 25 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 5 +  _li) += (f_flux - b_flux);
 
 _term =  nf1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 25 +  _li ,25 +  _li)  += _term;
 _MATELM1( 80 +  _li ,25 +  _li)  += _term;
 _MATELM1( 5 +  _li ,25 +  _li)  -= _term;
 _term =  nf1 * _zdsqvol * CB_f_ca [ _li] ;
 _MATELM1( 25 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 5 +  _li ,80 +  _li)  -= _term;
 _term =  nf2 * _zdsqvol ;
 _MATELM1( 25 +  _li ,5 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,5 +  _li)  -= _term;
 _MATELM1( 5 +  _li ,5 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + CB_ca_s [ _li ] <-> CB_ca_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 f_flux =  ns1 * _zdsqvol * CB_ca_s [ _li] * ca [ _li] ;
 b_flux =  ns2 * _zdsqvol * CB_ca_ca [ _li] ;
 _RHS1( 15 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 5 +  _li) += (f_flux - b_flux);
 
 _term =  ns1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 15 +  _li ,15 +  _li)  += _term;
 _MATELM1( 80 +  _li ,15 +  _li)  += _term;
 _MATELM1( 5 +  _li ,15 +  _li)  -= _term;
 _term =  ns1 * _zdsqvol * CB_ca_s [ _li] ;
 _MATELM1( 15 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 5 +  _li ,80 +  _li)  -= _term;
 _term =  ns2 * _zdsqvol ;
 _MATELM1( 15 +  _li ,5 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,5 +  _li)  -= _term;
 _MATELM1( 5 +  _li ,5 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + iCB [ _li ] <-> iCB_ca_s [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 f_flux =  nf1 * _zdsqvol * iCB [ _li] * ca [ _li] ;
 b_flux =  nf2 * _zdsqvol * iCB_ca_s [ _li] ;
 _RHS1( 120 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 100 +  _li) += (f_flux - b_flux);
 
 _term =  nf1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 120 +  _li ,120 +  _li)  += _term;
 _MATELM1( 80 +  _li ,120 +  _li)  += _term;
 _MATELM1( 100 +  _li ,120 +  _li)  -= _term;
 _term =  nf1 * _zdsqvol * iCB [ _li] ;
 _MATELM1( 120 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 100 +  _li ,80 +  _li)  -= _term;
 _term =  nf2 * _zdsqvol ;
 _MATELM1( 120 +  _li ,100 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,100 +  _li)  -= _term;
 _MATELM1( 100 +  _li ,100 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + iCB [ _li ] <-> iCB_f_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 f_flux =  ns1 * _zdsqvol * iCB [ _li] * ca [ _li] ;
 b_flux =  ns2 * _zdsqvol * iCB_f_ca [ _li] ;
 _RHS1( 120 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 110 +  _li) += (f_flux - b_flux);
 
 _term =  ns1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 120 +  _li ,120 +  _li)  += _term;
 _MATELM1( 80 +  _li ,120 +  _li)  += _term;
 _MATELM1( 110 +  _li ,120 +  _li)  -= _term;
 _term =  ns1 * _zdsqvol * iCB [ _li] ;
 _MATELM1( 120 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 110 +  _li ,80 +  _li)  -= _term;
 _term =  ns2 * _zdsqvol ;
 _MATELM1( 120 +  _li ,110 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,110 +  _li)  -= _term;
 _MATELM1( 110 +  _li ,110 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + iCB_f_ca [ _li ] <-> iCB_ca_ca [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 f_flux =  nf1 * _zdsqvol * iCB_f_ca [ _li] * ca [ _li] ;
 b_flux =  nf2 * _zdsqvol * iCB_ca_ca [ _li] ;
 _RHS1( 110 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 90 +  _li) += (f_flux - b_flux);
 
 _term =  nf1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 110 +  _li ,110 +  _li)  += _term;
 _MATELM1( 80 +  _li ,110 +  _li)  += _term;
 _MATELM1( 90 +  _li ,110 +  _li)  -= _term;
 _term =  nf1 * _zdsqvol * iCB_f_ca [ _li] ;
 _MATELM1( 110 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 90 +  _li ,80 +  _li)  -= _term;
 _term =  nf2 * _zdsqvol ;
 _MATELM1( 110 +  _li ,90 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,90 +  _li)  -= _term;
 _MATELM1( 90 +  _li ,90 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + iCB_ca_s [ _li ] <-> iCB_ca_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 f_flux =  ns1 * _zdsqvol * iCB_ca_s [ _li] * ca [ _li] ;
 b_flux =  ns2 * _zdsqvol * iCB_ca_ca [ _li] ;
 _RHS1( 100 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 90 +  _li) += (f_flux - b_flux);
 
 _term =  ns1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 100 +  _li ,100 +  _li)  += _term;
 _MATELM1( 80 +  _li ,100 +  _li)  += _term;
 _MATELM1( 90 +  _li ,100 +  _li)  -= _term;
 _term =  ns1 * _zdsqvol * iCB_ca_s [ _li] ;
 _MATELM1( 100 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 90 +  _li ,80 +  _li)  -= _term;
 _term =  ns2 * _zdsqvol ;
 _MATELM1( 100 +  _li ,90 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,90 +  _li)  -= _term;
 _MATELM1( 90 +  _li ,90 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + PV [ _li ] <-> PV_ca [ _li ] ( m1 * _zdsqvol , m2 * _zdsqvol )*/
 f_flux =  m1 * _zdsqvol * PV [ _li] * ca [ _li] ;
 b_flux =  m2 * _zdsqvol * PV_ca [ _li] ;
 _RHS1( 70 +  _li) -= (f_flux - b_flux);
 _RHS1( 80 +  _li) -= (f_flux - b_flux);
 _RHS1( 60 +  _li) += (f_flux - b_flux);
 
 _term =  m1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 70 +  _li ,70 +  _li)  += _term;
 _MATELM1( 80 +  _li ,70 +  _li)  += _term;
 _MATELM1( 60 +  _li ,70 +  _li)  -= _term;
 _term =  m1 * _zdsqvol * PV [ _li] ;
 _MATELM1( 70 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 60 +  _li ,80 +  _li)  -= _term;
 _term =  m2 * _zdsqvol ;
 _MATELM1( 70 +  _li ,60 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,60 +  _li)  -= _term;
 _MATELM1( 60 +  _li ,60 +  _li)  += _term;
 /*REACTION*/
  /* ~ mg [ _li ] + PV [ _li ] <-> PV_mg [ _li ] ( p1 * _zdsqvol , p2 * _zdsqvol )*/
 f_flux =  p1 * _zdsqvol * PV [ _li] * mg [ _li] ;
 b_flux =  p2 * _zdsqvol * PV_mg [ _li] ;
 _RHS1( 70 +  _li) -= (f_flux - b_flux);
 _RHS1( 130 +  _li) -= (f_flux - b_flux);
 _RHS1( 50 +  _li) += (f_flux - b_flux);
 
 _term =  p1 * _zdsqvol * mg [ _li] ;
 _MATELM1( 70 +  _li ,70 +  _li)  += _term;
 _MATELM1( 130 +  _li ,70 +  _li)  += _term;
 _MATELM1( 50 +  _li ,70 +  _li)  -= _term;
 _term =  p1 * _zdsqvol * PV [ _li] ;
 _MATELM1( 70 +  _li ,130 +  _li)  += _term;
 _MATELM1( 130 +  _li ,130 +  _li)  += _term;
 _MATELM1( 50 +  _li ,130 +  _li)  -= _term;
 _term =  p2 * _zdsqvol ;
 _MATELM1( 70 +  _li ,50 +  _li)  -= _term;
 _MATELM1( 130 +  _li ,50 +  _li)  -= _term;
 _MATELM1( 50 +  _li ,50 +  _li)  += _term;
 /*REACTION*/
  } }
   rates ( _threadargscomma_ v , cai ) ;
   /* ~ ca [ 0 ] + C0 <-> C1 ( c01 , c10 )*/
 f_flux =  c01 * C0 * ca [ 0] ;
 b_flux =  c10 * C1 ;
 _RHS1( 4) -= (f_flux - b_flux);
 _RHS1( 80 +  0) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  c01 * ca [ 0] ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 80 +  0 ,4)  += _term;
 _MATELM1( 3 ,4)  -= _term;
 _term =  c01 * C0 ;
 _MATELM1( 4 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 3 ,80 +  0)  -= _term;
 _term =  c10 ;
 _MATELM1( 4 ,3)  -= _term;
 _MATELM1( 80 +  0 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + C1 <-> C2 ( c12 , c21 )*/
 f_flux =  c12 * C1 * ca [ 0] ;
 b_flux =  c21 * C2 ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 80 +  0) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  c12 * ca [ 0] ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 80 +  0 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  c12 * C1 ;
 _MATELM1( 3 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 2 ,80 +  0)  -= _term;
 _term =  c21 ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 80 +  0 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + C2 <-> C3 ( c23 , c32 )*/
 f_flux =  c23 * C2 * ca [ 0] ;
 b_flux =  c32 * C3 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 80 +  0) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  c23 * ca [ 0] ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 80 +  0 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  c23 * C2 ;
 _MATELM1( 2 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 1 ,80 +  0)  -= _term;
 _term =  c32 ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 80 +  0 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + C3 <-> C4 ( c34 , c43 )*/
 f_flux =  c34 * C3 * ca [ 0] ;
 b_flux =  c43 * C4 ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 80 +  0) -= (f_flux - b_flux);
 _RHS1( 0) += (f_flux - b_flux);
 
 _term =  c34 * ca [ 0] ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 80 +  0 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  c34 * C3 ;
 _MATELM1( 1 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 0 ,80 +  0)  -= _term;
 _term =  c43 ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 80 +  0 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + O0 <-> O1 ( o01 , o10 )*/
 f_flux =  o01 * O0 * ca [ 0] ;
 b_flux =  o10 * O1 ;
 _RHS1( 49) -= (f_flux - b_flux);
 _RHS1( 80 +  0) -= (f_flux - b_flux);
 _RHS1( 48) += (f_flux - b_flux);
 
 _term =  o01 * ca [ 0] ;
 _MATELM1( 49 ,49)  += _term;
 _MATELM1( 80 +  0 ,49)  += _term;
 _MATELM1( 48 ,49)  -= _term;
 _term =  o01 * O0 ;
 _MATELM1( 49 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 48 ,80 +  0)  -= _term;
 _term =  o10 ;
 _MATELM1( 49 ,48)  -= _term;
 _MATELM1( 80 +  0 ,48)  -= _term;
 _MATELM1( 48 ,48)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + O1 <-> O2 ( o12 , o21 )*/
 f_flux =  o12 * O1 * ca [ 0] ;
 b_flux =  o21 * O2 ;
 _RHS1( 48) -= (f_flux - b_flux);
 _RHS1( 80 +  0) -= (f_flux - b_flux);
 _RHS1( 47) += (f_flux - b_flux);
 
 _term =  o12 * ca [ 0] ;
 _MATELM1( 48 ,48)  += _term;
 _MATELM1( 80 +  0 ,48)  += _term;
 _MATELM1( 47 ,48)  -= _term;
 _term =  o12 * O1 ;
 _MATELM1( 48 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 47 ,80 +  0)  -= _term;
 _term =  o21 ;
 _MATELM1( 48 ,47)  -= _term;
 _MATELM1( 80 +  0 ,47)  -= _term;
 _MATELM1( 47 ,47)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + O2 <-> O3 ( o23 , o32 )*/
 f_flux =  o23 * O2 * ca [ 0] ;
 b_flux =  o32 * O3 ;
 _RHS1( 47) -= (f_flux - b_flux);
 _RHS1( 80 +  0) -= (f_flux - b_flux);
 _RHS1( 46) += (f_flux - b_flux);
 
 _term =  o23 * ca [ 0] ;
 _MATELM1( 47 ,47)  += _term;
 _MATELM1( 80 +  0 ,47)  += _term;
 _MATELM1( 46 ,47)  -= _term;
 _term =  o23 * O2 ;
 _MATELM1( 47 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 46 ,80 +  0)  -= _term;
 _term =  o32 ;
 _MATELM1( 47 ,46)  -= _term;
 _MATELM1( 80 +  0 ,46)  -= _term;
 _MATELM1( 46 ,46)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + O3 <-> O4 ( o34 , o43 )*/
 f_flux =  o34 * O3 * ca [ 0] ;
 b_flux =  o43 * O4 ;
 _RHS1( 46) -= (f_flux - b_flux);
 _RHS1( 80 +  0) -= (f_flux - b_flux);
 _RHS1( 45) += (f_flux - b_flux);
 
 _term =  o34 * ca [ 0] ;
 _MATELM1( 46 ,46)  += _term;
 _MATELM1( 80 +  0 ,46)  += _term;
 _MATELM1( 45 ,46)  -= _term;
 _term =  o34 * O3 ;
 _MATELM1( 46 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 45 ,80 +  0)  -= _term;
 _term =  o43 ;
 _MATELM1( 46 ,45)  -= _term;
 _MATELM1( 80 +  0 ,45)  -= _term;
 _MATELM1( 45 ,45)  += _term;
 /*REACTION*/
  /* ~ C0 <-> O0 ( f0 , b0 )*/
 f_flux =  f0 * C0 ;
 b_flux =  b0 * O0 ;
 _RHS1( 4) -= (f_flux - b_flux);
 _RHS1( 49) += (f_flux - b_flux);
 
 _term =  f0 ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 49 ,4)  -= _term;
 _term =  b0 ;
 _MATELM1( 4 ,49)  -= _term;
 _MATELM1( 49 ,49)  += _term;
 /*REACTION*/
  /* ~ C1 <-> O1 ( f1 , b1 )*/
 f_flux =  f1 * C1 ;
 b_flux =  b1 * O1 ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 48) += (f_flux - b_flux);
 
 _term =  f1 ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 48 ,3)  -= _term;
 _term =  b1 ;
 _MATELM1( 3 ,48)  -= _term;
 _MATELM1( 48 ,48)  += _term;
 /*REACTION*/
  /* ~ C2 <-> O2 ( f2 , b2 )*/
 f_flux =  f2 * C2 ;
 b_flux =  b2 * O2 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 47) += (f_flux - b_flux);
 
 _term =  f2 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 47 ,2)  -= _term;
 _term =  b2 ;
 _MATELM1( 2 ,47)  -= _term;
 _MATELM1( 47 ,47)  += _term;
 /*REACTION*/
  /* ~ C3 <-> O3 ( f3 , b3 )*/
 f_flux =  f3 * C3 ;
 b_flux =  b3 * O3 ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 46) += (f_flux - b_flux);
 
 _term =  f3 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 46 ,1)  -= _term;
 _term =  b3 ;
 _MATELM1( 1 ,46)  -= _term;
 _MATELM1( 46 ,46)  += _term;
 /*REACTION*/
  /* ~ C4 <-> O4 ( f4 , b4 )*/
 f_flux =  f4 * C4 ;
 b_flux =  b4 * O4 ;
 _RHS1( 0) -= (f_flux - b_flux);
 _RHS1( 45) += (f_flux - b_flux);
 
 _term =  f4 ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 45 ,0)  -= _term;
 _term =  b4 ;
 _MATELM1( 0 ,45)  -= _term;
 _MATELM1( 45 ,45)  += _term;
 /*REACTION*/
  cai = ca [ 0 ] ;
   ca1 = ca [ 1 ] ;
   ca2 = ca [ 2 ] ;
   ca3 = ca [ 3 ] ;
   ca4 = ca [ 4 ] ;
   ca5 = ca [ 5 ] ;
   ca9 = ca [ 9 ] ;
   mgi = mg [ 0 ] ;
     } return _reset;
 }
 
double ssCB (  double _lkdf , double _lkds ) {
   double _lssCB;
 _lssCB = CBnull / ( 1.0 + kdf ( _threadargs_ ) + kds ( _threadargs_ ) + ( kdf ( _threadargs_ ) * kds ( _threadargs_ ) ) ) ;
   
return _lssCB;
 }
 
static void _hoc_ssCB(void) {
  double _r;
   _r =  ssCB (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double ssCBfast (  double _lkdf , double _lkds ) {
   double _lssCBfast;
 _lssCBfast = ( CBnull * kds ( _threadargs_ ) ) / ( 1.0 + kdf ( _threadargs_ ) + kds ( _threadargs_ ) + ( kdf ( _threadargs_ ) * kds ( _threadargs_ ) ) ) ;
   
return _lssCBfast;
 }
 
static void _hoc_ssCBfast(void) {
  double _r;
   _r =  ssCBfast (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double ssCBslow (  double _lkdf , double _lkds ) {
   double _lssCBslow;
 _lssCBslow = ( CBnull * kdf ( _threadargs_ ) ) / ( 1.0 + kdf ( _threadargs_ ) + kds ( _threadargs_ ) + ( kdf ( _threadargs_ ) * kds ( _threadargs_ ) ) ) ;
   
return _lssCBslow;
 }
 
static void _hoc_ssCBslow(void) {
  double _r;
   _r =  ssCBslow (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double ssCBca (  double _lkdf , double _lkds ) {
   double _lssCBca;
 _lssCBca = ( CBnull * kdf ( _threadargs_ ) * kds ( _threadargs_ ) ) / ( 1.0 + kdf ( _threadargs_ ) + kds ( _threadargs_ ) + ( kdf ( _threadargs_ ) * kds ( _threadargs_ ) ) ) ;
   
return _lssCBca;
 }
 
static void _hoc_ssCBca(void) {
  double _r;
   _r =  ssCBca (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double kdf (  ) {
   double _lkdf;
 _lkdf = ( cainull * nf1 ) / nf2 ;
   
return _lkdf;
 }
 
static void _hoc_kdf(void) {
  double _r;
   _r =  kdf (  );
 hoc_retpushx(_r);
}
 
double kds (  ) {
   double _lkds;
 _lkds = ( cainull * ns1 ) / ns2 ;
   
return _lkds;
 }
 
static void _hoc_kds(void) {
  double _r;
   _r =  kds (  );
 hoc_retpushx(_r);
}
 
double kdc (  ) {
   double _lkdc;
 _lkdc = ( cainull * m1 ) / m2 ;
   
return _lkdc;
 }
 
static void _hoc_kdc(void) {
  double _r;
   _r =  kdc (  );
 hoc_retpushx(_r);
}
 
double kdm (  ) {
   double _lkdm;
 _lkdm = ( mginull * p1 ) / p2 ;
   
return _lkdm;
 }
 
static void _hoc_kdm(void) {
  double _r;
   _r =  kdm (  );
 hoc_retpushx(_r);
}
 
double ssPV (  double _lkdc , double _lkdm ) {
   double _lssPV;
 _lssPV = PVnull / ( 1.0 + kdc ( _threadargs_ ) + kdm ( _threadargs_ ) ) ;
   
return _lssPV;
 }
 
static void _hoc_ssPV(void) {
  double _r;
   _r =  ssPV (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double ssPVca (  double _lkdc , double _lkdm ) {
   double _lssPVca;
 _lssPVca = ( PVnull * _lkdc ) / ( 1.0 + _lkdc + _lkdm ) ;
   
return _lssPVca;
 }
 
static void _hoc_ssPVca(void) {
  double _r;
   _r =  ssPVca (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double ssPVmg (  double _lkdc , double _lkdm ) {
   double _lssPVmg;
 _lssPVmg = ( PVnull * kdm ( _threadargs_ ) ) / ( 1.0 + kdc ( _threadargs_ ) + kdm ( _threadargs_ ) ) ;
   
return _lssPVmg;
 }
 
static void _hoc_ssPVmg(void) {
  double _r;
   _r =  ssPVmg (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double u (  double _lx , double _lth ) {
   double _lu;
 if ( _lx < _lth ) {
     _lu = 1.0 ;
     }
   else {
     _lu = 0.0 ;
     }
   
return _lu;
 }
 
static void _hoc_u(void) {
  double _r;
   _r =  u (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv , double _lca ) {
   double _lqt , _lalpha , _lbeta ;
 _lqt = pow( q10 , ( ( celsius - 25.0 ) / 10.0 ) ) ;
   c01 = 4.0 * k1 * onoffrate * _lqt ;
   c12 = 3.0 * k1 * onoffrate * _lqt ;
   c23 = 2.0 * k1 * onoffrate * _lqt ;
   c34 = 1.0 * k1 * onoffrate * _lqt ;
   o01 = 4.0 * k1 * onoffrate * _lqt ;
   o12 = 3.0 * k1 * onoffrate * _lqt ;
   o23 = 2.0 * k1 * onoffrate * _lqt ;
   o34 = 1.0 * k1 * onoffrate * _lqt ;
   c10 = 1.0 * Kc * k1 * onoffrate * _lqt ;
   c21 = 2.0 * Kc * k1 * onoffrate * _lqt ;
   c32 = 3.0 * Kc * k1 * onoffrate * _lqt ;
   c43 = 4.0 * Kc * k1 * onoffrate * _lqt ;
   o10 = 1.0 * Ko * k1 * onoffrate * _lqt ;
   o21 = 2.0 * Ko * k1 * onoffrate * _lqt ;
   o32 = 3.0 * Ko * k1 * onoffrate * _lqt ;
   o43 = 4.0 * Ko * k1 * onoffrate * _lqt ;
   _lalpha = exp ( Qo * FARADAY * 10.0 * _lv / R / ( 273.15 + celsius ) ) ;
   _lbeta = exp ( Qc * FARADAY * 10.0 * _lv / R / ( 273.15 + celsius ) ) ;
   f0 = pf0 * _lalpha * _lqt ;
   f1 = pf1 * _lalpha * _lqt ;
   f2 = pf2 * _lalpha * _lqt ;
   f3 = pf3 * _lalpha * _lqt ;
   f4 = pf4 * _lalpha * _lqt ;
   b0 = pb0 * _lbeta * _lqt ;
   b1 = pb1 * _lbeta * _lqt ;
   b2 = pb2 * _lbeta * _lqt ;
   b3 = pb3 * _lbeta * _lqt ;
   b4 = pb4 * _lbeta * _lqt ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<140;_i++) _p[_dlist1[_i]] = 0.0;}
 /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
   ca mg CB CB_f_ca CB_ca_s CB_ca_ca iCB iCB_f_ca iCB_ca_s iCB_ca_ca PV PV_ca PV_mg }
 */
 /* COMPARTMENT ( 1e10 ) * parea {
   pump pumpca }
 */
 /* ~ ca [ 0 ] < < ( - ica * PI * diam / ( 2.0 * FARADAY ) )*/
 f_flux = b_flux = 0.;
 Dca [ 0] += (b_flux =   ( - ica * PI * diam / ( 2.0 * FARADAY ) ) );
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 10 - 1 ; _li ++ ) {
   /* ~ ca [ _li ] < < ( - beta * vmax * vrat [ _li ] * ca [ _li ] / ( ca [ _li ] + kpmp2 / kpmp1 ) )*/
 f_flux = b_flux = 0.;
 Dca [ _li] += (b_flux =   ( - beta * vmax * vrat [ _li ] * ca [ _li ] / ( ca [ _li ] + kpmp2 / kpmp1 ) ) );
 /*FLUX*/
  } }
 {int  _li ;for ( _li = 0 ; _li <= 10 - 2 ; _li ++ ) {
   /* ~ ca [ _li ] <-> ca [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 f_flux =  DCa * _zfrat [ _li + 1 ] * ca [ _li] ;
 b_flux =  DCa * _zfrat [ _li + 1 ] * ca [ _li + 1] ;
 Dca [ _li] -= (f_flux - b_flux);
 Dca [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ mg [ _li ] <-> mg [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 f_flux =  DCa * _zfrat [ _li + 1 ] * mg [ _li] ;
 b_flux =  DCa * _zfrat [ _li + 1 ] * mg [ _li + 1] ;
 Dmg [ _li] -= (f_flux - b_flux);
 Dmg [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ CB [ _li ] <-> CB [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 f_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB [ _li] ;
 b_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB [ _li + 1] ;
 DCB [ _li] -= (f_flux - b_flux);
 DCB [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ CB_f_ca [ _li ] <-> CB_f_ca [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 f_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_f_ca [ _li] ;
 b_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_f_ca [ _li + 1] ;
 DCB_f_ca [ _li] -= (f_flux - b_flux);
 DCB_f_ca [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ CB_ca_s [ _li ] <-> CB_ca_s [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 f_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_ca_s [ _li] ;
 b_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_ca_s [ _li + 1] ;
 DCB_ca_s [ _li] -= (f_flux - b_flux);
 DCB_ca_s [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ CB_ca_ca [ _li ] <-> CB_ca_ca [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 f_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_ca_ca [ _li] ;
 b_flux =  Dcbd1 * _zfrat [ _li + 1 ] * CB_ca_ca [ _li + 1] ;
 DCB_ca_ca [ _li] -= (f_flux - b_flux);
 DCB_ca_ca [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ PV [ _li ] <-> PV [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 f_flux =  Dpar * _zfrat [ _li + 1 ] * PV [ _li] ;
 b_flux =  Dpar * _zfrat [ _li + 1 ] * PV [ _li + 1] ;
 DPV [ _li] -= (f_flux - b_flux);
 DPV [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ PV_ca [ _li ] <-> PV_ca [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 f_flux =  Dpar * _zfrat [ _li + 1 ] * PV_ca [ _li] ;
 b_flux =  Dpar * _zfrat [ _li + 1 ] * PV_ca [ _li + 1] ;
 DPV_ca [ _li] -= (f_flux - b_flux);
 DPV_ca [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ PV_mg [ _li ] <-> PV_mg [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 f_flux =  Dpar * _zfrat [ _li + 1 ] * PV_mg [ _li] ;
 b_flux =  Dpar * _zfrat [ _li + 1 ] * PV_mg [ _li + 1] ;
 DPV_mg [ _li] -= (f_flux - b_flux);
 DPV_mg [ _li + 1] += (f_flux - b_flux);
 
 /*REACTION*/
  } }
 _zdsq = diam * diam ;
 {int  _li ;for ( _li = 0 ; _li <= 10 - 1 ; _li ++ ) {
   _zdsqvol = _zdsq * vrat [ _li ] ;
   /* ~ ca [ _li ] + CB [ _li ] <-> CB_ca_s [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 f_flux =  nf1 * _zdsqvol * CB [ _li] * ca [ _li] ;
 b_flux =  nf2 * _zdsqvol * CB_ca_s [ _li] ;
 DCB [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DCB_ca_s [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ _li ] + CB [ _li ] <-> CB_f_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 f_flux =  ns1 * _zdsqvol * CB [ _li] * ca [ _li] ;
 b_flux =  ns2 * _zdsqvol * CB_f_ca [ _li] ;
 DCB [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DCB_f_ca [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ _li ] + CB_f_ca [ _li ] <-> CB_ca_ca [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 f_flux =  nf1 * _zdsqvol * CB_f_ca [ _li] * ca [ _li] ;
 b_flux =  nf2 * _zdsqvol * CB_ca_ca [ _li] ;
 DCB_f_ca [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DCB_ca_ca [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ _li ] + CB_ca_s [ _li ] <-> CB_ca_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 f_flux =  ns1 * _zdsqvol * CB_ca_s [ _li] * ca [ _li] ;
 b_flux =  ns2 * _zdsqvol * CB_ca_ca [ _li] ;
 DCB_ca_s [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DCB_ca_ca [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ _li ] + iCB [ _li ] <-> iCB_ca_s [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 f_flux =  nf1 * _zdsqvol * iCB [ _li] * ca [ _li] ;
 b_flux =  nf2 * _zdsqvol * iCB_ca_s [ _li] ;
 DiCB [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DiCB_ca_s [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ _li ] + iCB [ _li ] <-> iCB_f_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 f_flux =  ns1 * _zdsqvol * iCB [ _li] * ca [ _li] ;
 b_flux =  ns2 * _zdsqvol * iCB_f_ca [ _li] ;
 DiCB [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DiCB_f_ca [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ _li ] + iCB_f_ca [ _li ] <-> iCB_ca_ca [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 f_flux =  nf1 * _zdsqvol * iCB_f_ca [ _li] * ca [ _li] ;
 b_flux =  nf2 * _zdsqvol * iCB_ca_ca [ _li] ;
 DiCB_f_ca [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DiCB_ca_ca [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ _li ] + iCB_ca_s [ _li ] <-> iCB_ca_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 f_flux =  ns1 * _zdsqvol * iCB_ca_s [ _li] * ca [ _li] ;
 b_flux =  ns2 * _zdsqvol * iCB_ca_ca [ _li] ;
 DiCB_ca_s [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DiCB_ca_ca [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ _li ] + PV [ _li ] <-> PV_ca [ _li ] ( m1 * _zdsqvol , m2 * _zdsqvol )*/
 f_flux =  m1 * _zdsqvol * PV [ _li] * ca [ _li] ;
 b_flux =  m2 * _zdsqvol * PV_ca [ _li] ;
 DPV [ _li] -= (f_flux - b_flux);
 Dca [ _li] -= (f_flux - b_flux);
 DPV_ca [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ mg [ _li ] + PV [ _li ] <-> PV_mg [ _li ] ( p1 * _zdsqvol , p2 * _zdsqvol )*/
 f_flux =  p1 * _zdsqvol * PV [ _li] * mg [ _li] ;
 b_flux =  p2 * _zdsqvol * PV_mg [ _li] ;
 DPV [ _li] -= (f_flux - b_flux);
 Dmg [ _li] -= (f_flux - b_flux);
 DPV_mg [ _li] += (f_flux - b_flux);
 
 /*REACTION*/
  } }
 rates ( _threadargscomma_ v , cai ) ;
 /* ~ ca [ 0 ] + C0 <-> C1 ( c01 , c10 )*/
 f_flux =  c01 * C0 * ca [ 0] ;
 b_flux =  c10 * C1 ;
 DC0 -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 DC1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ 0 ] + C1 <-> C2 ( c12 , c21 )*/
 f_flux =  c12 * C1 * ca [ 0] ;
 b_flux =  c21 * C2 ;
 DC1 -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 DC2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ 0 ] + C2 <-> C3 ( c23 , c32 )*/
 f_flux =  c23 * C2 * ca [ 0] ;
 b_flux =  c32 * C3 ;
 DC2 -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 DC3 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ 0 ] + C3 <-> C4 ( c34 , c43 )*/
 f_flux =  c34 * C3 * ca [ 0] ;
 b_flux =  c43 * C4 ;
 DC3 -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 DC4 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ 0 ] + O0 <-> O1 ( o01 , o10 )*/
 f_flux =  o01 * O0 * ca [ 0] ;
 b_flux =  o10 * O1 ;
 DO0 -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 DO1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ 0 ] + O1 <-> O2 ( o12 , o21 )*/
 f_flux =  o12 * O1 * ca [ 0] ;
 b_flux =  o21 * O2 ;
 DO1 -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 DO2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ 0 ] + O2 <-> O3 ( o23 , o32 )*/
 f_flux =  o23 * O2 * ca [ 0] ;
 b_flux =  o32 * O3 ;
 DO2 -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 DO3 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ca [ 0 ] + O3 <-> O4 ( o34 , o43 )*/
 f_flux =  o34 * O3 * ca [ 0] ;
 b_flux =  o43 * O4 ;
 DO3 -= (f_flux - b_flux);
 Dca [ 0] -= (f_flux - b_flux);
 DO4 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C0 <-> O0 ( f0 , b0 )*/
 f_flux =  f0 * C0 ;
 b_flux =  b0 * O0 ;
 DC0 -= (f_flux - b_flux);
 DO0 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C1 <-> O1 ( f1 , b1 )*/
 f_flux =  f1 * C1 ;
 b_flux =  b1 * O1 ;
 DC1 -= (f_flux - b_flux);
 DO1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C2 <-> O2 ( f2 , b2 )*/
 f_flux =  f2 * C2 ;
 b_flux =  b2 * O2 ;
 DC2 -= (f_flux - b_flux);
 DO2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C3 <-> O3 ( f3 , b3 )*/
 f_flux =  f3 * C3 ;
 b_flux =  b3 * O3 ;
 DC3 -= (f_flux - b_flux);
 DO3 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C4 <-> O4 ( f4 , b4 )*/
 f_flux =  f4 * C4 ;
 b_flux =  b4 * O4 ;
 DC4 -= (f_flux - b_flux);
 DO4 += (f_flux - b_flux);
 
 /*REACTION*/
  cai = ca [ 0 ] ;
 ca1 = ca [ 1 ] ;
 ca2 = ca [ 2 ] ;
 ca3 = ca [ 3 ] ;
 ca4 = ca [ 4 ] ;
 ca5 = ca [ 5 ] ;
 ca9 = ca [ 9 ] ;
 mgi = mg [ 0 ] ;
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 5]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 15]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 25]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 35]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 50]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 60]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 70]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 80]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 90]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 100]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 110]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 120]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
 for (_i=0; _i < 10; _i++) { _p[_dlist1[_i + 130]] /= ( diam * diam * vrat [ ((int) _i ) ]);}
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1() {_reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<140;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 5) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 5, _i + 5) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 15) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 15, _i + 15) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 25) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 25, _i + 25) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 35) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 35, _i + 35) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 50) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 50, _i + 50) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 60) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 60, _i + 60) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 70) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 70, _i + 70) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 80) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 80, _i + 80) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 90) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 90, _i + 90) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 100) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 100, _i + 100) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 110) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 110, _i + 110) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 120) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 120, _i + 120) *= ( diam * diam * vrat [ ((int) _i ) ]);  } 
for (_i=0; _i < 10; _i++) {
  	_RHS1(_i + 130) *= ( diam * diam * vrat [ ((int) _i ) ]) ;
_MATELM1(_i + 130, _i + 130) *= ( diam * diam * vrat [ ((int) _i ) ]);  } }
 /* COMPARTMENT _li , diam * diam * vrat [ ((int) _i ) ] {
 ca mg CB CB_f_ca CB_ca_s CB_ca_ca iCB iCB_f_ca iCB_ca_s iCB_ca_ca PV PV_ca PV_mg }
 */
 /* COMPARTMENT ( 1e10 ) * parea {
 pump pumpca }
 */
 /* ~ ca [ 0 ] < < ( - ica * PI * diam / ( 2.0 * FARADAY ) )*/
 /*FLUX*/
  {int  _li ;for ( _li = 0 ; _li <= 10 - 1 ; _li ++ ) {
 /* ~ ca [ _li ] < < ( - beta * vmax * vrat [ _li ] * ca [ _li ] / ( ca [ _li ] + kpmp2 / kpmp1 ) )*/
 /*FLUX*/
  } }
 {int  _li ;for ( _li = 0 ; _li <= 10 - 2 ; _li ++ ) {
 /* ~ ca [ _li ] <-> ca [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li + 1 ,80 +  _li)  -= _term;
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 80 +  _li ,80 +  _li + 1)  -= _term;
 _MATELM1( 80 +  _li + 1 ,80 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ mg [ _li ] <-> mg [ _li + 1 ] ( DCa * _zfrat [ _li + 1 ] , DCa * _zfrat [ _li + 1 ] )*/
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 130 +  _li ,130 +  _li)  += _term;
 _MATELM1( 130 +  _li + 1 ,130 +  _li)  -= _term;
 _term =  DCa * _zfrat [ _li + 1 ] ;
 _MATELM1( 130 +  _li ,130 +  _li + 1)  -= _term;
 _MATELM1( 130 +  _li + 1 ,130 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ CB [ _li ] <-> CB [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 35 +  _li ,35 +  _li)  += _term;
 _MATELM1( 35 +  _li + 1 ,35 +  _li)  -= _term;
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 35 +  _li ,35 +  _li + 1)  -= _term;
 _MATELM1( 35 +  _li + 1 ,35 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ CB_f_ca [ _li ] <-> CB_f_ca [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 25 +  _li ,25 +  _li)  += _term;
 _MATELM1( 25 +  _li + 1 ,25 +  _li)  -= _term;
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 25 +  _li ,25 +  _li + 1)  -= _term;
 _MATELM1( 25 +  _li + 1 ,25 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ CB_ca_s [ _li ] <-> CB_ca_s [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 15 +  _li ,15 +  _li)  += _term;
 _MATELM1( 15 +  _li + 1 ,15 +  _li)  -= _term;
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 15 +  _li ,15 +  _li + 1)  -= _term;
 _MATELM1( 15 +  _li + 1 ,15 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ CB_ca_ca [ _li ] <-> CB_ca_ca [ _li + 1 ] ( Dcbd1 * _zfrat [ _li + 1 ] , Dcbd1 * _zfrat [ _li + 1 ] )*/
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 5 +  _li ,5 +  _li)  += _term;
 _MATELM1( 5 +  _li + 1 ,5 +  _li)  -= _term;
 _term =  Dcbd1 * _zfrat [ _li + 1 ] ;
 _MATELM1( 5 +  _li ,5 +  _li + 1)  -= _term;
 _MATELM1( 5 +  _li + 1 ,5 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ PV [ _li ] <-> PV [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 70 +  _li ,70 +  _li)  += _term;
 _MATELM1( 70 +  _li + 1 ,70 +  _li)  -= _term;
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 70 +  _li ,70 +  _li + 1)  -= _term;
 _MATELM1( 70 +  _li + 1 ,70 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ PV_ca [ _li ] <-> PV_ca [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 60 +  _li ,60 +  _li)  += _term;
 _MATELM1( 60 +  _li + 1 ,60 +  _li)  -= _term;
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 60 +  _li ,60 +  _li + 1)  -= _term;
 _MATELM1( 60 +  _li + 1 ,60 +  _li + 1)  += _term;
 /*REACTION*/
  /* ~ PV_mg [ _li ] <-> PV_mg [ _li + 1 ] ( Dpar * _zfrat [ _li + 1 ] , Dpar * _zfrat [ _li + 1 ] )*/
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 50 +  _li ,50 +  _li)  += _term;
 _MATELM1( 50 +  _li + 1 ,50 +  _li)  -= _term;
 _term =  Dpar * _zfrat [ _li + 1 ] ;
 _MATELM1( 50 +  _li ,50 +  _li + 1)  -= _term;
 _MATELM1( 50 +  _li + 1 ,50 +  _li + 1)  += _term;
 /*REACTION*/
  } }
 _zdsq = diam * diam ;
 {int  _li ;for ( _li = 0 ; _li <= 10 - 1 ; _li ++ ) {
 _zdsqvol = _zdsq * vrat [ _li ] ;
 /* ~ ca [ _li ] + CB [ _li ] <-> CB_ca_s [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 _term =  nf1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 35 +  _li ,35 +  _li)  += _term;
 _MATELM1( 80 +  _li ,35 +  _li)  += _term;
 _MATELM1( 15 +  _li ,35 +  _li)  -= _term;
 _term =  nf1 * _zdsqvol * CB [ _li] ;
 _MATELM1( 35 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 15 +  _li ,80 +  _li)  -= _term;
 _term =  nf2 * _zdsqvol ;
 _MATELM1( 35 +  _li ,15 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,15 +  _li)  -= _term;
 _MATELM1( 15 +  _li ,15 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + CB [ _li ] <-> CB_f_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 _term =  ns1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 35 +  _li ,35 +  _li)  += _term;
 _MATELM1( 80 +  _li ,35 +  _li)  += _term;
 _MATELM1( 25 +  _li ,35 +  _li)  -= _term;
 _term =  ns1 * _zdsqvol * CB [ _li] ;
 _MATELM1( 35 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 25 +  _li ,80 +  _li)  -= _term;
 _term =  ns2 * _zdsqvol ;
 _MATELM1( 35 +  _li ,25 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,25 +  _li)  -= _term;
 _MATELM1( 25 +  _li ,25 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + CB_f_ca [ _li ] <-> CB_ca_ca [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 _term =  nf1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 25 +  _li ,25 +  _li)  += _term;
 _MATELM1( 80 +  _li ,25 +  _li)  += _term;
 _MATELM1( 5 +  _li ,25 +  _li)  -= _term;
 _term =  nf1 * _zdsqvol * CB_f_ca [ _li] ;
 _MATELM1( 25 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 5 +  _li ,80 +  _li)  -= _term;
 _term =  nf2 * _zdsqvol ;
 _MATELM1( 25 +  _li ,5 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,5 +  _li)  -= _term;
 _MATELM1( 5 +  _li ,5 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + CB_ca_s [ _li ] <-> CB_ca_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 _term =  ns1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 15 +  _li ,15 +  _li)  += _term;
 _MATELM1( 80 +  _li ,15 +  _li)  += _term;
 _MATELM1( 5 +  _li ,15 +  _li)  -= _term;
 _term =  ns1 * _zdsqvol * CB_ca_s [ _li] ;
 _MATELM1( 15 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 5 +  _li ,80 +  _li)  -= _term;
 _term =  ns2 * _zdsqvol ;
 _MATELM1( 15 +  _li ,5 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,5 +  _li)  -= _term;
 _MATELM1( 5 +  _li ,5 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + iCB [ _li ] <-> iCB_ca_s [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 _term =  nf1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 120 +  _li ,120 +  _li)  += _term;
 _MATELM1( 80 +  _li ,120 +  _li)  += _term;
 _MATELM1( 100 +  _li ,120 +  _li)  -= _term;
 _term =  nf1 * _zdsqvol * iCB [ _li] ;
 _MATELM1( 120 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 100 +  _li ,80 +  _li)  -= _term;
 _term =  nf2 * _zdsqvol ;
 _MATELM1( 120 +  _li ,100 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,100 +  _li)  -= _term;
 _MATELM1( 100 +  _li ,100 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + iCB [ _li ] <-> iCB_f_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 _term =  ns1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 120 +  _li ,120 +  _li)  += _term;
 _MATELM1( 80 +  _li ,120 +  _li)  += _term;
 _MATELM1( 110 +  _li ,120 +  _li)  -= _term;
 _term =  ns1 * _zdsqvol * iCB [ _li] ;
 _MATELM1( 120 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 110 +  _li ,80 +  _li)  -= _term;
 _term =  ns2 * _zdsqvol ;
 _MATELM1( 120 +  _li ,110 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,110 +  _li)  -= _term;
 _MATELM1( 110 +  _li ,110 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + iCB_f_ca [ _li ] <-> iCB_ca_ca [ _li ] ( nf1 * _zdsqvol , nf2 * _zdsqvol )*/
 _term =  nf1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 110 +  _li ,110 +  _li)  += _term;
 _MATELM1( 80 +  _li ,110 +  _li)  += _term;
 _MATELM1( 90 +  _li ,110 +  _li)  -= _term;
 _term =  nf1 * _zdsqvol * iCB_f_ca [ _li] ;
 _MATELM1( 110 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 90 +  _li ,80 +  _li)  -= _term;
 _term =  nf2 * _zdsqvol ;
 _MATELM1( 110 +  _li ,90 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,90 +  _li)  -= _term;
 _MATELM1( 90 +  _li ,90 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + iCB_ca_s [ _li ] <-> iCB_ca_ca [ _li ] ( ns1 * _zdsqvol , ns2 * _zdsqvol )*/
 _term =  ns1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 100 +  _li ,100 +  _li)  += _term;
 _MATELM1( 80 +  _li ,100 +  _li)  += _term;
 _MATELM1( 90 +  _li ,100 +  _li)  -= _term;
 _term =  ns1 * _zdsqvol * iCB_ca_s [ _li] ;
 _MATELM1( 100 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 90 +  _li ,80 +  _li)  -= _term;
 _term =  ns2 * _zdsqvol ;
 _MATELM1( 100 +  _li ,90 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,90 +  _li)  -= _term;
 _MATELM1( 90 +  _li ,90 +  _li)  += _term;
 /*REACTION*/
  /* ~ ca [ _li ] + PV [ _li ] <-> PV_ca [ _li ] ( m1 * _zdsqvol , m2 * _zdsqvol )*/
 _term =  m1 * _zdsqvol * ca [ _li] ;
 _MATELM1( 70 +  _li ,70 +  _li)  += _term;
 _MATELM1( 80 +  _li ,70 +  _li)  += _term;
 _MATELM1( 60 +  _li ,70 +  _li)  -= _term;
 _term =  m1 * _zdsqvol * PV [ _li] ;
 _MATELM1( 70 +  _li ,80 +  _li)  += _term;
 _MATELM1( 80 +  _li ,80 +  _li)  += _term;
 _MATELM1( 60 +  _li ,80 +  _li)  -= _term;
 _term =  m2 * _zdsqvol ;
 _MATELM1( 70 +  _li ,60 +  _li)  -= _term;
 _MATELM1( 80 +  _li ,60 +  _li)  -= _term;
 _MATELM1( 60 +  _li ,60 +  _li)  += _term;
 /*REACTION*/
  /* ~ mg [ _li ] + PV [ _li ] <-> PV_mg [ _li ] ( p1 * _zdsqvol , p2 * _zdsqvol )*/
 _term =  p1 * _zdsqvol * mg [ _li] ;
 _MATELM1( 70 +  _li ,70 +  _li)  += _term;
 _MATELM1( 130 +  _li ,70 +  _li)  += _term;
 _MATELM1( 50 +  _li ,70 +  _li)  -= _term;
 _term =  p1 * _zdsqvol * PV [ _li] ;
 _MATELM1( 70 +  _li ,130 +  _li)  += _term;
 _MATELM1( 130 +  _li ,130 +  _li)  += _term;
 _MATELM1( 50 +  _li ,130 +  _li)  -= _term;
 _term =  p2 * _zdsqvol ;
 _MATELM1( 70 +  _li ,50 +  _li)  -= _term;
 _MATELM1( 130 +  _li ,50 +  _li)  -= _term;
 _MATELM1( 50 +  _li ,50 +  _li)  += _term;
 /*REACTION*/
  } }
 rates ( _threadargscomma_ v , cai ) ;
 /* ~ ca [ 0 ] + C0 <-> C1 ( c01 , c10 )*/
 _term =  c01 * ca [ 0] ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 80 +  0 ,4)  += _term;
 _MATELM1( 3 ,4)  -= _term;
 _term =  c01 * C0 ;
 _MATELM1( 4 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 3 ,80 +  0)  -= _term;
 _term =  c10 ;
 _MATELM1( 4 ,3)  -= _term;
 _MATELM1( 80 +  0 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + C1 <-> C2 ( c12 , c21 )*/
 _term =  c12 * ca [ 0] ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 80 +  0 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  c12 * C1 ;
 _MATELM1( 3 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 2 ,80 +  0)  -= _term;
 _term =  c21 ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 80 +  0 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + C2 <-> C3 ( c23 , c32 )*/
 _term =  c23 * ca [ 0] ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 80 +  0 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  c23 * C2 ;
 _MATELM1( 2 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 1 ,80 +  0)  -= _term;
 _term =  c32 ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 80 +  0 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + C3 <-> C4 ( c34 , c43 )*/
 _term =  c34 * ca [ 0] ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 80 +  0 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  c34 * C3 ;
 _MATELM1( 1 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 0 ,80 +  0)  -= _term;
 _term =  c43 ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 80 +  0 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + O0 <-> O1 ( o01 , o10 )*/
 _term =  o01 * ca [ 0] ;
 _MATELM1( 49 ,49)  += _term;
 _MATELM1( 80 +  0 ,49)  += _term;
 _MATELM1( 48 ,49)  -= _term;
 _term =  o01 * O0 ;
 _MATELM1( 49 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 48 ,80 +  0)  -= _term;
 _term =  o10 ;
 _MATELM1( 49 ,48)  -= _term;
 _MATELM1( 80 +  0 ,48)  -= _term;
 _MATELM1( 48 ,48)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + O1 <-> O2 ( o12 , o21 )*/
 _term =  o12 * ca [ 0] ;
 _MATELM1( 48 ,48)  += _term;
 _MATELM1( 80 +  0 ,48)  += _term;
 _MATELM1( 47 ,48)  -= _term;
 _term =  o12 * O1 ;
 _MATELM1( 48 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 47 ,80 +  0)  -= _term;
 _term =  o21 ;
 _MATELM1( 48 ,47)  -= _term;
 _MATELM1( 80 +  0 ,47)  -= _term;
 _MATELM1( 47 ,47)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + O2 <-> O3 ( o23 , o32 )*/
 _term =  o23 * ca [ 0] ;
 _MATELM1( 47 ,47)  += _term;
 _MATELM1( 80 +  0 ,47)  += _term;
 _MATELM1( 46 ,47)  -= _term;
 _term =  o23 * O2 ;
 _MATELM1( 47 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 46 ,80 +  0)  -= _term;
 _term =  o32 ;
 _MATELM1( 47 ,46)  -= _term;
 _MATELM1( 80 +  0 ,46)  -= _term;
 _MATELM1( 46 ,46)  += _term;
 /*REACTION*/
  /* ~ ca [ 0 ] + O3 <-> O4 ( o34 , o43 )*/
 _term =  o34 * ca [ 0] ;
 _MATELM1( 46 ,46)  += _term;
 _MATELM1( 80 +  0 ,46)  += _term;
 _MATELM1( 45 ,46)  -= _term;
 _term =  o34 * O3 ;
 _MATELM1( 46 ,80 +  0)  += _term;
 _MATELM1( 80 +  0 ,80 +  0)  += _term;
 _MATELM1( 45 ,80 +  0)  -= _term;
 _term =  o43 ;
 _MATELM1( 46 ,45)  -= _term;
 _MATELM1( 80 +  0 ,45)  -= _term;
 _MATELM1( 45 ,45)  += _term;
 /*REACTION*/
  /* ~ C0 <-> O0 ( f0 , b0 )*/
 _term =  f0 ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 49 ,4)  -= _term;
 _term =  b0 ;
 _MATELM1( 4 ,49)  -= _term;
 _MATELM1( 49 ,49)  += _term;
 /*REACTION*/
  /* ~ C1 <-> O1 ( f1 , b1 )*/
 _term =  f1 ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 48 ,3)  -= _term;
 _term =  b1 ;
 _MATELM1( 3 ,48)  -= _term;
 _MATELM1( 48 ,48)  += _term;
 /*REACTION*/
  /* ~ C2 <-> O2 ( f2 , b2 )*/
 _term =  f2 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 47 ,2)  -= _term;
 _term =  b2 ;
 _MATELM1( 2 ,47)  -= _term;
 _MATELM1( 47 ,47)  += _term;
 /*REACTION*/
  /* ~ C3 <-> O3 ( f3 , b3 )*/
 _term =  f3 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 46 ,1)  -= _term;
 _term =  b3 ;
 _MATELM1( 1 ,46)  -= _term;
 _MATELM1( 46 ,46)  += _term;
 /*REACTION*/
  /* ~ C4 <-> O4 ( f4 , b4 )*/
 _term =  f4 ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 45 ,0)  -= _term;
 _term =  b4 ;
 _MATELM1( 0 ,45)  -= _term;
 _MATELM1( 45 ,45)  += _term;
 /*REACTION*/
  cai = ca [ 0 ] ;
 ca1 = ca [ 1 ] ;
 ca2 = ca [ 2 ] ;
 ca3 = ca [ 3 ] ;
 ca4 = ca [ 4 ] ;
 ca5 = ca [ 5 ] ;
 ca9 = ca [ 9 ] ;
 mgi = mg [ 0 ] ;
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 140;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cao = _ion_cao;
  cai = _ion_cai;
  ica = _ion_ica;
  cai = _ion_cai;
     _ode_spec1 ();
  _ion_cai = cai;
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 140; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 static void _ode_synonym(int _cnt, double** _pp, Datum** _ppd) { 
 	int _i; 
	for (_i=0; _i < _cnt; ++_i) {_p = _pp[_i]; _ppvar = _ppd[_i];
 _ion_cai =  ca [ 0 ] ;
 }}
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse(&_cvsparseobj1, 140, _dlist1, _p, _ode_matsol1, &_coef1);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cao = _ion_cao;
  cai = _ion_cai;
  ica = _ion_ica;
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  BK_conc = BK_conc0;
  BK_ro = BK_ro0;
  C4 = C40;
  C3 = C30;
  C2 = C20;
  C1 = C10;
  C0 = C00;
 for (_i=0; _i<10; _i++) CB_ca_ca[_i] = CB_ca_ca0;
 for (_i=0; _i<10; _i++) CB_ca_s[_i] = CB_ca_s0;
 for (_i=0; _i<10; _i++) CB_f_ca[_i] = CB_f_ca0;
 for (_i=0; _i<10; _i++) CB[_i] = CB0;
  O4 = O40;
  O3 = O30;
  O2 = O20;
  O1 = O10;
  O0 = O00;
 for (_i=0; _i<10; _i++) PV_mg[_i] = PV_mg0;
 for (_i=0; _i<10; _i++) PV_ca[_i] = PV_ca0;
 for (_i=0; _i<10; _i++) PV[_i] = PV0;
 for (_i=0; _i<10; _i++) ca[_i] = ca0;
  dr_two = dr_two0;
 for (_i=0; _i<10; _i++) iCB_ca_ca[_i] = iCB_ca_ca0;
 for (_i=0; _i<10; _i++) iCB_ca_s[_i] = iCB_ca_s0;
 for (_i=0; _i<10; _i++) iCB_f_ca[_i] = iCB_f_ca0;
 for (_i=0; _i<10; _i++) iCB[_i] = iCB0;
 for (_i=0; _i<10; _i++) mg[_i] = mg0;
  pumpca = pumpca0;
  pump = pump0;
  rad_inner = rad_inner0;
  rad_outer = rad_outer0;
  vol_shell = vol_shell0;
 {
   if ( _zfactors_done  == 0.0 ) {
     _zfactors_done = 1.0 ;
     factors ( _threadargs_ ) ;
     }
   {int  _li ;for ( _li = 0 ; _li <= 10 - 1 ; _li ++ ) {
     ca [ _li ] = cainull ;
     mg [ _li ] = mginull ;
     CB [ _li ] = 0.8 * ssCB ( _threadargscomma_ kdf ( _threadargs_ ) , kds ( _threadargs_ ) ) ;
     CB_f_ca [ _li ] = 0.8 * ssCBfast ( _threadargscomma_ kdf ( _threadargs_ ) , kds ( _threadargs_ ) ) ;
     CB_ca_s [ _li ] = 0.8 * ssCBslow ( _threadargscomma_ kdf ( _threadargs_ ) , kds ( _threadargs_ ) ) ;
     CB_ca_ca [ _li ] = 0.8 * ssCBca ( _threadargscomma_ kdf ( _threadargs_ ) , kds ( _threadargs_ ) ) ;
     iCB [ _li ] = 0.2 * ssCB ( _threadargscomma_ kdf ( _threadargs_ ) , kds ( _threadargs_ ) ) ;
     iCB_f_ca [ _li ] = 0.2 * ssCBfast ( _threadargscomma_ kdf ( _threadargs_ ) , kds ( _threadargs_ ) ) ;
     iCB_ca_s [ _li ] = 0.2 * ssCBslow ( _threadargscomma_ kdf ( _threadargs_ ) , kds ( _threadargs_ ) ) ;
     iCB_ca_ca [ _li ] = 0.2 * ssCBca ( _threadargscomma_ kdf ( _threadargs_ ) , kds ( _threadargs_ ) ) ;
     PV [ _li ] = ssPV ( _threadargscomma_ kdc ( _threadargs_ ) , kdm ( _threadargs_ ) ) ;
     PV_ca [ _li ] = ssPVca ( _threadargscomma_ kdc ( _threadargs_ ) , kdm ( _threadargs_ ) ) ;
     PV_mg [ _li ] = ssPVmg ( _threadargscomma_ kdc ( _threadargs_ ) , kdm ( _threadargs_ ) ) ;
     } }
   parea = PI * diam ;
   ica = 0.0 ;
   ica_pmp = 0.0 ;
   pump = TotalPump ;
   pumpca = 0.0 ;
   BK_ro = 1e4 * 1e12 * BK_g / BK_sing_chan_g ;
   rad_outer = ( diam / 2.0 ) * 1e-6 ;
   dr_two = 0.2e-6 ;
   rad_inner = rad_outer - dr_two ;
   vol_shell = 1e3 * PI * ( rad_outer * rad_outer - rad_inner * rad_inner ) ;
   BK_conc = 1e3 * ( PI * BK_ro * 2.0 * rad_outer ) / ( N_A * vol_shell ) ;
   C0 = BK_conc * 0.98 ;
   C1 = BK_conc * 0.02 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cao = _ion_cao;
  cai = _ion_cai;
  ica = _ion_ica;
  cai = _ion_cai;
 initmodel();
  _ion_cai = cai;
  nrn_wrote_conc(_ca_sym, (&(_ion_cai)) - 1, _style_ca);
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cao = _ion_cao;
  cai = _ion_cai;
  ica = _ion_ica;
  cai = _ion_cai;
 { error = sparse(&_sparseobj1, 140, _slist1, _dlist1, _p, &t, dt, state,&_coef1, _linmat1);
 if(error){fprintf(stderr,"at line 226 in file cdp_AIS.mod:\n:	ica = ica_pmp\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 140; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } {
   }
  _ion_cai = cai;
}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = C4_columnindex;  _dlist1[0] = DC4_columnindex;
 _slist1[1] = C3_columnindex;  _dlist1[1] = DC3_columnindex;
 _slist1[2] = C2_columnindex;  _dlist1[2] = DC2_columnindex;
 _slist1[3] = C1_columnindex;  _dlist1[3] = DC1_columnindex;
 _slist1[4] = C0_columnindex;  _dlist1[4] = DC0_columnindex;
 for(_i=0;_i<10;_i++){_slist1[5+_i] = CB_ca_ca_columnindex + _i;  _dlist1[5+_i] = DCB_ca_ca_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[15+_i] = CB_ca_s_columnindex + _i;  _dlist1[15+_i] = DCB_ca_s_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[25+_i] = CB_f_ca_columnindex + _i;  _dlist1[25+_i] = DCB_f_ca_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[35+_i] = CB_columnindex + _i;  _dlist1[35+_i] = DCB_columnindex + _i;}
 _slist1[45] = O4_columnindex;  _dlist1[45] = DO4_columnindex;
 _slist1[46] = O3_columnindex;  _dlist1[46] = DO3_columnindex;
 _slist1[47] = O2_columnindex;  _dlist1[47] = DO2_columnindex;
 _slist1[48] = O1_columnindex;  _dlist1[48] = DO1_columnindex;
 _slist1[49] = O0_columnindex;  _dlist1[49] = DO0_columnindex;
 for(_i=0;_i<10;_i++){_slist1[50+_i] = PV_mg_columnindex + _i;  _dlist1[50+_i] = DPV_mg_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[60+_i] = PV_ca_columnindex + _i;  _dlist1[60+_i] = DPV_ca_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[70+_i] = PV_columnindex + _i;  _dlist1[70+_i] = DPV_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[80+_i] = ca_columnindex + _i;  _dlist1[80+_i] = Dca_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[90+_i] = iCB_ca_ca_columnindex + _i;  _dlist1[90+_i] = DiCB_ca_ca_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[100+_i] = iCB_ca_s_columnindex + _i;  _dlist1[100+_i] = DiCB_ca_s_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[110+_i] = iCB_f_ca_columnindex + _i;  _dlist1[110+_i] = DiCB_f_ca_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[120+_i] = iCB_columnindex + _i;  _dlist1[120+_i] = DiCB_columnindex + _i;}
 for(_i=0;_i<10;_i++){_slist1[130+_i] = mg_columnindex + _i;  _dlist1[130+_i] = Dmg_columnindex + _i;}
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/iain/GIT/BK_unitary_g/NEURON/simple_model_50pS/mod/cdp_AIS.mod";
static const char* nmodl_file_text = 
  ": Calcium ion accumulation with radial and longitudinal diffusion and pump\n"
  "\n"
  "NEURON {\n"
  "  SUFFIX cdpAIS\n"
  "  USEION ca READ cao, cai, ica WRITE cai\n"
  "  RANGE ica_pmp,ca1,ca2,ca3,ca4,ca5,ca9\n"
  ":RANGE pump_0  \n"
  "GLOBAL vrat, TotalPump\n"
  "    : vrat must be GLOBAL--see INITIAL block\n"
  "    : however TotalBuffer and TotalPump may be RANGE\n"
  ":    THREADSAFE\n"
  "}\n"
  "\n"
  "DEFINE Nannuli 10\n"
  "\n"
  "UNITS {\n"
  "	(mol)   = (1)\n"
  "	(molar) = (1/liter)\n"
  "	(mM)    = (millimolar)\n"
  "	(um)    = (micron)\n"
  "	(mA)    = (milliamp)\n"
  "	FARADAY = (faraday)  (10000 coulomb)\n"
  "	PI      = (pi)       (1)\n"
  "\n"
  "	: BK stuff\n"
  "	(mV) = (millivolt)\n"
  "    (S) = (siemens)\n"
  "    R = (k-mole) (joule/degC)\n"
  "}\n"
  "\n"
  "CONSTANT {\n"
  "    q10 = 2\n"
  "    : Avogardos constant.\n"
  "    N_A = 6.02214076e23 \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    v\n"
  "	celsius =34     (degC)\n"
  "        \n"
  "	:cainull =2.5e-4 (mM)\n"
  "	cainull = 45e-6 (mM)\n"
  "        mginull =.59    (mM)\n"
  "\n"
  "        DCa     = .233  (um2/ms)\n"
  "	Dbtc 	= 0.007 (um2/ms)\n"
  "       Ddmnpe = 0.08	(um2/ms)\n"
  "	\n"
  "	Dcbd1   = .028  (um2/ms)\n"
  "        Dcbd2   = 0     (um2/ms)\n"
  "        Dpar    = .043  (um2/ms)\n"
  "\n"
  ":	values for benzothiazole coumarin (BTC)\n"
  ":	BTCnull = 0	(mM)\n"
  ":	b1 = 5.33	(/ms mM)\n"
  ":	b2 = 0.08	(/ms)\n"
  "\n"
  ":	values for caged compound DMNPE-4\n"
  ":	DMNPEnull = 0	(mM)\n"
  ":	c1 = 5.63	(/ms mM)\n"
  ":	c2 = 0.107e-3	(/ms)\n"
  "\n"
  ":       values for Calbindin (2 high and 2 low affinity binding sites)\n"
  "\n"
  ":        CBnull=	.16             (mM)\n"
  "         CBnull=	.08             (mM)       \n"
  "        nf1   =43.5           (/ms mM)\n"
  "        nf2   =3.58e-2        (/ms)\n"
  "        ns1   =5.5            (/ms mM)\n"
  "        ns2   =0.26e-2        (/ms)\n"
  "\n"
  ":       values for Parvalbumin\n"
  "\n"
  ":        PVnull  = 0.08           (mM)\n"
  "        PVnull  = 0.04           (mM)\n"
  "        m1    = 1.07e2        (/ms mM)\n"
  "        m2    = 9.5e-4                (/ms)\n"
  "        p1    = 0.8           (/ms mM)\n"
  "        p2    = 2.5e-2                (/ms)        \n"
  "    \n"
  "  	kpmp1    = 3e3       (/mM-ms)\n"
  "  	kpmp2    = 1.75e1   (/ms)\n"
  "  	kpmp3    = 7.255e1  (/ms)\n"
  "  : to eliminate pump, set TotalPump to 0 in hoc\n"
  "	TotalPump = 1e-15\n"
  ":	TotalPump = 1e-8\n"
  "	beta  = 1(1)           :introducing beta to take care of other ER mechanisms(SERCA and leak channel density)\n"
  "    vmax =0.1\n"
  "    :Kp = 0.27e-3 (mM)\n"
  "    Kp = 2.7e-3 (mM)\n"
  "\n"
  "\n"
  "	: BK stuff\n"
  "    \n"
  "    Qo = 0.73\n"
  "    Qc = -0.58\n"
  "    \n"
  "    k1 = 1.0e3 (/mM)\n"
  "    onoffrate = 1 (/ms)\n"
  "    \n"
  "    L0 = 1576\n"
  "    Kc = 11.917e-3 (mM)\n"
  "    Ko = 1.065e-3 (mM)\n"
  "    \n"
  "    pf0 = 5.5e-3  (/ms)\n"
  "    pf1 = 8e-3  (/ms)\n"
  "    pf2 = 2e-3   (/ms)\n"
  "    pf3 = 884e-3  (/ms)\n"
  "    pf4 = 900e-3  (/ms)\n"
  "    \n"
  "    pb0 = 8669e-3 (/ms)\n"
  "    pb1 = 1127e-3 (/ms)\n"
  "    pb2 = 25.2e-3  (/ms)\n"
  "    pb3 = 1013e-3  (/ms)\n"
  "    pb4 = 125.7e-3  (/ms)\n"
  "\n"
  "    BK_sing_chan_g = 50\n"
  "    BK_g = 24 : X1 = 6, X2 = 12, X3 = 18, X4 = 24  \n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	diam      (um)\n"
  "	ica       (mA/cm2)\n"
  "	ica_pmp   (mA/cm2)\n"
  "	ibg     :background calcium current\n"
  ":	ica_pmp_last   (mA/cm2)\n"
  "	parea     (um)     : pump area per unit length\n"
  "	cai       (mM)\n"
  "	ca1\n"
  "	ca2\n"
  "	ca3\n"
  "	ca4\n"
  "	ca5\n"
  "	ca9\n"
  "	mgi	(mM)	\n"
  "	vrat[Nannuli]  (1) : dimensionless\n"
  "                     : numeric value of vrat[i] equals the volume \n"
  "                     : of annulus i of a 1um diameter cylinder\n"
  "                     : multiply by diam^2 to get volume per um length\n"
  "\n"
  "\n"
  "	: BK stuff\n"
  "\n"
  "    c01    (/ms)\n"
  "    c12    (/ms)\n"
  "    c23    (/ms)\n"
  "    c34    (/ms)\n"
  "    o01    (/ms)\n"
  "    o12    (/ms)\n"
  "    o23    (/ms)\n"
  "    o34    (/ms)\n"
  "    f0     (/ms)\n"
  "    f1     (/ms)\n"
  "    f2     (/ms)\n"
  "    f3     (/ms)\n"
  "    f4     (/ms)\n"
  "\n"
  "    c10    (/ms)\n"
  "    c21    (/ms)\n"
  "    c32    (/ms)\n"
  "    c43    (/ms)\n"
  "    o10    (/ms)\n"
  "    o21    (/ms)\n"
  "    o32    (/ms)\n"
  "    o43    (/ms)\n"
  "    b0     (/ms)\n"
  "    b1     (/ms)\n"
  "    b2     (/ms)\n"
  "    b3     (/ms)\n"
  "    b4     (/ms)\n"
  "\n"
  "}\n"
  "\n"
  "CONSTANT { cao = 2	(mM) }\n"
  "\n"
  "STATE {\n"
  "	: ca[0] is equivalent to cai\n"
  "	: ca[] are very small, so specify absolute tolerance\n"
  "	: let it be ~1.5 - 2 orders of magnitude smaller than baseline level\n"
  "	ca[Nannuli]		(mM)\n"
  "	mg[Nannuli]		(mM)	<1e-7>\n"
  "\n"
  "        CB[Nannuli]		(mM)\n"
  "        CB_f_ca[Nannuli]	(mM)\n"
  "        CB_ca_s[Nannuli]	(mM)\n"
  "        CB_ca_ca[Nannuli]	(mM)\n"
  "\n"
  "        iCB[Nannuli]		(mM)\n"
  "        iCB_f_ca[Nannuli]	(mM)\n"
  "        iCB_ca_s[Nannuli]	(mM)\n"
  "        iCB_ca_ca[Nannuli]	(mM)\n"
  "\n"
  "        PV[Nannuli]		(mM)\n"
  "        PV_ca[Nannuli]		(mM)\n"
  "        PV_mg[Nannuli]		(mM)\n"
  "	\n"
  "	pump			(mol/cm2) <1e-15>\n"
  "	pumpca			(mol/cm2) <1e-15>\n"
  "\n"
  "\n"
  "	C0 \n"
  "    C1 \n"
  "    C2 \n"
  "    C3 \n"
  "    C4 \n"
  "    O0 \n"
  "    O1 \n"
  "    O2 \n"
  "    O3 \n"
  "    O4 \n"
  "\n"
  "    : Number of BK channels per square micron of membrane\n"
  "    BK_ro \n"
  "    : Effective concentration of BK channels in the outer shell\n"
  "    BK_conc \n"
  "\n"
  "    dr_two\n"
  "    rad_outer\n"
  "    rad_inner\n"
  "    vol_shell\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD sparse\n"
  ":	ica_pmp_last = ica_pmp\n"
  ":	ica = ica_pmp\n"
  "}\n"
  "\n"
  "LOCAL factors_done\n"
  "\n"
  "INITIAL {\n"
  "	if (factors_done == 0) {  : flag becomes 1 in the first segment\n"
  "		factors_done = 1       :   all subsequent segments will have\n"
  "		factors()              :   vrat = 0 unless vrat is GLOBAL\n"
  "	}\n"
  "	FROM i=0 TO Nannuli-1 {\n"
  "		ca[i] = cainull\n"
  "		mg[i] = mginull\n"
  "\n"
  "\n"
  "		CB[i] = 0.8*ssCB( kdf(), kds())   \n"
  "	        CB_f_ca[i] = 0.8*ssCBfast( kdf(), kds())\n"
  "       	 	CB_ca_s[i] = 0.8*ssCBslow( kdf(), kds())\n"
  "        	CB_ca_ca[i] = 0.8*ssCBca( kdf(), kds())\n"
  "\n"
  "        	iCB[i] = 0.2*ssCB( kdf(), kds())\n"
  "        	iCB_f_ca[i] = 0.2*ssCBfast( kdf(), kds())\n"
  "        	iCB_ca_s[i] = 0.2*ssCBslow( kdf(), kds())\n"
  "        	iCB_ca_ca[i] = 0.2*ssCBca(kdf(), kds())\n"
  "\n"
  "        	PV[i] = ssPV( kdc(), kdm())\n"
  "        	PV_ca[i] = ssPVca(kdc(), kdm())\n"
  "        	PV_mg[i] = ssPVmg(kdc(), kdm())\n"
  "\n"
  "		\n"
  "	}\n"
  "  	parea = PI*diam\n"
  "	ica = 0\n"
  "	ica_pmp = 0\n"
  ":	ica_pmp_last = 0\n"
  "	pump = TotalPump\n"
  "	pumpca = 0\n"
  "\n"
  "	BK_ro = 1e4*1e12*BK_g/BK_sing_chan_g  : / m2 \n"
  "	rad_outer = (diam/2)*1e-6 : m\n"
  "	dr_two = 0.2e-6 : m\n"
  "	rad_inner = rad_outer-dr_two : m\n"
  "	vol_shell = 1e3 * PI * (rad_outer*rad_outer - rad_inner*rad_inner) : L / unit length\n"
  "\n"
  "	BK_conc = 1e3*(PI * BK_ro * 2 * rad_outer) / (N_A * vol_shell) : mM			  \n"
  "	\n"
  "	: BK channel almost entirely in C0 state initially at background Ca \n"
  "	C0=BK_conc*0.98\n"
  "	C1=BK_conc*0.02\n"
  "	: The rest are insignificant initially \n"
  "}\n"
  "\n"
  "LOCAL radii[Nannuli]\n"
  "LOCAL frat[Nannuli]  : scales the rate constants for model geometry\n"
  "\n"
  "PROCEDURE factors() {\n"
  "	LOCAL r, dr2, dr3\n"
  "	r = 1/2                : starts at edge (half diam)\n"
  "	dr2 = 0.2/diam	       : full thickness of outermost annulus, here 100nm\n"
  "\n"
  "	dr3 = (r-dr2)/(Nannuli-1)	:other shells thickness\n"
  "                         : half thickness of all other annuli\n"
  "\n"
  "        radii[0] = r\n"
  "	radii[1] = r - dr2\n"
  "        FROM i=2 TO Nannuli-1 {\n"
  "		radii[i] = radii[i-1]- dr3\n"
  "	printf(\"%f\\n\",radii[i])\n"
  "	}\n"
  "\n"
  "	vrat[0] = 0\n"
  "	frat[0] = 2*r\n"
  "	\n"
  "	FROM i=0 TO Nannuli-2 {\n"
  "		vrat[i] = PI*((radii[i]*radii[i])-(radii[i+1]*radii[i+1]))\n"
  "  	}\n"
  "	vrat[Nannuli-1] = PI*radii[Nannuli-1]*radii[Nannuli-1]\n"
  "	FROM i=1 TO Nannuli-1 {\n"
  "		if (i==1) {\n"
  "			frat[i] = 2*PI*radii[i]/(dr2+(dr3/2))\n"
  "		} else if (i>1&&i<(Nannuli-1)) { \n"
  "			frat[i] = 2*PI*radii[i]/dr3\n"
  "		} else if (i==(Nannuli-1)) {\n"
  "			frat[i] = 2*PI*radii[i]/((dr3/2)+radii[i])\n"
  "		}\n"
  "	}\n"
  "}\n"
  "\n"
  "LOCAL dsq, dsqvol  : \n"
  "                   :   or use in COMPARTMENT statement\n"
  "\n"
  "KINETIC state {\n"
  "  COMPARTMENT i, diam*diam*vrat[i] {ca mg CB CB_f_ca CB_ca_s CB_ca_ca iCB iCB_f_ca iCB_ca_s iCB_ca_ca PV PV_ca PV_mg}\n"
  "  COMPARTMENT (1e10)*parea {pump pumpca}\n"
  "	:pump\n"
  ":	~ ca[0] + pump <-> pumpca  (kpmp1*parea*(1e10), kpmp2*parea*(1e10))\n"
  ":	~ pumpca <-> pump   (kpmp3*parea*(1e10), 0)\n"
  ":  	CONSERVE pump + pumpca = TotalPump * parea * (1e10)\n"
  "\n"
  ":	ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea\n"
  "\n"
  "	: all currents except pump\n"
  "	: ica is Ca efflux\n"
  "	~ ca[0] << (-ica*PI*diam/(2*FARADAY))\n"
  "\n"
  "	:RADIAL DIFFUSION OF ca, mg and mobile buffers\n"
  "\n"
  "    FROM i=0 TO Nannuli-1 {\n"
  "    ~ ca[i] << (-beta*vmax*vrat[i]*ca[i] / (ca[i] + kpmp2/kpmp1))\n"
  "\n"
  "    }\n"
  "    \n"
  "	FROM i=0 TO Nannuli-2 {\n"
  "		~ ca[i] <-> ca[i+1]	(DCa*frat[i+1], DCa*frat[i+1])\n"
  "		~ mg[i] <-> mg[i+1]	(DCa*frat[i+1], DCa*frat[i+1])\n"
  "		~ CB[i] <-> CB[i+1]	(Dcbd1*frat[i+1], Dcbd1*frat[i+1])\n"
  "		~ CB_f_ca[i] <-> CB_f_ca[i+1]	(Dcbd1*frat[i+1], Dcbd1*frat[i+1])\n"
  "		~ CB_ca_s[i] <-> CB_ca_s[i+1]	(Dcbd1*frat[i+1], Dcbd1*frat[i+1])\n"
  "		~ CB_ca_ca[i] <-> CB_ca_ca[i+1]	(Dcbd1*frat[i+1], Dcbd1*frat[i+1])\n"
  "		~ PV[i] <-> PV[i+1]	(Dpar*frat[i+1], Dpar*frat[i+1])\n"
  "		~ PV_ca[i] <-> PV_ca[i+1]	(Dpar*frat[i+1], Dpar*frat[i+1])\n"
  "		~ PV_mg[i] <-> PV_mg[i+1] 	(Dpar*frat[i+1], Dpar*frat[i+1])\n"
  "	}\n"
  "	dsq = diam*diam\n"
  "	FROM i=0 TO Nannuli-1 {\n"
  "		dsqvol = dsq*vrat[i]\n"
  "		:Calbindin	\n"
  "		~ ca[i] + CB[i] <-> CB_ca_s[i] (nf1*dsqvol, nf2*dsqvol)\n"
  "	       	~ ca[i] + CB[i] <-> CB_f_ca[i] (ns1*dsqvol, ns2*dsqvol)\n"
  "        	~ ca[i] + CB_f_ca[i] <-> CB_ca_ca[i] (nf1*dsqvol, nf2*dsqvol)\n"
  "        	~ ca[i] + CB_ca_s[i] <-> CB_ca_ca[i] (ns1*dsqvol, ns2*dsqvol)\n"
  "\n"
  "        	~ ca[i] + iCB[i] <-> iCB_ca_s[i] (nf1*dsqvol, nf2*dsqvol)\n"
  "        	~ ca[i] + iCB[i] <-> iCB_f_ca[i] (ns1*dsqvol, ns2*dsqvol)\n"
  "        	~ ca[i] + iCB_f_ca[i] <-> iCB_ca_ca[i] (nf1*dsqvol, nf2*dsqvol)\n"
  "        	~ ca[i] + iCB_ca_s[i] <-> iCB_ca_ca[i] (ns1*dsqvol, ns2*dsqvol)\n"
  "\n"
  "\n"
  "		:Paravalbumin\n"
  "        	~ ca[i] + PV[i] <-> PV_ca[i] (m1*dsqvol, m2*dsqvol)\n"
  "        	~ mg[i] + PV[i] <-> PV_mg[i] (p1*dsqvol, p2*dsqvol)\n"
  "	}\n"
  "\n"
  "\n"
  "	rates(v, cai)\n"
  "    ~ ca[0] + C0 <-> C1      (c01,c10)\n"
  "    ~ ca[0] + C1 <-> C2      (c12,c21)\n"
  "    ~ ca[0] + C2 <-> C3      (c23,c32)\n"
  "    ~ ca[0] + C3 <-> C4      (c34,c43)\n"
  "    ~ ca[0] + O0 <-> O1      (o01,o10)\n"
  "    ~ ca[0] + O1 <-> O2      (o12,o21)\n"
  "    ~ ca[0] + O2 <-> O3      (o23,o32)\n"
  "    ~ ca[0] + O3 <-> O4      (o34,o43)\n"
  "    ~ C0 <-> O0      (f0, b0)\n"
  "    ~ C1 <-> O1      (f1, b1)\n"
  "    ~ C2 <-> O2      (f2, b2)\n"
  "    ~ C3 <-> O3      (f3, b3)\n"
  "    ~ C4 <-> O4      (f4, b4)\n"
  "\n"
  "\n"
  "  	cai = ca[0]\n"
  "  	ca1 = ca[1]\n"
  "  	ca2 = ca[2]\n"
  "  	ca3 = ca[3]\n"
  "  	ca4 = ca[4]\n"
  "  	ca5 = ca[5]\n"
  "  	ca9 = ca[9]\n"
  ":  	if (ca[0]<cainull){: keep minimum\n"
  ":	   cai=cainull }\n"
  "	mgi = mg[0]\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION ssCB( kdf(), kds()) (mM) {\n"
  "	ssCB = CBnull/(1+kdf()+kds()+(kdf()*kds()))\n"
  "}\n"
  "FUNCTION ssCBfast( kdf(), kds()) (mM) {\n"
  "	ssCBfast = (CBnull*kds())/(1+kdf()+kds()+(kdf()*kds()))\n"
  "}\n"
  "FUNCTION ssCBslow( kdf(), kds()) (mM) {\n"
  "	ssCBslow = (CBnull*kdf())/(1+kdf()+kds()+(kdf()*kds()))\n"
  "}\n"
  "FUNCTION ssCBca(kdf(), kds()) (mM) {\n"
  "	ssCBca = (CBnull*kdf()*kds())/(1+kdf()+kds()+(kdf()*kds()))\n"
  "}\n"
  "FUNCTION kdf() (1) {\n"
  "	kdf = (cainull*nf1)/nf2\n"
  "}\n"
  "FUNCTION kds() (1) {\n"
  "	kds = (cainull*ns1)/ns2\n"
  "}\n"
  "FUNCTION kdc() (1) {\n"
  "	kdc = (cainull*m1)/m2\n"
  "}\n"
  "FUNCTION kdm() (1) {\n"
  "	kdm = (mginull*p1)/p2\n"
  "}\n"
  "FUNCTION ssPV( kdc(), kdm()) (mM) {\n"
  "	ssPV = PVnull/(1+kdc()+kdm())\n"
  "}\n"
  "FUNCTION ssPVca( kdc(), kdm()) (mM) {\n"
  "	ssPVca = (PVnull*kdc)/(1+kdc+kdm)\n"
  "}\n"
  "FUNCTION ssPVmg( kdc(), kdm()) (mM) {\n"
  "	ssPVmg = (PVnull*kdm())/(1+kdc()+kdm())\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION u (x, th) {\n"
  "  	if (x<th) {\n"
  "    		u = 1\n"
  "  	} else {\n"
  "    		u = 0\n"
  "  	}\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE rates(v(mV), ca (mM)) { \n"
  "    LOCAL qt, alpha, beta\n"
  "    \n"
  "    qt = q10^((celsius-25 (degC))/10 (degC))\n"
  "    \n"
  "    c01 = 4*k1*onoffrate*qt\n"
  "    c12 = 3*k1*onoffrate*qt\n"
  "    c23 = 2*k1*onoffrate*qt\n"
  "    c34 = 1*k1*onoffrate*qt\n"
  "    o01 = 4*k1*onoffrate*qt\n"
  "    o12 = 3*k1*onoffrate*qt\n"
  "    o23 = 2*k1*onoffrate*qt\n"
  "    o34 = 1*k1*onoffrate*qt\n"
  "    \n"
  "    c10 = 1*Kc*k1*onoffrate*qt\n"
  "    c21 = 2*Kc*k1*onoffrate*qt\n"
  "    c32 = 3*Kc*k1*onoffrate*qt\n"
  "    c43 = 4*Kc*k1*onoffrate*qt\n"
  "    o10 = 1*Ko*k1*onoffrate*qt\n"
  "    o21 = 2*Ko*k1*onoffrate*qt\n"
  "    o32 = 3*Ko*k1*onoffrate*qt\n"
  "    o43 = 4*Ko*k1*onoffrate*qt\n"
  "    \n"
  "    alpha = exp(Qo*FARADAY*10*v/R/(273.15 + celsius))\n"
  "    beta  = exp(Qc*FARADAY*10*v/R/(273.15 + celsius))\n"
  "    \n"
  "    f0  = pf0*alpha*qt\n"
  "    f1  = pf1*alpha*qt\n"
  "    f2  = pf2*alpha*qt\n"
  "    f3  = pf3*alpha*qt\n"
  "    f4  = pf4*alpha*qt\n"
  "    \n"
  "    b0  = pb0*beta*qt\n"
  "    b1  = pb1*beta*qt\n"
  "    b2  = pb2*beta*qt\n"
  "    b3  = pb3*beta*qt\n"
  "    b4  = pb4*beta*qt\n"
  "\n"
  "}\n"
  ;
#endif
