/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__CaT3_1
#define _nrn_initial _nrn_initial__CaT3_1
#define nrn_cur _nrn_cur__CaT3_1
#define _nrn_current _nrn_current__CaT3_1
#define nrn_jacob _nrn_jacob__CaT3_1
#define nrn_state _nrn_state__CaT3_1
#define _net_receive _net_receive__CaT3_1 
#define castate castate__CaT3_1 
#define evaluate_fct evaluate_fct__CaT3_1 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg(int);
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define pcabar _p[0]
#define pcabar_columnindex 0
#define iCa _p[1]
#define iCa_columnindex 1
#define g _p[2]
#define g_columnindex 2
#define minf _p[3]
#define minf_columnindex 3
#define taum _p[4]
#define taum_columnindex 4
#define hinf _p[5]
#define hinf_columnindex 5
#define tauh _p[6]
#define tauh_columnindex 6
#define m _p[7]
#define m_columnindex 7
#define h _p[8]
#define h_columnindex 8
#define cai _p[9]
#define cai_columnindex 9
#define cao _p[10]
#define cao_columnindex 10
#define Dm _p[11]
#define Dm_columnindex 11
#define Dh _p[12]
#define Dh_columnindex 12
#define qt _p[13]
#define qt_columnindex 13
#define T _p[14]
#define T_columnindex 14
#define E _p[15]
#define E_columnindex 15
#define zeta _p[16]
#define zeta_columnindex 16
#define v _p[17]
#define v_columnindex 17
#define _g _p[18]
#define _g_columnindex 18
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_iCa	*_ppvar[2]._pval
#define _ion_diCadv	*_ppvar[3]._pval
 
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_evaluate_fct(void);
 static void _hoc_ghk2(void);
 static void _hoc_kelvinfkt(void);
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
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_CaT3_1", _hoc_setdata,
 "evaluate_fct_CaT3_1", _hoc_evaluate_fct,
 "ghk2_CaT3_1", _hoc_ghk2,
 "kelvinfkt_CaT3_1", _hoc_kelvinfkt,
 0, 0
};
#define ghk2 ghk2_CaT3_1
#define kelvinfkt kelvinfkt_CaT3_1
 extern double ghk2( _threadargsprotocomma_ double , double , double , double );
 extern double kelvinfkt( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define A_tau_h A_tau_h_CaT3_1
 double A_tau_h = 0.17746;
#define A_tau_m A_tau_m_CaT3_1
 double A_tau_m = -2.3199;
#define B_tau_h B_tau_h_CaT3_1
 double B_tau_h = 0.13402;
#define B_tau_m B_tau_m_CaT3_1
 double B_tau_m = 2.5712;
#define C_tau_h C_tau_h_CaT3_1
 double C_tau_h = 0.0076;
#define C_tau_m C_tau_m_CaT3_1
 double C_tau_m = 1.2757;
#define eca eca_CaT3_1
 double eca = 0;
#define k_tau_h2 k_tau_h2_CaT3_1
 double k_tau_h2 = -5.5845;
#define k_tau_h1 k_tau_h1_CaT3_1
 double k_tau_h1 = 6.2692;
#define k_tau_m2 k_tau_m2_CaT3_1
 double k_tau_m2 = 9.6306;
#define k_tau_m1 k_tau_m1_CaT3_1
 double k_tau_m1 = 30.655;
#define k_h_inf k_h_inf_CaT3_1
 double k_h_inf = 6.4635;
#define k_m_inf k_m_inf_CaT3_1
 double k_m_inf = -4.7056;
#define v0_tau_h2 v0_tau_h2_CaT3_1
 double v0_tau_h2 = -101.436;
#define v0_tau_h1 v0_tau_h1_CaT3_1
 double v0_tau_h1 = -58.535;
#define v0_tau_m2 v0_tau_m2_CaT3_1
 double v0_tau_m2 = -28.386;
#define v0_tau_m1 v0_tau_m1_CaT3_1
 double v0_tau_m1 = -48.048;
#define vshift vshift_CaT3_1
 double vshift = 0;
#define v0_h_inf v0_h_inf_CaT3_1
 double v0_h_inf = -75.118;
#define v0_m_inf v0_m_inf_CaT3_1
 double v0_m_inf = -42.206;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "eca_CaT3_1", "mV",
 "v0_m_inf_CaT3_1", "mV",
 "v0_h_inf_CaT3_1", "mV",
 "k_m_inf_CaT3_1", "mV",
 "k_h_inf_CaT3_1", "mV",
 "v0_tau_m1_CaT3_1", "mV",
 "k_tau_m1_CaT3_1", "mV",
 "v0_tau_m2_CaT3_1", "mV",
 "k_tau_m2_CaT3_1", "mV",
 "v0_tau_h1_CaT3_1", "mV",
 "k_tau_h1_CaT3_1", "mV",
 "k_tau_h2_CaT3_1", "mV",
 "pcabar_CaT3_1", "cm/s",
 "iCa_CaT3_1", "mA/cm2",
 "g_CaT3_1", "coulombs/cm3",
 "taum_CaT3_1", "ms",
 "tauh_CaT3_1", "ms",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "eca_CaT3_1", &eca_CaT3_1,
 "v0_m_inf_CaT3_1", &v0_m_inf_CaT3_1,
 "v0_h_inf_CaT3_1", &v0_h_inf_CaT3_1,
 "vshift_CaT3_1", &vshift_CaT3_1,
 "k_m_inf_CaT3_1", &k_m_inf_CaT3_1,
 "k_h_inf_CaT3_1", &k_h_inf_CaT3_1,
 "C_tau_m_CaT3_1", &C_tau_m_CaT3_1,
 "A_tau_m_CaT3_1", &A_tau_m_CaT3_1,
 "v0_tau_m1_CaT3_1", &v0_tau_m1_CaT3_1,
 "k_tau_m1_CaT3_1", &k_tau_m1_CaT3_1,
 "B_tau_m_CaT3_1", &B_tau_m_CaT3_1,
 "v0_tau_m2_CaT3_1", &v0_tau_m2_CaT3_1,
 "k_tau_m2_CaT3_1", &k_tau_m2_CaT3_1,
 "C_tau_h_CaT3_1", &C_tau_h_CaT3_1,
 "A_tau_h_CaT3_1", &A_tau_h_CaT3_1,
 "v0_tau_h1_CaT3_1", &v0_tau_h1_CaT3_1,
 "k_tau_h1_CaT3_1", &k_tau_h1_CaT3_1,
 "B_tau_h_CaT3_1", &B_tau_h_CaT3_1,
 "v0_tau_h2_CaT3_1", &v0_tau_h2_CaT3_1,
 "k_tau_h2_CaT3_1", &k_tau_h2_CaT3_1,
 0,0
};
 static DoubVec hoc_vdoub[] = {
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
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"CaT3_1",
 "pcabar_CaT3_1",
 0,
 "iCa_CaT3_1",
 "g_CaT3_1",
 "minf_CaT3_1",
 "taum_CaT3_1",
 "hinf_CaT3_1",
 "tauh_CaT3_1",
 0,
 "m_CaT3_1",
 "h_CaT3_1",
 0,
 0};
 static Symbol* _ca_sym;
 static Symbol* _Ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 19, _prop);
 	/*initialize range parameters*/
 	pcabar = 0.00025;
 	_prop->param = _p;
 	_prop->param_size = 19;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 prop_ion = need_memb(_Ca_sym);
 	_ppvar[2]._pval = &prop_ion->param[3]; /* iCa */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_diCadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _CaT_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	ion_reg("Ca", 2.0);
 	_ca_sym = hoc_lookup("ca_ion");
 	_Ca_sym = hoc_lookup("Ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 19, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "Ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "Ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 CaT3_1 /Users/iain/GIT/BK_unitary_g/NEURON/simple_model_250pS/mod/CaT.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double q10 = 1.0913;
 static double F = 9.6485e4;
 static double R = 8.3145;
static int _reset;
static char *modelname = "Low threshold calcium current Cerebellum Purkinje Cell Model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int evaluate_fct(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int castate(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   evaluate_fct ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / taum ;
   Dh = ( hinf - h ) / tauh ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 evaluate_fct ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taum )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tauh )) ;
  return 0;
}
 /*END CVODE*/
 static int castate (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   evaluate_fct ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taum)))*(- ( ( ( minf ) ) / taum ) / ( ( ( ( - 1.0 ) ) ) / taum ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tauh)))*(- ( ( ( hinf ) ) / tauh ) / ( ( ( ( - 1.0 ) ) ) / tauh ) - h) ;
   }
  return 0;
}
 
double ghk2 ( _threadargsprotocomma_ double _lv , double _lci , double _lco , double _lz ) {
   double _lghk2;
 E = ( 1e-3 ) * _lv ;
   zeta = ( _lz * F * E ) / ( R * T ) ;
   if ( fabs ( 1.0 - exp ( - zeta ) ) < 1e-6 ) {
     _lghk2 = ( 1e-6 ) * ( _lz * F ) * ( _lci - _lco * exp ( - zeta ) ) * ( 1.0 + zeta / 2.0 ) ;
     }
   else {
     _lghk2 = ( 1e-6 ) * ( _lz * zeta * F ) * ( _lci - _lco * exp ( - zeta ) ) / ( 1.0 - exp ( - zeta ) ) ;
     }
   
return _lghk2;
 }
 
static void _hoc_ghk2(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  ghk2 ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
static int  evaluate_fct ( _threadargsprotocomma_ double _lv ) {
   minf = 1.0 / pow( ( 1.0 + exp ( ( _lv - v0_m_inf - vshift ) / k_m_inf ) ) , ( 1.0 / 3.0 ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv - v0_h_inf - vshift ) / k_h_inf ) ) ;
   taum = 1.0 / ( C_tau_m + A_tau_m / ( 1.0 + exp ( ( v0_tau_m1 - _lv - vshift ) / k_tau_m1 ) ) + B_tau_m / ( 1.0 + exp ( ( v0_tau_m2 - _lv - vshift ) / k_tau_m2 ) ) ) / qt ;
   tauh = 1.0 / ( C_tau_h + A_tau_h / ( 1.0 + exp ( ( v0_tau_h1 - _lv - vshift ) / k_tau_h1 ) ) + B_tau_h / ( 1.0 + exp ( ( v0_tau_h2 - _lv - vshift ) / k_tau_h2 ) ) ) / qt ;
   g = ghk2 ( _threadargscomma_ _lv - vshift , cai , cao , 2.0 ) ;
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 evaluate_fct ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double kelvinfkt ( _threadargsprotocomma_ double _lt ) {
   double _lkelvinfkt;
 _lkelvinfkt = 273.19 + _lt ;
   
return _lkelvinfkt;
 }
 
static void _hoc_kelvinfkt(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  kelvinfkt ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_Ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_Ca_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   T = kelvinfkt ( _threadargscomma_ celsius ) ;
   qt = pow( q10 , ( ( celsius - 32.0 ) / 10.0 ) ) ;
   evaluate_fct ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   iCa = ( 1e3 ) * pcabar * m * m * m * h * g ;
   }
 _current += iCa;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _diCa;
  _diCa = iCa;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_diCadv += (_diCa - iCa)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_iCa += iCa ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
  cai = _ion_cai;
  cao = _ion_cao;
 {   castate(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = m_columnindex;  _dlist1[0] = Dm_columnindex;
 _slist1[1] = h_columnindex;  _dlist1[1] = Dh_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/iain/GIT/BK_unitary_g/NEURON/simple_model_250pS/mod/CaT.mod";
static const char* nmodl_file_text = 
  "TITLE Low threshold calcium current Cerebellum Purkinje Cell Model\n"
  "\n"
  "COMMENT\n"
  "\n"
  "Q10 is estimated from this work, Temperature dependence of T-type Calcium channel gating, NEUROSCIENCE\n"
  "written by Yunliang Zang according to the data provided by Stephane Diudone, compared with the summarised data from stephane,\n"
  "T type calcium channels has two gates. so the activation curve was refitted.\n"
  "The junction potential is -6.6 mV\n"
  "It does not work even changing it back to cai\n"
  "April 16th, 2015\n"
  "This version does not contribute to the calcium concentration and BK together with SK. \n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "        SUFFIX CaT3_1\n"
  ":        USEION ca READ cai, cao WRITE ica VALENCE 2\n"
  ":        NONSPECIFIC_CURRENT i\n"
  "	USEION ca READ cai,cao\n"
  "	USEION Ca WRITE iCa VALENCE 2\n"
  "        RANGE g, pcabar, minf, taum, hinf, tauh\n"
  "    	RANGE iCa, m ,h\n"
  ":    THREADSAFE\n"
  "    }\n"
  "\n"
  "UNITS {\n"
  "        (molar) = (1/liter)\n"
  "        (mV) =  (millivolt)\n"
  "        (mA) =  (milliamp)\n"
  "        (mM) =  (millimolar)\n"
  "\n"
  "}\n"
  "\n"
  "CONSTANT {\n"
  "    q10 = 1.0913        :estimate from Iftinca\n"
  "	F = 9.6485e4 (coulombs)\n"
  "	R = 8.3145 (joule/kelvin)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "        v               (mV)\n"
  "        celsius (degC)\n"
  "        eca (mV)\n"
  "	pcabar  = 2.5e-4 (cm/s)\n"
  "        cai = 1e-4  (mM)           : adjusted for eca=120 mV\n"
  "	cao = 2  (mM)\n"
  "	\n"
  "	v0_m_inf = -42.206 (mV)\n"
  "	v0_h_inf = -75.118 (mV)\n"
  "    vshift = 0.0 :-6.6			:liquid junction potential\n"
  "\n"
  "	k_m_inf = -4.7056 (mV)\n"
  "	k_h_inf = 6.4635  (mV)\n"
  "	\n"
  "	C_tau_m = 1.2757\n"
  "	A_tau_m = -2.3199\n"
  "	v0_tau_m1 = -48.048 (mV)\n"
  "	k_tau_m1 = 30.655 (mV)\n"
  "	B_tau_m = 2.5712\n"
  "	v0_tau_m2 = -28.386 (mV)\n"
  "	k_tau_m2 = 9.6306 (mV)\n"
  "	\n"
  "	C_tau_h = 0.0076\n"
  "	A_tau_h = 0.17746\n"
  "	v0_tau_h1 = -58.535 (mV)\n"
  "	k_tau_h1 = 6.2692 (mV)\n"
  "	B_tau_h = 0.13402\n"
  "	v0_tau_h2=-101.436\n"
  "	k_tau_h2 = -5.5845 (mV)\n"
  "}\n"
  "    \n"
  "\n"
  "STATE {\n"
  "        m h\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "        iCa     (mA/cm2)\n"
  "	g        (coulombs/cm3) \n"
  "        minf\n"
  "        taum   (ms)\n"
  "        hinf\n"
  "        tauh   (ms)\n"
  "        qt\n"
  "	T (kelvin)\n"
  "	E (volt)\n"
  "	zeta\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE castate METHOD cnexp \n"
  "\n"
  "       iCa = (1e3) *pcabar*m*m *m*h * g\n"
  "\n"
  "}\n"
  "\n"
  "DERIVATIVE castate {\n"
  "        evaluate_fct(v)\n"
  "\n"
  "        m' = (minf - m) / taum\n"
  "        h' = (hinf - h) / tauh\n"
  "}\n"
  "\n"
  "FUNCTION ghk2( v (mV), ci (mM), co (mM), z )  (coulombs/cm3) {\n"
  "    E = (1e-3) * v\n"
  "      zeta = (z*F*E)/(R*T)\n"
  "\n"
  "\n"
  "    if ( fabs(1-exp(-zeta)) < 1e-6 ) {\n"
  "        ghk2 = (1e-6) * (z*F) * (ci - co*exp(-zeta)) * (1 + zeta/2)\n"
  "    } else {\n"
  "        ghk2 = (1e-6) * (z*zeta*F) * (ci - co*exp(-zeta)) / (1-exp(-zeta))\n"
  "    }\n"
  "}\n"
  "\n"
  "\n"
  "UNITSOFF\n"
  "INITIAL {\n"
  "	\n"
  "	T = kelvinfkt (celsius)\n"
  "	    qt = q10^((celsius-32 (degC))/10 (degC))\n"
  "        evaluate_fct(v)\n"
  "        m = minf\n"
  "        h = hinf\n"
  "}\n"
  "\n"
  "PROCEDURE evaluate_fct(v(mV)) { \n"
  "\n"
  "        minf = 1.0 / ( 1 + exp((v  - v0_m_inf-vshift)/k_m_inf) )^(1/3)\n"
  "        \n"
  "        hinf = 1.0 / ( 1 + exp((v - v0_h_inf-vshift)/k_h_inf) )\n"
  "\n"
  "	taum = 1/( C_tau_m + A_tau_m / (1+exp((v0_tau_m1-v-vshift)/ k_tau_m1))+ B_tau_m/ (1+exp((v0_tau_m2-v-vshift)/k_tau_m2)))/qt\n"
  "\n"
  "	tauh = 1/( C_tau_h + A_tau_h / (1+exp((v0_tau_h1-v-vshift)/ k_tau_h1))+ B_tau_h/ (1+exp((v0_tau_h2-v-vshift)/k_tau_h2)))/qt\n"
  "\n"
  "	g = ghk2(v-vshift, cai, cao, 2)\n"
  "}\n"
  "\n"
  "FUNCTION kelvinfkt( t (degC) )  (kelvin) {\n"
  "    kelvinfkt = 273.19 + t\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
