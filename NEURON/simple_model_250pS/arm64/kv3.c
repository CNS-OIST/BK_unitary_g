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
 
#define nrn_init _nrn_init__Kv3
#define _nrn_initial _nrn_initial__Kv3
#define nrn_cur _nrn_cur__Kv3
#define _nrn_current _nrn_current__Kv3
#define nrn_jacob _nrn_jacob__Kv3
#define nrn_state _nrn_state__Kv3
#define _net_receive _net_receive__Kv3 
#define rateConst rateConst__Kv3 
#define state state__Kv3 
 
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
#define vshift _p[0]
#define vshift_columnindex 0
#define gbar _p[1]
#define gbar_columnindex 1
#define ik _p[2]
#define ik_columnindex 2
#define g _p[3]
#define g_columnindex 3
#define n _p[4]
#define n_columnindex 4
#define ek _p[5]
#define ek_columnindex 5
#define qt _p[6]
#define qt_columnindex 6
#define alpha _p[7]
#define alpha_columnindex 7
#define beta _p[8]
#define beta_columnindex 8
#define Dn _p[9]
#define Dn_columnindex 9
#define _g _p[10]
#define _g_columnindex 10
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
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
 static void _hoc_alphaFkt(void);
 static void _hoc_betaFkt(void);
 static void _hoc_rateConst(void);
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
 "setdata_Kv3", _hoc_setdata,
 "alphaFkt_Kv3", _hoc_alphaFkt,
 "betaFkt_Kv3", _hoc_betaFkt,
 "rateConst_Kv3", _hoc_rateConst,
 0, 0
};
#define alphaFkt alphaFkt_Kv3
#define betaFkt betaFkt_Kv3
 extern double alphaFkt( double );
 extern double betaFkt( double );
 /* declare global and static user variables */
#define ninf ninf_Kv3
 double ninf = 0;
#define tau tau_Kv3
 double tau = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gbar_Kv3", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ninf_Kv3", "1",
 "tau_Kv3", "ms",
 "vshift_Kv3", "mV",
 "gbar_Kv3", "S/cm2",
 "ik_Kv3", "mA/cm2",
 "g_Kv3", "S/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double n0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ninf_Kv3", &ninf_Kv3,
 "tau_Kv3", &tau_Kv3,
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
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Kv3",
 "vshift_Kv3",
 "gbar_Kv3",
 0,
 "ik_Kv3",
 "g_Kv3",
 0,
 "n_Kv3",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	vshift = 4;
 	gbar = 0.005;
 	_prop->param = _p;
 	_prop->param_size = 11;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
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

 void _kv3_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 11, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Kv3 /Users/iain/GIT/BK_unitary_g/NEURON/simple_model_250pS/mod/kv3.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double e0 = 1.60217646e-19;
 static double q10 = 2.7;
 static double ca = 0.22;
 static double cva = 16;
 static double cka = -26.5;
 static double cb = 0.22;
 static double cvb = 16;
 static double ckb = 26.5;
 static double zn = 1.9196;
static int _reset;
static char *modelname = "Voltage-gated potassium channel from Kv3 subunits";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rateConst(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rateConst ( _threadargscomma_ v ) ;
   Dn = alpha * ( 1.0 - n ) - beta * n ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rateConst ( _threadargscomma_ v ) ;
 Dn = Dn  / (1. - dt*( ( alpha )*( ( ( - 1.0 ) ) ) - ( beta )*( 1.0 ) )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
   rateConst ( _threadargscomma_ v ) ;
    n = n + (1. - exp(dt*(( alpha )*( ( ( - 1.0 ) ) ) - ( beta )*( 1.0 ))))*(- ( ( alpha )*( ( 1.0 ) ) ) / ( ( alpha )*( ( ( - 1.0 ) ) ) - ( beta )*( 1.0 ) ) - n) ;
   }
  return 0;
}
 
static int  rateConst (  double _lv ) {
   alpha = qt * alphaFkt ( _threadargscomma_ _lv ) ;
   beta = qt * betaFkt ( _threadargscomma_ _lv ) ;
   ninf = alpha / ( alpha + beta ) ;
   tau = 1.0 / ( alpha + beta ) ;
    return 0; }
 
static void _hoc_rateConst(void) {
  double _r;
   _r = 1.;
 rateConst (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alphaFkt (  double _lv ) {
   double _lalphaFkt;
 _lalphaFkt = ca * exp ( - ( _lv + cva + vshift ) / cka ) ;
   
return _lalphaFkt;
 }
 
static void _hoc_alphaFkt(void) {
  double _r;
   _r =  alphaFkt (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double betaFkt (  double _lv ) {
   double _lbetaFkt;
 _lbetaFkt = cb * exp ( - ( _lv + cvb + vshift ) / ckb ) ;
   
return _lbetaFkt;
 }
 
static void _hoc_betaFkt(void) {
  double _r;
   _r =  betaFkt (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
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
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  n = n0;
 {
   qt = pow( q10 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   rateConst ( _threadargscomma_ v ) ;
   n = ninf ;
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
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gbar * pow( n , 4.0 ) ;
   ik = g * ( v - ek ) ;
   }
 _current += ik;

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
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
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
  ek = _ion_ek;
 { error =  state();
 if(error){fprintf(stderr,"at line 100 in file kv3.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = n_columnindex;  _dlist1[0] = Dn_columnindex;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/iain/GIT/BK_unitary_g/NEURON/simple_model_250pS/mod/kv3.mod";
static const char* nmodl_file_text = 
  "TITLE Voltage-gated potassium channel from Kv3 subunits\n"
  "COMMENT\n"
  "Voltage-gated potassium channel with high threshold and fast activation/deactivation kinetics\n"
  "\n"
  "KINETIC SCHEME: Hodgkin-Huxley (n^4)\n"
  "n'= alpha * (1-n) - betha * n\n"
  "g(v) = gbar * n^4 * ( v-ek )\n"
  "\n"
  "The rate constants of activation (alpha) and deactivation (beta) were approximated by:\n"
  "\n"
  "alpha(v) = ca * exp(-(v+cva)/cka)\n"
  "beta(v) = cb * exp(-(v+cvb)/ckb)\n"
  "\n"
  "Parameters can, cvan, ckan, cbn, cvbn, ckbn are given in the CONSTANT block.\n"
  "Values derive from least-square fits to experimental data of G/Gmax(v) and taun(v) in Martina et al. J Neurophysiol 97:563-571, 2007\n"
  "Model includes a calculation of Kv gating current\n"
  "\n"
  "Reference: Akemann et al., Biophys. J. (2009) 96: 3959-3976\n"
  "Notice that there is another set of data related with Kv3 by McKay and Turner European Journal of Neuroscience, Vol. 20, pp. 729\n"
  "\n"
  "\n"
  "739, 2004\n"
  "in that paper, the activation threshold of Kv3 is much lower.\n"
  "Laboratory for Neuronal Circuit Dynamics\n"
  "RIKEN Brain Science Institute, Wako City, Japan\n"
  "http://www.neurodynamics.brain.riken.jp\n"
  "\n"
  "Date of Implementation: April 2007\n"
  "Contact: akemann@brain.riken.jp\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX Kv3\n"
  "	USEION k READ ek WRITE ik\n"
  "	RANGE gbar, g, ik,vshift\n"
  "	GLOBAL ninf, tau\n"
  ":    THREADSAFE\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mV) = (millivolt)\n"
  "	(mA) = (milliamp)\n"
  "	(nA) = (nanoamp)\n"
  "	(pA) = (picoamp)\n"
  "	(S)  = (siemens)\n"
  "	(mS) = (millisiemens)\n"
  "	(nS) = (nanosiemens)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)		\n"
  "}\n"
  "\n"
  "CONSTANT {\n"
  "	e0 = 1.60217646e-19 (coulombs)\n"
  "	q10 = 2.7\n"
  "\n"
  "	ca = 0.22 (1/ms)\n"
  "	cva = 16 (mV)\n"
  "	cka = -26.5 (mV)\n"
  "	cb = 0.22 (1/ms)\n"
  "	cvb = 16 (mV)\n"
  "	ckb = 26.5 (mV)\n"
  "	\n"
  "	zn = 1.9196 (1)		: valence of n-gate\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	vshift = 4 (mV)\n"
  "	gbar = 0.005 (S/cm2)   <0,1e9>\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	celsius (degC)\n"
  "	v (mV)\n"
  "	\n"
  "	ik (mA/cm2)\n"
  " \n"
  "	ek (mV)\n"
  "	g (S/cm2)\n"
  "	qt (1)\n"
  "\n"
  "	ninf (1)\n"
  "	tau (ms)\n"
  "	alpha (1/ms)\n"
  "	beta (1/ms)\n"
  "}\n"
  "\n"
  "STATE { n }\n"
  "\n"
  "INITIAL {\n"
  "	qt = q10^((celsius-22 (degC))/10 (degC))\n"
  "	rateConst(v)\n"
  "	n = ninf\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "      g = gbar * n^4 \n"
  "	ik = g * (v - ek)\n"
  "\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	rateConst(v)\n"
  "	n' = alpha * (1-n) - beta * n\n"
  "}\n"
  "\n"
  "PROCEDURE rateConst(v (mV)) {\n"
  "	alpha = qt * alphaFkt(v)\n"
  "	beta = qt * betaFkt(v)\n"
  "	ninf = alpha / (alpha + beta) \n"
  "	tau = 1 / (alpha + beta)\n"
  "}\n"
  "\n"
  "FUNCTION alphaFkt(v (mV)) (1/ms) {\n"
  "	alphaFkt = ca * exp(-(v+cva+vshift)/cka)\n"
  "}\n"
  "\n"
  "FUNCTION betaFkt(v (mV)) (1/ms) {\n"
  "	betaFkt = cb * exp(-(v+cvb+vshift)/ckb)\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
