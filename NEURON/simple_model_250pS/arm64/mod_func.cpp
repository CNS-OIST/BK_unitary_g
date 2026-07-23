#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _CaP_reg(void);
extern void _CaT_reg(void);
extern void _cdp_AIS_reg(void);
extern void _ih_reg(void);
extern void _kv3_reg(void);
extern void _mslo_reg(void);
extern void _nap_reg(void);
extern void _narsg_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"mod/CaP.mod\"");
    fprintf(stderr, " \"mod/CaT.mod\"");
    fprintf(stderr, " \"mod/cdp_AIS.mod\"");
    fprintf(stderr, " \"mod/ih.mod\"");
    fprintf(stderr, " \"mod/kv3.mod\"");
    fprintf(stderr, " \"mod/mslo.mod\"");
    fprintf(stderr, " \"mod/nap.mod\"");
    fprintf(stderr, " \"mod/narsg.mod\"");
    fprintf(stderr, "\n");
  }
  _CaP_reg();
  _CaT_reg();
  _cdp_AIS_reg();
  _ih_reg();
  _kv3_reg();
  _mslo_reg();
  _nap_reg();
  _narsg_reg();
}

#if defined(__cplusplus)
}
#endif
