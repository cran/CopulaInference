#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .C calls */
extern void estdep(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void taucopula(void *, void *, void *, void *);
extern void ppla(void *, void *, void *, void *,  void *);
extern void dpla(void *, void *, void *, void *, void *);
extern void hpla(void *, void *, void *, void *, void *, void *);
extern void taupla(void *, void *);
extern void prepare_data(void *, void *, void *, void *, void *, void *);
extern void bi_emp_cdf(void *, void *, void *, void *, void *, void *, void *, void *);
extern void emp_cdf(void *, void *, void *, void *, void *);
extern void est_dep(void *, void *, void *, void *, void *);
static const R_CMethodDef CEntries[] = {
  {"estdep",           (DL_FUNC) &estdep,            12},
  {"taucopula", (DL_FUNC) &taucopula, 4},
  { "ppla", (DL_FUNC) &ppla, 5},
  { "dpla", (DL_FUNC) &dpla, 5},
  { "hpla", (DL_FUNC) &hpla, 6},
  { "taupla", (DL_FUNC) &taupla, 2},
  { "prepare_data", (DL_FUNC) &prepare_data, 6},
  { "bi_emp_cdf", (DL_FUNC) &bi_emp_cdf, 8},
  { "emp_cdf", (DL_FUNC) &emp_cdf, 5},
  { "est_dep", (DL_FUNC) &est_dep, 5},
  {NULL, NULL, 0}
};

void R_init_CopulaInference(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
