#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void acfangl(void *, void *, void *, void *, void *, void *, void *, void *);
extern void acfdist(void *, void *, void *, void *, void *, void *, void *, void *);
extern void discretrajr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void fipatir(void *, void *, void *, void *, void *, void *, void *);
extern void optcutr(void *, void *, void *, void *, void *, void *);
extern void partrajr(void *, void *, void *, void *, void *, void *, void *);
extern void permutR2n(void *, void *, void *, void *, void *, void *);
extern void prepquart(void *, void *, void *, void *, void *, void *, void *);
extern void runsltr(void *, void *, void *, void *);
extern void testindepangl(void *, void *, void *, void *, void *, void *, void *);
extern void testindepdist(void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP contrastM(SEXP, SEXP, SEXP, SEXP);
extern SEXP dynprog(SEXP, SEXP);
extern SEXP findpathc(SEXP, SEXP, SEXP);
extern SEXP RasterPas(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP redistime(SEXP, SEXP, SEXP);
extern SEXP residtime(SEXP, SEXP, SEXP);
extern SEXP simulmod(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simulmodmv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"acfangl",       (DL_FUNC) &acfangl,        8},
    {"acfdist",       (DL_FUNC) &acfdist,        8},
    {"discretrajr",   (DL_FUNC) &discretrajr,   13},
    {"fipatir",       (DL_FUNC) &fipatir,        7},
    {"optcutr",       (DL_FUNC) &optcutr,        6},
    {"partrajr",      (DL_FUNC) &partrajr,       7},
    {"permutR2n",     (DL_FUNC) &permutR2n,      6},
    {"prepquart",     (DL_FUNC) &prepquart,      7},
    {"runsltr",       (DL_FUNC) &runsltr,        4},
    {"testindepangl", (DL_FUNC) &testindepangl,  7},
    {"testindepdist", (DL_FUNC) &testindepdist,  7},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"contrastM",  (DL_FUNC) &contrastM,   4},
    {"dynprog",    (DL_FUNC) &dynprog,     2},
    {"findpathc",   (DL_FUNC) &findpathc,    3},
    {"RasterPas",  (DL_FUNC) &RasterPas,   5},
    {"redistime",  (DL_FUNC) &redistime,   3},
    {"residtime",  (DL_FUNC) &residtime,   3},
    {"simulmod",   (DL_FUNC) &simulmod,    5},
    {"simulmodmv", (DL_FUNC) &simulmodmv, 13},
    {NULL, NULL, 0}
};

void R_init_adehabitatLT(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
