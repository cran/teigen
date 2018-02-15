#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void mycovwt(void *, void *, void *, void *, void *, void *, void *);
extern void mymaha(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"mycovwt", (DL_FUNC) &mycovwt, 7},
    {"mymaha",  (DL_FUNC) &mymaha,  7},
    {NULL, NULL, 0}
};

void R_init_teigen(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
