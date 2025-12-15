/*  macOS, otherwise  CBLAS/LAPACKE  */
#if defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif
