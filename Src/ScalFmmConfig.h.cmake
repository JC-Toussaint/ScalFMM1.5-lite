// See LICENCE file at project root

#ifndef SSCALFMMCONFIG_H
#define SSCALFMMCONFIG_H

///////////////////////////////////////////////////////
// Debug
///////////////////////////////////////////////////////

#cmakedefine SCALFMM_USE_LOG
#cmakedefine SCALFMM_USE_ASSERT
#cmakedefine SCALFMM_USE_MEM_STATS

///////////////////////////////////////////////////////
// Blas / Lapack (required by the Chebyshev kernel)
///////////////////////////////////////////////////////

#cmakedefine SCALFMM_USE_BLAS
// Fortran Mangling
#cmakedefine SCALFMM_BLAS_ADD_
#cmakedefine SCALFMM_BLAS_UPCASE
#cmakedefine SCALFMM_BLAS_NOCHANGE

///////////////////////////////////////////////////////
// Vectorized P2P
///////////////////////////////////////////////////////

#cmakedefine SCALFMM_USE_SSE
#cmakedefine SCALFMM_USE_AVX

#ifdef __INTEL_COMPILER
#pragma warning (disable : 858 )
#pragma warning (disable : 2326 )
#endif

///////////////////////////////////////////////////////
// Flags and libs used to compile
///////////////////////////////////////////////////////
#include <string>
const std::string SCALFMMCompileFlags("@SCALFMM_COMPILE_FLAGS@");
const std::string SCALFMMCompileLibs("@SCALFMM_COMPILE_LIBS@");

#endif // SSCALFMMCONFIG_H
