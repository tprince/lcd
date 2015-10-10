/* loops.f -- translated by f2c (version 19991025).
*/

//#define restrict
#ifdef _OPENMP
#include <omp.h>
#endif
extern "C"{
#include "g2c.h"
#include <stdlib.h>
}
#include <cmath>
#include <numeric>
#include <functional>
#include <vector>
#include <iterator>
#include <algorithm>

using namespace std;

/* Common Block Declarations */

extern "C" struct {
  real array[1000*1000];
}
cdata_;

#define cdata_1 cdata_

/* Table of constant values */

static real c_b3 = 1.f;
static integer c__1 = 1;
static real c_b393 = 0.f;
  extern "C"/* Subroutine */ int init_ (integer *, integer *, real *, real *,
	 real *, real *, real *, real *, real *, real *, const char *, ftnlen),
    forttime_ (real *), set1d_ (integer *, real *, real *, integer *);
  extern "C" /* Subroutine */ int check_ (real *, integer *, integer *, real *,
      const char *, ftnlen), dummy_ (integer *, integer *, real *, real *,
       real *, real *, real *, real *, real *, real *, real *);
  extern "C" real cs1d_ (integer *, real *);
  extern "C" real cs2d_ (integer *, real *);
  extern "C" /* Subroutine */ int s471s_ (void);
#if !defined __INTEL_COMPILER
#if defined __GNUC__
#define max(x,y) fmaxf(x,y)
#define min(x,y) fminf(x,y)
#else
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))
#endif
#endif

/* *********************************************************************** */
/*                TEST SUITE FOR VECTORIZING COMPILERS                  * */
/*                        (File 2 of 2)                                 * */
/*                                                                      * */
/*  Version:   2.0                                                      * */
/*  Date:      3/14/88                                                  * */
/*  Authors:   Original loops from a variety of                         * */
/*             sources. Collection and synthesis by                     * */
/*                                                                      * */
/*             David Callahan  -  Tera Computer                         * */
/*             Jack Dongarra   -  University of Tennessee               * */
/*             David Levine    -  Argonne National Laboratory           * */
/* *********************************************************************** */
/*  Version:   3.0                                                      * */
/*  Date:      1/4/91                                                   * */
/*  Authors:   David Levine    -  Executable version                    * */
/* *********************************************************************** */
/*                         ---DESCRIPTION---                            * */
/*                                                                      * */
/*  This test consists of a variety of  loops that represent different  * */
/*  constructs  intended   to  test  the   analysis  capabilities of a  * */
/*  vectorizing  compiler.  Each loop is  executed with vector lengths  * */
/*  of 10, 100, and 1000.   Also included  are several simple  control  * */
/*  loops  intended  to  provide  a  baseline  measure  for  comparing  * */
/*  compiler performance on the more complicated loops.                 * */
/*                                                                      * */
/*  The  output from a  run  of the test  consists of seven columns of  * */
/*  data:                                                               * */
/*     Loop:        The name of the loop.                               * */
/*     VL:          The vector length the loop was run at.              * */
/*     Seconds:     The time in seconds to run the loop.                * */
/*     Checksum:    The checksum calculated when running the test.      * */
/*     PreComputed: The precomputed checksum (64-bit IEEE arithmetic).  * */
/*     Residual:    A measure of the accuracy of the calculated         * */
/*                  checksum versus the precomputed checksum.           * */
/*     No.:         The number of the loop in the test suite.           * */
/*                                                                      * */
/*  The  residual  calculation  is  intended  as  a  check  that   the  * */
/*  computation  was  done  correctly  and  that 64-bit arithmetic was  * */
/*  used.   Small   residuals    from    non-IEEE    arithmetic    and  * */
/*  nonassociativity  of   some calculations  are   acceptable.  Large  * */
/*  residuals  from   incorrect  computations or  the  use   of 32-bit  * */
/*  arithmetic are not acceptable.                                      * */
/*                                                                      * */
/*  The test  output  itself  does not report   any  results;  it just  * */
/*  contains data.  Absolute  measures  such as Mflops and  total time  * */
/*  used  are  not   appropriate    metrics  for  this  test.   Proper  * */
/*  interpretation of the results involves correlating the output from  * */
/*  scalar and vector runs  and the  loops which  have been vectorized  * */
/*  with the speedup achieved at different vector lengths.              * */
/*                                                                      * */
/*  These loops  are intended only  as  a partial test of the analysis  * */
/*  capabilities of a vectorizing compiler (and, by necessity,  a test  * */
/*  of the speed and  features  of the underlying   vector  hardware).  * */
/*  These loops  are  by no means  a  complete  test  of a vectorizing  * */
/*  compiler and should not be interpreted as such.                     * */
/*                                                                      * */
/* *********************************************************************** */
/*                           ---DIRECTIONS---                           * */
/*                                                                      * */
/*  To  run this  test,  you will  need  to  supply  a  function named  * */
/*  second() that returns user CPU time.                                * */
/*                                                                      * */
/*  This test is distributed as two separate files, one containing the  * */
/*  driver  and  one containing the loops.   These  two files MUST  be  * */
/*  compiled separately.                                                * */
/*                                                                      * */
/*  Results must  be supplied from  both scalar and vector  runs using  * */
/*  the following rules for compilation:                                * */
/*                                                                      * */
/*     Compilation   of the  driver  file must  not  use any  compiler  * */
/*     optimizations (e.g., vectorization, function  inlining,  global  * */
/*     optimizations,...).   This   file   also must  not  be analyzed  * */
/*     interprocedurally to  gather information useful  in  optimizing  * */
/*     the test loops.                                                  * */
/*                                                                      * */
/*     The file containing the  loops must be compiled twice--once for  * */
/*     a scalar run and once for a vector run.                          * */
/*                                                                      * */
/*        For the scalar  run, global (scalar) optimizations should be  * */
/*        used.                                                         * */
/*                                                                      * */
/*        For  the  vector run,  in  addition   to  the  same   global  * */
/*        optimizations specified  in the scalar   run,  vectorization  * */
/*        and--if available--automatic  call generation to   optimized  * */
/*        library  routines,  function inlining,  and  interprocedural  * */
/*        analysis should be  used.  Note again that function inlining  * */
/*        and interprocedural  analysis  must  not be  used to  gather  * */
/*        information  about any of the  program  units  in the driver  * */
/*        program.                                                      * */
/*                                                                      * */
/*     No changes  may  be made  to   the   source code.   No compiler  * */
/*     directives may be used, nor may  a file  be split into separate  * */
/*     program units.  (The exception is  filling  in  the information  * */
/*     requested in subroutine "info" as described below.)              * */
/*                                                                      * */
/*     All files must be compiled to use 64-bit arithmetic.             * */
/*                                                                      * */
/*     The  outer  timing  loop  is  included  only   to increase  the  * */
/*     granularity of the calculation.  It should not be vectorizable.  * */
/*     If it is found to be so, please notify the authors.              * */
/*                                                                      * */
/*  All runs  must be  made  on a standalone  system  to minimize  any  * */
/*  external effects.                                                   * */
/*                                                                      * */
/*  On virtual memory computers,  runs should be  made with a physical  * */
/*  memory and working-set  size  large enough  that  any  performance  * */
/*  degradation from page  faults is negligible.   Also,  the  timings  * */
/*  should be repeatable  and you  must  ensure  that timing anomalies  * */
/*  resulting from paging effects are not present.                      * */
/*                                                                      * */
/*  You should edit subroutine "info"   (the  last subroutine  in  the  * */
/*  driver program) with  information specific to your  runs, so  that  * */
/*  the test output will be annotated automatically.                    * */
/*                                                                      * */
/*  Please return the following three files in an electronic format:    * */
/*                                                                      * */
/*  1. Test output from a scalar run.                                   * */
/*  2. Test output from a vector run.                                   * */
/*  3. Compiler output listing (source echo, diagnostics, and messages) * */
/*     showing which loops have been vectorized.                        * */
/*                                                                      * */
/*  The preferred media  for receipt, in order  of preference, are (1)  * */
/*  electronic mail, (2) 9-track  magnetic or  cartridge tape in  Unix  * */
/*  tar  format, (3) 5" IBM PC/DOS  floppy   diskette, or  (4) 9-track  * */
/*  magnetic  tape in  ascii  format,   80 characters per card,  fixed  * */
/*  records, 40 records per block, 1600bpi.  Please return to           * */
/*                                                                      * */
/*  David Levine       		                                       * */
/*  Mathematics and Computer Science Division                           * */
/*  Argonne National Laboratory                                         * */
/*  Argonne, Illinois 60439                                             * */
/*  levine@mcs.anl.gov                                                  * */
/* *********************************************************************** */
/* %1.1 */
/* Subroutine */ extern "C" int
s114_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer i__,j;
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  real t1, t2;
  integer nl;
  real chksum;


/*     linear dependence testing */
/*     transpose vectorization */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s114 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = (*ntimes << 1) / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma omp parallel for if(i__2 > 101)
      for (j = 2; j <= i__2; ++j) {
	  int i__3 = j - 1;
#if defined _OPENMP && _OPENMP >= 201307 && ! __MIC__
#pragma omp simd
#endif
	  for (int i__ = 1; i__ <= i__3; ++i__)
	      aa[i__ + j * aa_dim1] = aa[j + i__ * aa_dim1] + bb[i__ + j *
								 bb_dim1];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n << 1);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = (*ntimes / *n << 1) * ((*n * *n - *n) / 2);
  check_ (&chksum, &i__1, n, &t2, "s114 ", (ftnlen) 5);
  return 0;
}				/* s114_ */

/* %1.2 */
/* Subroutine */ extern "C" int
s125_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j, k;
  real t1, t2;
  integer nl;
  real chksum;


/*     induction variable recognition */
/*     induction variable in two loops; collapsing possible */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s125 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      k = 0;
      i__2 = *n;
	  i__3 = *n;
#pragma omp parallel for if(i__2 > 101)
      for (j = 1; j <= i__2; ++j) {
	  int k = (j-1)*i__3;
	  for (int i__ = 1; i__ <= i__3; ++i__)
	      cdata_1.array[k++] = aa[i__ + j * aa_dim1] + bb[i__ + j *
					    bb_dim1] * cc[i__ + j * cc_dim1];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  i__1 = *n * *n;
  chksum = cs1d_ (&i__1, cdata_1.array);
  i__1 = *ntimes / *n * *n * *n;
  check_ (&chksum, &i__1, n, &t2, "s125 ", (ftnlen) 5);
  return 0;
}				/* s125_ */

/* %1.2 */
/* Subroutine */ extern "C" int
s126_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * __restrict bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j, k;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs2d_ (integer *, real *);


/*     induction variable recognition */
/*     induction variable in two loops; recurrence in inner loop */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s126 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
#if defined _OPENMP && _OPENMP < 201307 || !defined __INTEL_COMPILER && __AVX__
      for (int j = 2; j <= i__3; ++j){
	  int k = j - 2 - i__3;
	  for (int i__ = 1; i__ <= i__2; ++i__)
	      bb[i__ + j * bb_dim1] = bb[i__ + j * bb_dim1 - bb_dim1] +
		cdata_1.array[k+= i__3] * cc[i__ + j * cc_dim1];
#else
#pragma omp parallel for simd if(i__2 > 103)
      for (i__ = 1; i__ <= i__2; ++i__) {
	  int k = i__ * i__3 - i__3 + 1;
	  for (int j = 2; j <= i__3; ++j)
	      bb[i__ + j * bb_dim1] = bb[i__ + j * bb_dim1 - bb_dim1] +
		cdata_1.array[k++ - 1] * cc[i__ + j * cc_dim1];
#endif
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &bb[bb_offset]);
  i__1 = *ntimes / *n * *n * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s126 ", (ftnlen) 5);
  return 0;
}				/* s126_ */

/* %1.3 */
/* Subroutine */ extern "C" int
s132_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__,
       real * d__, real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j, k, m;
  real t1, t2;
  integer nl;
  real chksum;


/*     global data flow analysis */
/*     loop with multiple dimension ambiguous subscripts */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  m = 1;
  j = m;
  k = m + 1;
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s132 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd safelen(32)
#endif
      for (i__ = 2; i__ <= i__2; ++i__)
	  aa[i__ + j * aa_dim1] = aa[i__ - 1 + k * aa_dim1] + b[i__] * c__[k];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes * *n - 1;
  check_ (&chksum, &i__1, n, &t2, "s132 ", (ftnlen) 5);
  return 0;
}				/* s132_ */

/* %1.4 */
/* Subroutine */ extern "C" int
s141_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real *__restrict bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j, k;
  real t1, t2;
  integer nl;
  real chksum;


/*     nonlinear dependence testing */
/*     walk a row in a symmetric packed array */
/*     element a(i,j) for (j>i) stored in location j*(j-1)/2+i */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s141 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma omp parallel for if(i__2 > 103)
      for (j = 1; j <= i__2; ++j) {
	  int k = j * (j - 1) / 2;
	  transform(&cdata_1.array[k],&cdata_1.array[k+j],
	    &bb[1+j*bb_dim1],&cdata_1.array[k],plus<float>());
//	  for (int i__ = 1; i__ <= j; ++i__)
//	      cdata_1.array[i__ + k - 1] += bb[i__ + j * bb_dim1];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  i__1 = *n * *n;
  chksum = cs1d_ (&i__1, cdata_1.array);
  i__1 = *ntimes / *n * *n * *n;
  check_ (&chksum, &i__1, n, &t2, "s141 ", (ftnlen) 5);
  return 0;
}				/* s141_ */

/* %1.6 */
/* Subroutine */ extern "C" int
s161_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b,
       real * __restrict c__, real * d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control flow */
/*     tests for recognition of loop independent dependences */
/*     between statements in mutually exclusive regions. */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s161 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
#if defined __AVX2__ && OPENMP >= 201307 || defined __MIC__
#pragma omp simd
#endif
      for (i__ = 1; i__ <= i__2; ++i__)
	  if (b[i__] >= 0.f)
	      a[i__] = c__[i__] + d__[i__] * e[i__];
	  else
	      c__[i__ + 1] = a[i__] + d__[i__] * d__[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &c__[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s161 ", (ftnlen) 5);
  return 0;
}				/* s161_ */

/* %1.6 */
/* Subroutine */ extern "C" int
s162_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b,
       real * c__, real * d__, real * e, real * aa, real * bb,
       real * cc, integer * k) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control flow */
/*     deriving assertions */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s162 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      if (*k > 0) {
	  i__2 = *n - 1;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
	  for (i__ = 1; i__ <= i__2; ++i__)
	      a[i__] = a[i__ + *k] + b[i__] * c__[i__];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s162 ", (ftnlen) 5);
  return 0;
}				/* s162_ */

/* %1.7 */
/* Subroutine */ extern "C" int
s171_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc, integer * inc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     symbolics */
/*     symbolic dependence tests */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s171 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
      for (i__ = 1; i__ <= i__2; ++i__)
	  a[i__ * *inc] += b[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s171 ", (ftnlen) 5);
  return 0;
}				/* s171_ */

/* %1.7 */
/* Subroutine */ extern "C" int
s172_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc, integer * n1,
       integer * n3) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     symbolics */
/*     vectorizable if n3 .ne. 0 */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s172 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      i__3 = *n3;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
      for (i__ = *n1; i__ <= i__2; i__ += i__3)
	  a[i__] += b[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s172 ", (ftnlen) 5);
  return 0;
}				/* s172_ */

/* %1.7 */
/* Subroutine */ extern "C" int
s175_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc, integer * inc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     symbolics */
/*     symbolic dependence tests */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s175 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - *inc;	//should be *n - *inc
      i__3 = *inc;
//      if(i__3 == 1)
//	transform(&a[1],&a[i__2]+1,&b[1],&a[2],plus<float>());
//      else
#if defined _OPENMP && _OPENMP >= 201307 && (defined __MIC__ || defined __AVX2__)
#pragma omp simd
#endif
      for (i__ = 1; i__ <= i__2; i__ += i__3)
	  a[i__] = a[i__ + i__3] + b[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s175 ", (ftnlen) 5);
  return 0;
}				/* s175_ */

/* %1.7 */
/* Subroutine */ extern "C" int
s176_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b,
       real * c__, real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j, m;
  real t1, t2;
  integer nl;
  real chksum;


/*     symbolics */
/*     convolution */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  m = *n / 2;
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s176 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = (*ntimes << 2) / *n;
	i__3 = i__2 = m;
  for (nl = 1; nl <= i__1; ++nl) {
#if defined __INTEL_COMPILER
      // set up class in reverse order; count time spent
      vector<float> Cr(m);
      reverse_copy(&c__[1],&c__[i__3]+1,Cr.begin());
#pragma omp parallel for if(i__3 > 103)
    for (i__ = 1; i__ <= i__3; ++i__) 
	a[i__] += inner_product(Cr.begin(),Cr.end(),&b[i__],0.f);
#else
	for (i__ = 1; i__ <= i__3; ++i__) 
	    for (int j = 1; j <= i__2; ++j) 
		a[j] += b[j+i__2-i__] * c__[i__];
#endif
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n << 2);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes / *n << 2) * (*n / 2) * (*n / 2);
  check_ (&chksum, &i__1, n, &t2, "s176 ", (ftnlen) 5);
  return 0;
}				/* s176_ */


/* ********************************************************** */
/*                                                         * */
/*                      VECTORIZATION                      * */
/*                                                         * */
/* ********************************************************** */
/* %2.1 */
/* Subroutine */ extern "C" int
s211_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * __restrict b,
       real * c__, real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     statement reordering */
/*     statement reordering allows vectorization */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s211 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
// pragma used as portable way to imply nofusion
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
      for (i__ = 2; i__ <= i__2; ++i__)
	  b[i__] = b[i__ + 1] - e[i__] * d__[i__];
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
      for (i__ = 2; i__ <= i__2; ++i__)
	  a[i__] = b[i__ - 1] + c__[i__] * d__[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 2);
  check_ (&chksum, &i__1, n, &t2, "s211 ", (ftnlen) 5);
  return 0;
}				/* s211_ */

/* %2.3 */
/* Subroutine */ extern "C" int
s232_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop interchange */
/*     interchanging of triangular loops */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s232 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = (*ntimes << 1) / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma omp parallel for schedule(guided) if(i__2 > 101)
      for (j = 2; j <= i__2; ++j) {
	  int i__3 = j;
	  for (int i__ = 2; i__ <= i__3; ++i__)
	      aa[i__ + j * aa_dim1] = aa[i__ - 1 + j * aa_dim1] * aa[i__ -
				     1 + j * aa_dim1] + bb[i__ + j * bb_dim1];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n << 1);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = (*ntimes / *n << 1) * ((*n * *n - *n) / 2);
  check_ (&chksum, &i__1, n, &t2, "s232 ", (ftnlen) 5);
  return 0;
}				/* s232_ */

/* %2.3 */
/* Subroutine */ extern "C" int
s233_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * __restrict aa, real * __restrict bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop interchange */
/*     interchanging with one of two inner loops */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s233 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
#pragma omp parallel if(i__2 > 53)
      {
#pragma omp for nowait
      for (int j = 2; j <= i__3; ++j)
	  for (int i__ = 2; i__ <= i__2; ++i__)
	      bb[i__ + j * bb_dim1] = bb[i__ - 1 + j*bb_dim1] +
		    cc[i__ + j * cc_dim1];

#if defined _OPENMP && _OPENMP >= 201307 && (defined __INTEL_COMPILER || !defined __AVX__)
#pragma omp for simd
	  for (int i__ = 2; i__ <= i__2; ++i__)
      for (int j = 2; j <= i__3; ++j)
#else
#pragma omp single
      for (int j = 2; j <= i__3; ++j)
	  for (int i__ = 2; i__ <= i__2; ++i__)
#endif
	      aa[i__ + j * aa_dim1] = aa[i__ + (j - 1) * aa_dim1] + cc[i__
							       + j * cc_dim1];
      }
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
      }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &aa[aa_offset]) + cs2d_ (n, &bb[bb_offset]);
  i__1 = *ntimes / *n * (*n - 1) * ((*n << 1) - 2);
  check_ (&chksum, &i__1, n, &t2, "s233 ", (ftnlen) 5);
  return 0;
}				/* s233_ */

/* %2.3 */
/* Subroutine */ extern "C" int
s234_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop interchange */
/*     if loop to do loop, interchanging with if loop necessary */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s234 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
	int i__2 = *n;
#if defined _OPENMP && _OPENMP >= 201307 && (defined __INTEL_COMPILER || !defined __AVX__)
#pragma omp parallel for simd if(*n > 53)
	for(int i__ = 1; i__<= i__2; ++i__)
	    for(int j = 2; j <= i__2; ++j)
#else
	    for(int j = 2; j <= i__2; ++j)
	for(int i__ = 1; i__<= i__2; ++i__)
#endif
	    aa[i__ + j * aa_dim1] = aa[i__ + (j - 1) * aa_dim1] +
		bb[i__ + (j - 1) * bb_dim1] * cc[i__ + (j - 1) * cc_dim1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes / *n * *n * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s234 ", (ftnlen) 5);
  return 0;
}				/* s234_ */

/* %2.3 */
/* Subroutine */ extern "C" int
s235_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b, real * c__,
       real * d__, real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop interchanging */
/*     imperfectly nested loops */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s235 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	  a[i__] += b[i__] * c__[i__];
#if defined _OPENMP && _OPENMP >= 201307 && (defined __INTEL_COMPILER || !defined __AVX__)
#pragma omp parallel for simd if(i__2 > 53)
      for (i__ = 1; i__ <= i__2; ++i__)
	  for (int j = 2; j <= i__3; ++j)
#else
	  for (int j = 2; j <= i__3; ++j)
      for (int i__ = 1; i__ <= i__2; ++i__)
#endif
	      aa[i__ + j * aa_dim1] = aa[i__ + (j - 1) * aa_dim1] + bb[i__
						       + j * bb_dim1] * a[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &aa[aa_offset]) + cs1d_ (n, &a[1]);
  i__1 = *ntimes / *n * *n * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s235 ", (ftnlen) 5);
  return 0;
}				/* s235_ */

/* %2.5 */
/* Subroutine */ extern "C" int
s255_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  integer im1, im2;
  real x,y;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     scalar and array expansion */
/*     carry around variables, 2 levels */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s255 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      x = b[*n];
      y = b[*n - 1];
      i__2 = *n;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
      for (i__ = 1; i__ <= i__2; ++i__) {
	  a[i__] = (b[i__] + x + y) * .333f;
	  y = x;
	  x = b[i__];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s255 ", (ftnlen) 5);
  return 0;
}				/* s255_ */

/* %2.5 */
/* Subroutine */ extern "C" int
s257_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__,
       real * d__, real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     scalar and array expansion */
/*     array expansion */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s257 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
#pragma omp parallel if(i__2 > 103)
      {
#pragma omp single
      for (i__ = 2; i__ <= i__2; ++i__) {
	  a[i__] = aa[i__ + i__3 * aa_dim1] - a[i__ - 1];
	  aa[i__ + i__3 * aa_dim1] = a[i__] + bb[i__ + i__3 * bb_dim1];
	}
#pragma omp for nowait
      for (j = 1; j < i__3; ++j){
	  float tmp= a[1];
	  for (int i__ = 2; i__ <= i__2; ++i__) {
	      tmp = aa[i__ + j * aa_dim1] - tmp;
	      aa[i__ + j * aa_dim1] = tmp + bb[i__ + j * bb_dim1];
	    }
	}
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs1d_ (n, &a[1]) + cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes / *n * (*n - 1) * *n;
  check_ (&chksum, &i__1, n, &t2, "s257 ", (ftnlen) 5);
  return 0;
}				/* s257_ */

/* %2.7 */
/* Subroutine */ extern "C" int
s273_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * __restrict b,
       real * __restrict c__, real * d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control flow */
/*     simple loop with dependent conditional */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s273 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
      for (i__ = 1; i__ <= i__2; ++i__) {
	  float tmp= d__[i__] * e[i__];
	  c__[i__] += (a[i__] += tmp) * d__[i__];
	  b[i__] += (a[i__] < 0.f)? tmp: 0.f;
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]) + cs1d_ (n, &c__[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s273 ", (ftnlen) 5);
  return 0;
}				/* s273_ */


/* %2.7 */
/* Subroutine */ extern "C" int
s275_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     control flow */
/*     if around inner loop, interchanging needed */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s275 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
	for (j = 2; j <= i__3; ++j)
//OK if i__2 <= aa_dim1 or aa_dim1 > 64
#if defined __AVX2__ && OPENMP >= 201307 || defined __MIC__
#pragma omp simd safelen(32)
#endif
	    for (i__ = 2; i__ <= i__2; ++i__)
		aa[i__ + j * aa_dim1] = (aa[i__ + aa_dim1] > 0.f)? 
		    aa[i__ + (j - 1) * aa_dim1] + bb[ i__ + j * bb_dim1] *
		    cc[i__ + j * cc_dim1]: aa[i__ + j * aa_dim1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes / *n * (*n - 1) * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s275 ", (ftnlen) 5);
  return 0;
}				/* s275_ */

/* %2.8 */
/* Subroutine */ extern "C" int
s281_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * __restrict b,
       real * c__, real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  real x;
  integer nl;
  real chksum;


/*     crossing thresholds */
/*     index set splitting */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s281 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      i__ = 1;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
    for (i__= i__; i__ <= (i__2+1)/2; ++i__)
	a[i__] = (b[i__] = a[i__2 - i__ + 1] + b[i__] * c__[i__])- 1.f;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
    for (i__= (i__2+3)/2; i__ <= i__2; ++i__)
	a[i__] = (b[i__] = a[i__2 - i__ + 1] + b[i__] * c__[i__])- 1.f;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s281 ", (ftnlen) 5);
  return 0;
}				/* s281_ */

/* %2.10 */
/* Subroutine */ extern "C" int
s2101_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     diagonals */
/*     main diagonal calculation */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s2101", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#if _OPENMP >= 201307
#if defined __INTEL_COMPILER || !defined __AVX__
#pragma omp parallel for simd if(i__2 > 103)
#endif
#else
#pragma omp parallel for if(i__2 > 103)
#endif
      for (i__ = 1; i__ <= i__2; ++i__)
	  aa[i__ + i__ * aa_dim1] += bb[i__ + i__ * bb_dim1] * cc[i__ + i__
								  * cc_dim1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s2101", (ftnlen) 5);
  return 0;
}				/* s2101_ */

/* %2.12 */
/* Subroutine */ extern "C" int
s2102_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     diagonals */
/*     identity matrix, best results vectorize both inner and outer loops */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s2102", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
#pragma omp parallel for if(i__2 > 129)
	for (j = 1; j <= i__2; ++j) 
	    for(int i = 1; i <= i__2; ++i)
		aa[i + j * aa_dim1] = i==j;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes / *n * *n * *n;
  check_ (&chksum, &i__1, n, &t2, "s2102", (ftnlen) 5);
  return 0;
}				/* s2102_ */

/* %3.1 */
/* Subroutine */ extern "C" int
s318_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc, integer * __restrict inc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;
  real r__1;

  /* Local variables */
  integer i__, k;
  integer index;
  real t1, t2;
  integer nl;
  real chksum, max__;


/*     reductions */
/*     isamax, max absolute value, increments not equal to 1 */


  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s318 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      k = 1;
      index = 1;
      max__ = ABS (a[1]);
      i__2 = *n;
#if defined _OPENMP 
#if _OPENMP >= 201307 && !defined __INTEL_COMPILER
#pragma omp simd lastprivate(index)
#endif
#endif
// icc fixed as of above build date; gcc PR60117
// leaving out reduction misleads compiler (breaks later icc)
// zero trip case has undefined result
      for (int i__ = 2; i__ <= i__2; ++i__)
	  if ((r__1 = a[k += *inc], ABS (r__1)) > max__) {
	      index = i__;
	      max__ = ABS (r__1);
	    }
      chksum = max__ + (real) index;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = max__ + (real) index;
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s318 ", (ftnlen) 5);
  return 0;
}				/* s318_ */

/* %3.2 */
/* Subroutine */ extern "C" int
s322_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b,
       real * c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     recurrences */
/*     second order linear recurrence */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s322 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      float tmp = a[2];
      i__2 = *n;
      for (i__ = 3; i__ <= i__2; ++i__)
	  a[i__] = tmp = (a[i__] + a[i__ - 2] * c__[i__]) + tmp * b[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 2);
  check_ (&chksum, &i__1, n, &t2, "s322 ", (ftnlen) 5);
  return 0;
}				/* s322_ */

/* %3.5 */
/* Subroutine */ extern "C" int
s352_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum, dot;


/*     loop rerolling */
/*     unrolled dot product */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s352 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes * 5;
  for (nl = 1; nl <= i__1; ++nl) {
      dot = 0.f;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; i__ += 5)
	  dot += a[i__] * b[i__] + a[i__ + 1] * b[i__ + 1] + a[i__ + 2]
	    * b[i__ + 2] + a[i__ + 3] * b[i__ + 3] + a[i__ + 4] * b[i__ + 4];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes * 5);
  chksum = dot;
  i__1 = *ntimes * 5 * (*n / 5);
  check_ (&chksum, &i__1, n, &t2, "s352 ", (ftnlen) 5);
  return 0;
}				/* s352_ */

/* %3.1 */
/* Subroutine */ extern "C" int
s3110_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * __restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;
  integer xindex, yindex;
  real max__;


/*     reductions */
/*     if to max with index reduction, 2 dimensions */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s3110", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      max__ = aa[aa_dim1 + 1];
      xindex = yindex = 1;
      i__2 = i__3 = *n;
#if defined _OPENMP && _OPENMP < 201307
#pragma omp parallel for if(i__2 > 103)
#else
#pragma omp parallel for if(i__2 > 103) reduction(max: max__) lastprivate(xindex,yindex)
#endif
      for (j = 1; j <= i__2; ++j) {
	  float *maxj=max_element(&aa[1+j*aa_dim1],&aa[i__3+1+j*aa_dim1]);
#if defined _OPENMP && _OPENMP < 201307
#pragma omp critical
#endif
	    if(*maxj > max__ || *maxj == max__ && j < yindex){
		max__= *maxj;
		xindex=maxj- &aa[1+j*aa_dim1] + 1;
		yindex=j;
		}
	}
      chksum = max__ + (real) xindex + (real) yindex;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = max__ + (real) xindex + (real) yindex;
  i__1 = *ntimes / *n * *n * *n;
  check_ (&chksum, &i__1, n, &t2, "s3110", (ftnlen) 5);
  return 0;
}				/* s3110_ */

/* %4.4 */
/* Subroutine */ extern "C" int
s442_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b,
       real * c__, real * d__, real * e, real * aa,
       real * bb, real * cc, integer * indx) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     non-logical if's */
/*     computed goto */

  /* Parameter adjustments */
  --indx;
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s442 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#if defined __INTEL_COMPILER
#pragma omp parallel for num_threads(2) if(i__2 > 103)
#endif
      for (i__ = 1; i__ <= i__2; ++i__)
	  switch (indx[i__]) {
	    case 1: a[i__] += b[i__] * b[i__];
	      continue;
	    case 2: a[i__] += c__[i__] * c__[i__];
	      continue;
	    case 3: a[i__] += d__[i__] * d__[i__];
	      continue;
	    case 4: a[i__] += e[i__] * e[i__];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s442 ", (ftnlen) 5);
  return 0;
}				/* s442_ */

/* %4.5 */
#include <math.h>
/* Subroutine */ extern "C" int
s451_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b,
       real * c__, real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     intrinsic functions */
/*     intrinsics */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s451 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#if _OPENMP < 201307 || !defined __INTEL_COMPILER && defined __AVX__
#pragma omp parallel for num_threads(4) if(i__2 > 103)
#else
#pragma omp parallel for simd num_threads(4) if(i__2 > 103)
#endif
      for (i__ = 1; i__ <= i__2; ++i__)
	  a[i__] = sinf (b[i__]) + cosf (c__[i__]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s451 ", (ftnlen) 5);
  return 0;
}				/* s451_ */

/* %4.9 */
/* Subroutine */ extern "C" int
s491_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * __restrict a, real * b,
       real * c__, real * d__, real * e, real * aa,
       real * bb, real * cc, integer * ip) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     vector semantics */
/*     indirect addressing on lhs, store in sequence */

  /* Parameter adjustments */
  --ip;
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s491 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#if defined _OPENMP && _OPENMP >= 201307
#pragma omp simd
#endif
      for (i__ = 1; i__ <= i__2; ++i__)
	  a[ip[i__]] = b[i__] + c__[i__] * d__[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s491 ", (ftnlen) 5);
  return 0;
}				/* s491_ */

/* %5.1 */
/* Subroutine */ extern "C" int
vpvpv_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * __restrict a, real * b,
	real * c__, real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     vector plus vector plus vector */

  /* Parameter adjustments */
  cc_dim1 = *ld;
  cc_offset = 1 + cc_dim1 * 1;
  cc -= cc_offset;
  bb_dim1 = *ld;
  bb_offset = 1 + bb_dim1 * 1;
  bb -= bb_offset;
  aa_dim1 = *ld;
  aa_offset = 1 + aa_dim1 * 1;
  aa -= aa_offset;
  --e;
  --d__;
  --c__;
  --b;
  --a;

  /* Function Body */
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "vpvpv", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	  a[i__] += b[i__] + c__[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vpvpv", (ftnlen) 5);
  return 0;
}				/* vpvpv_ */

