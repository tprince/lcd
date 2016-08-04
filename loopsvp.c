/* loops.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include "g2c.h"
#include <math.h>

/* Common Block Declarations */

struct {
  real array[1000000];
}
cdata_;

#define cdata_1 cdata_

/* Table of constant values */

static real c_b3 = 1.f;
static integer c__1 = 1;
static real c_b393 = 0.f;
  extern /* Subroutine */ int init_ (integer *, integer *, real *, real *,
	 real *, real *, real *, real *, real *, real *, char *, ftnlen),
    forttime_ (real *), set1d_ (integer *, real *, real *, integer *);
  extern /* Subroutine */ int check_ (real *, integer *, integer *, real *,
	  char *, ftnlen), dummy_ (integer *, integer *, real *, real *,
	   real *, real *, real *, real *, real *, real *, real *);
  extern real cs1d_ (integer *, real *);
  extern real cs2d_ (integer *, real *);

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
/* Subroutine */ int
#if ! __INTEL_COMPILER && __GNUC__
__attribute__((optimize("no-tree-vectorize"))) 
#endif
s111_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  extern /* Subroutine */ int init_ (integer *, integer *, real *, real *,
				     real *, real *, real *, real *, real *,
				     real *, char *, ftnlen),
    forttime_ (real *);
  integer i__;
  extern /* Subroutine */ int check_ (real *, integer *, integer *, real *,
				      char *, ftnlen), dummy_ (integer *,
							       integer *,
							       real *, real *,
							       real *, real *,
							       real *, real *,
							       real *, real *,
							       real *);
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     linear dependence testing */
/*     no dependence - vectorizable */

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
	 &bb[bb_offset], &cc[cc_offset], "s111 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes << 1;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
//vectorization should trigger "seems inefficient"
#if ! __KNC__ && _OPENMP >= 201307
#pragma omp for simd safelen(1)
#endif
      for (i__ = 2; i__ <= i__2; i__ += 2)
	  a[i__] = a[i__ - 1] + b[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes << 1);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes << 1) * (*n / 2);
  check_ (&chksum, &i__1, n, &t2, "s111 ", (ftnlen) 5);
  return 0;
}				/* s111_ */

/* %1.1 */
/* Subroutine */ int
s114_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real *  bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs2d_ (integer *, real *);


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
#pragma omp parallel for if(i__2 > 103)
      for (j = 2; j <= i__2; ++j) {
	  int i__3 = j - 1;
#pragma unroll(0)
#if !  __KNC__ && _OPENMP >= 201307
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

/* %1.1 */
/* Subroutine */ int
s119_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real *  bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     linear dependence testing */
/*     no dependence - vectorizable */

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
	 &bb[bb_offset], &cc[cc_offset], "s119 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
      for (j = 2; j <= i__2; ++j) {
//OK if i__3 <= aa_dim1 or aa_dim1 >= 64
#if !  __KNC__ && _OPENMP >= 201307
#pragma omp simd safelen(32)
#endif
	    for (i__ = 2; i__ <= i__3; ++i__) 
		aa[i__ + j * aa_dim1] = aa[i__ - 1 + (j - 1) * aa_dim1] + bb[
			i__ + j * bb_dim1];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes / *n * (*n - 1) * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s119 ", (ftnlen) 5);
  return 0;
}				/* s119_ */

/* %1.2 */
/* Subroutine */ int
s122_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc, integer * n1,
       integer * n3) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j, k;
  real t1, t2;
  integer nl;
  real chksum;


/*     induction variable recognition */
/*     variable lower and upper bound, and stride */

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
	 &bb[bb_offset], &cc[cc_offset], "s122 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
    int nk = *n - (*n - *n1)/ *n3;
      j = 1;
      k = 0;
      i__2 = *n;
      i__3 = *n3;
#if _OPENMP >= 201307
#pragma omp simd
#endif
#if !  __INTEL_COMPILER 
// original loop direction
      for (i__ = *n1; i__ <= i__2; i__ += i__3)
	  a[i__] += b[*n - (k += j) + 1];
#else
// reversed for Intel Xeon compilers
      for (i__ = i__2; i__ >= *n1; i__ -= i__3)
	  a[i__] += b[nk++];
#endif
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s122 ", (ftnlen) 5);
  return 0;
}				/* s122_ */

/* %1.2 */
/* Subroutine */ int
s125_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real *  aa, real *  bb, real *  cc) {
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
      i__2 = i__3 = *n;
#pragma omp parallel for if(i__2 > 103)
      for (j = 1; j <= i__2; ++j) {
	  int k = (j-1)*i__3;
#pragma vector nontemporal
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

#if defined __AVX__
#include <immintrin.h>
#elif defined __SSE2__
#include <xmmintrin.h>
#endif
/* %1.2 */
/* Subroutine */ int
s126_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * restrict bb, real *  cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j, k;
  real t1, t2;
  integer nl;
  real chksum;


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
#if defined __AVX2__
// actually this code will run (slower) on SNB, maybe OK on IVB
#pragma omp parallel for if(i__2 > 103)
      for (i__ = 1; i__ <= i__2; i__ += 8) {
	  int k = i__ * i__3 - i__3;
	  __m256 tmp = _mm256_loadu_ps(&bb[i__ + bb_dim1]);
	  for (int j = 2; j <= i__3; ++j){
	      __m256 tmp1 = _mm256_set_ps(cdata_1.array[k+7*i__3],
		  cdata_1.array[k+6*i__3],cdata_1.array[k+5*i__3],
		  cdata_1.array[k+4*i__3],cdata_1.array[k+3*i__3],
		  cdata_1.array[k+2*i__3],cdata_1.array[k+1*i__3],
		  cdata_1.array[k+0*i__3]);
	      tmp=_mm256_add_ps(tmp,_mm256_mul_ps(tmp1,
	       _mm256_loadu_ps(&cc[i__ + j * cc_dim1])));
	      // this might break if 32-byte alignment isn't supported
	      _mm256_store_ps(&bb[i__ + j * bb_dim1],tmp);
	      ++k;
	      }
	  }
#else
#if __SSE3__ && !  __INTEL_COMPILER
#pragma omp parallel for if(i__2 > 103)
      for (i__ = 1; i__ <= i__2; i__ += 4) {
	  int k = i__ * i__3 - i__3;
	  __m128 tmp = _mm_loadu_ps(&bb[i__ + bb_dim1]);
	  for (int j = 2; j <= i__3; ++j){
	      __m128 tmp1 = _mm_set_ps(cdata_1.array[k+3*i__3],
		  cdata_1.array[k+2*i__3],cdata_1.array[k+1*i__3],
		  cdata_1.array[k+0*i__3]);
	      __m128 tmp2 = _mm_loadu_ps(&cc[i__ + j * cc_dim1]);
	      tmp=_mm_add_ps(tmp,_mm_mul_ps(tmp1,tmp2));
	      _mm_store_ps(&bb[i__ + j * bb_dim1],tmp);
	      ++k;
	      }
	  }
#ifdef __AVX__
_mm256_zeroupper();	// icpc uses AVX-128 but not g++
#endif
#else
#ifndef __SUNPRO_CC
#warning "SSE3/4/AVX unseen, dropping to C source"
#endif
#if _OPENMP && _OPENMP < 201307 || __AVX__ && !  __INTEL_COMPILER
#pragma omp parallel for if(i__2 > 103)
#else
#pragma omp parallel for simd if(i__2 > 103)
#endif
      for (i__ = 1; i__ <= i__2; ++i__ ) {
	  int k = i__ * i__3 - i__3;
	  for (int j = 2; j <= i__3; ++j)
	      bb[i__ + j * bb_dim1] = bb[i__ + (j - 1) * bb_dim1] +
		cdata_1.array[k++] * cc[i__ + j * cc_dim1];
	}
#endif
#endif
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
/* Subroutine */ int
s132_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real *  b, real *  c__,
       real * d__, real * e, real * restrict aa, real * bb, real * cc) {
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
#if  _OPENMP >= 201307
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
/* Subroutine */ int
s141_ (integer * ntimes, integer * ld, integer * n, real *
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
#pragma omp parallel if(i__2 > 103)
      {
      float avgchunk=(i__2+1)*i__2/2;
      int nt;
      if(i__2 > 103){
#if defined _OPENMP
	  nt=omp_get_num_threads();
// fair number of array elements per thread
	  avgchunk=avgchunk/nt;
#else
	  nt=1;
#endif
}	else
	  nt=1;
#pragma omp for
      for(int m= 1;m <= nt; ++m){
	  int jbeg=sqrtf(.25f+2*(m-1)*avgchunk)+1;
// closest approximation to targeted chunk size
	  int jend=sqrtf(.25f+2*m*avgchunk);
	  for (int j = jbeg; j <= jend; ++j) {
	      int k = j * (j - 1) / 2;
	      for (int i__ = 1; i__ <= j; ++i__)
		  cdata_1.array[i__ + k - 1] += bb[i__ + j * bb_dim1];
	    }
	}
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
/* Subroutine */ int
s162_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b,
       real *  c__, real * d__, real * e, real * aa, real * bb,
       real * cc, integer *  k) {
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
	  i__2 = *n - 1;	// should be *n - *k
// the point of this benchmark is to see the compiler using the conditional
// to resolve overlap direction (if not done at run time)
#pragma omp simd
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
/* Subroutine */ int
s173_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, k;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     symbolics */
/*     expression in loop bounds and subscripts */

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
  k = *n / 2;
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s173 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes << 1;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = k;
//clearly no overlap, so ivdep should not be needed
#pragma vector unaligned
#pragma omp simd
      for (i__ = 1; i__ <= i__2; ++i__)
	  a[i__ + k] = a[i__] + b[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes << 1);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes << 1) * (*n / 2);
  check_ (&chksum, &i__1, n, &t2, "s173 ", (ftnlen) 5);
  return 0;
}				/* s173_ */


/* %1.7 */
/* Subroutine */ int
s175_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc, integer * inc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


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
      int len;
      i__2 = *n - *inc;
      i__3 = *inc;
#if __MIC__ || __AVX2__
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

/* %2.3 */
/* Subroutine */ int
s232_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real *  bb, real * cc) {
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
#pragma omp parallel if(i__2 > 103)
      {
      float avgchunk=(i__2+1)*i__2/2;
      int nt;
      if(i__2 > 103){
#if defined _OPENMP
	  nt=omp_get_num_threads();
// fair number of array elements per thread
	  avgchunk=avgchunk/nt;
#else
	  nt=1;
#endif
}	else
	  nt=1;
#pragma omp for
      for(int m= 1;m <= nt; ++m){
	  int jbeg=sqrtf(.25f+2*(m-1)*avgchunk)+2;
// closest approximation to targeted chunk size
	  int jend=sqrtf(.25f+2*m*avgchunk)+1;
	  int j;
	  for (j = jbeg; j < jend; j += 2) {
	      float tmp = aa[1 + j*aa_dim1],tmp2 = aa[1 + (j+1)*aa_dim1];
	      for (int i__ = 2; i__ <= j; ++i__){
		  aa[i__ + j * aa_dim1] = tmp = tmp * tmp
		    + bb[i__ + j * bb_dim1];
		  aa[i__ + j * aa_dim1+aa_dim1] = tmp2 = tmp2 * tmp2
		    + bb[i__ + j * bb_dim1+bb_dim1];
		  }
	      aa[j + 1 + j * aa_dim1+aa_dim1] = tmp2 * tmp2
		+ bb[j + 1 + j * bb_dim1+bb_dim1];
	    }
	  if(j == jend)
	      for (int i__ = 2; i__ <= j; ++i__)
		  aa[i__ + j * aa_dim1] = aa[i__ + j * aa_dim1-1] * 
		    aa[i__ + j * aa_dim1-1] + bb[i__ + j * bb_dim1];
	}
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
/* Subroutine */ int
s233_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real *  cc) {
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
#pragma novector
      for (int j = 2; j <= i__3; ++j)
	  for (int i__ = 2; i__ <= i__2; ++i__)
	      bb[i__ + j * bb_dim1] = bb[i__ - 1 + j*bb_dim1] +
		    cc[i__ + j * cc_dim1];

#if __MIC__
#pragma omp for simd
	  for (i__ = 2; i__ <= i__2; ++i__)
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

/* %2.4 */
/* Subroutine */ int
s244_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real *  c__, real *  d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     node splitting */
/*     false dependence cycle breaking */

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
	 &bb[bb_offset], &cc[cc_offset], "s244 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
#pragma vector aligned
#pragma omp simd // aligned(a,b,c__,d__)
      for (i__ = 1; i__ <= i__2; ++i__) {
	  a[i__] = b[i__] + c__[i__] * d__[i__];
	  b[i__] += c__[i__];
	}
	a[i__] = b[i__ - 1] + a[i__] * d__[i__ - 1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s244 ", (ftnlen) 5);
  return 0;
}				/* s244_ */

/* %2.5 */
/* Subroutine */ int
s251_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b,
       real *  c__, real *  d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real s;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     scalar and array expansion */
/*     scalar expansion */

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
	 &bb[bb_offset], &cc[cc_offset], "s251 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma vector aligned
#pragma omp simd // aligned(a,b,c__,d__)
      for (i__ = 1; i__ <= i__2; ++i__) {
	  s = b[i__] + c__[i__] * d__[i__];
	  a[i__] = s * s;
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s251 ", (ftnlen) 5);
  return 0;
}				/* s251_ */


/* %2.5 */
/* Subroutine */ int
s257_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *  a, real * b, real * c__,
       real * d__, real * e, real * restrict aa, real *  bb,
       real * cc) {
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
#pragma novector
#pragma unroll(2)
      for (i__ = 2; i__ <= i__2; ++i__) {
	  a[i__] = aa[i__ + i__3 * aa_dim1] - a[i__ - 1];
	  aa[i__ + i__3 * aa_dim1] = a[i__] + bb[i__ + i__3 * bb_dim1];
	}
#pragma omp for nowait
      for (j = 1; j < i__3; ++j){
	  float tmp= a[1];
#pragma unroll(2)
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

/* %2.5 */
/* Subroutine */ int
s258_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *  a, real * restrict b,
       real *  c__, real *  d__, real * restrict e,
       real *  aa, real * bb, real * cc)
{
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real s,sv[*n];
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     scalar and array expansion */
/*     wrap-around scalar under an if */

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
	 &bb[bb_offset], &cc[cc_offset], "s258 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      s = 0.f;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	  sv[i__ - 1] = s = (a[i__] > 0.f)? d__[i__] * d__[i__]: s;
#pragma omp simd
      for (i__ = 1; i__ <= i__2; ++i__) {
	  b[i__] = sv[i__ - 1] * c__[i__] + d__[i__];
	  e[i__] = sv[i__ - 1] * aa[i__ + aa_dim1]+ aa[i__ + aa_dim1];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &b[1]) + cs1d_ (n, &e[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s258 ", (ftnlen) 5);
  return 0;
}				/* s258_ */

/* %2.7 */
/* Subroutine */ int
s271_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b,
       real *  c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control flow */
/*     loop with singularity handling */

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
	 &bb[bb_offset], &cc[cc_offset], "s271 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma vector aligned
#pragma omp simd // aligned(a,b,c__)
      for (i__ = 1; i__ <= i__2; ++i__)
	  a[i__] += (b[i__] > 0.f ? b[i__] : 0.f) * c__[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s271 ", (ftnlen) 5);
  return 0;
}				/* s271_ */

/* %2.7 */
/* Subroutine */ int
s275_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real *  bb, real *  cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs2d_ (integer *, real *);


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
#if __MIC__ || __AVX2__
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
/* Subroutine */ int
s281_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real *  c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
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
#pragma omp simd
    for (i__= 1; i__ <= (i__2+1)/2; ++i__)
	a[i__] = (b[i__] = a[i__2 - i__ + 1] + b[i__] * c__[i__])- 1.f;
#pragma omp simd
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
/* Subroutine */ int
s2101_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * restrict aa, real *  bb, real *  cc) {
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
      aa_dim1++,bb_dim1++,cc_dim1++;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#if _OPENMP >= 201307 && (__INTEL_COMPILER || !  __AVX__)
#pragma omp parallel for simd if(i__2 > 103)
#endif
      for (i__ = 1; i__ <= i__2; ++i__)
	  aa[i__ * aa_dim1] += bb[i__ * bb_dim1] * cc[i__ * cc_dim1];
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
/* Subroutine */ int
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
	for (j = 1; j <= i__2; ++j) {
	    for(int i = 1; i <= i__2; ++i)
		aa[i + j * aa_dim1] = i==j;
	}
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
/* Subroutine */ int
s318_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc, integer * inc) {
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
#if _OPENMP >= 201307
#pragma omp simd lastprivate(index) reduction(max: max__)
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

/* %3.1 */
/* Subroutine */ int
s3110_ (integer * ntimes, integer * ld, integer * n, real *
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
#if _OPENMP && _OPENMP < 201307 || __KNC__
#pragma omp parallel for if(i__2 > 103)
#else
#pragma omp parallel for if(i__2 > 103) reduction(max: max__) lastprivate(xindex,yindex)
#endif
      for (j = 1; j <= i__2; ++j) {
	  int indxj=0;
	  float maxj=max__;
#if __INTEL_COMPILER
// this is risky (and breaks on MIC with compiler 15.0 beta) but fast when works
#pragma simd firstprivate(maxj,indxj) lastprivate(maxj,indxj)
// gcc omp simd reduction lastprivate shows a small loss on AVX
#elif _OPENMP >= 201307 
#pragma omp simd reduction(max: maxj) lastprivate(indxj)
#endif
	  for (int i__ = 1; i__ <= i__3; ++i__)
	      if (aa[i__ + j * aa_dim1] > maxj){
		  maxj = aa[i__ + j * aa_dim1];
		  indxj = i__;
		  }
#if _OPENMP && _OPENMP < 201307 || __KNC__
#pragma omp critical
#endif
	    if(maxj > max__) {
		max__= maxj;
		xindex=indxj;
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

/* %3.2 */
/* Subroutine */ int
s322_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b,
       real *  c__, real * d__, real * e, real * aa, real * bb,
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

/* %3.2 */
/* Subroutine */ int
s323_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real *  c__, real *  d__, real *  e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     recurrences */
/*     coupled recurrence */

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
	 &bb[bb_offset], &cc[cc_offset], "s323 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      float tmp = b[1];
      i__2 = *n;
      for (i__ = 2; i__ <= i__2; ++i__) {
	  a[i__] = tmp + c__[i__] * d__[i__];
	  b[i__] =  tmp += c__[i__] * e[i__]+ c__[i__] *d__[i__];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s323 ", (ftnlen) 5);
  return 0;
}				/* s323_ */

/* %3.5 */
/* Subroutine */ int
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
      // 50% gain for observing ordering on Core i7
#if !defined __KNC__
#pragma novector
#endif
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

/* %4.1 */
/* Subroutine */ int
s413_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real *  c__, real *  d__, real *  e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop recognition */
/*     if loop to do loop, code on both sides of increment */

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
	 &bb[bb_offset], &cc[cc_offset], "s413 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__ = 1;
      b[i__] += d__[i__] * e[i__];
#pragma omp simd
      for (i__ = 2; i__ < *n; ++i__) {
	  b[i__] += d__[i__] * e[i__];
	  a[i__] = c__[i__] + d__[i__] * e[i__];
	}
      a[i__] = c__[i__] + d__[i__] * e[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s413 ", (ftnlen) 5);
  return 0;
}				/* s413_ */

/* %4.2 */
/* Subroutine */ int
s422_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *  a, real * b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
#define x ((real *)&cdata_1 + 4)
  real t1, t2;
  integer nl;
  real chksum;


/*     storage classes and equivalencing */
/*     common and equivalence statement */
/*     anti-dependence, threshold of 4 */

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
  set1d_ (n, x, &c_b393, &c__1);
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s422 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
// the point of the benchmark is for the compiler to analyze the source
// and see that vectorization with source/destination overlap is safe.
// gcc 4.7 does this.  ifort does this.  icc/icl 12.x do not.
#pragma omp simd
      for (i__ = 1; i__ <= i__2; ++i__)
	  x[i__ - 1] = cdata_1.array[i__ + 7] + a[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, x);
  i__1 = *ntimes * (*n - 8);
  check_ (&chksum, &i__1, n, &t2, "s422 ", (ftnlen) 5);
  return 0;
}				/* s422_ */

#undef x


/* %4.2 */
/* Subroutine */ int
s423_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *  a, real * b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
#define x ((real *)&cdata_1 + 63)
  real t1, t2;
  integer nl;
  real chksum;


/*     storage classes and equivalencing */
/*     common and equivalenced variables - with anti-dependence */

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
  set1d_ (n, x, &c_b3, &c__1);
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s423 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
// vectorlength assertion nearly OK even if loop were reversed
#pragma omp simd safelen(32)
      for (i__ = 1; i__ <= i__2; ++i__)
	  cdata_1.array[i__] = x[i__ - 1] + a[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, cdata_1.array);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s423 ", (ftnlen) 5);
  return 0;
}				/* s423_ */

#undef x


/* %4.2 */
/* Subroutine */ int
s424_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *  a, real * b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
#define x ((real *)&cdata_1 + 63)
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);

/*     storage classes and equivalencing */
/*     common and equivalenced variables - overlap */
/*     vectorizeable in strips of 64 or less */

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
  set1d_ (n, x, &c_b393, &c__1);
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s424 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
// vectorlength assertion is OK by inspection
#pragma omp simd safelen(32)
      for (i__ = 1; i__ <= i__2; ++i__)
	  x[i__] = cdata_1.array[i__ - 1] + a[i__];
      dummy_ (ld, n, x, &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, x);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s424 ", (ftnlen) 5);
  return 0;
}				/* s424_ */

#undef x

/* %4.4 */
/* Subroutine */ int
s441_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b,
       real *  c__, real *  d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     non-logical if's */
/*     arithmetic if translated automatically to if..else */

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
	 &bb[bb_offset], &cc[cc_offset], "s441 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma vector aligned
      for (i__ = 1; i__ <= i__2; ++i__)
	a[i__] +=(d__[i__] <= 0.f?b[i__]:c__[i__]) *
		(d__[i__]==0.f?b[i__]:c__[i__]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s441 ", (ftnlen) 5);
  return 0;
}				/* s441_ */
/* %4.4 */
/* Subroutine */ int
s442_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b,
       real *  c__, real *  d__, real *  e, real * aa,
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
#ifdef _OPENMP
    int nt= 2;
#endif
      i__2 = *n;
#if defined __INTEL_COMPILER
#pragma omp parallel for num_threads(nt) if(i__2 > 103)
#if !defined __KNC__
#pragma novector
#endif
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
/* Subroutine */ int
s451_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b,
       real *  c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


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
#ifdef _OPENMP
    int maxt=omp_get_max_threads();
    int nt= maxt>6? 6:maxt;
#endif
      i__2 = *n;
#if _OPENMP < 201307 || (! __INTEL_COMPILER) &&  __AVX__
#pragma omp parallel for num_threads(nt) if(i__2 > 103)
#else
#pragma omp parallel for simd num_threads(nt) if(i__2 > 103)
#pragma vector nontemporal
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
/* Subroutine */ int
s491_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real *  b,
       real *  c__, real *  d__, real * e, real * aa,
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
#pragma omp simd
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
/* Subroutine */ int
vpvpv_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real *  b,
	real *  c__, real * d__, real * e, real * aa, real * bb,
	real * cc) {
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

