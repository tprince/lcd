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
#include <cilk/cilk.h>

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
s111_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
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
#ifndef __MIC__
#pragma novector
#pragma unroll(4)
#endif
      a[2:i__2/2:2] = a[1:i__2/2:2] + b[2:i__2:2];
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
#if defined __SSE2__
#include <xmmintrin.h>
#endif
#if defined __AVX__
#include <immintrin.h>
#endif
/* Subroutine */ int
s112_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     linear dependence testing */
/*     loop reversal */

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
	 &bb[bb_offset], &cc[cc_offset], "s112 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      // loop must vectorize backwards on account of data overlap
#if 1
#ifndef __MIC__
#pragma novector
#endif
	  a[*n:*n-1:-1]= a[*n-1:*n-1:-1] + b[*n-1:*n-1:-1];
#else
      for (i__ = *n - 1; i__ >= 1; --i__)
	  a[i__ + 1] = a[i__] + b[i__];
#endif
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s112 ", (ftnlen) 5);
  return 0;
}				/* s112_ */

/* %1.1 */
/* Subroutine */ int
s113_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     linear dependence testing */
/*     a(i)=a(1) but no actual dependence cycle */

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
	 &bb[bb_offset], &cc[cc_offset], "s113 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[2:*n-1] = a[1] + b[2:*n-1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s113 ", (ftnlen) 5);
  return 0;
}				/* s113_ */

/* %1.1 */
/* Subroutine */ int
s114_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * cc) {
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
      cilk_for (int j = 2; j <= i__2; ++j)
	  aa[1+j*aa_dim1:j-1]=aa[j+aa_dim1:j-1:aa_dim1]+bb[1+j*aa_dim1:j-1];
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
s115_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *restrict a, real * b, real * c__, real * d__,
       real * e, real *restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     linear dependence testing */
/*     triangular saxpy loop */

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
	 &bb[bb_offset], &cc[cc_offset], "s115 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = (*ntimes << 1) / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
	for (j = 2; j <= i__2; ++j) {
	    i__3 = j-1;
	    a[j] -= __sec_reduce_add(aa[1+j*aa_dim1:i__3]*a[1:i__3]);
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n << 1);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes / *n << 1) * ((*n * *n - *n) / 2);
  check_ (&chksum, &i__1, n, &t2, "s115 ", (ftnlen) 5);
  return 0;
}				/* s115_ */

/* %1.1 */
/* Subroutine */ int
s116_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     linear dependence testing */

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
	 &bb[bb_offset], &cc[cc_offset], "s116 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes * 5;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 5;
      for (i__ = 1; i__ <= i__2; i__ += 5) {
	  a[i__] *= a[i__ + 1];
	  a[i__ + 1] *= a[i__ + 2];
	  a[i__ + 2] *= a[i__ + 3];
	  a[i__ + 3] *= a[i__ + 4];
	  a[i__ + 4] *= a[i__ + 5];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes * 5);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * 5 * (*n / 5);
  check_ (&chksum, &i__1, n, &t2, "s116 ", (ftnlen) 5);
  return 0;
}				/* s116_ */

/* %1.1 */
/* Subroutine */ int
s118_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * b, real * c__,
       real * d__, real * e, real * aa, real * restrict bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     linear dependence testing */
/*     potential dot product recursion */

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
	 &bb[bb_offset], &cc[cc_offset], "s118 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = (*ntimes << 1) / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 2; i__ <= i__2; ++i__) {
	  float sum = 0;
	  i__3 = i__ - 1;
	  a[i__] += __sec_reduce_add(bb[i__+i__3*bb_dim1:i__3:-bb_dim1]*a[i__-i__3:i__3]);
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n << 1);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes / *n << 1) * ((*n * *n - *n) / 2);
  check_ (&chksum, &i__1, n, &t2, "s118 ", (ftnlen) 5);
  return 0;
}				/* s118_ */

/* %1.1 */
/* Subroutine */ int
s119_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * cc) {
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
      i__2 = *n;
      for (j = 2; j <= i__2; ++j) {
	  i__3 = *n;
	aa[2+j*aa_dim1:i__3-1]=aa[1+(j-1)*aa_dim1:i__3-1]+bb[2+j*aa_dim1:i__3-1];
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
s121_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     induction variable recognition */
/*     loop with possible ambiguity because of scalar store */

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
	 &bb[bb_offset], &cc[cc_offset], "s121 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
      a[1:i__2]=a[2:i__2]+b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s121 ", (ftnlen) 5);
  return 0;
}				/* s121_ */

/* %1.2 */
/* Subroutine */ int
s122_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
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
      j = 1;
      k = 0;
      i__2 = *n;
      i__3 = *n3;
      a[*n1:(i__2-*n1)/i__3+1:i__3] += b[*n:(i__2-*n1)/i__3+1: -j];
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
s123_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     induction variable recognition */
/*     induction variable under an if */

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
	 &bb[bb_offset], &cc[cc_offset], "s123 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes << 1;
  for (nl = 1; nl <= i__1; ++nl) {
      j = 0;
      i__2 = *n / 2;
      for (i__ = 1; i__ <= i__2; ++i__) {
	  a[++j] = b[i__] + d__[i__] * e[i__];
	  if (c__[i__] > 0.f)
	      a[++j] = c__[i__] + d__[i__] * e[i__];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes << 1);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes << 1) * (*n / 2);
  check_ (&chksum, &i__1, n, &t2, "s123 ", (ftnlen) 5);
  return 0;
}				/* s123_ */

/* %1.2 */
/* Subroutine */ int
s124_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     induction variable recognition */
/*     induction variable under both sides of if (same value) */

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
	 &bb[bb_offset], &cc[cc_offset], "s124 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes << 1;
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(e+1,64);
  __assume_aligned(a+1,64);
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n / 2;
	a[1:i__2] = (b[1:i__2] > 0.f?
	    b[1:i__2]:c__[1:i__2]) + d__[1:i__2] * e[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes << 1);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes << 1) * (*n / 2);
  check_ (&chksum, &i__1, n, &t2, "s124 ", (ftnlen) 5);
  return 0;
}				/* s124_ */

/* %1.2 */
/* Subroutine */ int
s125_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * restrict cc) {
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
      cilk_for (int j = 1; j <= i__2; ++j) {
	  int k = (j-1)*i__3;
	  cdata_1.array[k:i__3]=aa[1+j*aa_dim1:i__3]+bb[1+j*bb_dim1:i__3]*cc[1+j*cc_dim1:i__3];
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
/* Subroutine */ int
s126_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * restrict bb, real * restrict cc) {
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
      cilk_for (int i__ = 1; i__ <= i__2; i__ += 8) {
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
	      // this will break if 32-byte alignment isn't supported
	      _mm256_store_ps(&bb[i__ + j * bb_dim1],tmp);
	      ++k;
	      }
	  }
#else
#if defined __SSE2__
      cilk_for (int i__ = 1; i__ <= i__2; i__ += 4) {
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
#else
      cilk_for (int i__ = 1; i__ <= i__2; ++i__ ) {
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

/* %1.2 */
/* Subroutine */ int
s127_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     induction variable recognition */
/*     induction variable with multiple increments */

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
	 &bb[bb_offset], &cc[cc_offset], "s127 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes << 1;
  for (nl = 1; nl <= i__1; ++nl) {
      j = 0;
      i__2 = *n / 2;
      a[1:i__2:2] = b[1:i__2] + c__[1:i__2] * d__[1:i__2];
      a[2:i__2:2] = b[1:i__2] + d__[1:i__2] * e[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes << 1);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes << 1) * (*n / 2);
  check_ (&chksum, &i__1, n, &t2, "s127 ", (ftnlen) 5);
  return 0;
}				/* s127_ */

/* %1.2 */
/* Subroutine */ int
s128_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j, k;
  real t1, t2;
  integer nl;
  real chksum;


/*     induction variables */
/*     coupled induction variables */

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
	 &bb[bb_offset], &cc[cc_offset], "s128 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes << 1;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n / 2;
#ifndef __MIC__
#pragma novector
#endif
      b[1:i__2:2] = (a[1:i__2] = b[1:i__2:2] - d__[1:i__2]) + c__[1:i__2:2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes << 1);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = (*ntimes << 1) * (*n / 2);
  check_ (&chksum, &i__1, n, &t2, "s128 ", (ftnlen) 5);
  return 0;
}				/* s128_ */

/* %1.3 */
/* Subroutine */ int
s131_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, m;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     global data flow analysis */
/*     forward substitution */

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
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s131 ", (ftnlen) 5);
  if (a[1] > 0.f)
      a[1] = b[1];
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - m;
      a[1:i__2] = a[m+1:i__2] + b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s131 ", (ftnlen) 5);
  return 0;
}				/* s131_ */

/* %1.3 */
/* Subroutine */ int
s132_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * restrict b, real * restrict c__,
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
	aa[2+j*aa_dim1:i__2-1] = aa[1+k*aa_dim1:i__2-1] + b[2:i__2-1] * c__[k];
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
       real * e, real * aa, real *restrict bb, real * cc) {
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
      cilk_for(int j = 1; j <= *n; ++j) {
	  int k = j * (j - 1) / 2;
	    cdata_1.array[k:j] += bb[1 + j * bb_dim1:j];
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

/* %1.5 */
/* Subroutine */ int
s151_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1;

  /* Local variables */
  real t1, t2;
  integer nl;
  real chksum;
  extern /* Subroutine */ int s151s_ (real *, real *, integer *, integer );


/*     interprocedural data flow analysis */
/*     passing parameter information into a subroutine */

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
	 &bb[bb_offset], &cc[cc_offset], "s151 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      s151s_ (&a[1], &b[1], n, 1);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s151 ", (ftnlen) 5);
  return 0;
}				/* s151_ */

/* Subroutine */ int
s151s_ (real * restrict a, real * restrict b, integer * n, integer  m) {
  /* System generated locals */
  integer i__1;

  /* Local variables */
  integer i__;

  /* Parameter adjustments */
  --b;
  --a;

  /* Function Body */
  i__1 = *n - m;
  a[1:i__1] = a[1 + m:i__1] + b[1:i__1];
  return 0;
}				/* s151s_ */

/* %1.5 */
/* Subroutine */ int
s152_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);
  extern /* Subroutine */ int s152s_ (real *, real *, real *, integer *);


/*     interprocedural data flow analysis */
/*     collecting information from a subroutine */

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
	 &bb[bb_offset], &cc[cc_offset], "s152 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  __assume_aligned(b+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(e+1,64);
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma omp simd
      for (i__ = 1; i__ <= i__2; ++i__) {
	  b[i__] = d__[i__] * e[i__];
	  s152s_ (&a[1], &b[1], &c__[1], &i__);
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s152 ", (ftnlen) 5);
  return 0;
}				/* s152_ */

/* Subroutine */ int
s152s_ (real * restrict a, real * restrict b, real * restrict c__,
	integer * i__) {
  /* Parameter adjustments */
  --c__;
  --b;
  --a;

  /* Function Body */
  a[*i__] += b[*i__] * c__[*i__];
  return 0;
}				/* s152s_ */

/* %1.6 */
/* Subroutine */ int
s161_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(e+1,64);
  __assume_aligned(a+1,64);
      a[1:i__2] = (b[1:i__2] >= 0.f)? c__[1:i__2] + d__[1:i__2] * e[1:i__2]: a[1:i__2];
      c__[2:i__2] = (b[1:i__2] >= 0.f)? c__[2:i__2] : d__[1:i__2] * d__[1:i__2]+ a[1:i__2];
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
/* Subroutine */ int
s162_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
       real * cc, integer * restrict k) {
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
	  i__2 = *n - *k;
	  a[1:i__2] = a[1 + *k:i__2] + b[1:i__2] * c__[1:i__2];
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
s171_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc, integer * inc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

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
	 &bb[bb_offset], &cc[cc_offset], "s171 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  __assume_aligned(b+1,64);
  __assume_aligned(a+1,64);
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
	a[*inc:i__2:*inc] += b[1:i__2];
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
/* Subroutine */ int
s172_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
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
  extern real cs1d_ (integer *, real *);


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
      a[*n1:(i__2-*n1)/i__3+1:i__3] += b[*n1:(i__2-*n1)/i__3+1:i__3];
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
      i__2 = *n / 2;
	a[1 + i__2:i__2] = a[1:i__2] + b[1:i__2];
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
s174_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *restrict a, real *restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     symbolics */
/*     loop with subscript that may seem ambiguous */

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
	 &bb[bb_offset], &cc[cc_offset], "s174 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes << 1;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n / 2;
	a[1:i__2] = a[1 + i__2:i__2] + b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes << 1);
  chksum = cs1d_ (n, &a[1]);
  i__1 = (*ntimes << 1) * (*n / 2);
  check_ (&chksum, &i__1, n, &t2, "s174 ", (ftnlen) 5);
  return 0;
}				/* s174_ */

/* %1.7 */
/* Subroutine */ int
s175_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
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
      i__2 = *n - 1;
      i__3 = *inc;
      a[1:i__2/i__3:i__3] = a[1 + i__3:i__2/i__3:i__3] +
		 b[1:i__2/i__3:i__3];
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
/* Subroutine */ int
s176_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
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
  for (nl = 1; nl <= i__1; ++nl) {
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(a+1,64);
	if(m > 201)cilk_for (int i__ = 1; i__ <= m; ++i__) 
	    a[i__] += __sec_reduce_add(b[i__:m]*c__[m:m:-1]);
	else for (int i__ = 1; i__ <= m; ++i__){
	    float tmp = c__[i__];
	    a[1:m] += b[1+m-i__:m]*tmp;
	    }
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
/* Subroutine */ int
s211_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
       real * bb, real * cc) {
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
      i__2 = *n - 3;
      // optimized for Core i7 and MIC
      a[2] = b[1] + c__[2] * d__[2];
      a[3:i__2] = ( b[2:i__2] = b[3:i__2] - e[2:i__2] * d__[2:i__2] )
	  + c__[3:i__2] * d__[3:i__2];
      i__2 += 2;
      b[i__2] = b[i__2 + 1] - e[i__2] * d__[i__2];
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

/* %2.1 */
/* Subroutine */ int
s212_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     statement reordering */
/*     dependency needing temporary */

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
	 &bb[bb_offset], &cc[cc_offset], "s212 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
      b[1:i__2] += a[2:i__2] * d__[1:i__2];
      a[1:i__2] *= c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s212 ", (ftnlen) 5);
  return 0;
}				/* s212_ */

/* %2.2 */
/* Subroutine */ int
s221_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop distribution */
/*     loop that is partially recursive */

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
	 &bb[bb_offset], &cc[cc_offset], "s221 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      float tmp = b[1];
      i__2 = *n;
      a[2:i__2-1] += c__[2:i__2-1] * d__[2:i__2-1];
      for (i__ = 2; i__ <= i__2; ++i__)
	  b[i__] = (a[i__] + d__[i__]) + b[i__ - 1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s221 ", (ftnlen) 5);
  return 0;
}				/* s221_ */

/* %2.2 */
/* Subroutine */ int
s222_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop distribution */
/*     partial loop vectorization, recurrence in middle */

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
	 &bb[bb_offset], &cc[cc_offset], "s222 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[2:i__2-1] += b[2:i__2-1] * c__[2:i__2-1];
//#pragma nofusion	//fusion breaks it here
      for (i__ = 2; i__ <= i__2; ++i__)
	  b[i__] = b[i__ - 1] * b[i__ - 1] * a[i__];
//#pragma nofusion
      a[2:i__2-1] -= b[2:i__2-1] * c__[2:i__2-1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s222 ", (ftnlen) 5);
  return 0;
}				/* s222_ */

/* %2.3 */
/* Subroutine */ int
s231_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs2d_ (integer *, real *);


/*     loop interchange */
/*     loop with multiple dimension recursion */

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
	 &bb[bb_offset], &cc[cc_offset], "s231 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
    for (int j = 2; j <= i__3; ++j)
	aa[1 + j * aa_dim1:i__2] = aa[1 + (j - 1) * aa_dim1:i__2] + bb[1 
		    + j * bb_dim1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes / *n * *n * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s231 ", (ftnlen) 5);
  return 0;
}				/* s231_ */

/* %2.3 */
/* Subroutine */ int
s232_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * cc) {
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
	  cilk_for(int j = 1; j < *n; j += 2) {
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
	  if(1 & *n)
	      for (int i__ = 2; i__ <= *n; ++i__)
		  aa[i__ + *n * aa_dim1] = aa[i__ + *n * aa_dim1-1] * 
		    aa[i__ + *n * aa_dim1-1] + bb[i__ + *n * bb_dim1];
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
       real * e, real * restrict aa, real * restrict bb, real * restrict cc) {
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
      for (j = 2; j <= i__3; ++j)
	  aa[2 + j * aa_dim1:i__2-1] = aa[2 + (j - 1) * aa_dim1:i__2-1] + cc[2
						       + j * cc_dim1:i__2-1];
      cilk_for (int j = 2; j <= i__3; ++j)
	  for (int i__ = 2; i__ <= i__2; ++i__) 
	      bb[i__ + j * bb_dim1] = bb[i__ - 1 + j*bb_dim1] +
		    cc[i__ + j * cc_dim1];
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
/* Subroutine */ int
s234_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * restrict cc) {
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
	for(int j = 2; j <= *n; ++j)
	    aa[1 + j * aa_dim1:*n] = aa[1 + (j - 1) * aa_dim1:*n] +
		bb[1 + (j - 1) * bb_dim1:*n] * cc[1 + (j - 1) * cc_dim1:*n];
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
/* Subroutine */ int
s235_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * b, real * c__,
       real * d__, real * e, real * restrict aa, real * restrict bb,
       real * cc) {
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
	a[1:i__2] += b[1:i__2] * c__[1:i__2];
      for (int j = 2; j <= i__3; ++j)
	  aa[1 + j * aa_dim1:i__2] = aa[1 + (j - 1) * aa_dim1:i__2] + bb[1
					       + j * bb_dim1:i__2] * a[1:i__2];
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

/* %2.4 */
/* Subroutine */ int
s241_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
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
/*     preloading necessary to allow vectorization */

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
	 &bb[bb_offset], &cc[cc_offset], "s241 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(a+1,64);
      for (i__ = 1; i__ <= i__2; ++i__) {
	  float tmp = a[i__ + 1] * d__[i__];
	  b[i__] = (a[i__] = b[i__] * c__[i__] * d__[i__]) * tmp;
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s241 ", (ftnlen) 5);
  return 0;
}				/* s241_ */

/* %2.4 */
/* Subroutine */ int
s242_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
       real * bb, real * cc, real * restrict s1, real * restrict s2) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     node splitting */

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
	 &bb[bb_offset], &cc[cc_offset], "s242 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
    float tmp[i__2 = *n];
      // optimized for Core i7
      tmp[1:i__2 - 1] = *s1 + *s2 + b[2:i__2-1] + c__[2:i__2-1] + d__[2:i__2-1];
      for (i__ = 2; i__ <= i__2; ++i__)
	  a[i__] = tmp[i__ - 1] + a[i__ - 1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s242 ", (ftnlen) 5);
  return 0;
}				/* s242_ */

/* %2.4 */
/* Subroutine */ int
s243_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
	 &bb[bb_offset], &cc[cc_offset], "s243 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
	a[1:i__2] = ( b[1:i__2] += (c__[1:i__2]+e[1:i__2]) * d__[1:i__2] ) +
	    a[2:i__2] * d__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s243 ", (ftnlen) 5);
  return 0;
}				/* s243_ */

/* %2.4 */
/* Subroutine */ int
s244_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
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
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(a+1,64);
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
      a[1:i__2] = b[1:i__2] + c__[1:i__2] * d__[1:i__2];
      b[1:i__2] += c__[1:i__2];
	a[*n] = b[*n - 1] + a[*n] * d__[*n - 1];
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
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
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
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(a+1,64);
      a[1:i__2] =(b[1:i__2] + c__[1:i__2] * d__[1:i__2])*
	  (b[1:i__2] + c__[1:i__2] * d__[1:i__2]);
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
s252_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real s, t;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     scalar and array expansion */
/*     loop with ambiguous scalar temporary */

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
	 &bb[bb_offset], &cc[cc_offset], "s252 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1]=b[1]*c__[1];
      a[2:i__2-1]= b[2:i__2-1]*c__[2:i__2-1]+b[1:i__2-1]*c__[1:i__2-1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s252 ", (ftnlen) 5);
  return 0;
}				/* s252_ */

/* %2.5 */
/* Subroutine */ int
s253_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
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


/*     scalar and array expansion */
/*     scalar expansion, assigned under if */

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
	 &bb[bb_offset], &cc[cc_offset], "s253 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
#ifdef __linux__
  int tmp[i__2 = *n] __attribute__ ((aligned(64)));
#else
  __declspec(align(64)) int tmp[i__2 = *n];
#endif
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(a+1,64);
      tmp[0:i__2]= a[1:i__2] > b[1:i__2];
      c__[1:i__2] += tmp[0:i__2]? a[1:i__2] -= 
	(tmp[0:i__2-1]? b[1:i__2]:0.f) * d__[1:i__2]:0.f;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &c__[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s253 ", (ftnlen) 5);
  return 0;
}				/* s253_ */

/* %2.5 */
/* Subroutine */ int
s254_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *  a, real *  b, real * c__,
//       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real x;
  real t1, t2;
  integer nl;
  real chksum;


/*     scalar and array expansion */
/*     carry around variable */

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
	 &bb[bb_offset], &cc[cc_offset], "s254 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      x = b[*n];
      i__2 = *n;
#if 1
	  i__=1; a[i__] = (b[i__] + x) * .5f;
      a[2:i__2-1] = (b[2:i__2-1]+b[1:i__2-1])*.5f;
#else
      a[1:i__2] = (b[1:i__2]+__sec_rotate(b[1:i__2],1))*.5f;
#endif
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s254 ", (ftnlen) 5);
  return 0;
}				/* s254_ */

/* %2.5 */
/* Subroutine */ int
s255_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real x, y;
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
#if 1
#pragma simd firstprivate(x,y)
      for (i__ = 1; i__ <= i__2; ++i__) {
	  a[i__] = (b[i__] + x + y) * .333f;
	  y = x;
	  x = b[i__];
	}
#else
      a[1:i__2]=(b[1:i__2]+__sec_rotate(b[1:i__2],1)+
	__sec_rotate(b[1:i__2],2))*.333f;
#endif
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
/* Subroutine */ int
s256_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * b, real * c__,
       real * d__, real * e, real * restrict aa, real * restrict bb,
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
	 &bb[bb_offset], &cc[cc_offset], "s256 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
      for (j = 2; j <= i__3; ++j){
	  aa[1 + j * aa_dim1:i__2] += - a[j-1] + bb[1 + j * bb_dim1:i__2];
	  a[j]= aa[i__2 + j * aa_dim1] - bb[i__2 + j * bb_dim1];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs1d_ (n, &a[1]) + cs2d_ (n, &aa[aa_offset]);
  i__1 = *ntimes / *n * *n * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s256 ", (ftnlen) 5);
  return 0;
}				/* s256_ */

/* %2.5 */
/* Subroutine */ int
s257_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * b, real * c__,
       real * d__, real * e, real * restrict aa, real * restrict bb,
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
      cilk_for (int j = 1; j < i__3; ++j){
	  float tmp= a[1];
	  for (int i__ = 2; i__ <= i__2; ++i__) {
	      tmp = aa[i__ + j * aa_dim1] - tmp;
	      aa[i__ + j * aa_dim1] = tmp + bb[i__ + j * bb_dim1];
	    }
	}
      for (i__ = 2; i__ <= i__2; ++i__) {
	  a[i__] = aa[i__ + i__3 * aa_dim1] - a[i__ - 1];
	  aa[i__ + i__3 * aa_dim1] = a[i__] + bb[i__ + i__3 * bb_dim1];
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
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e,
       real * restrict aa, real * bb, real * cc)
{
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real s[*n],ss;
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
      i__2 = *n;
      ss = 0;
      for (i__ = 1; i__ <= i__2; ++i__)
	  s[i__-1] = ss = (a[i__] > 0.f)? d__[i__] * d__[i__]: ss;
#pragma nofusion
      b[1:i__2] = s[0:i__2] * c__[1:i__2] + d__[1:i__2];
      e[1:i__2] = s[0:i__2] * aa[1 + aa_dim1:i__2 + aa_dim1] 
	+ aa[1 + aa_dim1:i__2 + aa_dim1];
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

/* %2.6 */
/* Subroutine */ int
s261_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     scalar renaming */

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
	 &bb[bb_offset], &cc[cc_offset], "s261 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      // optimized for Core i7
//      a[2:i__2-1] += b[2:i__2-1] + c__[1:i__2-1];
//      #pragma nofusion
//      c__[2:i__2-1]*= d__[2:i__2-1];
      a[2] += b[2] + c__[1];
      a[3:i__2-2] += b[3:i__2-2] + (c__[2:i__2-2]*= d__[2:i__2-2]);
      c__[i__2] *= d__[i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &c__[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s261 ", (ftnlen) 5);
  return 0;
}				/* s261_ */

/* %2.7 */
/* Subroutine */ int
s271_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(a+1,64);
      a[1:i__2] += (b[1:i__2] > 0.f ? b[1:i__2] : 0.f) * c__[1:i__2];
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
s272_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
       real * bb, real * cc, real * restrict t) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control flow */
/*     loop with independent conditional */

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
	 &bb[bb_offset], &cc[cc_offset], "s272 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma vector aligned
      a[1:i__2] += (e[1:i__2] >= *t? c__[1:i__2]: 0.f) * d__[1:i__2];
      b[1:i__2] += (e[1:i__2] >= *t? c__[1:i__2]: 0.f) * c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s272 ", (ftnlen) 5);
  return 0;
}				/* s272_ */

/* %2.7 */
/* Subroutine */ int
s273_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
  __assume_aligned(e+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(a+1,64);
      c__[1:i__2] += (a[1:i__2] += d__[1:i__2] * e[1:i__2]) * d__[1:i__2];
      b[1:i__2] += (a[1:i__2] < 0.f)? d__[1:i__2] * e[1:i__2]: 0.f;
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
/* Subroutine */ int
s274_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
/*     complex loop with dependent conditional */

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
	 &bb[bb_offset], &cc[cc_offset], "s274 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
  __assume_aligned(e+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(b+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(a+1,64);
      a[1:i__2] = c__[1:i__2] + e[1:i__2] * d__[1:i__2];
      b[1:i__2] += a[1:i__2] > 0.f ? a[1:i__2] : 0.f;
      a[1:i__2] = a[1:i__2] > 0.f ? a[1:i__2] : e[1:i__2] * d__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s274 ", (ftnlen) 5);
  return 0;
}				/* s274_ */

/* %2.7 */
/* Subroutine */ int
s275_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * restrict cc) {
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
	for (int j = 2; j <= i__3; ++j)
	    aa[2 + j * aa_dim1:i__2-1] = (aa[2 + aa_dim1:i__2-1] > 0.f)? 
		aa[2 + (j - 1) * aa_dim1:i__2-1] + bb[ 2 + j * bb_dim1:i__2-1] *
		cc[2 + j * cc_dim1:i__2-1]: aa[2 + j * aa_dim1:i__2-1];
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

/* %2.7 */
/* Subroutine */ int
s276_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  integer mid;
  extern real cs1d_ (integer *, real *);


/*     control flow */
/*     if test using loop index */

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
	 &bb[bb_offset], &cc[cc_offset], "s276 ", (ftnlen) 5);
  forttime_ (&t1);
  mid = *n / 2;
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
	a[1:mid-1] += b[1:mid-1] * c__[1:mid-1];
	a[mid:i__2-mid+1] += b[mid:i__2-mid+1] * d__[mid:i__2-mid+1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s276 ", (ftnlen) 5);
  return 0;
}				/* s276_ */

/* %2.7 */
/* Subroutine */ int
s277_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
/*     test for dependences arising from guard variable computation. */

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
	 &bb[bb_offset], &cc[cc_offset], "s277 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
      for (i__ = 1; i__ <= i__2; ++i__)
	  if (a[i__] < 0.f) {
	      b[i__ + 1] = c__[i__] + d__[i__] * e[i__];
	      if (b[i__] < 0.f)
		  a[i__] += c__[i__] * d__[i__];
	    }
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s277 ", (ftnlen) 5);
  return 0;
}				/* s277_ */

/* %2.7 */
/* Subroutine */ int
s278_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
/*     if/goto to block if-then-else */

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
	 &bb[bb_offset], &cc[cc_offset], "s278 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma vector aligned
      for (i__ = 1; i__ <= i__2; ++i__)
	  if (a[i__] <= 0.f)
	      b[i__] = -b[i__] + d__[i__] * e[i__];
	  else
	      c__[i__] = -c__[i__] + d__[i__] * e[i__];
      a[1:i__2] = b[1:i__2] + c__[1:i__2] * d__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]) + cs1d_ (n, &c__[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s278 ", (ftnlen) 5);
  return 0;
}				/* s278_ */

/* %2.7 */
/* Subroutine */ int
s279_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
/*     vector if/gotos */

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
	 &bb[bb_offset], &cc[cc_offset], "s279 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
#pragma unroll(4)
#pragma vector aligned
      for (i__ = 1; i__ <= i__2; ++i__) 
	  if (a[i__] > 0.f)
	      c__[i__] = -c__[i__] + e[i__] * e[i__];
	  else
	      if ((b[i__] = -b[i__] + d__[i__] * d__[i__])> a[i__])
		  c__[i__] += d__[i__] * e[i__];
      a[1:i__2] = b[1:i__2] + c__[1:i__2] * d__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]) + cs1d_ (n, &c__[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s279 ", (ftnlen) 5);
  return 0;
}				/* s279_ */

/* %2.7 */
/* Subroutine */ int
s2710_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * restrict d__, real * restrict e,
	real * aa, real * bb, real * cc, real * restrict x) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control flow */
/*     scalar and vector ifs */

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
	 &bb[bb_offset], &cc[cc_offset], "s2710", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      int ng10 = *n > 10;
      int xg0 = *x > 0;
      i__2 = *n;
#pragma vector aligned
#if defined __MIC__
#pragma simd
#endif
      for (i__ = 1; i__ <= i__2; ++i__)
	  if (a[i__] > b[i__]) {
	      a[i__] += b[i__] * d__[i__];
	      c__[i__] = ng10 ? c__[i__]+ d__[i__] * d__[i__]:
		  d__[i__] * e[i__] + 1.f;
	    }
	  else {
	      b[i__] = a[i__] + e[i__] * e[i__];
	      c__[i__] = xg0 ? a[i__] + d__[i__] * d__[i__]:
		  c__[i__] + e[i__] * e[i__];
	    }
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]) + cs1d_ (n, &b[1]) + cs1d_ (n, &c__[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s2710", (ftnlen) 5);
  return 0;
}				/* s2710_ */

/* %2.7 */
/* Subroutine */ int
s2711_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
/*     semantic if removal */

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
	 &bb[bb_offset], &cc[cc_offset], "s2711", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
  __assume_aligned(c__+1,64);
  __assume_aligned(b+1,64);
  __assume_aligned(a+1,64);
      i__2 = *n;
      a[1:i__2] += (b[1:i__2] != 0.f?b[1:i__2]:0.f) * c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s2711", (ftnlen) 5);
  return 0;
}				/* s2711_ */

/* %2.7 */
/* Subroutine */ int
s2712_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
/*     if to elemental min */

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
	 &bb[bb_offset], &cc[cc_offset], "s2712", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
  __assume_aligned(c__+1,64);
  __assume_aligned(b+1,64);
  __assume_aligned(a+1,64);
      a[1:i__2] += (a[1:i__2] > b[1:i__2]?b[1:i__2]:0.f) * c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s2712", (ftnlen) 5);
  return 0;
}				/* s2712_ */

/* %2.8 */
/* Subroutine */ int
s281_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real x;
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
	a[1:(i__2+1)/2] = ( b[1:(i__2+1)/2] = a[i__2:(i__2+1)/2:-1] +
	    b[1:(i__2+1)/2] * c__[1:(i__2+1)/2] ) - 1.f;
	a[(i__2+3)/2:i__2/2] = ( b[(i__2+3)/2:i__2/2] = a[i__2/2:i__2/2:-1] +
	    b[(i__2+3)/2:i__2/2] * c__[(i__2+3)/2:i__2/2] ) - 1.f;
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

/* %2.9 */
/* Subroutine */ int
s291_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  integer im1;


/*     loop peeling */
/*     wrap around variable, 1 level */

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
	 &bb[bb_offset], &cc[cc_offset], "s291 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      im1 = *n;
      i__2 = *n;
      i__ = 1;
	  a[i__] = (b[i__] + b[im1]) * .5f;
      a[2:i__2] = (b[2:i__2]+b[1:i__2 -1])*.5f;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s291 ", (ftnlen) 5);
  return 0;
}				/* s291_ */

/* %2.9 */
/* Subroutine */ int
s292_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  real x, y;


/*     loop peeling */
/*     wrap around variable, 2 levels */

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
	 &bb[bb_offset], &cc[cc_offset], "s292 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      x = b[*n];
      y = b[*n - 1];
#pragma simd firstprivate(x,y)
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
  check_ (&chksum, &i__1, n, &t2, "s292 ", (ftnlen) 5);
  return 0;
}				/* s292_ */

/* %2.9 */
/* Subroutine */ int
s293_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop peeling */
/*     a(i)=a(1) with actual dependence cycle, loop is vectorizable */

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
	 &bb[bb_offset], &cc[cc_offset], "s293 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      float tmp=a[1];
      i__2 = *n;
#pragma unroll(4)
      a[1:i__2] = tmp;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s293 ", (ftnlen) 5);
  return 0;
}				/* s293_ */

/* %2.10 */
/* Subroutine */ int
s2101_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * restrict aa, real * restrict bb, real * restrict cc) {
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
      aa[aa_dim1:i__2:aa_dim1] += bb[bb_dim1:i__2:bb_dim1] * cc[cc_dim1:i__2:cc_dim1];
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
	cilk_for (int j = 1; j <= i__2; ++j) 
	    aa[1 + j * aa_dim1:i__3] = (__sec_implicit_index(0)+1)==j;
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

/* %2.11 */
/* Subroutine */ int
s2111_ (integer * ntimes, integer * ld, integer * n, real *
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


/*     wavefronts */

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
	 &bb[bb_offset], &cc[cc_offset], "s2111", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = i__3 = *n;
      for (j = 2; j <= i__2; ++j) {
	  float tmp = aa[1 + j*aa_dim1];
	  for (i__ = 2; i__ <= i__3; ++i__)
	      aa[i__ + j * aa_dim1] = tmp += aa[i__ + (j - 1) * aa_dim1];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum = cs2d_ (n, &aa[aa_offset]);
  if (chksum == 0.f)
      chksum = 3.f;
  i__1 = *ntimes / *n * (*n - 1) * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s2111", (ftnlen) 5);
  return 0;
}				/* s2111_ */


/* ********************************************************** */
/*                                                         * */
/*                   IDIOM RECOGNITION                     * */
/*                                                         * */
/* ********************************************************** */
/* %3.1 */
/* Subroutine */ int
s311_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum, sum;


/*     reductions */
/*     sum reduction */

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
	 &bb[bb_offset], &cc[cc_offset], "s311 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      sum = __sec_reduce_add(a[1:i__2]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = sum;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s311 ", (ftnlen) 5);
  return 0;
}				/* s311_ */

/* %3.1 */
/* Subroutine */ int
s312_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  real prod;
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     reductions */
/*     product reduction */

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
	 &bb[bb_offset], &cc[cc_offset], "s312 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      prod = 1.f;
      i__2 = *n;
      prod = __sec_reduce_mul(a[1:i__2]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = prod;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s312 ", (ftnlen) 5);
  return 0;
}				/* s312_ */

/* %3.1 */
/* Subroutine */ int
s313_ (integer * ntimes, integer * ld, integer * n, real *
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


/*     reductions */
/*     dot product */

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
	 &bb[bb_offset], &cc[cc_offset], "s313 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      dot = __sec_reduce_add(a[1:i__2] * b[1:i__2]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = dot;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s313 ", (ftnlen) 5);
  return 0;
}				/* s313_ */

/* %3.1 */
/* Subroutine */ int
s314_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real x;
  real t1, t2;
  integer nl;
  real chksum;


/*     reductions */
/*     if to max reduction */

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
	 &bb[bb_offset], &cc[cc_offset], "s314 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      x = __sec_reduce_max(a[1:i__2]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = x;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s314 ", (ftnlen) 5);
  return 0;
}				/* s314_ */

/* %3.1 */
/* Subroutine */ int
s315_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real x;
  integer index;
  real t1, t2;
  integer nl;
  real chksum;


/*     reductions */
/*     if to max with index reduction, 1 dimension */

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
	 &bb[bb_offset], &cc[cc_offset], "s315 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      index = __sec_reduce_max_ind(a[1:i__2]);
      x = a[index];
      chksum = x + (real) index;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = x + (real) index;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s315 ", (ftnlen) 5);
  return 0;
}				/* s315_ */

/* %3.1 */
/* Subroutine */ int
s316_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real x;
  real t1, t2;
  integer nl;
  real chksum;


/*     reductions */
/*     if to min reduction */

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
	 &bb[bb_offset], &cc[cc_offset], "s316 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      x = __sec_reduce_min(a[1:i__2]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = x;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s316 ", (ftnlen) 5);
  return 0;
}				/* s316_ */

/* %3.1 */
/* Subroutine */ int
s317_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real q;
  real t1, t2;
  integer nl;
  real chksum;


/*     reductions */
/*     product reduction, vectorize with */
/*     1. scalar expansion of factor, and product reduction */
/*     2. closed form solution: q = factor**n */

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
	 &bb[bb_offset], &cc[cc_offset], "s317 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      q = 1.f;
      i__2 = *n;
#pragma simd reduction(*: q)
      for (i__ = 1; i__ <= i__2; ++i__)
	  q *= .99999f;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = q;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s317 ", (ftnlen) 5);
  return 0;
}				/* s317_ */

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
      i__2 = *n;
      index = __sec_reduce_max_ind(ABS(a[1:i__2:*inc]));
      max__ = ABS(a[index]);
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
s319_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum, sum;


/*     reductions */
/*     coupled reductions */

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
	 &bb[bb_offset], &cc[cc_offset], "s319 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
  __assume_aligned(d__+1,64);
  __assume_aligned(e+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(b+1,64);
  __assume_aligned(a+1,64);
      sum = __sec_reduce_add(
	   (a[1:i__2]=c__[1:i__2]+d__[1:i__2]) +
	    (b[1:i__2]=c__[1:i__2]+e[1:i__2]) );
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = sum;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s319 ", (ftnlen) 5);
  return 0;
}				/* s319_ */

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
#ifdef __MIC__
#pragma omp parallel for if(i__2 > 103)
#else
#pragma omp parallel for if(i__2 > 103) reduction(max: max__) lastprivate(xindex,yindex)
#endif
      for (j = 1; j <= i__2; ++j) {
	  int indxj;
	  float maxj;
	  indxj= __sec_reduce_max_ind(aa[1 + j * aa_dim1:i__3]);
	  maxj= aa[indxj];
#ifdef __MIC__
#pragma omp critical
#endif
	    if(maxj > max__) {
		max__= maxj;
		xindex=indxj - j * aa_dim1;
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

/* %3.1 */
/* Subroutine */ int
s3111_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum, sum;


/*     reductions */
/*     conditional sum reduction */

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
	 &bb[bb_offset], &cc[cc_offset], "s3111", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      sum = __sec_reduce_add(a[1:i__2] > 0.f ? a[1:i__2] : 0.f);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &sum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = sum;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s3111", (ftnlen) 5);
  return 0;
}				/* s3111_ */

/* %3.1 */
/* Subroutine */ int
s3112_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
	real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum, sum;


/*     reductions */
/*     sum reduction saving running sums */

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
	 &bb[bb_offset], &cc[cc_offset], "s3112", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      sum = 0.f;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	  b[i__] = sum += a[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &b[1]) + sum;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s3112", (ftnlen) 5);
  return 0;
}				/* s3112_ */

/* %3.1 */
/* Subroutine */ int
s3113_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;
  real r__1;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum, max__;


/*     reductions */
/*     maximum of absolute value */

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
	 &bb[bb_offset], &cc[cc_offset], "s3113", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      max__ = __sec_reduce_max(abs (a[1:i__2]));
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = max__;
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s3113", (ftnlen) 5);
  return 0;
}				/* s3113_ */

/* %3.2 */
/* Subroutine */ int
s321_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     recurrences */
/*     first order linear recurrence */

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
	 &bb[bb_offset], &cc[cc_offset], "s321 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 2; i__ <= i__2; ++i__)
	  a[i__] += a[i__ - 1] * b[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s321 ", (ftnlen) 5);
  return 0;
}				/* s321_ */

/* %3.2 */
/* Subroutine */ int
s322_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
      i__2 = *n;
#pragma unroll(4)
      for (i__ = 2; i__ <= i__2; ++i__) {
	  a[i__] = b[i__ - 1] + c__[i__] * d__[i__];
	  b[i__] = b[i__ - 1] + c__[i__] * (e[i__]+d__[i__]);
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

/* %3.3 */
/* Subroutine */ int
s331_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     search loops */
/*     if to last-1 */

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
	 &bb[bb_offset], &cc[cc_offset], "s331 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
	for (i__ = i__2; i__ >= 1; --i__)
	    if (a[i__] < 0.f)
		break;
	j = i__ >= 1? i__: -1;
//      j = *n - __sec_reduce_max_ind((a[*n:*n:-1] < 0.f) ? 1.f : 0.f) + 1;
      chksum = (real) j;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = (real) j;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s331 ", (ftnlen) 5);
  return 0;
}				/* s331_ */

/* %3.3 */
/* Subroutine */ int
s332_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc, real * t) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  integer index;
  real t1, t2;
  integer nl;
  real chksum;


/*     search loops */
/*     first value greater than threshold */

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
	 &bb[bb_offset], &cc[cc_offset], "s332 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
#if defined __AVX2__ || defined __MIC__
      int tmp[i__2 = *n];
#else
      float tmp[i__2 = *n];
#endif
      tmp[0:i__2]= a[1:i__2] > *t;
      index = __sec_reduce_max_ind(tmp[0:i__2]) + 1;
      if(index > i__2) index = -1;
      chksum = a[index] + (real) index;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &chksum);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = a[index] + (real) index;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s332 ", (ftnlen) 5);
  return 0;
}				/* s332_ */

/* %3.4 */
/* Subroutine */ int
s341_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     packing */
/*     pack positive values */

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
	 &bb[bb_offset], &cc[cc_offset], "s341 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      j = 0;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	  if (b[i__] > 0.f)
	      a[++j] = b[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s341 ", (ftnlen) 5);
  return 0;
}				/* s341_ */

/* %3.4 */
/* Subroutine */ int
s342_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     packing */
/*     unpacking */

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
	 &bb[bb_offset], &cc[cc_offset], "s342 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      j = 0;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	  if (a[i__] > 0.f)
	      a[i__] = b[++j];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s342 ", (ftnlen) 5);
  return 0;
}				/* s342_ */

/* %3.4 */
/* Subroutine */ int
s343_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__, j, k;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     packing */
/*     pack 2-d array into one dimension */

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
	 &bb[bb_offset], &cc[cc_offset], "s343 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
      k = 0;
      i__2 = i__3 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	  for (j = 1; j <= i__3; ++j)
	      if (bb[i__ + j * bb_dim1] > 0.f)
		  cdata_1.array[++k - 1] = aa[i__ + j * aa_dim1];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  i__1 = *n * *n;
  chksum = cs1d_ (&i__1, cdata_1.array);
  i__1 = *ntimes / *n * *n * *n;
  check_ (&chksum, &i__1, n, &t2, "s343 ", (ftnlen) 5);
  return 0;
}				/* s343_ */

/* %3.5 */
/* Subroutine */ int
s351_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real alpha;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     loop rerolling */
/*     unrolled saxpy */

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
	 &bb[bb_offset], &cc[cc_offset], "s351 ", (ftnlen) 5);
  forttime_ (&t1);
  alpha = c__[1];
  i__1 = *ntimes * 5;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; i__ += 5) {
	  a[i__] += alpha * b[i__];
	  a[i__ + 1] += alpha * b[i__ + 1];
	  a[i__ + 2] += alpha * b[i__ + 2];
	  a[i__ + 3] += alpha * b[i__ + 3];
	  a[i__ + 4] += alpha * b[i__ + 4];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes * 5);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * 5 * (*n / 5);
  check_ (&chksum, &i__1, n, &t2, "s351 ", (ftnlen) 5);
  return 0;
}				/* s351_ */

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

/* %3.5 */
/* Subroutine */ int
s353_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc, integer * ip) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real alpha;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     loop rerolling */
/*     unrolled sparse saxpy */

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
	 &bb[bb_offset], &cc[cc_offset], "s353 ", (ftnlen) 5);
  forttime_ (&t1);
  alpha = c__[1];
  i__1 = *ntimes * 5;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; i__ += 5) {
	  a[i__] += alpha * b[ip[i__]];
	  a[i__ + 1] += alpha * b[ip[i__ + 1]];
	  a[i__ + 2] += alpha * b[ip[i__ + 2]];
	  a[i__ + 3] += alpha * b[ip[i__ + 3]];
	  a[i__ + 4] += alpha * b[ip[i__ + 4]];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes * 5);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * 5 * (*n / 5);
  check_ (&chksum, &i__1, n, &t2, "s353 ", (ftnlen) 5);
  return 0;
}				/* s353_ */


/* ********************************************************** */
/*                                                         * */
/*                 LANGUAGE COMPLETENESS                   * */
/*                                                         * */
/* ********************************************************** */
/* %4.1 */
/* Subroutine */ int
s411_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
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


/*     loop recognition */
/*     if loop to do loop, zero trip */

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
	 &bb[bb_offset], &cc[cc_offset], "s411 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
  __assume_aligned(c__+1,64);
  __assume_aligned(b+1,64);
  __assume_aligned(a+1,64);
      i__2 = *n;
      a[1:i__2] += b[1:i__2] * c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s411 ", (ftnlen) 5);
  return 0;
}				/* s411_ */

/* %4.1 */
/* Subroutine */ int
s412_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real *restrict a, real *restrict b,
       real *restrict c__, real * d__,
       real * e, real * aa, real * bb, real * cc, integer * inc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2, i__3;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop recognition */
/*     if loop with variable increment */

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
	 &bb[bb_offset], &cc[cc_offset], "s412 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      i__3 = *inc;
      for (i__ = 0; (i__ += i__3) <= i__2;)
	  a[i__] += b[i__] * c__[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s412 ", (ftnlen) 5);
  return 0;
}				/* s412_ */

/* %4.1 */
/* Subroutine */ int
s413_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
      b[1:*n-1] += d__[1:*n-1] * e[1:*n-1];
      a[2:*n-1] = c__[2:*n-1] + d__[2:*n-1] * e[2:*n-1];
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

/* %4.1 */
/* Subroutine */ int
s414_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * restrict aa, real * restrict bb, real * restrict cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, j;
  real t1, t2;
  integer nl;
  real chksum;


/*     loop recognition */
/*     if loop to do loop, interchanging with do necessary */

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
	 &bb[bb_offset], &cc[cc_offset], "s414 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes / *n;
  for (nl = 1; nl <= i__1; ++nl) {
	i__2 = *n;
	for (j = 2; j <= i__2; ++j) {
	i__ = 1;
	while (i__ <= *n) {
	    aa[i__ + j * aa_dim1] = aa[i__ + (j - 1) * aa_dim1] + bb[i__ + (j 
		    - 1) * bb_dim1] * cc[i__ + (j - 1) * cc_dim1];
	    ++i__;
	    }
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes / *n);
  chksum =
    cs2d_ (n, &aa[aa_offset]) + cs2d_ (n, &bb[bb_offset]) + cs2d_ (n,
								   &cc
								   [cc_offset]);
  i__1 = *ntimes / *n * *n * (*n - 2);
  check_ (&chksum, &i__1, n, &t2, "s414 ", (ftnlen) 5);
  return 0;
}				/* s414_ */

/* %4.1 */
/* Subroutine */ int
s415_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
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


/*     loop recognition */
/*     while loop */

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
	 &bb[bb_offset], &cc[cc_offset], "s415 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
#if defined __AVX2__ || defined __MIC__
      int tmp[i__2 = *n];
#else
      float tmp[i__2 = *n];
#endif
      tmp[0:i__2]= a[1:i__2] < 0.f;
      i__ = __sec_reduce_max_ind(tmp[0:i__2]);
      a[1:i__] += b[1:i__] * c__[1:i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s415 ", (ftnlen) 5);
  return 0;
}				/* s415_ */

#include <stdio.h>
/* %4.2 */
/* Subroutine */ int
s421_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;
  static real equiv_0[1000];

  /* Local variables */
  integer i__;
#define x (equiv_0)
#define y (equiv_0)
  real t1, t2;
  integer nl;
  real chksum;


/*     storage classes and equivalencing */
/*     incorrect comment from original author:
**     equivalence- no overlap */

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
  set1d_ (n, y, &c_b3, &c__1);
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s421 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
      x[0:i__2] = y[1:i__2] + a[1:i__2];
      dummy_ (ld, n, x, &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, x);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s421 ", (ftnlen) 5);
  return 0;
}				/* s421_ */

#undef y
#undef x


/* %4.2 */
/* Subroutine */ int
s422_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * b, real * c__,
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
      x[0:i__2] = cdata_1.array[8:i__2] + a[1:i__2];
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
       ctime, real * dtime, real * restrict a, real * b, real * c__,
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
//#pragma omp simd safelen(16)
      cdata_1.array[1:i__2] = x[0:i__2] + a[1:i__2];
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
       ctime, real * dtime, real * restrict a, real * b, real * c__,
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
      x[1:i__2] = cdata_1.array[0:i__2-1] + a[1:i__2];
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


/* %4.3 */
/* Subroutine */ int
s431_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     parameters */
/*     parameter statement */

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
	 &bb[bb_offset], &cc[cc_offset], "s431 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] += b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s431 ", (ftnlen) 5);
  return 0;
}				/* s431_ */

/* %4.3 */
/* Subroutine */ int
s432_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* Initialized data */

  static integer k1 = 1;
  static integer k2 = 2;

  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, k;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     parameters */
/*     data statement */

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
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[4 * aa_offset / 4],
	 &bb[4 * bb_offset / 4], &cc[4 * cc_offset / 4], "s432 ", (ftnlen) 5);
  k = (k1 << 1) - k2;
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] = a[1+k:i__2] + b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s432 ", (ftnlen) 5);
  return 0;
}				/* s432_ */

/* %4.4 */
/* Subroutine */ int
s441_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
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
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(a+1,64);
	a[1:i__2] +=(d__[1:i__2] <= 0.f?b[1:i__2]:c__[1:i__2]) *
		(d__[1:i__2]==0.f?b[1:i__2]:c__[1:i__2]);
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
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e, real * aa,
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
      for (int i__ = 1; i__ <= i__2; ++i__)
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

/* %4.4 */
/* Subroutine */ int
s443_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
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
	 &bb[bb_offset], &cc[cc_offset], "s443 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
  __assume_aligned(b+1,64);
  __assume_aligned(c__+1,64);
  __assume_aligned(d__+1,64);
  __assume_aligned(a+1,64);
      a[1:i__2] += b[1:i__2] * (d__[1:i__2] <= 0.f?c__[1:i__2]:b[1:i__2]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s443 ", (ftnlen) 5);
  return 0;
}				/* s443_ */

/* %4.5 */
#include <math.h>
/* Subroutine */ int
s451_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
      i__2 = *n;
//      cilk_for (int i__ = 1; i__ <= i__2; ++i__)
	  a[1:i__2] = sinf (b[1:i__2]) + cosf (c__[1:i__2]);
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

/* %4.5 */
/* Subroutine */ int
s452_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
/*     seq function */

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
	 &bb[bb_offset], &cc[cc_offset], "s452 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] = b[1:i__2] + c__[1:i__2] * (__sec_implicit_index(0)+1);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s452 ", (ftnlen) 5);
  return 0;
}				/* s452_ */

/* %4.5 */
/* Subroutine */ int
s453_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
       real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real s;
  real t1, t2;
  integer nl;
  real chksum;


/*     intrinsic functions */
/*     seq function */

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
	 &bb[bb_offset], &cc[cc_offset], "s453 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
	a[1:i__2] = (__sec_implicit_index(0)+1)*2 * b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s453 ", (ftnlen) 5);
  return 0;
}				/* s453_ */

/* %4.7 */
/* Subroutine */ int
s471_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * a, real * b, real * c__, real * d__,
       real * e, real * aa, real * bb, real * cc, real * c471s) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, m;
  real x[1000];
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);
  extern /* Subroutine */ int s471s_ (void);


/*     call statements */

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
  m = *n;
  set1d_ (n, x, &c_b393, &c__1);
  init_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	 &bb[bb_offset], &cc[cc_offset], "s471 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = m;
      cilk_for (int i__ = 1; i__ <= i__2; ++i__) {
	  x[i__ - 1] = b[i__] + d__[i__] * d__[i__];
	  s471s_ ();
	  b[i__] = c__[i__] + d__[i__] * e[i__];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes) - *n * *ntimes * *c471s;
  chksum = cs1d_ (n, x) + cs1d_ (n, &b[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s471 ", (ftnlen) 5);
  return 0;
}				/* s471_ */

/* %4.8 */
/* Subroutine */ int
s481_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
       real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     non-local goto's */
/*     stop statement */

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
	 &bb[bb_offset], &cc[cc_offset], "s481 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
	  if (d__[i__] < 0.f)
	      exit (EXIT_FAILURE);
	  a[i__] += b[i__] * c__[i__];
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s481 ", (ftnlen) 5);
  return 0;
}				/* s481_ */

/* %4.8 */
/* Subroutine */ int
s482_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * d__, real * e, real * aa, real * bb,
       real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     non-local goto's */
/*     other loop exit with code before exit */

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
	 &bb[bb_offset], &cc[cc_offset], "s482 ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
	  a[i__] += b[i__] * c__[i__];
	  if (c__[i__] > b[i__])
	      return 0;
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s482 ", (ftnlen) 5);
  return 0;
}				/* s482_ */

/* %4.9 */
/* Subroutine */ int
s491_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * e, real * aa,
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
      a[ip[1:i__2]] = b[1:i__2] + c__[1:i__2] * d__[1:i__2];
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

/* %4.11 */
/* Subroutine */ int
s4112_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
	real * d__, real * e, real * aa, real * bb, real * cc, integer * ip,
	real * restrict s) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     indirect addressing */
/*     sparse saxpy */

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
	 &bb[bb_offset], &cc[cc_offset], "s4112", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] += b[ip[1:i__2]] * *s;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s4112", (ftnlen) 5);
  return 0;
}				/* s4112_ */

/* %4.11 */
/* Subroutine */ int
s4113_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * d__, real * e, real * aa, real * bb,
	real * cc, integer * ip) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  extern real cs1d_ (integer *, real *);


/*     indirect addressing */
/*     indirect addressing on rhs and lhs */

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
	 &bb[bb_offset], &cc[cc_offset], "s4113", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[ip[1:i__2]] = b[ip[1:i__2]] + c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s4113", (ftnlen) 5);
  return 0;
}				/* s4113_ */

/* %4.11 */
/* Subroutine */ int
s4114_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * restrict d__, real * e, real * aa,
	real * bb, real * cc, integer * ip, integer * n1) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__, k;
  real t1, t2;
  integer nl;
  real chksum;


/*     indirect addressing */
/*     mix indirect addressing with variable lower and upper bounds */

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
	 &bb[bb_offset], &cc[cc_offset], "s4114", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] = b[1:i__2] + c__[i__2 - ip[1:i__2] + 1] * d__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s4114", (ftnlen) 5);
  return 0;
}				/* s4114_ */

/* %4.11 */
/* Subroutine */ int
s4115_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * aa, real * bb, real * cc, integer * ip) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum, sum;


/*     indirect addressing */
/*     sparse dot product */

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
	 &bb[bb_offset], &cc[cc_offset], "s4115", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      sum = __sec_reduce_add(a[1:i__2] * b[ip[1:i__2]]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = sum;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s4115", (ftnlen) 5);
  return 0;
}				/* s4115_ */

/* %4.11 */
/* Subroutine */ int
s4116_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * aa, real * bb, real * cc, integer * ip, integer * j,
	integer * inc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;
  integer off;
  real sum;


/*     indirect addressing */
/*     more complicated sparse sdot */

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
	 &bb[bb_offset], &cc[cc_offset], "s4116", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n - 1;
      sum = __sec_reduce_add(a[*inc + 1:i__2] * aa[ip[1:i__2] + *j * aa_dim1]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = sum;
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s4116", (ftnlen) 5);
  return 0;
}				/* s4116_ */

/* %4.11 */
/* Subroutine */ int
s4117_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * restrict d__, real * e, real * aa,
	real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     indirect addressing */
/*     seq function */

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
	 &bb[bb_offset], &cc[cc_offset], "s4117", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 2; i__ <= i__2; ++i__)
	  a[i__] = b[i__] + c__[i__ / 2] * d__[i__];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * (*n - 1);
  check_ (&chksum, &i__1, n, &t2, "s4117", (ftnlen) 5);
  return 0;
}				/* s4117_ */

/* %4.12 */
/* Subroutine */ int
s4121_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * d__, real * e, real * aa, real * bb,
	real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     statement functions */
/*     elementwise multiplication */

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
	 &bb[bb_offset], &cc[cc_offset], "s4121", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] += b[1:i__2] * c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "s4121", (ftnlen) 5);
  return 0;
}				/* s4121_ */

#include <string.h>
/* %5.1 */
/* Subroutine */ int
va_ (integer * ntimes, integer * ld, integer * n, real *
     ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
     real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     vector assignment */

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
	 &bb[bb_offset], &cc[cc_offset], "va   ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
//      memcpy(&a[1],&b[1],i__2*sizeof(a[1]));
	  a[1:*n] = b[1:*n];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "va   ", (ftnlen) 5);
  return 0;
}				/* va_ */

/* %5.1 */
/* Subroutine */ int
vag_ (integer * ntimes, integer * ld, integer * n, real *
      ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
      real * d__, real * e, real * aa, real * bb, real * cc, integer * ip) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     vector assignment, gather */

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
	 &bb[bb_offset], &cc[cc_offset], "vag  ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] = b[ip[1:i__2]];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vag  ", (ftnlen) 5);
  return 0;
}				/* vag_ */

/* %5.1 */
/* Subroutine */ int
vas_ (integer * ntimes, integer * ld, integer * n, real *
      ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
      real * d__, real * e, real * aa, real * bb, real * cc, integer * ip) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     vector assignment, scatter */

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
	 &bb[bb_offset], &cc[cc_offset], "vas  ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[ip[1:i__2]] = b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vas  ", (ftnlen) 5);
  return 0;
}				/* vas_ */

/* %5.1 */
/* Subroutine */ int
vif_ (integer * ntimes, integer * ld, integer * n, real *
      ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
      real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     vector if */

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
	 &bb[bb_offset], &cc[cc_offset], "vif  ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] = (b[1:i__2] > 0.f)? b[1:i__2]: a[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vif  ", (ftnlen) 5);
  return 0;
}				/* vif_ */

/* %5.1 */
/* Subroutine */ int
vpv_ (integer * ntimes, integer * ld, integer * n, real *
      ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
      real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     vector plus vector */

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
	 &bb[bb_offset], &cc[cc_offset], "vpv  ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] += b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vpv  ", (ftnlen) 5);
  return 0;
}				/* vpv_ */

/* %5.1 */
/* Subroutine */ int
vtv_ (integer * ntimes, integer * ld, integer * n, real *
      ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
      real * d__, real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     vector times vector */

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
	 &bb[bb_offset], &cc[cc_offset], "vtv  ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] *= b[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vtv  ", (ftnlen) 5);
  return 0;
}				/* vtv_ */

/* %5.1 */
/* Subroutine */ int
vpvtv_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
/*     vector plus vector times vector */

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
	 &bb[bb_offset], &cc[cc_offset], "vpvtv", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] += b[1:i__2] * c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vpvtv", (ftnlen) 5);
  return 0;
}				/* vpvtv_ */

/* %5.1 */
/* Subroutine */ int
vpvts_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b, real * c__,
	real * d__, real * e, real * aa, real * bb, real * cc,
	real * restrict s) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     vector plus vector times scalar */

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
	 &bb[bb_offset], &cc[cc_offset], "vpvts", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] += b[1:i__2] * *s;
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vpvts", (ftnlen) 5);
  return 0;
}				/* vpvts_ */

/* %5.1 */
/* Subroutine */ int
vpvpv_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
      a[1:i__2] += b[1:i__2] + c__[1:i__2];
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

/* %5.1 */
/* Subroutine */ int
vtvtv_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * restrict a, real * restrict b,
	real * restrict c__, real * d__, real * e, real * aa, real * bb,
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
/*     vector times vector times vector */

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
	 &bb[bb_offset], &cc[cc_offset], "vtvtv", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      a[1:i__2] *= b[1:i__2] * c__[1:i__2];
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, &a[1]);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vtvtv", (ftnlen) 5);
  return 0;
}				/* vtvtv_ */

/* %5.1 */
/* Subroutine */ int
vsumr_ (integer * ntimes, integer * ld, integer * n, real *
	ctime, real * dtime, real * a, real * b, real * c__, real * d__,
	real * e, real * aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real t1, t2;
  integer nl;
  real chksum, sum;


/*     control loops */
/*     vector sum reduction */

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
	 &bb[bb_offset], &cc[cc_offset], "vsumr", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      sum = __sec_reduce_add(a[1:i__2]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = sum;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vsumr", (ftnlen) 5);
  return 0;
}				/* vsumr_ */

/* %5.1 */
/* Subroutine */ int
vdotr_ (integer * ntimes, integer * ld, integer * n, real *
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


/*     control loops */
/*     vector dot product reduction */

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
	 &bb[bb_offset], &cc[cc_offset], "vdotr", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      dot = __sec_reduce_add(a[1:i__2] * b[1:i__2]);
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = dot;
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vdotr", (ftnlen) 5);
  return 0;
}				/* vdotr_ */

/* %5.1 */
/* Subroutine */ int
vbor_ (integer * ntimes, integer * ld, integer * n, real *
       ctime, real * dtime, real * restrict a, real * restrict b,
       real * restrict c__, real * restrict d__, real * restrict e,
       real * restrict aa, real * bb, real * cc) {
  /* System generated locals */
  integer aa_dim1, aa_offset, bb_dim1, bb_offset, cc_dim1, cc_offset, i__1,
    i__2;

  /* Local variables */
  integer i__;
  real s[1000], a1, b1, c1, d1, e1, f1;
  real t1, t2;
  integer nl;
  real chksum;


/*     control loops */
/*     basic operations rates, isolate arithmetic from memory traffic */
/*     all combinations of three, 59 flops for 6 loads and 1 store. */

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
	 &bb[bb_offset], &cc[cc_offset], "vbor ", (ftnlen) 5);
  forttime_ (&t1);
  i__1 = *ntimes;
  for (nl = 1; nl <= i__1; ++nl) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
	  a1 = a[i__];
	  b1 = b[i__];
	  c1 = c__[i__];
	  d1 = d__[i__];
	  e1 = e[i__];
	  f1 = aa[i__ + aa_dim1];
	  a1 *= ((d1 *(e1 + f1)+ e1 * f1) + c1 *((e1 + f1) + d1))
	      + b1 *(((e1 + f1) + d1) + c1);
	  b1 *= (d1 *(e1 + f1)+ e1 * f1) + c1 *((e1 + f1) + d1);
	  c1 *= d1 *(e1 + f1)+ e1 * f1;
	  d1 *= e1 * f1;
	  s[i__ - 1] = a1 * b1 * c1 * d1;
	}
      dummy_ (ld, n, &a[1], &b[1], &c__[1], &d__[1], &e[1], &aa[aa_offset],
	      &bb[bb_offset], &cc[cc_offset], &c_b3);
    }
  forttime_ (&t2);
  t2 = t2 - t1 - *ctime - *dtime * (real) (*ntimes);
  chksum = cs1d_ (n, s);
  i__1 = *ntimes * *n;
  check_ (&chksum, &i__1, n, &t2, "vbor ", (ftnlen) 5);
  return 0;
}				/* vbor_ */
