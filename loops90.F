C***********************************************************************
C                TEST SUITE FOR VECTORIZING COMPILERS                  *
C                        (File 2 of 2)                                 *
C                                                                      *
C  Version:   2.0                                                      *
C  Date:      3/14/88                                                  *
C  Authors:   Original loops from a variety of                         *
C             sources. Collection and synthesis by                     *
C                                                                      *
C             David Callahan  -  Tera Computer                         *
C             Jack Dongarra   -  University of Tennessee               *
C             David Levine    -  Argonne National Laboratory           *
C***********************************************************************
C  Version:   3.0                                                      *
C  Date:      1/4/91                                                   *
C  Authors:   David Levine    -  Executable version                    *
C***********************************************************************
C                         ---DESCRIPTION---                            *
C                                                                      *
C  This test consists of a variety of  loops that represent different  *
C  constructs  intended   to  test  the   analysis  capabilities of a  *
C  vectorizing  compiler.  Each loop is  executed with vector lengths  *
C  of 10, 100, and 1000.   Also included  are several simple  control  *
C  loops  intended  to  provide  a  baseline  measure  for  comparing  *
C  compiler performance on the more complicated loops.                 *
C                                                                      *
C  The  output from a  run  of the test  consists of seven columns of  *
C  data:                                                               *
C     Loop:        The name of the loop.                               *
C     VL:          The vector length the loop was run at.              *
C     Seconds:     The time in seconds to run the loop.                *
C     Checksum:    The checksum calculated when running the test.      *
C     PreComputed: The precomputed checksum (64-bit IEEE arithmetic).  *
C     Residual:    A measure of the accuracy of the calculated         *
C                  checksum versus the precomputed checksum.           *
C     No.:         The number of the loop in the test suite.           *
C                                                                      *
C  The  residual  calculation  is  intended  as  a  check  that   the  *
C  computation  was  done  correctly  and  that 64-bit arithmetic was  *
C  used.   Small   residuals    from    non-IEEE    arithmetic    and  *
C  nonassociativity  of   some calculations  are   acceptable.  Large  *
C  residuals  from   incorrect  computations or  the  use   of 32-bit  *
C  arithmetic are not acceptable.                                      *
C                                                                      *
C  The test  output  itself  does not report   any  results;  it just  *
C  contains data.  Absolute  measures  such as Mflops and  total time  *
C  used  are  not   appropriate    metrics  for  this  test.   Proper  *
C  interpretation of the results involves correlating the output from  *
C  scalar and vector runs  and the  loops which  have been vectorized  *
C  with the speedup achieved at different vector lengths.              *
C                                                                      *
C  These loops  are intended only  as  a partial test of the analysis  *
C  capabilities of a vectorizing compiler (and, by necessity,  a test  *
C  of the speed and  features  of the underlying   vector  hardware).  *
C  These loops  are  by no means  a  complete  test  of a vectorizing  *
C  compiler and should not be interpreted as such.                     *
C                                                                      *
C***********************************************************************
C                           ---DIRECTIONS---                           *
C                                                                      *
C  To  run this  test,  you will  need  to  supply  a  function named  *
C  second() that returns user CPU time.                                *
C                                                                      *
C  This test is distributed as two separate files, one containing the  *
C  driver  and  one containing the loops.   These  two files MUST  be  *
C  compiled separately.                                                *
C                                                                      *
C  Results must  be supplied from  both scalar and vector  runs using  *
C  the following rules for compilation:                                *
C                                                                      *
C     Compilation   of the  driver  file must  not  use any  compiler  *
C     optimizations (e.g., vectorization, function  inlining,  global  *
C     optimizations,...).   This   file   also must  not  be analyzed  *
C     interprocedurally to  gather information useful  in  optimizing  *
C     the test loops.                                                  *
C                                                                      *
C     The file containing the  loops must be compiled twice--once for  *
C     a scalar run and once for a vector run.                          *
C                                                                      *
C        For the scalar  run, global (scalar) optimizations should be  *
C        used.                                                         *
C                                                                      *
C        For  the  vector run,  in  addition   to  the  same   global  *
C        optimizations specified  in the scalar   run,  vectorization  *
C        and--if available--automatic  call generation to   optimized  *
C        library  routines,  function inlining,  and  interprocedural  *
C        analysis should be  used.  Note again that function inlining  *
C        and interprocedural  analysis  must  not be  used to  gather  *
C        information  about any of the  program  units  in the driver  *
C        program.                                                      *
C                                                                      *
C     No changes  may  be made  to   the   source code.   No compiler  *
C     directives may be used, nor may  a file  be split into separate  *
C     program units.  (The exception is  filling  in  the information  *
C     requested in subroutine "info" as described below.)              *
C                                                                      *
C     All files must be compiled to use 64-bit arithmetic.             *
C                                                                      *
C     The  outer  timing  loop  is  included  only   to increase  the  *
C     granularity of the calculation.  It should not be vectorizable.  *
C     If it is found to be so, please notify the authors.              *
C                                                                      *
C  All runs  must be  made  on a standalone  system  to minimize  any  *
C  external effects.                                                   *
C                                                                      *
C  On virtual memory computers,  runs should be  made with a physical  *
C  memory and working-set  size  large enough  that  any  performance  *
C  degradation from page  faults is negligible.   Also,  the  timings  *
C  should be repeatable  and you  must  ensure  that timing anomalies  *
C  resulting from paging effects are not present.                      *
C                                                                      *
C  You should edit subroutine "info"   (the  last subroutine  in  the  *
C  driver program) with  information specific to your  runs, so  that  *
C  the test output will be annotated automatically.                    *
C                                                                      *
C  Please return the following three files in an electronic format:    *
C                                                                      *
C  1. Test output from a scalar run.                                   *
C  2. Test output from a vector run.                                   *
C  3. Compiler output listing (source echo, diagnostics, and messages) *
C     showing which loops have been vectorized.                        *
C                                                                      *
C  The preferred media  for receipt, in order  of preference, are (1)  *
C  electronic mail, (2) 9-track  magnetic or  cartridge tape in  Unix  *
C  tar  format, (3) 5 IBM PC/DOS  floppy   diskette, or  (4) 9-track  *
C  magnetic  tape in  ascii  format,   80 characters per card,  fixed  *
C  records, 40 records per block, 1600bpi.  Please return to           *
C                                                                      *
C  David Levine       		                                       *
C  Mathematics and Computer Science Division                           *
C  Argonne National Laboratory                                         *
C  Argonne, Illinois 60439                                             *
C  levine@mcs.anl.gov                                                  *
C***********************************************************************
C%1.1
      subroutine s111(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     linear dependence testing
C     no dependence - vectorizable
C
      integer ntimes,ld,n,i,nl
      real a(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
#if __INTEL_COMPILER
      real :: b(n)
#else
      real,volatile :: b(n)     ! stop vectorization
#endif
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s111 ')
      call forttime(t1)
      do nl= 1,2*ntimes
	a(2:n:2)= a(1:n-1:2)+b(2:n:2)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,2*ntimes*(n/2),n,t2,'s111 ')
      return
      end
C%1.1
      subroutine s112(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     linear dependence testing
C     loop reversal
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s112 ')
      call forttime(t1)
      do nl= 1,ntimes
	a(n:2:-1)= a(n-1:1:-1)+b(n-1:1:-1)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s112 ')
      return
      end
C%1.1
      subroutine s113(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     linear dependence testing
C     a(i)=a(1) but no actual dependence cycle
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s113 ')
      call forttime(t1)
      do nl= 1,ntimes
	a(2:n)= a(1)+b(2:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s113 ')
      return
      end
C%1.1
      subroutine s114(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     linear dependence testing
C     transpose vectorization
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s114 ')
      call forttime(t1)
      do nl= 1,2*ntimes/n
!$omp parallel do if(n>103)
	  do j= 2,n
	      aa(1:j-1,j)= aa(j,1:j-1)+bb(1:j-1,j)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*(ntimes/n)))
      chksum= cs2d(n,aa)
      call check(chksum,2*(ntimes/n)*((n*n-n)/2),n,t2,'s114 ')
      return
      end
C%1.1
      subroutine s115(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     linear dependence testing
C     triangular saxpy loop
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s115 ')
      call forttime(t1)
      do nl= 1,2*ntimes/n
	  do j= 2,n
	      a(j)= a(j)-dot_product(aa(:j-1,j),a(:j-1))
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*(ntimes/n)))
      chksum= cs1d(n,a)
      call check(chksum,2*(ntimes/n)*((n*n-n)/2),n,t2,'s115 ')
      return
      end
C%1.1
      subroutine s116(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     linear dependence testing
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s116 ')
      call forttime(t1)
      do nl= 1,5*ntimes
	  a(:n-1)= a(2:n)*a(:n-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(5*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,5*ntimes*(n/5),n,t2,'s116 ')
      return
      end
C%1.1
      subroutine s118(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     linear dependence testing
C     potential dot product recursion
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s118 ')
      call forttime(t1)
      do nl= 1,2*ntimes/n
	  do i= 2,n
	      a(i)= a(i)+dot_product(bb(i,i-1:1:-1),a(:i-1))
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*(ntimes/n)))
      chksum= cs1d(n,a)
      call check(chksum,2*(ntimes/n)*((n*n-n)/2),n,t2,'s118 ')
      return
      end
C%1.1
      subroutine s119(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     linear dependence testing
C     no dependence - vectorizable
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s119 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	  do j= 2,n
	      aa(2:n,j)= aa(:n-1,j-1)+bb(2:n,j)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)
      call check(chksum,(ntimes/n)*(n-1)*(n-1),n,t2,'s119 ')
      return
      end
C%1.2
      subroutine s121(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     induction variable recognition
C     loop with possible ambiguity because of scalar store
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s121 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(:n-1)= a(2:n)+b(:n-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s121 ')
      return
      end
C%1.2
      subroutine s122(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,n1,n3)
      use lcd_mod
C
C     induction variable recognition
C     variable lower and upper bound, and stride
C
      integer ntimes,ld,n,i,nl,j,k,n1,n3
      real a(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
#if __INTEL_COMPILER
      real :: b(n)
#else
      real,volatile :: b(n)     ! stop vectorization
#endif
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s122 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(n1:n:n3)= a(n1:n:n3)+b(n:n-(n-n1)/n3:-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s122 ')
      return
      end
C%1.2
      subroutine s123(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     induction variable recognition
C     induction variable under an if
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s123 ')
      call forttime(t1)
      do nl= 1,2*ntimes
	  j= 0
	  do i= 1,n/2
	      j= j+1
	      a(j)= b(i)+d(i)*e(i)
	      if(c(i) >  0.)then
		  j= j+1
		  a(j)= c(i)+d(i)*e(i)
		endif
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,2*ntimes*(n/2),n,t2,'s123 ')
      return
      end
C%1.2
      subroutine s124(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     induction variable recognition
C     induction variable under both sides of if (same value)
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s124 ')
      call forttime(t1)
      do nl= 1,2*ntimes
!dir$ vector aligned
	a(:n/2)= merge(b(:n/2),c(:n/2),b(:n/2)>0)+d(:n/2)*e(:n/2)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,2*ntimes*(n/2),n,t2,'s124 ')
      return
      end
C%1.2
      subroutine s125(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     induction variable recognition
C     induction variable in two loops; collapsing possible
C
      integer ntimes,ld,n,i,nl,j,k,nn
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,array
      parameter(nn= 1000)
      common/cdata /array(nn*nn)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s125 ')
      call forttime(t1)
      do nl= 1,ntimes/n
!$omp parallel do if(n>101)
	  do j= 1,n
		array((j-1)*n+1:j*n)= aa(:n,j)+bb(:n,j)*cc(:n,j)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs1d(n*n,array)
      call check(chksum,(ntimes/n)*n*n,n,t2,'s125 ')
      return
      end
C%1.2
      subroutine s126(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     induction variable recognition
C     induction variable in two loops; recurrence in inner loop
C
      integer ntimes,ld,n,i,nl,j,k,nn
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d,array
      parameter(nn= 1000)
      common/cdata /array(nn*nn)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s126 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	  do j= 2,n-1,2
	      bb(:n,j)= bb(:n,j-1)+array(j-1:n*n+j-1:n)*cc(:n,j)
	      bb(:n,j+1)= bb(:n,j)+array(j:n*n+j:n)*cc(:n,j+1)
	    enddo
	  if(j==n) bb(:n,j)= bb(:n,j-1)+array(j-1:n*n+j-1:n)*cc(:n,j)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,bb)
      call check(chksum,(ntimes/n)*n*(n-1),n,t2,'s126 ')
      return
      end
C%1.2
      subroutine s127(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     induction variable recognition
C     induction variable with multiple increments
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s127 ')
      call forttime(t1)
      do nl= 1,2*ntimes
!dir$ ivdep
	forall(i=1:n/2)a(i*2-1:i*2)= b(i)+[c(i),e(i)]*d(i)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,2*ntimes*(n/2),n,t2,'s127 ')
      return
      end
C%1.2
      subroutine s128(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     induction variables
C     coupled induction variables
C
      integer ntimes,ld,n,i,nl,j,k
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s128 ')
      call forttime(t1)
      do nl= 1,2*ntimes
#if ! __MIC__ && ! __AVX512F__
!dir$ no vector
!dir$ unroll(4)
#endif
	  a(:n/2)= b(:n-1:2)-d(:n/2)
#if ! __MIC__ && ! __AVX512F__
!dir$ no vector
!dir$ unroll(4)
#endif
	  b(:n-1:2)= a(:n/2)+c(:n-1:2)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,2*ntimes*(n/2),n,t2,'s128 ')
      return
      end
C%1.3
      subroutine s131(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     global data flow analysis
C     forward substitution
C
      integer ntimes,ld,n,i,nl,m
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      m= 1
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s131 ')
      if(a(1) >  0) a(1)= b(1)
      call forttime(t1)
      do nl= 1,ntimes
	    a(:n-m)= a(m+1:n)+b(:n-m)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s131 ')
      return
      end
C%1.3
      subroutine s132(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     global data flow analysis
C     loop with multiple dimension ambiguous subscripts
C
      integer ntimes,ld,n,i,nl,j,k,m
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      m= 1
      j= m
      k= m+1
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s132 ')
      call forttime(t1)
      do nl= 1,ntimes
	aa(2:n,j)= aa(:n-1,k)+b(2:n)*c(k)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs2d(n,aa)
      call check(chksum,ntimes*n-1,n,t2,'s132 ')
      return
      end
C%1.4
      subroutine s141(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     nonlinear dependence testing
C     walk a row in a symmetric packed array
C     element a(i,j) for (j>i) stored in location j*(j-1)/2+i
C
#ifdef _OPENMP
      USE omp_lib
#endif
      integer ntimes,ld,n,i,nl,j,k,nn
      integer nt,m,jbeg,jend
      real avgchunk
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,array
      parameter(nn= 1000)
      common/cdata /array(nn*nn)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s141 ')
      call forttime(t1)
      do nl= 1,ntimes/n
!$omp parallel private(k,jbeg,jend,nt,avgchunk) if(n>103)
	  avgchunk=(n+1)*n/2
	  if(n>103)then
#if defined _OPENMP
	      nt=omp_get_num_threads()
! fair number of array elements per thread
	      avgchunk=avgchunk/nt
#else
	      nt=1
#endif
	    else
	      nt=1
	    endif
!$omp do
	  do m= 1,nt
	    jbeg=sqrt(.25+2*(m-1)*avgchunk)+1
! closest approximation to targeted chunk size
	      jend=sqrt(.25+2*m*avgchunk)
	      do j= jbeg,jend
		  k= j*(j-1)/2
		  array(k+1:k+j)= array(k+1:k+j)+bb(:j,j)
		enddo
	    enddo
!$omp end parallel
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs1d(n*n,array)
      call check(chksum,(ntimes/n)*n*n,n,t2,'s141 ')
      return
      end
C%1.5
      subroutine s151(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     interprocedural data flow analysis
C     passing parameter information into a subroutine
C
      integer ntimes,ld,n,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s151 ')
      call forttime(t1)
      do nl= 1,ntimes
	  call s151s(a,b,n,1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s151 ')
      return
      contains
      subroutine s151s(a,b,n,m)
      integer i,n,m
      real a(n),b(n)
	a(:n-m)= a(1+m:n)+b(:n-m)
      return
      end subroutine s151s
      end
C%1.5
      subroutine s152(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     interprocedural data flow analysis
C     collecting information from a subroutine
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s152 ')
      call forttime(t1)
      do nl= 1,ntimes
	  do i= 1,n
	      b(i)= d(i)*e(i)
	      call s152s(a,b,c,i)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s152 ')
      return
      contains
      subroutine s152s(a,b,c,i)
      integer i
      real a(*),b(*),c(*)
      a(i)= a(i)+b(i)*c(i)
      return
      end subroutine s152s
      end
C%1.6
      subroutine s161(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     tests for recognition of loop independent dependences
C     between statements in mutually exclusive regions.
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s161 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ assume_aligned a:64,d:64,b:64,e:64,c:64
!dir$ vector always
	where(b(:n-1) < 0) c(2:n)= a(:n-1)+d(:n-1)*d(:n-1)
!dir$ vector aligned
	where(.not. b(:n-1) < 0) a(:n-1)= c(:n-1)+d(:n-1)*e(:n-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,c)
      call check(chksum,ntimes*(n-1),n,t2,'s161 ')
      return
      end
C%1.6
      subroutine s162(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,k)
      use lcd_mod
C
C     control flow
C     deriving assertions
C
      integer ntimes,ld,n,i,nl,k
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s162 ')
      call forttime(t1)
      do nl= 1,ntimes
	  if(k > 0) a(:n-k)= a(1+k:n)+b(:n-k)*c(:n-k)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s162 ')
      return
      end
C%1.7
      subroutine s171(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,inc)
      use lcd_mod
C
C     symbolics
C     symbolic dependence tests
C
      integer ntimes,ld,n,i,nl,inc
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s171 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ assume_aligned a:64,b:64
	    a(inc:inc*n:inc)= a(inc:inc*n:inc)+b(:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s171 ')
      return
      end
C%1.7
      subroutine s172(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,n1,n3)
      use lcd_mod
C
C     symbolics
C     vectorizable if n3 .ne. 0
C
      integer ntimes,ld,n,i,nl,n1,n3
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s172 ')
      call forttime(t1)
      do nl= 1,ntimes
	a(n1:n:n3)= a(n1:n:n3)+b(n1:n:n3)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s172 ')
      return
      end
C%1.7
      subroutine s173(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     symbolics
C     expression in loop bounds and subscripts
C
      integer ntimes,ld,n,i,nl,k
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      k= n/2
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s173 ')
      call forttime(t1)
      do nl= 1,2*ntimes
!dir$ vector unaligned
	a(1+k:2*k)= a(:k)+b(:k)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,2*ntimes*(n/2),n,t2,'s173 ')
      return
      end
C%1.7
      subroutine s174(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     symbolics
C     loop with subscript that may seem ambiguous
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s174 ')
      call forttime(t1)
      do nl= 1,2*ntimes
	a(:n/2)= a(n/2+1:n)+b(:n/2)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,2*ntimes*(n/2),n,t2,'s174 ')
      return
      end
C%1.7
      subroutine s175(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,inc)
      use lcd_mod
C
C     symbolics
C     symbolic dependence tests
C
      integer ntimes,ld,n,i,nl,inc
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s175 ')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n-inc:inc)= a(inc+1:n:inc)+b(:n-inc:inc)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s175 ')
      return
      end
C%1.7
      subroutine s176(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     symbolics
C     convolution
C
      integer ntimes,ld,n,i,nl,j,m
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      m= n/2
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s176 ')
      call forttime(t1)
      do nl= 1,4*ntimes/n
!$omp parallel do if(n>103)
	do i= 1,m  
	  a(i)= a(i)+dot_product(b(i:i+m-1),c(m:1:-1))
	enddo
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(4*(ntimes/n)))
      chksum= cs1d(n,a)
      call check(chksum,4*(ntimes/n)*(n/2)*(n/2),n,t2,'s176 ')
      return
      end
C
C**********************************************************
C                                                         *
C                      VECTORIZATION                      *
C                                                         *
C**********************************************************
C%2.1
      subroutine s211(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     statement reordering
C     statement reordering allows vectorization
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s211 ')
      call forttime(t1)
      do nl= 1,ntimes
	  b(2:n-1)= b(3:n)-e(2:n-1)*d(2:n-1)
	  a(2:n-1)= b(1:n-2)+c(2:n-1)*d(2:n-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-2),n,t2,'s211 ')
      return
      end
C%2.1
      subroutine s212(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     statement reordering
C     dependency needing temporary
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s212 ')
      call forttime(t1)
      do nl= 1,ntimes
	  b(:n-1)= b(:n-1)+a(2:n)*d(:n-1)
	  a(:n-1)= a(:n-1)*c(:n-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s212 ')
      return
      end
C%2.2
      subroutine s221(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop distribution
C     loop that is partially recursive
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s221 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(2:n)= a(2:n)+c(2:n)*d(2:n)
	  do i= 2,n
	      b(i)= (d(i)+a(i))+b(i-1)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s221 ')
      return
      end
C%2.2
      subroutine s222(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop distribution
C     partial loop vectorization, recurrence in middle
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s222 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(2:n)= a(2:n)+b(2:n)*c(2:n)
	  do i= 2,n
	      b(i)= b(i-1)*b(i-1)*a(i)
	    enddo
          a(2:n)= a(2:n)-b(2:n)*c(2:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s222 ')
      return
      end
C%2.3
      subroutine s231(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop interchange
C     loop with multiple dimension recursion
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s231 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	do j= 2,n
	  aa(:n,j)= aa(:n,j-1)+bb(:n,j)
	enddo
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)
      call check(chksum,(ntimes/n)*n*(n-1),n,t2,'s231 ')
      return
      end
C%2.3
      subroutine s232(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop interchange
C     interchanging of triangular loops
C
#ifdef _OPENMP
      USE omp_lib
#endif
      integer ntimes,ld,n,i,nl,j
      integer nt,m,jbeg,jend
      real avgchunk
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s232 ')
      call forttime(t1)
      do nl= 1,2*ntimes/n
!$omp parallel private(nt,avgchunk,jbeg,jend) if(n>103)
	  avgchunk=(n-1)*n/2
	  if(n>103)then
#if defined _OPENMP
	      nt=omp_get_num_threads()
! fair number of array elements per thread
	      avgchunk=avgchunk/nt
#else
	      nt=1
#endif
	    else
	      nt=1
	    endif
!$omp do
	  do m= 1,nt
	      jbeg=sqrt(.25+2*(m-1)*avgchunk)+2
	      jend=sqrt(.25+2*m*avgchunk)+1
	      do j= jbeg,jend
		do i= 2,j
		  aa(i,j)= aa(i-1,j)*aa(i-1,j)+bb(i,j)
		  enddo
		enddo
	    enddo
!$omp end parallel
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(2*(ntimes/n)))
      chksum= cs2d(n,aa)
      call check(chksum,2*(ntimes/n)*((n*n-n)/2),n,t2,'s232 ')
      return
      end
C%2.3
      subroutine s233(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop interchange
C     interchanging with one of two inner loops
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s233 ')
      call forttime(t1)
      do nl= 1,ntimes/n
!$omp parallel if(n>103)
!$omp do
!dir$ no vector
	do j= 2,n
	  do i= 2,n
	    bb(i,j)= bb(i-1,j)+cc(i,j)
	    enddo
	  enddo
!$omp end do nowait
!$omp single
	do j= 2,n
	    aa(2:n,j)= aa(2:n,j-1)+cc(2:n,j)
	    enddo
!$omp end single
!$omp end parallel
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)+cs2d(n,bb)
      call check(chksum,(ntimes/n)*(n-1)*(2*n-2),n,t2,'s233 ')
      return
      end
C%2.3
      subroutine s234(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop interchange
C     if loop to do loop, interchanging with if loop necessary
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s234 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	  do j=2,n
	      aa(:n,j)= aa(:n,j-1)+bb(:n,j-1)*cc(:n,j-1)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)
      call check(chksum,(ntimes/n)*n*(n-1),n,t2,'s234 ')
      return
      end
C%2.3
      subroutine s235(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop interchanging
C     imperfectly nested loops
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s235 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	  a(:n)= a(:n)+b(:n)*c(:n)
	  do j= 2,n
	    aa(:n,j)= aa(:n,j-1)+bb(:n,j)*a(:n)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)+cs1d(n,a)
      call check(chksum,(ntimes/n)*n*(n-1),n,t2,'s235 ')
      return
      end
C%2.4
      subroutine s241(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     node splitting
C     preloading necessary to allow vectorization
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      real tmp
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s241 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ assume_aligned a:64,b:64,c:64,d:64
	  do i= 1,n-1
	      tmp= a(i+1)*d(i)
	      a(i)= b(i)*c(i)*d(i)
	      b(i)= a(i)*tmp
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s241 ')
      return
      end
C%2.4
      subroutine s242(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,s1,s2)
      use lcd_mod
C
C     node splitting
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,s1,s2
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s242 ')
      call forttime(t1)
      do nl= 1,ntimes
	  do i= 2,n
	    a(i)= ((s1+s2)+b(i)+c(i)+d(i))+a(i-1)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s242 ')
      return
      end
C%2.4
      subroutine s243(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     node splitting
C     false dependence cycle breaking
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s243 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  b(:n-1)= b(:n-1)+(c(:n-1)+e(:n-1))*d(:n-1)
	  a(:n-1)= b(:n-1)+a(2:n)*d(:n-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s243 ')
      return
      end
C%2.4
      subroutine s244(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     node splitting
C     false dependence cycle breaking
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s244 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:n-1)= b(:n-1)+c(:n-1)*d(:n-1)
!dir$ vector aligned
	  b(:n-1)= c(:n-1)+b(:n-1)
	a(n)= b(n-1)+a(n)*d(n-1)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s244 ')
      return
      end
C%2.5
      subroutine s251(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar and array expansion
C     scalar expansion
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,s
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s251 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:n)= (b(:n)+c(:n)*d(:n))**2
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s251 ')
      return
      end
C%2.5
      subroutine s252(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar and array expansion
C     loop with ambiguous scalar temporary
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,s,t
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s252 ')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)=b(:n)*c(:n)+eoshift(b(:n)*c(:n),1)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s252 ')
      return
      end
C%2.5
      subroutine s253(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar and array expansion
C     scalar expansion, assigned under if
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,s
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s253 ')
      call forttime(t1)
      do nl= 1,ntimes
#if __INTEL_COMPILER || __GNUC_ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)
!dir$ vector aligned
	  do concurrent( i= 1:n, a(i) > b(i))
	    a(i)= a(i)-b(i)*d(i)
	    c(i)= a(i)+c(i)
	  enddo
#else
	  forall ( i= 1:n, a(i) > b(i))
	    a(i)= a(i)-b(i)*d(i)
	    c(i)= a(i)+c(i)
	  endforall
#endif
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,c)
      call check(chksum,ntimes*n,n,t2,'s253 ')
      return
      end
C%2.5
      subroutine s254(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar and array expansion
C     carry around variable
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,x
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s254 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(:n)= (b(:n)+cshift(b(:n),1))*.5
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s254 ')
      return
      end
C%2.5
      subroutine s255(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar and array expansion
C     carry around variables, 2 levels
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,x,y
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s255 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(:n)= (b(:n)+cshift(b(:n),1)+cshift(b(:n),2))*.333
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s255 ')
      return
      end
C%2.5
      subroutine s256(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar and array expansion
C     array expansion
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s256 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	do j= 2,n
        !dir$ distribute point
C ifort is sensitive to operand order
	    aa(:n,j)= aa(:n,j)+bb(:n,j)-a(j-1)
	    a(j)=aa(n,j)-bb(n,j)
	  enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs1d(n,a)+cs2d(n,aa)
      call check(chksum,(ntimes/n)*n*(n-1),n,t2,'s256 ')
      return
      end
C%2.5
      subroutine s257(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar and array expansion
C     array expansion
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,cs2d
      real tmp
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s257 ')
      call forttime(t1)
      do nl= 1,ntimes/n
!$omp parallel if(n>103)
!$omp single
	  do i= 2,n
		a(i)= aa(i,n)-a(i-1)
		aa(i,n)= a(i)+bb(i,n)
	      enddo
!$omp end single
!$omp do private(tmp)
	do j= 1,n-1
	  tmp= a(1)
	  do i= 2,n
		tmp= aa(i,j)-tmp
		aa(i,j)= tmp+bb(i,j)
	      enddo
	    enddo
!$omp end do nowait
!$omp end parallel
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs1d(n,a)+cs2d(n,aa)
      call check(chksum,(ntimes/n)*(n-1)*n,n,t2,'s257 ')
      return
      end
C%2.5
      subroutine s258(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar and array expansion
C     wrap-around scalar under an if
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,s(n),ss
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s258 ')
      call forttime(t1)
      do nl= 1,ntimes
	  ss = 0
	  do i= 1,n
	      ss= merge(d(i)**2,ss,a(i) > 0.)
	      s(i) = ss
	    enddo
!dir$ no fusion
	  b(:n)= s(:n)*c(:n)+d(:n)
	  e(:n)= s(:n)*aa(:n,1)+aa(:n,1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,b)+cs1d(n,e)
      call check(chksum,ntimes*n,n,t2,'s258 ')
      return
      end
C%2.6
      subroutine s261(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     scalar renaming
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,t
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s261 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ assume_aligned a:64,b:64,c:64,d:64
	  c(2:n)= c(2:n)*d(2:n)
!dir$ nofusion
	  a(2:n)= a(2:n)+b(2:n)+c(:n-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,c)
      call check(chksum,ntimes*(n-1),n,t2,'s261 ')
      return
      end
C%2.7
      subroutine s271(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     loop with singularity handling
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s271 ')
      call forttime(t1)
      do nl= 1,ntimes
	    a(:n)= a(:n)+max(b(:n),0.)*c(:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s271 ')
      return
      end
C%2.7
      subroutine s272(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,t)
      use lcd_mod
C
C     control flow
C     loop with independent conditional
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,t
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s272 ')
      call forttime(t1)
      do nl= 1,ntimes
#if __INTEL_COMPILER || __GNUC_ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7)
!dir$ vector aligned
	  do concurrent( i= 1:n, e(i) >= t)
	    a(i)= a(i)+c(i)*d(i)
	    b(i)= b(i)+c(i)*c(i)
	    enddo
#else
	  forall ( i= 1:n, e(i) >= t)
	    a(i)= a(i)+c(i)*d(i)
	    b(i)= b(i)+c(i)*c(i)
	    endforall
#endif
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*n,n,t2,'s272 ')
      return
      end
C%2.7
      subroutine s273(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     simple loop with dependent conditional
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s273 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:n)= a(:n)+d(:n)*e(:n)
!dir$ vector aligned
	  c(:n)= c(:n)+a(:n)*d(:n)
!dir$ vector aligned
          b(:n)= b(:n)+merge(d(:n),0.,a(:n)<0)*e(:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)+cs1d(n,c)
      call check(chksum,ntimes*n,n,t2,'s273 ')
      return
      end
C%2.7
      subroutine s274(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     complex loop with dependent conditional
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s274 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:n)= c(:n)+e(:n)*d(:n)
!dir$ vector aligned
	  b(:n)= max(a(:n),0.)+b(:n)
!dir$ vector aligned
	  where(a(:n) <= 0.) a(:n)= d(:n)*e(:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*n,n,t2,'s274 ')
      return
      end
C%2.7
      subroutine s275(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     if around inner loop, interchanging needed
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s275 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	  do j= 2,n
	      where(aa(2:n,1) > 0.)
		aa(2:n,j)= aa(2:n,j-1)+bb(2:n,j)*cc(2:n,j)
	      endwhere
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)
      call check(chksum,(ntimes/n)*(n-1)*(n-1),n,t2,'s275 ')
      return
      end
C%2.7
      subroutine s276(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     if test using loop index
C
      integer ntimes,ld,n,i,nl,mid
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s276 ')
      call forttime(t1)
      mid= n/2
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:mid-1)= a(:mid-1)+b(:mid-1)*c(:mid-1)
	  a(mid:n)= a(mid:n)+b(mid:n)*d(mid:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s276 ')
      return
      end
C%2.7
      subroutine s277(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     test for dependences arising from guard variable computation.
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s277 ')
      call forttime(t1)
      do nl= 1,ntimes
! not valid for single do concurrent
!dir$ assume_aligned a:64,b:64,c:64,d:64,e:64
	  where(a(:n-1) < 0)b(2:n)= c(:n-1)+d(:n-1)*e(:n-1)
!dir$ nofusion
          a(:n-1)= a(:n-1)+merge(d(:n-1),0.,max(a(:n-1),b(:n-1)) < 0)
     &*c(:n-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s277 ')
      return
      end
C%2.7
      subroutine s278(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     if/goto to block if-then-else
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s278 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  where(a(:n) <= 0.) b(:n)= -b(:n)+d(:n)*e(:n)
!dir$ vector aligned
	  where(.not. a(:n) <= 0.) c(:n)= -c(:n)+d(:n)*e(:n)
!dir$ vector aligned
	  a(:n)= b(:n)+c(:n)*d(:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)+cs1d(n,c)
      call check(chksum,ntimes*n,n,t2,'s278 ')
      return
      end
C%2.7
      subroutine s279(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     vector if/gotos
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s279 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  where(a(:n) > 0) c(:n)= -c(:n)+e(:n)*e(:n)
!dir$ vector aligned
	  where(.not. a(:n) > 0)
		b(:n)= -b(:n)+d(:n)*d(:n)
		c(:n)= c(:n)+merge(d(:n)*e(:n),0.,b(:n) > a(:n)) 
	    endwhere
!dir$ vector aligned
	  a(:n)= b(:n)+c(:n)*d(:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)+cs1d(n,c)
      call check(chksum,ntimes*n,n,t2,'s279 ')
      return
      end
C%2.7
      subroutine s2710(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,x)
      use lcd_mod
C
C     control flow
C     scalar and vector ifs
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,x
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s2710')
      call forttime(t1)
      do nl= 1,ntimes
	where(a(:n) > b(:n))
	    a(:n)= a(:n)+b(:n)*d(:n)
	    c(:n)= merge(c(:n)+d(:n)*d(:n),1.+d(:n)*e(:n),n > 10)
	  elsewhere
	      b(:n)= a(:n)+e(:n)*e(:n)
	    c(:n)=merge(a(:n)+d(:n)*d(:n),c(:n)+e(:n)*e(:n),x>0.)
	  endwhere
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)+cs1d(n,c)
      call check(chksum,ntimes*n,n,t2,'s2710')
      return
      end
C%2.7
      subroutine s2711(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     semantic if removal
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s2711')
      call forttime(t1)
      do nl= 1,ntimes
	  a(:n)= a(:n)+merge(b(:n),0.,b(:n).ne.0.) *c(:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s2711')
      return
      end
C%2.7
      subroutine s2712(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control flow
C     if to elemental min
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s2712')
      call forttime(t1)
      do nl= 1,ntimes
	  a(:n)= a(:n)+merge(b(:n),0.,a(:n) > b(:n)) *c(:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s2712')
      return
      end
C%2.8
      subroutine s281(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     crossing thresholds
C     index set splitting
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,x
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s281 ')
      call forttime(t1)
      do nl= 1,ntimes
	  b(:(n+1)/2)= a(n:n/2+1:-1)+b(:(n+1)/2)*c(:(n+1)/2)
	  a(:(n+1)/2)= b(:(n+1)/2)-1.0
	  b((n+3)/2:n)= a(n/2:1:-1)+b((n+3)/2:n)*c((n+3)/2:n)
	  a((n+3)/2:n)= b((n+3)/2:n)-1.0
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*n,n,t2,'s281 ')
      return
      end
C%2.9
      subroutine s291(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop peeling
C     wrap around variable, 1 level
C
      integer ntimes,ld,n,i,nl,im1
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s291 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:n)= (b(:n)+cshift(b(:n),1))*.5
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s291 ')
      return
      end
C%2.9
      subroutine s292(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop peeling
C     wrap around variable, 2 levels
C
      integer ntimes,ld,n,i,nl,im1,im2
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s292 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(:n)= (b(:n)+cshift(b(:n),1)+cshift(b(:n),2))*.333
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s292 ')
      return
      end
C%2.9
      subroutine s293(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop peeling
C     a(i)=a(1) with actual dependence cycle, loop is vectorizable
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s293 ')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(1)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s293 ')
      return
      end
C%2.10
      subroutine s2101(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     diagonals
C     main diagonal calculation
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s2101')
      call forttime(t1)
      do nl= 1,ntimes
	  forall(i= 1:n) aa(i,i)= aa(i,i)+bb(i,i)*cc(i,i)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs2d(n,aa)
      call check(chksum,ntimes*n,n,t2,'s2101')
      return
      end
C%2.12
      subroutine s2102(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     diagonals
C     identity matrix, best results vectorize both inner and outer loops
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s2102')
      call forttime(t1)
      do nl= 1,ntimes/n
          forall(j=1:n,i= 1:n) aa(j,i)= merge(1.,0.,j==i)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)
      call check(chksum,(ntimes/n)*n*n,n,t2,'s2102')
      return
      end
C%2.11
      subroutine s2111(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     wavefronts
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s2111')
      call forttime(t1)
      do nl= 1,ntimes/n
	  do j= 2,n
	    do i= 2,n
	      aa(i,j)= aa(i-1,j)+aa(i,j-1)
	      enddo
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)
      if(chksum == 0.) chksum= 3.0
      call check(chksum,(ntimes/n)*(n-1)*(n-1),n,t2,'s2111')
      return
      end
C
C**********************************************************
C                                                         *
C                   IDIOM RECOGNITION                     *
C                                                         *
C**********************************************************
C%3.1
      subroutine s311(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     sum reduction
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,sumvar
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s311 ')
      call forttime(t1)
      do nl= 1,ntimes
	    sumvar= sum(a(:n))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,sumvar)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= sumvar
      call check(chksum,ntimes*n,n,t2,'s311 ')
      return
      end
C%3.1
      subroutine s312(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     product reduction
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,prod
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s312 ')
      call forttime(t1)
      do nl= 1,ntimes
	  prod= product(a(:n))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,prod)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= prod
      call check(chksum,ntimes*n,n,t2,'s312 ')
      return
      end
C%3.1
      subroutine s313(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     dot product
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,dot
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s313 ')
      call forttime(t1)
      do nl= 1,ntimes
	dot= dot_product(a(:n),b(:n))
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,dot)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= dot
      call check(chksum,ntimes*n,n,t2,'s313 ')
      return
      end
C%3.1
      subroutine s314(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     if to max reduction
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,x
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s314 ')
      call forttime(t1)
      do nl= 1,ntimes
	  x= maxval(a(:n))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,x)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= x
      call check(chksum,ntimes*n,n,t2,'s314 ')
      return
      end
C%3.1
      subroutine s315(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     if to max with index reduction, 1 dimension
C
      integer ntimes,ld,n,i,nl,index
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,x
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s315 ')
      call forttime(t1)
      do nl= 1,ntimes
	    index= maxloc(a(:n),dim=1)
	    x= a(index)
	  chksum= x+float(index)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,chksum)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= x+float(index)
      call check(chksum,ntimes*n,n,t2,'s315 ')
      return
      end
C%3.1
      subroutine s316(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     if to min reduction
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,x
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s316 ')
      call forttime(t1)
      do nl= 1,ntimes
	  x= minval(a(:n))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,x)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= x
      call check(chksum,ntimes*n,n,t2,'s316 ')
      return
      end
C%3.1
      subroutine s317(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     product reduction, vectorize with
C     1. scalar expansion of factor, and product reduction
C     2. closed form solution: q = factor**n
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,q,factor
      parameter(factor= .99999)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s317 ')
      call forttime(t1)
      do nl= 1,ntimes
	  q= 1.
	  do i= 1,n
	    q= factor*q
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,q)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= q
      call check(chksum,ntimes*n,n,t2,'s317 ')
      return
      end
C%3.1
      subroutine s318(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,inc)
      use lcd_mod
C
C     reductions
C     isamax, max absolute value, increments not equal to 1
C
C
      integer ntimes,ld,n,i,nl,inc,k,index
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,max
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s318 ')
      call forttime(t1)
      do nl= 1,ntimes
	  index= maxloc(abs(a(:1+(n-1)*inc:inc)),dim=1)
	  max= abs(a((index-1)*inc+1))
	  chksum= max+float(index)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,chksum)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= max+float(index)
      call check(chksum,ntimes*(n-1),n,t2,'s318 ')
      return
      end
C%3.1
      subroutine s319(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     coupled reductions
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,sumvar
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s319 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	      a(:n)= c(:n)+d(:n)
!dir$ vector aligned
	      b(:n)= c(:n)+e(:n)
	      sumvar=sum(a(:n)+b(:n))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,sumvar)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= sumvar
      call check(chksum,ntimes*n,n,t2,'s319 ')
      return
      end
C%3.1
      subroutine s3110(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     if to max with index reduction, 2 dimensions
C
      integer ntimes,ld,n,i,nl,j,xindex,yindex,ml(2)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,max
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s3110')
      call forttime(t1)
      do nl= 1,ntimes/n
!$omp parallel workshare if(n>103)
	  ml= maxloc(aa(:n,:n))
!$omp end parallel workshare
	  xindex= ml(1)
	  yindex= ml(2)
	  max= aa(xindex,yindex)
	  chksum= max+float(xindex)+float(yindex)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,chksum)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= max+float(xindex)+float(yindex)
      call check(chksum,(ntimes/n)*n*n,n,t2,'s3110')
      return
      end
C%3.1
      subroutine s3111(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     conditional sum reduction
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,sumvar
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s3111')
      call forttime(t1)
      do nl= 1,ntimes
	  sumvar= sum(max(a(:n),0.))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,sumvar)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= sumvar
      call check(chksum,ntimes*n,n,t2,'s3111')
      return
      end
C%3.1
      subroutine s3112(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     sum reduction saving running sums
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,sum
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s3112')
      call forttime(t1)
      do nl= 1,ntimes
	  sum= 0.
	  do i= 1,n
	      sum= sum+a(i)
	      b(i)= sum
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,sum)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,b)+sum
      call check(chksum,ntimes*n,n,t2,'s3112')
      return
      end
C%3.1
      subroutine s3113(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     reductions
C     maximum of absolute value
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,max
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s3113')
      call forttime(t1)
      do nl= 1,ntimes
	  max= maxval(abs(a(:n)))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,max)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= max
      call check(chksum,ntimes*(n-1),n,t2,'s3113')
      return
      end
C%3.2
      subroutine s321(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     recurrences
C     first order linear recurrence
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s321 ')
      call forttime(t1)
      do nl= 1,ntimes
	  do i= 2,n
	    a(i)= a(i)+a(i-1)*b(i)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s321 ')
      return
      end
C%3.2
      subroutine s322(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     recurrences
C     second order linear recurrence
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s322 ')
      call forttime(t1)
      do nl= 1,ntimes
	  do i= 3,n
	    a(i)= (a(i)+a(i-2)*c(i))+a(i-1)*b(i)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-2),n,t2,'s322 ')
      return
      end
C%3.2
      subroutine s323(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     recurrences
C     coupled recurrence
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s323 ')
      call forttime(t1)
      do nl= 1,ntimes
	  do i= 2,n
	      a(i)= b(i-1)+c(i)*d(i)
	      b(i)= b(i-1)+(c(i)*(e(i)+d(i)))
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s323 ')
      return
      end
C%3.3
      subroutine s331(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     search loops
C     if to last-1
C
      integer ntimes,ld,n,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s331 ')
      call forttime(t1)
      do nl= 1,ntimes
	  j= findloc(a(:n) < 0,.true.,.true.)
	  chksum= float(j)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,chksum)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= float(j)
      call check(chksum,ntimes*n,n,t2,'s331 ')
      return
      contains
      integer function findloc(a,val,back)result(i)
      logical a(:),val,back
      if(back)then
      do i=size(a),1,-1
        if(a(i).eqv.val)exit
        enddo
      else
      do i=1,size(a)
        if(a(i).eqv.val)exit
        enddo
        endif
        if(i>size(a))i=0
        return
      end function findloc
      end
C%3.3
      subroutine s332(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,t)
      use lcd_mod
C
C     search loops
C     first value greater than threshold
C
      integer ntimes,ld,n,nl,index
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,t,value
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s332 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  index= findloc(a(:n) > t,.true.)
	  value= a(index)
          chksum= value+float(index)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,chksum)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= value+float(index)
      call check(chksum,ntimes*n,n,t2,'s332 ')
      return
      contains
      integer function findloc(a,val)result(i)
      logical a(:),val
      do i=1,size(a)
        if(a(i).eqv.val)exit
        enddo
        if(i>size(a))i=0
        return
      end function findloc
      end
C%3.4
      subroutine s341(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     packing
C     pack positive values
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s341 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(:n)= pack(b(:n),b(:n)>0)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s341 ')
      return
      end
C%3.4
      subroutine s342(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     packing
C     unpacking
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s342 ')
      call forttime(t1)
      do nl= 1,ntimes
	  a(:n)= unpack(b,a(:n)>0.,a(:n))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s342 ')
      return
      end
C%3.4
      subroutine s343(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     packing
C     pack 2-d array into one dimension
C
      integer ntimes,ld,n,i,nl,j,k,nn
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,array
      parameter(nn= 1000)
      common/cdata /array(nn*nn)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s343 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	  k= 0
	  do i= 1,n
	      array(k+1:k+n)= pack(aa(i,:n),bb(i,:n) > 0)
	      k= k + count(bb(i,:n) > 0)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs1d(n*n,array)
      call check(chksum,(ntimes/n)*n*n,n,t2,'s343 ')
      return
      end
C%3.5
      subroutine s351(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop rerolling
C     unrolled saxpy
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,alpha,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s351 ')
      call forttime(t1)
      alpha= c(1)
      do nl= 1,5*ntimes
	  do i= 1,n,5
	      a(i)= a(i)+alpha*b(i)
	      a(i+1)= a(i+1)+alpha*b(i+1)
	      a(i+2)= a(i+2)+alpha*b(i+2)
	      a(i+3)= a(i+3)+alpha*b(i+3)
	      a(i+4)= a(i+4)+alpha*b(i+4)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(5*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,5*ntimes*(n/5),n,t2,'s351 ')
      return
      end
C%3.5
      subroutine s352(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop rerolling
C     unrolled dot product
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,dot
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s352 ')
      call forttime(t1)
      do nl= 1,5*ntimes
	  dot= 0.
#if ! __MIC__ && ! __AVX512F__
!dir$ no vector
!$omp simd safelen(1)
#endif
	  do i= 1,n,5
	    dot= dot+(a(i)*b(i)+a(i+1)*b(i+1)+a(i+2)*b(i+2)+
     &a(i+3)*b(i+3)+a(i+4)*b(i+4))
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,dot)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(5*ntimes))
      chksum= dot
      call check(chksum,5*ntimes*(n/5),n,t2,'s352 ')
      return
      end
C%3.5
      subroutine s353(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip)
      use lcd_mod
C
C     loop rerolling
C     unrolled sparse saxpy
C
      integer ntimes,ld,n,i,nl,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,alpha,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s353 ')
      call forttime(t1)
      alpha= c(1)
      do nl= 1,5*ntimes
	  a(:n)= a(:n)+alpha*b(ip(:n))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(5*ntimes))
      chksum= cs1d(n,a)
      call check(chksum,5*ntimes*(n/5),n,t2,'s353 ')
      return
      end
C
C**********************************************************
C                                                         *
C                 LANGUAGE COMPLETENESS                   *
C                                                         *
C**********************************************************
C%4.1
      subroutine s411(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop recognition
C     if loop to do loop, zero trip
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s411 ')
      call forttime(t1)
      do nl= 1,ntimes
	  i= 0
	  do
	      i= i+1
	      if(i >  n) exit
	      a(i)= a(i)+b(i)*c(i)
23493       enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s411 ')
      return
      end
C%4.1
      subroutine s412(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,inc)
      use lcd_mod
C
C     loop recognition
C     if loop with variable increment
C
      integer ntimes,ld,n,i,nl,inc
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s412 ')
      call forttime(t1)
      do nl= 1,ntimes
	  i= 0
	  do
	      i= i+inc
	      if(i > n) exit
	      a(i)= a(i)+b(i)*c(i)
23500       enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s412 ')
      return
      end
C%4.1
      subroutine s413(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop recognition
C     if loop to do loop, code on both sides of increment
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s413 ')
      call forttime(t1)
      do nl= 1,ntimes
	  b(:n-1)= b(:n-1)+d(:n-1)*e(:n-1)
!dir$ nofusion
	  a(2:n)= c(2:n)+d(2:n)*e(2:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)+cs1d(n,b)
      call check(chksum,ntimes*(n-1),n,t2,'s413 ')
      return
      end
C%4.1
      subroutine s414(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop recognition
C     if loop to do loop, interchanging with do necessary
C
      integer ntimes,ld,n,i,nl,j
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs2d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s414 ')
      call forttime(t1)
      do nl= 1,ntimes/n
	do j= 2,n
	  i= 1
	  dowhile(i <= n)
	    aa(i,j)= aa(i,j-1)+bb(i,j-1)*cc(i,j-1)
	    i= i+1
	    enddo
	 enddo
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes/n))
      chksum= cs2d(n,aa)+cs2d(n,bb)+cs2d(n,cc)
      call check(chksum,(ntimes/n)*n*(n-2),n,t2,'s414 ')
      return
      end
C%4.1
      subroutine s415(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     loop recognition
C     while loop
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      integer tmp(n)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s415 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  i= findloc(a(:n) < 0,.true.)
!dir$ vector aligned
	  a(:i-1)= a(:i-1)+b(:i-1)*c(:i-1)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s415 ')
      return
      contains
      integer function findloc(a,val)result(i)
      logical a(:),val
      do i=1,size(a)
        if(a(i).eqv.val)exit
        enddo
        if(i>size(a))i=0
        return
      end function findloc
      end
C%4.2
      subroutine s421(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     storage classes and equivalencing
C     equivalence- no overlap
C authors comment on overlap is false
C but is vectorizable by checking direction at compile or run time
C
      integer ntimes,ld,n,i,nl,nn
      parameter(nn= 1000)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      real x(nn),y(nn)
      equivalence(x(1),y(1))
      call set1d(n,x,0.0,1)
      call set1d(n,y,1.0,1)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s421 ')
      call forttime(t1)
      do nl= 1,ntimes
	x(:n-1)= y(2:n)+a(:n-1)
	call dummy(ld,n,x,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,x)
      call check(chksum,ntimes*(n-1),n,t2,'s421 ')
      return
      end
C%4.2
      subroutine s422(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     storage classes and equivalencing
C     common and equivalence statement
C     anti-dependence, threshold of 4
C
      integer ntimes,ld,n,i,nl,nn,vl
      parameter(nn= 1000,vl= 64)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      common/cdata /array(nn*nn)
      real x(nn),array
      equivalence(x(1),array(5))
      real t1,t2,chksum,ctime,dtime,cs1d
      call set1d(n,x,0.0,1)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s422 ')
      call forttime(t1)
      do nl= 1,ntimes
	x(:n)= array(9:n+8)+a(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,x)
      call check(chksum,ntimes*(n-8),n,t2,'s422 ')
      return
      end
C%4.2
      subroutine s423(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     storage classes and equivalencing
C     common and equivalenced variables - with anti-dependence
C
      integer ntimes,ld,n,i,nl,nn,vl
      parameter(nn= 1000,vl= 64)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real x(nn),array
      common/cdata /array(nn*nn)
      equivalence(array(vl),x(1))
      real t1,t2,chksum,ctime,dtime,cs1d
      call set1d(n,x,1.0,1)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s423 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ simd
	array(2:n)= x(:n-1)+a(:n-1)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,array)
      call check(chksum,ntimes*(n-1),n,t2,'s423 ')
      return
      end
C%4.2
      subroutine s424(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     storage classes and equivalencing
C     common and equivalenced variables - overlap
C     vectorizeable in strips of 64 or less
C
      integer ntimes,ld,n,i,nl,nn,vl
      parameter(nn= 1000,vl= 64)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real x(nn),array
      real t1,t2,chksum,ctime,dtime,cs1d
      common/cdata /array(nn*nn)
      equivalence(array(vl),x(1))
      call set1d(n,x,0.0,1)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s424 ')
      call forttime(t1)
      do nl= 1,ntimes
cdec$ loop count avg (1000)
	x(2:n)= array(:n-1)+a(:n-1)
	call dummy(ld,n,x,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,x)
      call check(chksum,ntimes*(n-1),n,t2,'s424 ')
      return
      end
C%4.3
      subroutine s431(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     parameters
C     parameter statement
C
      integer ntimes,ld,n,i,nl,k,k1,k2
      parameter(k1= 1,k2= 2,k= 2*k1-k2)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s431 ')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(k+1:n+k)+b(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s431 ')
      return
      end
C%4.3
      subroutine s432(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     parameters
C     data statement
C
      integer ntimes,ld,n,i,nl,k,k1,k2
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      data k1,k2/1,2/
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s432 ')
      k= 2*k1-k2
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(k+1:n+k)+b(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s432 ')
      return
      end
C%4.4
      subroutine s441(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     non-logical ifs
C     arithmetic if translated automatically to if..else
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s441 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:n)= a(:n)+merge(b(:n),c(:n),d(:n)<=0)*
     &		merge(c(:n),b(:n),d(:n)/=0)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s441 ')
      return
      end
C%4.4
      subroutine s442(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,indx)
      use lcd_mod
C
C     non-logical ifs
C     computed goto
C
#ifdef _OPENMP
      USE omp_lib
#endif
      integer ntimes,ld,n,i,nl,indx(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s442 ')
      call forttime(t1)
      do nl= 1,ntimes
#if defined __INTEL_COMPILER
!$omp parallel do num_threads(2) if(n>103)
!dir$ no vector
#endif
	  do i= 1,n
	      select case(indx(i))
   		case(1)
			a(i) = a(i) + b(i) * b(i)
   		case(2)
			a(i) = a(i) + c(i) * c(i)
   		case(3)
			a(i) = a(i) + d(i) * d(i)
   		case(4)
			a(i) = a(i) + e(i) * e(i)
		endselect
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s442 ')
      return
      end
C%4.4
      subroutine s443(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     non-logical ifs
C     arithmetic if translated automatically to if..else
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s443 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:n)= a(:n)+b(:n)*merge(c(:n),b(:n),d(:n) <= 0)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s443 ')
      return
      end
C%4.5
      subroutine s451(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     intrinsic functions
C     intrinsics
C
#ifdef _OPENMP
      USE omp_lib
#endif
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s451 ')
      call forttime(t1)
      do nl= 1,ntimes
!$omp parallel workshare num_threads(min(6,omp_get_max_threads()))
!$omp& if(n>103)
	a(:n)= sin(b(:n))+cos(c(:n))
!$omp end parallel workshare
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s451 ')
      return
      end
C%4.5
      subroutine s452(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     intrinsic functions
C     seq function
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s452 ')
      call forttime(t1)
      do nl= 1,ntimes
	  forall( i= 1:n) a(i)= b(i)+c(i)*i
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s452 ')
      return
      end
C%4.5
      subroutine s453(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     intrinsic functions
C     seq function
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,s
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s453 ')
      call forttime(t1)
      do nl= 1,ntimes
	  forall( i= 1:n) a(i)= i*2*b(i)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s453 ')
      return
      end
C%4.7
      subroutine s471(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,c471s)
      use lcd_mod
C
C     call statements
C
      integer ntimes,ld,n,i,nl,nn,m
      parameter(nn= 1000)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,c471s,x(nn)
      m= n
      call set1d(n,x,0.0,1)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s471 ')
      call forttime(t1)
      do nl= 1,ntimes
!$omp parallel do if(n>103)
	  do i= 1,m
	      x(i)= b(i)+d(i)*d(i)
	      call s471s
	      b(i)= c(i)+d(i)*e(i)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))-((n*ntimes)*c471s)
      chksum= cs1d(n,x)+cs1d(n,b)
      call check(chksum,ntimes*n,n,t2,'s471 ')
      return
      end
C%4.8
      subroutine s481(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     non-local gotos
C     stop statement
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s481 ')
      call forttime(t1)
      do nl= 1,ntimes
          i=findloc(d(:n) < 0.,.true.) 
          if(i<1)i= n +1
          a(:i-1)= a(:i-1)+b(:i-1)*c(:i-1)
          if(d(i) < 0.) stop 'stop 1'
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s481 ')
      return
      stop
      contains
      integer function findloc(a,val)result(i)
      logical a(:),val
      do i=1,size(a)
        if(a(i).eqv.val)exit
        enddo
        if(i>size(a))i=0
        return
      end function findloc
      end
C%4.8
      subroutine s482(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     non-local gotos
C     other loop exit with code before exit
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s482 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
          i=findloc(c(:n) > b(:n),.true.)
          if(i<1)i=n
          a(:i)= a(:i)+b(:i)*c(:i)
          if(c(i) > b(i))return
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s482 ')
      return
      contains
      integer function findloc(a,val)result(i)
      logical a(:),val
      do i=1,size(a)
        if(a(i).eqv.val)exit
        enddo
        if(i>size(a))i=0
        return
      end function findloc
      end
C%4.9
      subroutine s491(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip)
      use lcd_mod
C
C     vector semantics
C     indirect addressing on lhs, store in sequence
C
      integer ntimes,ld,n,i,nl,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s491 ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector always
	  do i= 1,n
	    a(ip(i))= b(i)+c(i)*d(i)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s491 ')
      return
      end
C%4.11
      subroutine s4112(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip,s)
      use lcd_mod
C
C     indirect addressing
C     sparse saxpy
C
      integer ntimes,ld,n,i,nl,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,s
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s4112')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(:n)+b(ip(:n))*s
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s4112')
      return
      end
C%4.11
      subroutine s4113(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip)
      use lcd_mod
C
C     indirect addressing
C     indirect addressing on rhs and lhs
C
      integer ntimes,ld,n,i,nl,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s4113')
      call forttime(t1)
      do nl= 1,ntimes
	  do i= 1,n
	    a(ip(i))= b(ip(i))+c(i)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s4113')
      return
      end
C%4.11
      subroutine s4114(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip,n1)
      use lcd_mod
C
C     indirect addressing
C     mix indirect addressing with variable lower and upper bounds
C
      integer ntimes,ld,n,i,nl,k,n1,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s4114')
      call forttime(t1)
      do nl= 1,ntimes
	  a(n1:n)= b(n1:n)+c(n-ip(n1:n)+1)*d(n1:n)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s4114')
      return
      end
C%4.11
      subroutine s4115(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip)
      use lcd_mod
C
C     indirect addressing
C     sparse dot product
C
      integer ntimes,ld,n,i,nl,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,sum
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s4115')
      call forttime(t1)
      do nl= 1,ntimes
	sum= dot_product(a(:n),b(ip(:n)))
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= sum
      call check(chksum,ntimes*n,n,t2,'s4115')
      return
      end
C%4.11
      subroutine s4116(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip,j,i
     &nc)
      use lcd_mod
C
C     indirect addressing
C     more complicated sparse sdot
C
      integer ntimes,ld,n,i,nl,j,off,inc,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,sum
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s4116')
      call forttime(t1)
      do nl= 1,ntimes
	  sum= dot_product(a(1+inc:n-1+inc),aa(ip(:n-1),j))
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= sum
      call check(chksum,ntimes*(n-1),n,t2,'s4116')
      return
      end
C%4.11
      subroutine s4117(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     indirect addressing
C     seq function
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s4117')
      call forttime(t1)
      do nl= 1,ntimes
	  forall(i= 2:n) a(i)= b(i)+c(i/2)*d(i)
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*(n-1),n,t2,'s4117')
      return
      end
C%4.12
      subroutine s4121(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     statement functions
C     elementwise multiplication
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,f,x,y
      f(x,y)= x*y
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'s4121')
      call forttime(t1)
      do nl= 1,ntimes
	  do i= 1,n
	    a(i)= a(i)+f(b(i),c(i))
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'s4121')
      return
      end
C%5.1
      subroutine va(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector assignment
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'va   ')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= b(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'va   ')
      return
      end
C%5.1
      subroutine vag(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip)
      use lcd_mod
C
C     control loops
C     vector assignment, gather
C
      integer ntimes,ld,n,i,nl,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vag  ')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= b(ip(:n))
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vag  ')
      return
      end
C%5.1
      subroutine vas(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,ip)
      use lcd_mod
C
C     control loops
C     vector assignment, scatter
C
      integer ntimes,ld,n,i,nl,ip(n)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vas  ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector always
	  do i= 1,n
	    a(ip(i))= b(i)
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vas  ')
      return
      end
C%5.1
      subroutine vif(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector if
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vif  ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  a(:n)= merge(b(:n),a(:n),b(:n) > 0.) 
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vif  ')
      return
      end
C%5.1
      subroutine vpv(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector plus vector
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vpv  ')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(:n)+b(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vpv  ')
      return
      end
C%5.1
      subroutine vtv(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector times vector
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vtv  ')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(:n)*b(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vtv  ')
      return
      end
C%5.1
      subroutine vpvtv(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector plus vector times vector
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vpvtv')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(:n)+b(:n)*c(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vpvtv')
      return
      end
C%5.1
      subroutine vpvts(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc,s)
      use lcd_mod
C
C     control loops
C     vector plus vector times scalar
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d,s
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vpvts')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(:n)+b(:n)*s
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vpvts')
      return
      end
C%5.1
      subroutine vpvpv(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector plus vector plus vector
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vpvpv')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= b(:n)+c(:n)+a(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vpvpv')
      return
      end
C%5.1
      subroutine vtvtv(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector times vector times vector
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vtvtv')
      call forttime(t1)
      do nl= 1,ntimes
	a(:n)= a(:n)*b(:n)*c(:n)
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,a)
      call check(chksum,ntimes*n,n,t2,'vtvtv')
      return
      end
C%5.1
      subroutine vsumr(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector sum reduction
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,sumvar
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vsumr')
      call forttime(t1)
      do nl= 1,ntimes
	sumvar= sum(a(:n))
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,sumvar)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= sumvar
      call check(chksum,ntimes*n,n,t2,'vsumr')
      return
      end
C%5.1
      subroutine vdotr(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     vector dot product reduction
C
      integer ntimes,ld,n,i,nl
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,dot
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vdotr')
      call forttime(t1)
      do nl= 1,ntimes
	dot= dot_product(a(:n),b(:n))
	call dummy(ld,n,a,b,c,d,e,aa,bb,cc,dot)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= dot
      call check(chksum,ntimes*n,n,t2,'vdotr')
      return
      end
C%5.1
      subroutine vbor(ntimes,ld,n,ctime,dtime,a,b,c,d,e,aa,bb,cc)
      use lcd_mod
C
C     control loops
C     basic operations rates, isolate arithmetic from memory traffic
C     all combinations of three, 59 flops for 6 loads and 1 store.
C	17 flops actual, with Fortran approved association
C
      integer ntimes,ld,n,i,nl,nn
      parameter(nn= 1000)
      real a(n),b(n),c(n),d(n),e(n),aa(ld,n),bb(ld,n),cc(ld,n)
      real t1,t2,chksum,ctime,dtime,cs1d
      real a1,b1,c1,d1,e1,f1,s(nn)
      call init(ld,n,a,b,c,d,e,aa,bb,cc,'vbor ')
      call forttime(t1)
      do nl= 1,ntimes
!dir$ vector aligned
	  do i= 1,n
	      a1= a(i)
	      b1= b(i)
	      c1= c(i)
	      d1= d(i)
	      e1= e(i)
	      f1= aa(i,1)
	      a1=a1*(((d1*(e1+f1)+e1*f1)+c1*((e1+f1)+d1))+
     &		b1*(((e1+f1)+d1)+c1))
	      b1= b1*((d1*(e1+f1)+e1*f1)+c1*((e1+f1)+d1))
	      c1= c1*(d1*(e1+f1)+e1*f1)
	      d1= (e1*f1)*d1
	      s(i)= a1*b1*c1*d1
	    enddo
	  call dummy(ld,n,a,b,c,d,e,aa,bb,cc,1.)
	enddo
      call forttime(t2)
      t2= t2-t1-ctime-(dtime*float(ntimes))
      chksum= cs1d(n,s)
      call check(chksum,ntimes*n,n,t2,'vbor ')
      return
      end
