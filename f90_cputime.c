/* SET CLOCK_RATE AND RECOMPILE FOR TARGET MACHINE */
#undef _WIN32
/* set not to use API calls even on Windows */
	#ifdef _WIN32
	#include <windows.h>
	#endif
unsigned long long int rdtsc( ) 
{
#ifdef _M_IA64
  
unsigned __int64 __getReg(int whichReg);
#pragma intrinsic(__getReg);
#define INL_REGID_APITC 3116
  
  return  __getReg(INL_REGID_APITC);
#elif defined(_WIN32)
	unsigned long long int qpc;
	(void)QueryPerformanceCounter((LARGE_INTEGER *)&qpc);
	return qpc;
#elif defined(__GNUC__)
#if defined i386
   long long a;
   asm volatile("rdtsc":"=A" (a));
   return a;
#elif defined __x86_64
   unsigned int _hi,_lo;
   asm volatile("rdtsc":"=a"(_lo),"=d"(_hi));
   return ((unsigned long long int)_hi << 32) | _lo;
#else
	unsigned long result;
/* gcc-IA64 version */
	__asm__ __volatile__("mov %0=ar.itc" : "=r"(result) :: "memory");
	while (__builtin_expect ((int) result == -1, 0))
		__asm__ __volatile__("mov %0=ar.itc" : "=r"(result) ::
"memory");
	return result;

#endif
#elif defined(_WIN64)
  _asm 
  {
   xor	rax,rax
   xor	eax,eax
   _emit 0x0f /* rdtsc */
   _emit 0x31
   shl	rdx,32
   or	rax,rdx
  }
#elif defined(_M_IX86)
  _asm 
  {
   _emit 0x0f /* rdtsc */
   _emit 0x31
  }
return; 
#else
#error "only supports IA64,IX86,GNUC"
#endif
}
 
double G77_etime_0 (float tarray[2])
{
   static int win32_platform = -1;
   double usertime, systime;
 
     {
       static double clock_per=1./(long long)CLOCK_RATE;
       static unsigned long long int old_count;
       unsigned long long count;
       if(!old_count){
	#ifdef _WIN32
	unsigned long long int qpf;
	if(QueryPerformanceFrequency((LARGE_INTEGER *)&qpf))
	    clock_per=1./(long long)qpf;
	#endif
	old_count=rdtsc();
	}
 
       count = rdtsc();
       tarray[0] = usertime = (long long)(count - old_count) * clock_per;
       tarray[1] = 0;
     }
   return usertime ;
 
}
void f90_cputime4(float *time){	// Intel Fortran 7 call
 	float tarray[2];
 	*time=(float)G77_etime_0 (tarray);
}
void forttime_(float *time){	// Intel Fortran 8
 	float tarray[2];
 	*time=(float)G77_etime_0 (tarray);
}
