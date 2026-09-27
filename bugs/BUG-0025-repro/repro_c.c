/* BUG-0025, the C side: lnsrch of csrc/internal/qnewtopt.c with a NaN
   objective. Build and run (from the repo root):
     gcc -Icsrc/internal bugs/BUG-0025-repro/repro_c.c csrc/internal/qnewtopt.c \
         csrc/internal/nlatools.c $(pkg-config --cflags --libs gsl) -lm -o /tmp/r
     timeout 10 /tmp/r; echo exit=$?          (124 = killed by the timeout) */
#include <stdio.h>
#include <math.h>
#include "fue.h"
#include "nlatools.h"

FILE *outputv = NULL;
real  macheps;

real lnsrch( int n, real *xk, real fk, real *gk, real *dk, real *xkp1,
             real *fkp1, real maxstep, real steptol, int *retcode,
             int *maxtaken, real (*func)(real *) );

static real f( real *x )                 /* NaN beyond x = 1 */
{
   if ( x[1] > 1.0 ) return NAN;
   return ( x[1] - 3.0 ) * ( x[1] - 3.0 );
}

int main( void )
{
   real xk[2] = { 0, 0.0 }, gk[2] = { 0, -6.0 }, dk[2] = { 0, 6.0 };
   real xkp1[2], fkp1, lambda;
   int  retcode, maxtaken;
   macheps = 2.2e-16;
   lambda = lnsrch( 1, xk, 9.0, gk, dk, xkp1, &fkp1, 100.0, 1e-7,
                    &retcode, &maxtaken, f );
   printf( "returned: lambda %g retcode %d x %g f %g\n", lambda, retcode,
           xkp1[1], fkp1 );
   return 0;
}
