
/*****************************************************************************/
/*  DRVMLEST.C                                                               */
/*  Driver module for exact maximum likelihood estimation.                   */
/*  Copyright (C) José Alberto Mauricio, 1995.                               */
/*                                                                           */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by     */
/*  the Free Software Foundation; either version 2 of the License, or       */
/*  (at your option) any later version.                                      */
/*                                                                           */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             */
/*  GNU General Public License for more details.                             */
/*                                                                           */
/*  You should have received a copy of the GNU General Public License        */
/*  along with this program; if not, write to the Free Software Foundation,  */
/*  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.     */
/*****************************************************************************/

#include "fue.h"            /* Header file (prototype declarations)          */
#include "nlatools.h"            /* Header file (prototype declarations)     */
extern real macheps;          /* Machine epsilon (global: declared in DRV.C) */
extern FILE *outputv;         /* Output file (global: declared in DRV.C)     */

/*****************************************************************************/

real pi10x, pi20x;
struct Tvarma varmax;
void (*castx)( real *, struct Tvarma *,int *, int, int );

/*****************************************************************************/
/*****************************************************************************/

void est( void (*cast)( real *, struct Tvarma *, int *, int, int ),
          int npar, real *par, real *dev, real **cov, int maxits, int nrits,
          real grtol, real sptol, real xitol, int chkma,
          real **a, real *sigma2, real *logelf, int *ifault )

/*****************************************************************************/
/*                                                                           */
/* cast ->   : casting routine [drv.c].                                      */
/* npar ->   : n� of parameters to estimate.                                 */
/* par <->   : -> initial estimate (guess); <- final estimate (solution).    */
/* dev <-    : estimated standard deviations at final estimate.              */
/* cov <-    : estimated covariance matrix at final estimate.                */
/* maxits -> : maximum number of iterations.                                 */
/* nrits ->  : report optimizer progress every [nrits] iterations.           */
/* grtol ->  : stopping tolerance 1: relative gradient value.                */
/* sptol ->  : stopping tolerance 2: relative parameter change.              */
/* xitol ->  : tolerance for changing to quick recursions when computing lf. */
/* chkma ->  : flag for invertibility checking.                              */
/* a <-      : residual vector at final estimate.                            */
/* sigma2 <- : final estimate of "concentrated" residual variance.           */
/* logelf <- : exact log-likelihood at final estimate.                       */
/* ifault <- : [0]: successful return;                                       */
/*             [1]: bad initial estimate: Q is not positive definite;        */
/*             [2]: bad initial estimate: AR has at least one unit root;     */
/*             [3]: bad initial estimate: AR is strictly non-stationary;     */
/*             [4]: bad initial estimate: MA is strictly non-invertible;     */
/*             [5]: bad initial estimate: unknown numerical problems (rare); */
/*             [6]: bad initial estimate: see cast().                        */
/*                                                                           */
/*****************************************************************************/

{
   const real LOG2PI = 1.837877066;

/* Declaration of local variables:                                           */

   int  i, j;
   real pi1, pi2, pi3, **mtmp, *vtmp;

/* Declaration of local functions:                                           */

   void raxopt( real (*func)( real [] ), real *fk, int n, real *xk, real **b,
                int maxits, int nrits, real gradtol, real steptol );
   void fdhess( real (*func)(real *), int n, real *x, real f,
                real eta, real **H );
   real objcfunc( real *x );


/* [1]: First computation of the log-likelihood: initialize pi10x & pi20x:   */

   *ifault = 0;                               /* Initialize fault indicator. */

   varmax.xitol = xitol;                      /* Estimation method.          */
   varmax.chkma = chkma;                      /* Check for invertibility.    */

   (*cast)( par, &varmax, ifault, 1, 0 );     /* Allocate VARMA structure.   */

   if ( *ifault > 0 ) return;                 /* Bad initial estimates (1).  */

   flikam( varmax.m, varmax.n, varmax.p, varmax.q, varmax.mu, varmax.phi,
           varmax.theta, varmax.qq, varmax.w, 1.0, varmax.xitol, varmax.chkma,
           TRUE, varmax.a, &pi10x, &pi20x, &pi3, ifault );

   if ( *ifault > 0 ) return;                 /* Bad initial estimates (2):  */
                                              /* ifault = 1-2-3-4-5.         */

/* [2]: Maximization of the specified objective function:                    */

   mtmp  = matrix( 1, npar, 1, npar );
   vtmp  = vector( 1, npar );
   castx = cast;                              /* Assign casting routine.     */

   raxopt( objcfunc, &pi1, npar, par, mtmp, maxits, nrits, grtol, sptol );

/* This is an alternative way of computing the second derivative matrix:     */

/* fdhess( objcfunc, npar, par, pi1, macheps, mtmp );                        */
/* choldcp( mtmp, npar, &pi2, &pi3, ifault );                                */

/* [3]: Sample estimation of the variance-covariance matrix:                 */

   for ( i = 1; i <= npar; i++ )
       {
       for ( j = 1; j <= npar; j++ ) vtmp[j] = 0.0;
       vtmp[i] = 1.0;
       cholsol( mtmp, npar, vtmp );
       for ( j = 1; j <= npar; j++ )
           cov[j][i] = ( 2.0 * pi1 * vtmp[j] ) / varmax.n;
       dev[i] = sqrt( cov[i][i] );
       }

   free_vector( vtmp, 1, npar );
   free_matrix( mtmp, 1, npar, 1, npar );


/* [4]: Last computation of the log-likelihood: residuals and sigma2 (note   */
/*      that this last calculation is performed through elf and not          */
/*      through flikam):                                                     */

   (*cast)( par, &varmax, ifault, 0, 0 );

   elf( varmax.m, varmax.n, varmax.p, varmax.q, varmax.mu, varmax.phi,
        varmax.theta, varmax.qq, varmax.w, 1.0, varmax.xitol, varmax.chkma,
        TRUE, a, &pi1, &pi2, &pi3, ifault );

   *logelf = -0.5 * varmax.m * varmax.n * ( LOG2PI - log( varmax.m )
             - log( varmax.n ) + 1.0 )
             - 0.5 * varmax.n * ( varmax.m * log( pi1 ) + log( pi2 ) );
   *sigma2 = pi1 / (varmax.n * varmax.m);

   (*cast)( par, &varmax, ifault, 0, 1 );     /* Deallocate VARMA structure. */
}

/*****************************************************************************/

real objcfunc( real *x )

{
   int  ifault;
   real pi1, pi2, pi3;

/* [1]: Put parameter vector into VARMA structure [since (*castx) receives   */
/*      ifault as a formal parameter, it is possible to check for model      */
/*      adequacy within (*castx), prior to the standard tests carried out    */
/*      within subroutine flikam]:                                           */

   ifault = 0;
   (*castx)( x, &varmax, &ifault, 0, 0 );
   if ( ifault > 0 )                            /* ifault = 6-7-8- ...       */
      return( 1.0 );

/* [2]: Compute objective function and return:                               */

   flikam( varmax.m, varmax.n, varmax.p, varmax.q, varmax.mu, varmax.phi,
           varmax.theta, varmax.qq, varmax.w, 1.0, varmax.xitol, varmax.chkma,
           FALSE, varmax.a, &pi1, &pi2, &pi3, &ifault );

   if ( ifault > 0 )                            /* ifault = 1-2-3-4-5.       */
      return( 1.0 );
   else
      return( pow( (pi1 / pi10x), varmax.m ) * (pi2 / pi20x) );
}

/*****************************************************************************/
