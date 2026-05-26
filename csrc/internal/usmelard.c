
/*****************************************************************************/
/*  USMELARD.C                                                               */
/*  Computation of the log-likelihood function of a scalar ARMA(p,q) model.  */
/*  Source: Melard, G. (1984) ALGORITHM AS 197, 104-114.                     */
/*  Copyright (C) José Alberto Mauricio, 1996.                               */
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
/*****************************************************************************/
/*****************************************************************************/

void flikam( int m, int n, int mp, int mq, real *mu, real ***p, real ***q,
             real **qq, real **w, real sigma2, real toler, int chkma, int atf,
             real **at, real *sumsq, real *fact, real *loglik, int *ifault )

{
   const real LOG2PI = 1.837877066;

/* Declaration of variables and function twacf():                            */

   int  i, j, jfrom, k, last, loop, mpp1, mqp1, mrp1, mr, mxpq, mxpqp1, nexti;
   real detman, detcar, r, epsil1, vw1, a, aor, vl1, alf, flj;
   real *vw, *vl, *vk, *e;

   void twacf( real ***p, real ***q, int mp, int mq, real *acf, int ma,
               real *cvli, int mxpqp1, real *alpha, int mxpq, int *ifault );

/* [0]: Initialization:                                                      */

   detman  = 1.0;
   detcar  = 0.0;
   *fact   = 0.0;                                /* Determinant (see [7]).   */
   *sumsq  = 0.0;                                /* Quadratic form.          */
   *loglik = 0.0;                                /* Log-likelihood value.    */
   *ifault =   0;                                /* Error flag.              */
   epsil1  = 1.0e-10;

   mxpq    = ( mp > mq ) ? mp : mq;              /* Max(p,q).                */
   mxpqp1  = mxpq + 1;                           /* Max(p,q) + 1.            */
   mqp1    = mq + 1;                             /* q + 1.                   */
   mpp1    = mp + 1;                             /* p + 1.                   */
   mr      = ( mp > mqp1 ) ? mp : mqp1;          /* Max(p,q+1).              */
   mrp1    = mr + 1;                             /* Max(p,q+1) + 1.          */
   nexti   = n + 1;

/* [1]: Check for strict non-invertibility through chekma:                   */

   if ( (mq > 0) && (chkma == 1) )
      {
      vw = vector( 1, m * mq );
      vl = vector( 1, m * mq );
      vk = vector( 1, m * mq );
      chekma( m, mq, q, vw, vl, vk, ifault );
      free_vector( vk, 1, m * mq );
      free_vector( vl, 1, m * mq );
      free_vector( vw, 1, m * mq );
      }

   if ( *ifault > 0 )                  /* Strict non-invertibility detected: */
      {                                /* set flag and return.               */
      *ifault = 4;
      return;
      }

/* [2]: Calculation of the autocovariance function of a process with unit    */
/*      innovation variance (vw) and the covariances between the variable    */
/*      and the lagged innovations (vl):                                     */

   vw = vector( 1, mrp1 );
   vl = vector( 1, mrp1 );
   vk = vector( 1, mr );

   twacf( p, q, mp, mq, vw, mxpqp1, vl, mxpqp1, vk, mxpq, ifault );

   if ( *ifault > 0 )               /* Autoregressive unit root(s) detected: */
      {                             /* set flag and return.                  */
      *ifault = 2;
      free_vector( vk, 1, mr );
      free_vector( vl, 1, mrp1 );
      free_vector( vw, 1, mrp1 );
      return;
      }

/* [3]: Computation of the first column of matrix P (vk):                    */

   vk[1] = vw[1];

   if ( mr > 1 ) for ( k = 2; k <= mr; k++ )
      {
      vk[k] = 0.0;
      if ( k <= mp )
         for ( j = k; j <= mp; j++ ) vk[k] += p[j][1][1] * vw[j+2-k];
      if ( k <= mqp1 )
         for ( j = k; j <= mqp1; j++ ) vk[k] -= q[j-1][1][1] * vl[j+1-k];
      }

/* [4]: Computation of the initial vectors L and K (vl and Vk):              */

   r = vk[1];
   vl[mr] = 0.0;
   for ( j = 1; j <= mr; j++ )
       {
       vw[j] = 0.0;
       if ( j != mr ) vl[j] = vk[j+1];
       if ( j <= mp ) vl[j] += p[j][1][1] * r;
       vk[j] = vl[j];
       }

/* [5]: Initialization and loop on time (main computations):                 */

   last        = mpp1 - mq;
   loop        = mp;
   jfrom       = mpp1;
   vw[mpp1]    = 0.0;
   vl[mxpqp1]  = 0.0;

   e = vector( 1, n );      /* This will store the standardized innovations. */

   for ( i = 1; i <= n; i++ )
       {

    /* Test for skipped updating: */

       if ( i != last ) goto m1;
       loop = ( mp < mq ) ? mp : mq;
       jfrom = loop + 1;

    /* Test for switching: */

       if ( mq <= 0 )
          {
          nexti = i;
          break;
          }
m1:    if ( r <= epsil1 )               /* Strict non-stationarity detected: */
          {                             /* set flag and return.              */
          *ifault = 3;
          free_vector( e, 1, n );
          free_vector( vk, 1, mr );
          free_vector( vl, 1, mrp1 );
          free_vector( vw, 1, mrp1 );
          return;
          }
       if ( (fabs( r - 1.0) < toler) && (i > mxpq) )
          {
          nexti = i;
          break;
          }

    /* Updating scalars: */

       detman *= r;
       while ( fabs( detman ) >= 1.0 )
          {
          detman *= 0.0625;
          detcar += 4.0;
          }
       while ( fabs( detman ) < 0.0625 )
          {
          detman *= 16.0;
          detcar -=  4.0;
          }
       vw1      = vw[1];
       a        = (w[1][i] - mu[1]) - vw1;
       at[1][i] = a;                             /* Innovation (residual).   */
       e[i]     = a / sqrt( r );                 /* Standardized innovation. */
       aor      = a / r;
       *sumsq  += a * aor;
       vl1      = vl[1];
       alf      = vl1 / r;
       r       -= alf * vl1;

    /* Updating vectors: */

       if ( loop != 0 ) for ( j = 1; j <= loop; j++ )
          {
          flj    = vl[j+1] + p[j][1][1] * vl1;
          vw[j]  = vw[j+1] + p[j][1][1] * vw1 + aor * vk[j];
          vl[j]  = flj - alf * vk[j];
          vk[j] -= alf * flj;
          }
       if ( jfrom <= mq ) for ( j = jfrom; j <= mq; j++ )
          {
          vw[j]  = vw[j+1] + aor * vk[j];
          vl[j]  = vl[j+1] - alf * vk[j];
          vk[j] -= alf * vl[j+1];
          }
       if ( jfrom <= mp ) for ( j = jfrom; j <= mp; j++ )
          vw[j] = vw[j+1] + p[j][1][1] * (w[1][i] - mu[1]);
       }

/* [6]: Quick recursions calculations:                                       */

   if ( nexti <= n )
      {
//    *ifault = -nexti;
      for ( i = nexti; i <= n; i++ ) e[i] = (w[1][i] - mu[1]);
      if ( mp != 0 ) for ( i = nexti; i <= n; i++ )
         for ( j = 1; j <= mp; j++ ) e[i] -= p[j][1][1] * (w[1][i-j] - mu[1]);
      if ( mq != 0 ) for ( i = nexti; i <= n; i++ )
         for ( j = 1; j <= mq; j++ ) e[i] += q[j][1][1] * e[i-j];
      for ( i = nexti; i <= n; i++ )
          *sumsq += e[i] * e[i];
      for ( i = nexti; i <= n; i++ ) at[1][i] = e[i];
      }

/* [7]: Return determinant = detman * pow( 2.0, detcar ) and log-likelihood  */
/*      (fact is returned as the appropriate factor in objective function):  */

   *fact   = log( detman * pow( 2.0, detcar ) ) / n;
   *fact   = exp( *fact );
   *loglik = -0.5 * (n * LOG2PI + log( *fact ) + *sumsq );

/* [8]: Deallocate workspace and return:                                     */

   free_vector( e, 1, n );
   free_vector( vk, 1, mr );
   free_vector( vl, 1, mrp1 );
   free_vector( vw, 1, mrp1 );
}

/*****************************************************************************/
/*****************************************************************************/

void twacf( real ***p, real ***q, int mp, int mq, real *acf, int ma,
            real *cvli, int mxpqp1, real *alpha, int mxpq, int *ifault )

{
   int  i, j, k, kc, j1, miim1p, mikp;
   real divv, epsil2;

/* [0]: Initialization and return if mp = mq = 0:                            */

   *ifault = 0;
   epsil2  = 1.0e-10;
   acf[1]  = 1.0;
   cvli[1] = 1.0;

   if ( ma == 1 ) return;
   for ( i = 2; i <= ma; i++ ) acf[i] = 0.0;
   if ( mxpqp1 == 1 ) return;
   for ( i = 2; i <= mxpqp1; i++ ) cvli[i] = 0.0;
   for ( k = 1; k <= mxpq; k++ ) alpha[k] = 0.0;

/* [1]: Computation of the ACF of the moving average part (acf):             */

   if ( mq != 0 ) for ( k = 1; k <= mq; k++ )
      {
      cvli[k+1] = -q[k][1][1];
      acf[k+1]  = -q[k][1][1];
      if ( mq != k )
         for ( j = 1; j <= mq-k; j++ ) acf[k+1] += q[j][1][1] * q[j+k][1][1];
      acf[1] += q[k][1][1] * q[k][1][1];
      }

/* [2]: Initialization of cvli = T.W.-S PHI; return if mp = 0:               */

   if ( mp == 0 ) return;

   for ( k = 1; k <= mp; k++ )
       {
       alpha[k] = p[k][1][1];
       cvli[k]  = p[k][1][1];
       }

/* [3]: Computation of T.W.-S ALPHA and DELTA; DELTA stored in acf which is  */
/*      gradually overwritten:                                               */

   for ( k = 1; k <= mxpq; k++ )
       {
       kc = mxpq - k;
       if ( kc < mp )
          {
          divv = 1.0 - alpha[kc+1] * alpha[kc+1];
          if ( divv <= epsil2 )
             {
             *ifault = 1;
             return;
             }
          if ( kc != 0 )
             for ( j = 1; j <= kc; j++ )
                 alpha[j] = (cvli[j] + alpha[kc+1] * cvli[kc+1-j]) / divv;
          }
       if ( kc < mq )
          {
          j1 = ( kc + 1 - mp > 1 ) ? kc + 1 - mp : 1;
          for ( j = j1; j <= kc; j++ ) acf[j+1] += acf[kc+2] * alpha[kc+1-j];
          }
       if ( kc < mp )
          for ( j = 1; j <= kc; j++ ) cvli[j] = alpha[j];
       }

/* [4]: Computation of T.W.-S NU (store in cvli and copy into acf):          */

   acf[1] *= 0.5;
   for ( k = 1; k <= mxpq; k++ ) if ( k <= mp )
       {
       divv = 1.0 - alpha[k] * alpha[k];
       for ( j = 1; j <= k+1; j++ )
           cvli[j] = (acf[j] + alpha[k] * acf[k+2-j]) / divv;
       for ( j = 1; j <= k+1; j++ ) acf[j] = cvli[j];
       }

/* [5]: Computation of acf, which is gradually overwritten:                  */

   for ( i = 1; i <= ma; i++ )
       {
       miim1p = ( i - 1 < mp ) ? i - 1 : mp;
       if ( miim1p != 0 )
          for ( j = 1; j <= miim1p; j++ ) acf[i] += p[j][1][1] * acf[i-j];
       }
   acf[1] *= 2.0;

/* [6]: Computation of cvli (return when mq = 0):                            */

   cvli[1] = 1.0;
   if ( mq <= 0 ) return;
   for ( k = 1; k <= mq; k++ )
       {
       cvli[k+1] = -q[k][1][1];
       if ( mp != 0 )
          {
          mikp = ( k < mp ) ? k : mp;
          for ( j = 1; j <= mikp; j++ ) cvli[k+1] += p[j][1][1] * cvli[k+1-j];
          }
        }
}

/*****************************************************************************/
