/*****************************************************************************/
/*  fdhess_se.c -- part of fue.                                              */
/*                                                                           */
/*  Standard errors from fdhess AT the optimum (fue BUG-0015).               */
/*                                                                           */
/*  raxopt leaves the Hessian ACCUMULATED by BFGS along the path: good to    */
/*  steer the search, not the curvature at the optimum. It depends on the    */
/*  path (two runs of one model gave SE(mu) 0.073 and 0.028) and a search    */
/*  that starts at the optimum never builds it. fdhess is the call Mauricio  */
/*  left commented out in drvmlest.c; it is the default since the study in   */
/*  drvarma-python docs/STUDY-standard-errors.md (exact GLS within 0.35%).   */
/*  Guards, as in drvarma and drtran (atsw-gui):                             */
/*  - a neighbour the objective refuses means the optimum is on the          */
/*    boundary: no unrestricted Hessian exists there;                        */
/*  - choldcp is a MODIFIED Cholesky that patches pivots, so the Hessian is  */
/*    first checked with a plain one;                                        */
/*  - either way the BFGS factor is kept and est_se_how says why; unless     */
/*    raxopt did not iterate: it starts at the identity, so there is no      */
/*    BFGS Hessian either, and dev/cov are NAN.                              */
/*                                                                           */
/*  Copyright (C) 2026 D.E. Guerrero. GPL-2.0-or-later, as the rest of fue.  */
/*****************************************************************************/

#include <math.h>
#include "fue.h"
#include "nlatools.h"
#include "fdhess_se.h"

extern real macheps;
extern int  qn_last_nit;          /* qnewtopt.c: iterations of the last raxopt */

int est_fdhess = 1;
int est_se_how = EST_SE_BFGS;
static int objc_rejects = 0;

void fdhess( real (*func)(real *), int n, real *x, real f, real eta, real **H );

const char *est_se_label( int how )
{
   switch ( how )
      {
      case EST_SE_FDHESS:   return "fdhess";
      case EST_SE_BOUNDARY: return "bfgs (fdhess: the optimum is on the boundary "
                                   "of the admissible region)";
      case EST_SE_NOTPD:    return "bfgs (fdhess: the Hessian is not positive definite)";
      case EST_SE_NONE_BOUNDARY:
         return "none (fdhess: the optimum is on the boundary of the admissible "
                "region; the search did not move, so it built no BFGS Hessian)";
      case EST_SE_NONE_NOTPD:
         return "none (fdhess: the Hessian is not positive definite; the search "
                "did not move, so it built no BFGS Hessian)";
      default:              return "bfgs";
      }
}

real objc_reject( void )
{
   objc_rejects++;
   return 1.0;
}

real objc_finite( real f )
{
   return isfinite( f ) ? f : objc_reject();
}

int est_se_bfgs( void )
{
   return est_se_how == EST_SE_BFGS || est_se_how == EST_SE_BOUNDARY ||
          est_se_how == EST_SE_NOTPD;
}

/* 1 if H (k x k, symmetric) is positive definite by a plain Cholesky.       */
static int strict_pd( real **H, int k )
{
   int i, j, l, ok = 1;
   real s, **L = matrix( 1, k, 1, k );
   for ( j = 1; j <= k && ok; j++ )
       {
       s = H[j][j];
       for ( l = 1; l < j; l++ ) s -= L[j][l] * L[j][l];
       if ( !(s > 0.0) || !isfinite( s ) ) { ok = 0; break; }
       L[j][j] = sqrt( s );
       for ( i = j + 1; i <= k; i++ )
           {
           s = H[i][j];
           for ( l = 1; l < j; l++ ) s -= L[i][l] * L[j][l];
           L[i][j] = s / L[j][j];
           }
       }
   free_matrix( L, 1, k, 1, k );
   return ok;
}

void fdhess_cov( real (*objective)( real * ), int npar, real *par, real f,
                 real **cov, real *dev, int n )
{
   int  i, j, pfault;
   real d1, d2, *v, **H;

   v = vector( 1, npar );  H = matrix( 1, npar, 1, npar );
   objc_rejects = 0;
   fdhess( objective, npar, par, f, macheps, H );

   if ( objc_rejects > 0 )
      est_se_how = EST_SE_BOUNDARY;
   else if ( !strict_pd( H, npar ) )
      est_se_how = EST_SE_NOTPD;
   else
      {
      pfault = 0;
      choldcp( H, npar, &d1, &d2, &pfault );
      if ( pfault ) est_se_how = EST_SE_NOTPD;
      else
         {
         est_se_how = EST_SE_FDHESS;
         for ( i = 1; i <= npar; i++ )
             {
             for ( j = 1; j <= npar; j++ ) v[j] = 0.0;
             v[i] = 1.0;
             cholsol( H, npar, v );
             for ( j = 1; j <= npar; j++ )
                 cov[j][i] = ( 2.0 * f * v[j] ) / n;
             dev[i] = sqrt( cov[i][i] );
             }
         }
      }
   if ( ( est_se_how == EST_SE_BOUNDARY || est_se_how == EST_SE_NOTPD ) &&
        qn_last_nit == 0 )
      {
      est_se_how = ( est_se_how == EST_SE_BOUNDARY ) ? EST_SE_NONE_BOUNDARY
                                                     : EST_SE_NONE_NOTPD;
      for ( i = 1; i <= npar; i++ )
          { dev[i] = NAN; for ( j = 1; j <= npar; j++ ) cov[i][j] = NAN; }
      }
   free_matrix( H, 1, npar, 1, npar );  free_vector( v, 1, npar );
}
