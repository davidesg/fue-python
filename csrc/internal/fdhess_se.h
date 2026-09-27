/*****************************************************************************/
/*  fdhess_se.h -- part of fue.                                              */
/*                                                                           */
/*  Standard errors from Mauricio's fdhess AT the optimum (fue BUG-0015).    */
/*  Kept out of drvmlest.c, which stays Mauricio's source with a handful of  */
/*  lines marked BUG-0015 (docs/PROVENANCE.md, section 2).                   */
/*                                                                           */
/*  Copyright (C) 2026 D.E. Guerrero. GPL-2.0-or-later, as the rest of fue.  */
/*****************************************************************************/

#ifndef FDHESS_SE_H
#define FDHESS_SE_H

/* Include after fue.h (which has no include guard): it needs `real`.     */

/* Which Hessian gave the standard errors (est_se_how):                      */
#define EST_SE_BFGS          0  /* the one BFGS accumulated (fdhess not asked)*/
#define EST_SE_FDHESS        1  /* fdhess at the optimum                      */
#define EST_SE_BOUNDARY      2  /* BFGS kept: the optimum is on the boundary  */
#define EST_SE_NOTPD         3  /* BFGS kept: fdhess not positive definite    */
#define EST_SE_NONE_BOUNDARY 4  /* neither: boundary, and BFGS never built    */
#define EST_SE_NONE_NOTPD    5  /* neither: not PD, and BFGS never built      */

extern int est_fdhess;          /* 1: fdhess (default); 0: the BFGS matrix    */
extern int est_se_how;          /* <- set by est()                           */

const char *est_se_label( int how );

/* cov = 2 F H^-1 / n from fdhess at par, guarded; sets est_se_how.         */
void fdhess_cov( real (*objective)( real * ), int npar, real *par, real f,
                 real **cov, real *dev, int n );

/* 1 when est() must compute cov/dev from the BFGS factor.                  */
int  est_se_bfgs( void );

/* objcfunc's inadmissible point (1.0), counted: that is the boundary test. */
real objc_reject( void );
/* A non-finite objective is inadmissible too (as in drvarma and drtran).   */
real objc_finite( real f );

#endif
