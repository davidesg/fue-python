
/***************************************************************************
 *   fue.h                                                                 *
 *   Exact maximum likelihood estimation of time series models:            *
 *   Copyright (C) 2009 by A.B Treadway & D.E Guerrero & J.A. Mauricio     *
 *   abtreadway@telefonica.net &  warriord@rocketmail.com                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#define TRUE    1
#define FALSE   0
#define OK      1
#define WRONG   0
#define MAXSTR 90
#define FREE_ARG char*

/*****************************************************************************/

/*typedef char * STRING;                   old header                        */
typedef double real;

/*****************************************************************************/

struct Tvarma                    /* Standard multiple ARMA model structure:  */
    {
    int  m; int n; int p; int q;
    real *mu; real ***phi; real ***theta;
    real **qq; real **w; real **a;
    real sigma2; real logelf; real xitol; int chkma;
    };

/*****************************************************************************/

struct Tseries                   /* Standard single time series structure:   */
    {
    char *name;                  /* Time series name (string).               */
    int  nobs;                   /* Number of observations.                  */
    int  freq;                   /* Frequency (observations per year).       */
    int  numbering;              /* only numbering option (binary)           */
    int  begtime;                /* Beginning period.                        */
    int  begyear;                /* Beginning year.                          */
    int  endtime;                /* Ending period.                           */
    int  endyear;                /* Ending year.                             */
    real mean;                   /* Sample mean.                             */
    real var;                    /* Sample variance.                         */
    real skew;                   /* Skewness coefficient.                    */
    real kurt;                   /* Kurtosis coefficient.                    */
    real jarquebera;             /* Jarque-Bera statistical                  */ 
    int  max;                    /* Maximum value.                           */
    int  min;                    /* Minimum value.                           */
    real refactor;                /* reescaling factor                        */
    real *data;                  /* Vector of time series observations.      */
    };

/*****************************************************************************/

struct Tusmodel          /* Seasonal US model with deterministic components: */
    {
    char *name;          /* Model name (string).                             */
    char *residuals;     /* Model residuals name                             */
    real boxlam;         /* Box-Cox parameter lambda: either 0.0 or 1.0.     */
    int  sper;           /* Seasonal period: either 1(A), 4(Q) or 12(M).     */
    int  nrdiff;         /* Number of regular differences.                   */
    int  nadiff;         /* Number of (complete) annual differences.         */
    int  *ifadf;         /* Individual factors of the annual difference.     */
    int  ornsop;         /* Order of the resulting non-stationary operator.  */
    real *rnsop;         /* Resulting non-stationary operator.               */
    real cbands;
    real mu;             /* Mean parameter.                                  */
    int  Imu;            /* Flag (0-1): estimation of mu.                    */
    real sigma2;         /* Estimated residual variance.                     */

 /* Section [1]: Deterministic variables (intervention + seasonal):          */

    int  NdetVar;        /* Number of deterministic variables (detvars)      */
    int  *Nomega;        /* Order of omega(B) for each detvar (order s).     */
    real **Omega;        /* Matrix of omegas: one row for each detvar.       */
    int  **Imega;        /* Matrix of flags (0-1): omegas to be estimated.   */
    int  *Ndelta;        /* Order of delta(B) for each detvar (order r).     */
    real **Delta;        /* Matrix of deltas: one row for each detvar.       */
    int  **Ielta;        /* Matrix of flags (0-1): deltas to be estimated.   */

 /* Section [2]: Standard regular-annual AR and MA factors:                  */

    int  NumAr1;         /* Number of regular AR factors.                    */
    int  *p1;            /* Number of phis for each factor (AR order).       */
    real **Ar1;          /* Matrix of phis: one row for each factor.         */
    int  **Ia1;          /* Matrix of flags (0-1): phis to be estimated.     */

    int  NumAr2;         /* Number of annual AR factors.                     */
    int  *p2;            /* Number of PHIS for each factor (AR order).       */
    real **Ar2;          /* Matrix of PHIS: one row for each factor.         */
    int  **Ia2;          /* Matrix of flags (0-1): PHIS to be estimated.     */

    int  NumMa1;         /* Number of regular MA factors.                    */
    int  *q1;            /* Number of thetas for each factor (MA order).     */
    real **Ma1;          /* Matrix of thetas: one row for each factor.       */
    int  **Im1;          /* Matrix of flags (0-1): thetas to be estimated.   */

    int  NumMa2;         /* Number of annual MA factors.                     */
    int  *q2;            /* Number of THETAS for each factor (MA order).     */
    real **Ma2;          /* Matrix of THETAS: one row for each factor.       */
    int  **Im2;          /* Matrix of flags (0-1): THETAS to be estimated.   */

 /* Section [3]: AR and MA factors of order 2 with fixed frequency:          */

    int  NumAr1f;        /* Number of regular AR factors (fixed frequency).  */
    real *pfre1;         /* Frequency for each factor.                       */
    real **Ar1f;         /* Matrix of phis: one row for each factor.         */
    int  *Ia1f;          /* Vector of flags (0-1): phis to be estimated.     */

    int  NumAr2f;        /* Number of annual AR factors (fixed frequency).   */
    real *pfre2;         /* Frequency for each factor.                       */
    real **Ar2f;         /* Matrix of PHIS: one row for each factor.         */
    int  *Ia2f;          /* Vector of flags (0-1): PHIS to be estimated.     */

    int  NumMa1f;        /* Number of regular MA factors (fixed frequency).  */
    real *qfre1;         /* Frequency for each factor.                       */
    real **Ma1f;         /* Matrix of thetas: one row for each factor.       */
    int  *Im1f;          /* Vector of flags (0-1): thetas to be estimated.   */

    int  NumMa2f;        /* Number of annual MA factors (fixed frequency).   */
    real *qfre2;         /* Frequency for each factor.                       */
    real **Ma2f;         /* Matrix of THETAS: one row for each factor.       */
    int  *Im2f;          /* Vector of flags (0-1): THETAS to be estimated.   */


/* DG 06/21/04 ***************************************************************/

 /* Section [4]: AR and MA seasonal factors of order S-1 :                   */
/*
    int  NumAr1S;        /* Number of seasonal AR factors.                   */
/*    real *pSre1;         /* Frequency for each factor.                       */
/*    real **Ar1S;         /* Matrix of phis: one row for each factor.         */
/*    int  *Ia1S;          /* Vector of flags (0-1): thetas to be estimated.   */

/*    int  NumMa1S;        /* Number of seasonal MA factors.                   */
/*    real *qSre1;         /* Frequency for each factor.                       */
/*    real **Ma1S;         /* Matrix of thetas: one row for each factor.       */
/*    int  *Im1S;          /* Vector of flags (0-1): thetas to be estimated.   */

/* END DG ********************************************************************/
    };
/*****************************************************************************/

void est( void (*cast)( real *, struct Tvarma *, int *, int, int ),
          int npar, real *par, real *dev, real **cov, int maxits, int nrits,
          real grtol, real sptol, real xitol, int chkma,
          real **a, real *sigma2, real *logelf, int *ifault );

void elf( int m, int n, int p, int q, real *mu, real ***phi, real ***theta,
          real **qq, real **w, real sigma2, real delta, int chkma, int atf,
          real **a, real *f1, real *f2, real *logelf, int *ifault );

void flikam( int m, int n, int mp, int mq, real *mu, real ***p, real ***q,
             real **qq, real **w, real sigma2, real toler, int chkma, int atf,
             real **at, real *sumsq, real *fact, real *loglik, int *ifault );

void marma( int k, int n, int p, int q, real *mu, real ***phi, real ***theta,
            real **qq, real **w, real sigma2, real xtol, int chkma, int atf,
            real **v, real *r1, real *r2, real *rlogl, int *ifault );

 void chekma( int m, int q, real ***theta, real *wr, real *wi,
             real *wmod, int *ifault );

/*void gsl_chekma( int m, int q, real ***theta, real *wr, real *wi,
             real *wmod, int *ifault );	 

/*****************************************************************************/


/*****************************************************************************/

void ObsToDate( int beg_per, int beg_sub, int obs_no, int freq,
                int *per, int *sub );
void DateToObs( int beg_per, int beg_sub, int per, int sub, int freq,
                int *obs_no );
real Mean( real *data, int nobs );
real Stdev( real *data, int nobs );
int  MaxVal( real *data, int nobs );
int  MinVal( real *data, int nobs );
real Skew( real *data, int nobs );
real Kurt( real *data, int nobs );
real JarqueBera( real Skew, real Kurt, int nobs);
void Acf( struct Tseries *ser, int lags, real *corr );
void Acf_dist( struct Tseries *ser, int lags );
void Acf_disttex( struct Tseries *ser, char *distputf );

void Pacf( int lags, real *pcorr );
real ChiTest( real *corr, int lags, int nobs );
void Pacf_New (double *ser , int nobs, double *pcorr, int lags);


void File_StatSer( struct Tseries *ser );
void File_StatSer_X11( struct Tseries *ser, int npar, int tsnobs, int tmornsop, int tsby, double cbands, char *outx11, double boxlam, double refactor );
void File_PlotSer( struct Tseries *ser );
void File_HistSer( struct Tseries *ser );
void File_CorrSer( struct Tseries *ser, int npar );

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double Acf_Pacf_Max ( struct Tseries *ser,  int lags );
void BoxCox  ( double *DataInput, double *DataOutput, double boxlam, double boxm, int nobs, double refactor, int geometric);

