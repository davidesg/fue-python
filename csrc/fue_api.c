/*
 * fue_api.c  —  Public FUE estimation API.
 *
 * Bridges fue_api.h (cffi-visible surface) and the internal estimation engine.
 * The mathematical core (elfvarma, usmelard, qnewtopt, nlatools, drvmlest) is
 * unchanged. This file contains:
 *
 *   populate_globals()  — FueModelSpec  →  Tm + Ts + DataMat
 *   count_npar_build_par() — count free params, fill initial par[]
 *   cast_us()           — extracted from fue.c:3645 (uses module-level globals)
 *   unscramble()        — extracted from fue.c:3939
 *   CalcNonsOp()        — extracted from fue.c:4303
 *   calcnu()            — extracted from fue.c:4489
 *   Cos2()              — extracted from fue.c:4290
 *   BoxCox_simple()     — local version (no geometric Jacobian)
 *   fue_estimate()      — main API entry point
 *
 * Copyright (C) 2009-2026 A.B Treadway, D.E. Guerrero & J.A. Mauricio
 * License: GPL v2 or later.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "fue_api.h"
#include "internal/fue.h"
#include "internal/nlatools.h"

/* ── Globals required by the internal estimation engine ─────────────────── */
/* cast_us() accesses these as module-level globals, matching fue.c layout.  */

struct Tusmodel Tm;
struct Tseries  Ts;
double        **DataMat;

real    macheps = 2.22e-16;
FILE   *outputv = NULL;       /* engine diagnostics suppressed in API mode  */

/* ── fue_strerror ───────────────────────────────────────────────────────── */

const char *fue_strerror(int ifault)
{
    switch (ifault) {
    case  0: return "OK";
    case  1: return "bad initial estimate: Q not positive definite";
    case  2: return "bad initial estimate: AR has a unit root";
    case  3: return "bad initial estimate: AR is non-stationary";
    case  4: return "bad initial estimate: MA is non-invertible";
    case  5: return "bad initial estimate: unknown numerical problem";
    case  6: return "bad initial estimate: rejected by cast routine";
    default: return "unknown error";
    }
}

/* ── fue_defaults ───────────────────────────────────────────────────────── */

void fue_defaults(FueModelSpec *spec)
{
    memset(spec, 0, sizeof(*spec));
    spec->boxlam      = 1.0;
    spec->refactor    = 1.0;
    spec->estimate_mu = 0;
    spec->maxits      = 500;
    spec->grtol       = pow(DBL_EPSILON, 1.1 / 3.0);
    spec->sptol       = pow(DBL_EPSILON, 2.0 / 3.0);
    spec->xitol       = 1e-3;
    spec->chkma       = 1;
    spec->eml         = 1;
}

/* ── Result allocation ──────────────────────────────────────────────────── */

static FueResult *alloc_result(int npar, int nresiduals)
{
    FueResult *r = calloc(1, sizeof(*r));
    if (!r) return NULL;
    r->npar       = npar;
    r->nresiduals = nresiduals;
    if (npar > 0) {
        r->params     = calloc((size_t)npar,       sizeof(double));
        r->std_errors = calloc((size_t)npar,       sizeof(double));
        r->cov_matrix = calloc((size_t)npar * npar, sizeof(double));
    }
    if (nresiduals > 0)
        r->residuals  = calloc((size_t)nresiduals, sizeof(double));
    return r;
}

void fue_result_free(FueResult *r)
{
    if (!r) return;
    free(r->params);
    free(r->std_errors);
    free(r->cov_matrix);
    free(r->residuals);
    free(r);
}

/* ── Helper: Cos2 ───────────────────────────────────────────────────────── */

static double Cos2(double freq, int sper)
{
    const double PI = 3.141592654;
    return cos(2.0 * PI * freq / sper);
}

/* ── calcnu: impulse response of ω(B)/δ(B) ─────────────────────────────── */
/* Extracted from fue.c:4489                                                 */

static void calcnu(double *omega, int s, double *delta, int r,
                   double *nu, int lags)
{
    int i, j;
    double sum1, sum2;

    nu[0] = omega[0];
    for (j = 1; j <= lags; j++) {
        sum1 = 0.0;
        if (r > 0)
            for (i = 1; i <= j; i++)
                if (i <= r) sum1 += delta[i] * nu[j - i];
        sum2 = 0.0;
        if (s > 0)
            if (j <= s) sum2 = omega[j];
        nu[j] = sum1 - sum2;
    }
}

/* ── CalcNonsOp: non-stationary differencing operator ──────────────────── */
/* Extracted from fue.c:4303                                                 */

static void CalcNonsOp(int sp, int d, int ds, int *ifds, int ord, double *op)
{
    double *pol1, *pol2, *pol3, *pol4;
    int i, j, k, pp, pp1, pp2;

    pp1  = d + ds * sp;
    pol1 = vector(0, pp1);
    pol2 = vector(0, pp1);
    for (i = 1; i <= pp1; i++) pol1[i] = 0.0;
    pol1[0] = -1.0;
    pp      = 0;

    if (ds > 0)
        for (k = 1; k <= ds; k++) {
            for (i = 0; i <= pp1; i++) pol2[i] = 0.0;
            for (j = 0; j <= pp + sp; j++)
                if      (j >= 0  && j <  sp) pol2[j] =  pol1[j];
                else if (j >= sp && j <= pp) pol2[j] =  pol1[j] - pol1[j - sp];
                else if (j >  pp && j <= pp + sp) pol2[j] = -pol1[j - sp];
            pp += sp;
            for (i = 0; i <= pp; i++) pol1[i] = pol2[i];
        }
    if (d > 0)
        for (k = 1; k <= d; k++) {
            for (i = 0; i <= pp1; i++) pol2[i] = 0.0;
            for (j = 0; j <= pp + 1; j++)
                if      (j >= 0 && j < 1)  pol2[j] =  pol1[j];
                else if (j >= 1 && j <= pp) pol2[j] =  pol1[j] - pol1[j - 1];
                else if (j >  pp && j <= pp + 1) pol2[j] = -pol1[j - 1];
            pp += 1;
            for (i = 0; i <= pp; i++) pol1[i] = pol2[i];
        }

    free_vector(pol2, 0, pp1);

    pp2  = ord - pp1;
    pol2 = vector(0, pp2 > 0 ? pp2 : 1);
    pol3 = vector(0, pp2 > 0 ? pp2 : 1);
    pol4 = vector(0, 2);
    for (i = 1; i <= pp2; i++) pol2[i] = 0.0;
    pol3[0] = -1.0;
    pp      =   0;

    if (ifds && pp2 > 0) {
        #define APPLY_FACTOR(deg) \
            do { for (i=1;i<=pp2;i++) pol2[i]=0.0; \
                 pol2[0]=pol3[0]=-1.0; \
                 for (i=0;i<=pp;i++) for (j=0;j<=(deg);j++) pol2[j+i]-=pol4[j]*pol3[i]; \
                 pp+=(deg); for (i=1;i<=pp;i++) pol3[i]=pol2[i]; } while(0)

        if (((sp == 12) && ifds[0]) || ((sp == 4) && ifds[0])) {
            pol4[0]=-1.0; pol4[1]=1.0;
            APPLY_FACTOR(1);
        }
        if ((sp == 12) && ifds[1]) {
            pol4[0]=-1.0; pol4[1]=sqrt(3.0); pol4[2]=-1.0;
            APPLY_FACTOR(2);
        }
        if ((sp == 12) && ifds[2]) {
            pol4[0]=-1.0; pol4[1]=1.0; pol4[2]=-1.0;
            APPLY_FACTOR(2);
        }
        if (((sp == 12) && ifds[3]) || ((sp == 4) && ifds[1])) {
            pol4[0]=-1.0; pol4[1]=0.0; pol4[2]=-1.0;
            APPLY_FACTOR(2);
        }
        if ((sp == 12) && ifds[4]) {
            pol4[0]=-1.0; pol4[1]=-1.0; pol4[2]=-1.0;
            APPLY_FACTOR(2);
        }
        if ((sp == 12) && ifds[5]) {
            pol4[0]=-1.0; pol4[1]=-sqrt(3.0); pol4[2]=-1.0;
            APPLY_FACTOR(2);
        }
        if (((sp == 12) && ifds[6]) || ((sp == 4) && ifds[2])) {
            pol4[0]=-1.0; pol4[1]=-1.0;
            APPLY_FACTOR(1);
        }
        #undef APPLY_FACTOR
    }
    for (i = 0; i <= (pp2 > 0 ? pp2 : 0); i++) pol2[i] = pol3[i];

    free_vector(pol4, 0, 2);
    free_vector(pol3, 0, pp2 > 0 ? pp2 : 1);

    for (i = 1; i <= ord; i++) op[i] = 0.0;
    op[0] = -1.0;
    for (i = 0; i <= pp1; i++)
        for (j = 0; j <= (pp2 > 0 ? pp2 : 0); j++)
            op[j + i] -= pol2[j] * pol1[i];

    free_vector(pol2, 0, pp2 > 0 ? pp2 : 1);
    free_vector(pol1, 0, pp1);
}

/* ── unscramble: expand factored AR/MA polynomials to single operator ───── */
/* Extracted from fue.c:3939 (inactive sections omitted).                    */

static void unscramble(struct Tusmodel *Model, int *OrderAr, int *OrderMa,
                       double *ArFactor, double *MaFactor, int *ifault)
{
    int i, j, k, pr, pa, qr, qa, p1, p2, q1, q2, itmp1, itmp2, pq;
    double *phir, *phia, *thetar, *thetaa, *tmp;

    pr = pa = qr = qa = pq = 0;
    p1 = p2 = q1 = q2 = 0;

    for (i = 1; i <= Model->NumAr1;  i++) pr += Model->p1[i];
    for (i = 1; i <= Model->NumAr1f; i++) pr += 2;
    for (i = 1; i <= Model->NumAr2;  i++) pa += Model->p2[i];
    for (i = 1; i <= Model->NumMa1;  i++) qr += Model->q1[i];
    for (i = 1; i <= Model->NumMa1f; i++) qr += 2;
    for (i = 1; i <= Model->NumMa2;  i++) qa += Model->q2[i];

    itmp1 = (pr > pa) ? pr : pa;
    itmp2 = (qr > qa) ? qr : qa;
    pq    = (itmp1 > itmp2) ? itmp1 : itmp2;

    phir   = vector(0, pr > 0 ? pr : 1);
    phia   = vector(0, pa > 0 ? pa : 1);
    thetar = vector(0, qr > 0 ? qr : 1);
    thetaa = vector(0, qa > 0 ? qa : 1);
    tmp    = vector(0, pq > 0 ? pq : 1);

    phir[0] = phia[0] = thetar[0] = thetaa[0] = tmp[0] = -1.0;

    /* [2]: Regular AR operator */
    for (k = 1; k <= Model->NumAr1; k++) {
        for (i = 1; i <= pr; i++) phir[i] = 0.0;
        phir[0] = tmp[0] = -1.0;
        for (i = 0; i <= p1; i++)
            for (j = 0; j <= Model->p1[k]; j++)
                phir[j + i] -= Model->Ar1[k][j] * tmp[i];
        p1 += Model->p1[k];
        for (i = 1; i <= p1; i++) tmp[i] = phir[i];
    }
    for (i = 0; i <= p1; i++) phir[i] = tmp[i];

    /* [3]: Regular AR(2) with fixed frequency */
    for (k = 1; k <= Model->NumAr1f; k++) {
        if (Model->Ar1f[k][2] > 0.0) { *ifault = 1; goto u1; }
        Model->Ar1f[k][1] = 2.0 * Cos2(Model->pfre1[k], Model->sper)
                            * sqrt(-Model->Ar1f[k][2]);
        for (i = 1; i <= pr; i++) phir[i] = 0.0;
        phir[0] = tmp[0] = -1.0;
        for (i = 0; i <= p1; i++)
            for (j = 0; j <= 2; j++)
                phir[j + i] -= Model->Ar1f[k][j] * tmp[i];
        p1 += 2;
        for (i = 1; i <= p1; i++) tmp[i] = phir[i];
    }
    for (i = 0; i <= p1; i++) phir[i] = tmp[i];

    /* [4]: Annual AR operator */
    for (k = 1; k <= Model->NumAr2; k++) {
        for (i = 1; i <= pa; i++) phia[i] = 0.0;
        phia[0] = tmp[0] = -1.0;
        for (i = 0; i <= p2; i++)
            for (j = 0; j <= Model->p2[k]; j++)
                phia[j + i] -= Model->Ar2[k][j] * tmp[i];
        p2 += Model->p2[k];
        for (i = 1; i <= p2; i++) tmp[i] = phia[i];
    }
    for (i = 0; i <= p2; i++) phia[i] = tmp[i];

    /* [6]: Regular MA operator */
    for (k = 1; k <= Model->NumMa1; k++) {
        if ((Model->q1[k] == 1) && (fabs(Model->Ma1[k][1]) > 1.0))
            Model->Ma1[k][1] = 1.0 / Model->Ma1[k][1];
        for (i = 1; i <= qr; i++) thetar[i] = 0.0;
        thetar[0] = tmp[0] = -1.0;
        for (i = 0; i <= q1; i++)
            for (j = 0; j <= Model->q1[k]; j++)
                thetar[j + i] -= Model->Ma1[k][j] * tmp[i];
        q1 += Model->q1[k];
        for (i = 1; i <= q1; i++) tmp[i] = thetar[i];
    }

    /* [7]: Regular MA(2) with fixed frequency */
    for (k = 1; k <= Model->NumMa1f; k++) {
        if (Model->Ma1f[k][2] < -1.0) Model->Ma1f[k][2] = 1.0 / Model->Ma1f[k][2];
        if (Model->Ma1f[k][2] > 0.0)  { *ifault = 1; goto u1; }
        Model->Ma1f[k][1] = 2.0 * Cos2(Model->qfre1[k], Model->sper)
                            * sqrt(-Model->Ma1f[k][2]);
        for (i = 1; i <= qr; i++) thetar[i] = 0.0;
        thetar[0] = tmp[0] = -1.0;
        for (i = 0; i <= q1; i++)
            for (j = 0; j <= 2; j++)
                thetar[j + i] -= Model->Ma1f[k][j] * tmp[i];
        q1 += 2;
        for (i = 1; i <= q1; i++) tmp[i] = thetar[i];
    }
    for (i = 0; i <= q1; i++) thetar[i] = tmp[i];

    /* [8]: Annual MA operator */
    for (k = 1; k <= Model->NumMa2; k++) {
        if ((Model->q2[k] == 1) && (fabs(Model->Ma2[k][1]) > 1.0))
            Model->Ma2[k][1] = 1.0 / Model->Ma2[k][1];
        for (i = 1; i <= qa; i++) thetaa[i] = 0.0;
        thetaa[0] = tmp[0] = -1.0;
        for (i = 0; i <= q2; i++)
            for (j = 0; j <= Model->q2[k]; j++)
                thetaa[j + i] -= Model->Ma2[k][j] * tmp[i];
        q2 += Model->q2[k];
        for (i = 1; i <= q2; i++) tmp[i] = thetaa[i];
    }
    for (i = 0; i <= q2; i++) thetaa[i] = tmp[i];

    /* [10]: Combine regular × annual into full operators */
    *OrderAr = pr + Model->sper * pa;
    *OrderMa = qr + Model->sper * qa;

    for (i = 1; i <= *OrderAr; i++) ArFactor[i] = 0.0;
    for (i = 1; i <= *OrderMa; i++) MaFactor[i] = 0.0;
    ArFactor[0] = MaFactor[0] = -1.0;

    for (i = 0; i <= pa; i++)
        for (j = 0; j <= pr; j++)
            ArFactor[j + i * Model->sper] -= phir[j] * phia[i];

    for (i = 0; i <= qa; i++)
        for (j = 0; j <= qr; j++)
            MaFactor[j + i * Model->sper] -= thetar[j] * thetaa[i];

u1:
    free_vector(tmp,    0, pq   > 0 ? pq   : 1);
    free_vector(thetaa, 0, qa   > 0 ? qa   : 1);
    free_vector(thetar, 0, qr   > 0 ? qr   : 1);
    free_vector(phia,   0, pa   > 0 ? pa   : 1);
    free_vector(phir,   0, pr   > 0 ? pr   : 1);
}

/* ── cast_us: translate Tm+Ts+DataMat into Tvarma for the engine ─────────  */
/* Extracted from fue.c:3645. Uses module-level Tm, Ts, DataMat globals.     */

static void cast_us(real *x, struct Tvarma *armax,
                    int *ifaultx, int firstx, int lastx)
{
    int i, j, k, itmp;
    int OrderAr, OrderMa, *NuLag;
    double *vtmp1, *vtmp2, **Nu, tmp1, tmp2;

    itmp = 0;

    /* [1]: Unpack parameter vector x into Tm */

    for (i = 1; i <= Tm.NdetVar; i++)
        for (j = 0; j <= Tm.Nomega[i]; j++)
            if (Tm.Imega[i][j] == 1) { itmp++; Tm.Omega[i][j] = x[itmp]; }

    for (i = 1; i <= Tm.NdetVar; i++)
        for (j = 1; j <= Tm.Ndelta[i]; j++)
            if (Tm.Ielta[i][j] == 1) { itmp++; Tm.Delta[i][j] = x[itmp]; }

    for (i = 1; i <= Tm.NumAr1; i++)
        for (j = 1; j <= Tm.p1[i]; j++)
            if (Tm.Ia1[i][j] == 1) { itmp++; Tm.Ar1[i][j] = x[itmp]; }

    for (i = 1; i <= Tm.NumAr2; i++)
        for (j = 1; j <= Tm.p2[i]; j++)
            if (Tm.Ia2[i][j] == 1) { itmp++; Tm.Ar2[i][j] = x[itmp]; }

    for (i = 1; i <= Tm.NumMa1; i++)
        for (j = 1; j <= Tm.q1[i]; j++)
            if (Tm.Im1[i][j] == 1) { itmp++; Tm.Ma1[i][j] = x[itmp]; }

    for (i = 1; i <= Tm.NumMa2; i++)
        for (j = 1; j <= Tm.q2[i]; j++)
            if (Tm.Im2[i][j] == 1) { itmp++; Tm.Ma2[i][j] = x[itmp]; }

    for (i = 1; i <= Tm.NumAr1f; i++)
        if (Tm.Ia1f[i] == 1) { itmp++; Tm.Ar1f[i][2] = x[itmp]; }

    for (i = 1; i <= Tm.NumMa1f; i++)
        if (Tm.Im1f[i] == 1) { itmp++; Tm.Ma1f[i][2] = x[itmp]; }

    if (Tm.Imu == 1) { itmp++; Tm.mu = x[itmp]; }

    /* [2]: Compute AR/MA orders */

    OrderAr = OrderMa = 0;
    for (i = 1; i <= Tm.NumAr1;  i++) OrderAr += Tm.p1[i];
    for (i = 1; i <= Tm.NumAr2;  i++) OrderAr += Tm.p2[i] * Tm.sper;
    for (i = 1; i <= Tm.NumAr1f; i++) OrderAr += 2;
    for (i = 1; i <= Tm.NumMa1;  i++) OrderMa += Tm.q1[i];
    for (i = 1; i <= Tm.NumMa2;  i++) OrderMa += Tm.q2[i] * Tm.sper;
    for (i = 1; i <= Tm.NumMa1f; i++) OrderMa += 2;

    armax->m = 1;
    armax->n = Ts.nobs - Tm.ornsop;
    armax->p = OrderAr;
    armax->q = OrderMa;

    /* [3]: First allocation of VARMA structure */

    if (firstx) {
        armax->mu    = vector(1, armax->m);
        armax->phi   = tensor(0, armax->p, 1, armax->m, 1, armax->m);
        armax->theta = tensor(0, armax->q, 1, armax->m, 1, armax->m);
        armax->qq    = matrix(1, armax->m, 1, armax->m);
        armax->w     = matrix(1, armax->m, 1, armax->n);
        armax->a     = matrix(1, armax->m, 1, armax->n);

        for (i = 1; i <= armax->m; i++) {
            armax->mu[i] = 0.0;
            for (j = 1; j <= armax->m; j++) {
                for (k = 0; k <= armax->p; k++) armax->phi[k][i][j]   = 0.0;
                for (k = 0; k <= armax->q; k++) armax->theta[k][i][j] = 0.0;
                armax->qq[i][j] = 0.0;
            }
            for (j = 1; j <= armax->n; j++) {
                armax->w[i][j] = 0.0;
                armax->a[i][j] = 0.0;
            }
        }
    }

    /* [4]: Expand AR/MA factors via unscramble */

    *ifaultx = 0;
    vtmp1 = vector(0, OrderAr);
    vtmp2 = vector(0, OrderMa);

    unscramble(&Tm, &OrderAr, &OrderMa, vtmp1, vtmp2, ifaultx);
    for (i = 1; i <= OrderAr; i++) armax->phi[i][1][1]   = vtmp1[i];
    for (i = 1; i <= OrderMa; i++) armax->theta[i][1][1] = vtmp2[i];

    free_vector(vtmp2, 0, OrderMa);
    free_vector(vtmp1, 0, OrderAr);

    armax->mu[1]    = Tm.mu;
    armax->qq[1][1] = 1.0;

    /* [5]: Compute differenced noise series */

    vtmp1 = vector(1, Ts.nobs);

    if (Tm.NdetVar > 0) {
        NuLag = ivector(1, Tm.NdetVar);
        Nu    = (double **)calloc((size_t)(Tm.NdetVar + 1), sizeof(double *));
        for (j = 1; j <= Tm.NdetVar; j++) {
            NuLag[j] = (Tm.Ndelta[j] == 0) ? Tm.Nomega[j] : 40;
            Nu[j]    = vector(0, NuLag[j]);
            calcnu(Tm.Omega[j], Tm.Nomega[j],
                   Tm.Delta[j], Tm.Ndelta[j],
                   Nu[j], NuLag[j]);
        }
        for (i = 1; i <= Ts.nobs; i++) {
            tmp1 = 0.0;
            for (j = 1; j <= Tm.NdetVar; j++) {
                tmp2 = 0.0;
                for (k = 0; k <= NuLag[j]; k++)
                    if (i - k >= 1) tmp2 += Nu[j][k] * DataMat[j][i - k];
                tmp1 += tmp2;
            }
            vtmp1[i] = DataMat[0][i] - tmp1;
        }
        for (j = Tm.NdetVar; j >= 1; j--) free_vector(Nu[j], 0, NuLag[j]);
        free(Nu);
        free_ivector(NuLag, 1, Tm.NdetVar);
    } else {
        for (i = 1; i <= Ts.nobs; i++) vtmp1[i] = DataMat[0][i];
    }

    for (j = Tm.ornsop + 1; j <= Ts.nobs; j++) {
        tmp1 = 0.0;
        for (i = 1; i <= Tm.ornsop; i++) tmp1 -= Tm.rnsop[i] * vtmp1[j - i];
        armax->w[1][j - Tm.ornsop] = vtmp1[j] + tmp1;
    }

    free_vector(vtmp1, 1, Ts.nobs);

    /* [6]: Last deallocation */

    if (lastx == 1) {
        free_matrix(armax->a,     1, armax->m, 1, armax->n);
        free_matrix(armax->w,     1, armax->m, 1, armax->n);
        free_matrix(armax->qq,    1, armax->m, 1, armax->m);
        free_tensor(armax->theta, 0, armax->q, 1, armax->m, 1, armax->m);
        free_tensor(armax->phi,   0, armax->p, 1, armax->m, 1, armax->m);
        free_vector(armax->mu,    1, armax->m);
    }
}

/* ── populate_globals: FueModelSpec → Tm + Ts + DataMat ─────────────────── */

static int populate_globals(const FueModelSpec *spec)
{
    int i, j, obs;

    memset(&Tm, 0, sizeof(Tm));
    memset(&Ts, 0, sizeof(Ts));

    /* ── Tseries ── */
    Ts.nobs    = spec->nobs;
    Ts.freq    = spec->sper;
    Ts.begyear = spec->begyear;
    Ts.begtime = spec->begtime;
    Ts.refactor = spec->refactor > 0.0 ? spec->refactor : 1.0;
    Ts.data    = vector(1, Ts.nobs);
    for (i = 1; i <= Ts.nobs; i++) Ts.data[i] = spec->data[i - 1];

    /* ── Tusmodel scalars ── */
    Tm.sper    = spec->sper > 0 ? spec->sper : 1;
    Tm.boxlam  = spec->boxlam;
    Tm.nrdiff  = spec->nrdiff;
    Tm.nadiff  = spec->nadiff;
    Tm.mu      = spec->mu0;
    Tm.Imu     = spec->estimate_mu ? 1 : 0;

    /* ornsop = order of non-stationary operator (regular + full seasonal diffs) */
    Tm.ornsop  = Tm.nrdiff + Tm.nadiff * Tm.sper;

    /* ifadf: individual annual difference factors (all zeros = full diff only) */
    Tm.ifadf   = (int *)calloc(8, sizeof(int));

    /* rnsop: non-stationary operator coefficients */
    if (Tm.ornsop > 0) {
        Tm.rnsop = vector(0, Tm.ornsop);
        CalcNonsOp(Tm.sper, Tm.nrdiff, Tm.nadiff, Tm.ifadf,
                   Tm.ornsop, Tm.rnsop);
    } else {
        Tm.rnsop = vector(0, 1);
        Tm.rnsop[0] = -1.0;
    }

    /* ── Interventions (deterministic variables) ── */
    Tm.NdetVar = spec->ninterventions;

    if (Tm.NdetVar > 0) {
        Tm.Nomega = ivector(1, Tm.NdetVar);
        Tm.Omega  = (double **)calloc((size_t)(Tm.NdetVar + 1), sizeof(double *));
        Tm.Imega  = (int   **)calloc((size_t)(Tm.NdetVar + 1), sizeof(int    *));
        Tm.Ndelta = ivector(1, Tm.NdetVar);
        Tm.Delta  = (double **)calloc((size_t)(Tm.NdetVar + 1), sizeof(double *));
        Tm.Ielta  = (int   **)calloc((size_t)(Tm.NdetVar + 1), sizeof(int    *));

        for (i = 1; i <= Tm.NdetVar; i++) {
            const FueIntervention *itv = &spec->interventions[i - 1];
            Tm.Nomega[i] = itv->nomega;
            Tm.Omega[i]  = vector(0, Tm.Nomega[i]);
            Tm.Imega[i]  = ivector(0, Tm.Nomega[i]);
            for (j = 0; j <= Tm.Nomega[i]; j++) {
                Tm.Omega[i][j] = itv->omega[j];
                Tm.Imega[i][j] = itv->omega_free[j] ? 1 : 0;
            }

            Tm.Ndelta[i] = itv->ndelta;
            if (Tm.Ndelta[i] > 0) {
                Tm.Delta[i] = vector(1, Tm.Ndelta[i]);
                Tm.Ielta[i] = ivector(1, Tm.Ndelta[i]);
                for (j = 1; j <= Tm.Ndelta[i]; j++) {
                    Tm.Delta[i][j] = itv->delta[j - 1];
                    Tm.Ielta[i][j] = itv->delta_free[j - 1] ? 1 : 0;
                }
            } else {
                Tm.Delta[i] = NULL;
                Tm.Ielta[i] = NULL;
            }
        }
    }

    /* ── AR factors (regular) ── */
    Tm.NumAr1 = spec->nar1;
    if (Tm.NumAr1 > 0) {
        Tm.p1  = ivector(1, Tm.NumAr1);
        Tm.Ar1 = (double **)calloc((size_t)(Tm.NumAr1 + 1), sizeof(double *));
        Tm.Ia1 = (int   **)calloc((size_t)(Tm.NumAr1 + 1), sizeof(int    *));
        for (i = 1; i <= Tm.NumAr1; i++) {
            const FueFactor *f = &spec->ar1[i - 1];
            Tm.p1[i]  = f->order;
            Tm.Ar1[i] = vector(0, f->order);
            Tm.Ia1[i] = ivector(0, f->order);
            for (j = 1; j <= f->order; j++) {
                Tm.Ar1[i][j] = f->coefs[j - 1];
                Tm.Ia1[i][j] = f->coef_free[j - 1] ? 1 : 0;
            }
        }
    }

    /* ── AR factors (seasonal) ── */
    Tm.NumAr2 = spec->nar2;
    if (Tm.NumAr2 > 0) {
        Tm.p2  = ivector(1, Tm.NumAr2);
        Tm.Ar2 = (double **)calloc((size_t)(Tm.NumAr2 + 1), sizeof(double *));
        Tm.Ia2 = (int   **)calloc((size_t)(Tm.NumAr2 + 1), sizeof(int    *));
        for (i = 1; i <= Tm.NumAr2; i++) {
            const FueFactor *f = &spec->ar2[i - 1];
            Tm.p2[i]  = f->order;
            Tm.Ar2[i] = vector(0, f->order);
            Tm.Ia2[i] = ivector(0, f->order);
            for (j = 1; j <= f->order; j++) {
                Tm.Ar2[i][j] = f->coefs[j - 1];
                Tm.Ia2[i][j] = f->coef_free[j - 1] ? 1 : 0;
            }
        }
    }

    /* ── MA factors (regular) ── */
    Tm.NumMa1 = spec->nma1;
    if (Tm.NumMa1 > 0) {
        Tm.q1  = ivector(1, Tm.NumMa1);
        Tm.Ma1 = (double **)calloc((size_t)(Tm.NumMa1 + 1), sizeof(double *));
        Tm.Im1 = (int   **)calloc((size_t)(Tm.NumMa1 + 1), sizeof(int    *));
        for (i = 1; i <= Tm.NumMa1; i++) {
            const FueFactor *f = &spec->ma1[i - 1];
            Tm.q1[i]  = f->order;
            Tm.Ma1[i] = vector(0, f->order);
            Tm.Im1[i] = ivector(0, f->order);
            for (j = 1; j <= f->order; j++) {
                Tm.Ma1[i][j] = f->coefs[j - 1];
                Tm.Im1[i][j] = f->coef_free[j - 1] ? 1 : 0;
            }
        }
    }

    /* ── MA factors (seasonal) ── */
    Tm.NumMa2 = spec->nma2;
    if (Tm.NumMa2 > 0) {
        Tm.q2  = ivector(1, Tm.NumMa2);
        Tm.Ma2 = (double **)calloc((size_t)(Tm.NumMa2 + 1), sizeof(double *));
        Tm.Im2 = (int   **)calloc((size_t)(Tm.NumMa2 + 1), sizeof(int    *));
        for (i = 1; i <= Tm.NumMa2; i++) {
            const FueFactor *f = &spec->ma2[i - 1];
            Tm.q2[i]  = f->order;
            Tm.Ma2[i] = vector(0, f->order);
            Tm.Im2[i] = ivector(0, f->order);
            for (j = 1; j <= f->order; j++) {
                Tm.Ma2[i][j] = f->coefs[j - 1];
                Tm.Im2[i][j] = f->coef_free[j - 1] ? 1 : 0;
            }
        }
    }

    /* Fixed-frequency factors not supported in API (disabled in fue.c too) */
    Tm.NumAr1f = Tm.NumAr2f = Tm.NumMa1f = Tm.NumMa2f = 0;

    /* Normalize polynomial 0th element to -1 (fue.c:862-869).
       unscramble() uses Ar?[i][0] in its convolution; leaving it at 0
       produces a different landscape in factor-phi space. */
    for (i = 1; i <= Tm.NumAr1; i++) Tm.Ar1[i][0] = -1.0;
    for (i = 1; i <= Tm.NumAr2; i++) Tm.Ar2[i][0] = -1.0;
    for (i = 1; i <= Tm.NumMa1; i++) Tm.Ma1[i][0] = -1.0;
    for (i = 1; i <= Tm.NumMa2; i++) Tm.Ma2[i][0] = -1.0;

    /* ── DataMat: series data + intervention indicators ── */

    DataMat = matrix(0, Tm.NdetVar, 1, Ts.nobs);

    /* DataMat[0]: Box-Cox transformed series */
    if (spec->boxlam == 0.0)
        for (i = 1; i <= Ts.nobs; i++)
            DataMat[0][i] = Ts.refactor * log(Ts.data[i]);
    else
        for (i = 1; i <= Ts.nobs; i++)
            DataMat[0][i] = Ts.refactor * Ts.data[i];

    /* DataMat[i]: intervention indicator series */
    for (i = 1; i <= Tm.NdetVar; i++) {
        const FueIntervention *itv = &spec->interventions[i - 1];
        obs = itv->obs_index + 1;  /* convert 0-based → 1-based */
        for (j = 1; j <= Ts.nobs; j++) DataMat[i][j] = 0.0;

        switch (itv->type) {
        case FUE_ITV_PULSE:
            if (obs >= 1 && obs <= Ts.nobs) DataMat[i][obs] = 1.0;
            break;
        case FUE_ITV_STEP:
            if (obs >= 1 && obs <= Ts.nobs)
                for (j = obs; j <= Ts.nobs; j++) DataMat[i][j] = 1.0;
            break;
        case FUE_ITV_RAMP:
            if (obs >= 1 && obs <= Ts.nobs)
                for (j = obs; j <= Ts.nobs; j++) DataMat[i][j] = j - obs + 1;
            break;
        case FUE_ITV_SEASONAL:
            /* obs_index is the period within the year (1-based after conversion) */
            for (j = 1; j <= Ts.nobs; j++)
                if (((j - Ts.begtime) % Ts.freq) + 1 == obs) DataMat[i][j] = 1.0;
            break;
        case FUE_ITV_COS: {
            const double PI = 3.141592654;
            double k = itv->harmonic;
            for (j = 1; j <= Ts.nobs; j++)
                DataMat[i][j] = cos(2.0 * PI * k / Ts.freq * j);
            break;
        }
        case FUE_ITV_SIN: {
            const double PI = 3.141592654;
            double k = itv->harmonic;
            for (j = 1; j <= Ts.nobs; j++)
                DataMat[i][j] = sin(2.0 * PI * k / Ts.freq * j);
            break;
        }
        case FUE_ITV_ALTER:
            for (j = 1; j <= Ts.nobs; j++)
                DataMat[i][j] = (j % 2 == 0) ? 1.0 : -1.0;
            break;
        }
    }

    return 0;
}

/* ── count_npar_build_par ────────────────────────────────────────────────── */
/* Replicated from fue.c:921-1050. Returns npar and fills par[1..npar].      */

static int count_npar_build_par(real *par)
{
    int i, j, n = 0;

    for (i = 1; i <= Tm.NdetVar; i++)
        for (j = 0; j <= Tm.Nomega[i]; j++)
            if (Tm.Imega[i][j] == 1) { n++; if (par) par[n] = Tm.Omega[i][j]; }

    for (i = 1; i <= Tm.NdetVar; i++)
        for (j = 1; j <= Tm.Ndelta[i]; j++)
            if (Tm.Ielta[i][j] == 1) { n++; if (par) par[n] = Tm.Delta[i][j]; }

    for (i = 1; i <= Tm.NumAr1; i++)
        for (j = 1; j <= Tm.p1[i]; j++)
            if (Tm.Ia1[i][j] == 1) { n++; if (par) par[n] = Tm.Ar1[i][j]; }

    for (i = 1; i <= Tm.NumAr2; i++)
        for (j = 1; j <= Tm.p2[i]; j++)
            if (Tm.Ia2[i][j] == 1) { n++; if (par) par[n] = Tm.Ar2[i][j]; }

    for (i = 1; i <= Tm.NumMa1; i++)
        for (j = 1; j <= Tm.q1[i]; j++)
            if (Tm.Im1[i][j] == 1) { n++; if (par) par[n] = Tm.Ma1[i][j]; }

    for (i = 1; i <= Tm.NumMa2; i++)
        for (j = 1; j <= Tm.q2[i]; j++)
            if (Tm.Im2[i][j] == 1) { n++; if (par) par[n] = Tm.Ma2[i][j]; }

    for (i = 1; i <= Tm.NumAr1f; i++)
        if (Tm.Ia1f[i] == 1) { n++; if (par) par[n] = Tm.Ar1f[i][2]; }

    for (i = 1; i <= Tm.NumMa1f; i++)
        if (Tm.Im1f[i] == 1) { n++; if (par) par[n] = Tm.Ma1f[i][2]; }

    if (Tm.Imu == 1) { n++; if (par) par[n] = Tm.mu; }

    return n;
}

/* ── fue_estimate ───────────────────────────────────────────────────────── */

FueResult *fue_estimate(const FueModelSpec *spec)
{
    real *par, *dev, **cov, **a;
    real  sigma2, logelf;
    int   npar, ifault, nresiduals;
    int   i, j;
    FueResult *result;

    if (!spec || !spec->data || spec->nobs <= 0) {
        result = calloc(1, sizeof(*result));
        if (result) result->ifault = -1;
        return result;
    }

    /* Redirect engine diagnostics to null sink (outputv=NULL causes crashes) */
    if (outputv == NULL) {
#ifdef _WIN32
        outputv = fopen("NUL", "w");
#else
        outputv = fopen("/dev/null", "w");
#endif
        if (!outputv) outputv = stderr;
    }

    /* [1]: Populate module-level globals (Tm, Ts, DataMat) */
    if (populate_globals(spec) != 0) return NULL;

    /* [2]: Count free parameters (first pass, no par vector) */
    npar = count_npar_build_par(NULL);
    nresiduals = Ts.nobs - Tm.ornsop;

    result = alloc_result(npar, nresiduals);
    if (!result) return NULL;

    if (npar == 0) {
        result->ifault = 0;
        return result;
    }

    /* [3]: Allocate optimizer workspace */
    par = vector(1, npar);
    dev = vector(1, npar);
    cov = matrix(1, npar, 1, npar);
    a   = matrix(1, 1, 1, nresiduals);

    /* [4]: Fill initial parameter vector */
    count_npar_build_par(par);

    /* [5]: Estimate
       xitol sign convention from fue.c:1087:
         exact ML (eml=1) → negative xitol forces exact Melard recursions in cxi()
         approx ML (eml=0) → positive xitol allows cheap approximation */
    {
        real abs_xitol = spec->xitol > 0 ? spec->xitol : 1e-3;
        real xitol_signed = spec->eml ? -abs_xitol : abs_xitol;
        est(cast_us, npar, par, dev, cov,
            spec->maxits > 0 ? spec->maxits : 200,
            10,
            spec->grtol  > 0 ? spec->grtol  : 1e-5,
            spec->sptol  > 0 ? spec->sptol  : 1e-7,
            xitol_signed,
            spec->chkma,
            a, &sigma2, &logelf, &ifault);
    }

    /* [6]: Pack results */
    result->ifault = ifault;
    result->sigma2 = sigma2;
    result->loglik = logelf;
    if (npar > 0 && nresiduals > 0) {
        result->aic = -2.0 * logelf + 2.0 * npar;
        result->bic = -2.0 * logelf + npar * log((double)nresiduals);
    }

    for (i = 1; i <= npar; i++) {
        result->params[i - 1]     = par[i];
        result->std_errors[i - 1] = dev[i];
        for (j = 1; j <= npar; j++)
            result->cov_matrix[(i-1)*npar + (j-1)] = cov[i][j];
    }
    for (i = 1; i <= nresiduals; i++)
        result->residuals[i - 1] = a[1][i];

    /* [7]: Free optimizer workspace */
    free_matrix(a,   1, 1,    1, nresiduals);
    free_matrix(cov, 1, npar, 1, npar);
    free_vector(dev, 1, npar);
    free_vector(par, 1, npar);

    return result;
}
