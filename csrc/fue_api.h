/*
 * fue_api.h  —  Public C API for the FUE estimation engine.
 *
 * This header is the only surface exposed to Python via cffi.
 * No function pointers, no nested variable-length arrays.
 *
 * Copyright (C) 2009-2026 A.B Treadway, D.E. Guerrero & J.A. Mauricio
 * License: GPL v2 or later.
 */

#ifndef FUE_API_H
#define FUE_API_H

#ifdef __cplusplus
extern "C" {
#endif

/* ── Hard limits (match Tusmodel capacity) ─────────────────────────────── */

#define FUE_MAX_DETVARS  64   /* max interventions                         */
#define FUE_MAX_FACTORS   8   /* max AR or MA factors per type             */
#define FUE_MAX_POLYORD  16   /* max polynomial order per factor           */

/* ── Intervention types ─────────────────────────────────────────────────── */

#define FUE_ITV_PULSE     0   /* isolated impulse (pulse indicator)        */
#define FUE_ITV_STEP      1   /* permanent level shift (step indicator)    */
#define FUE_ITV_RAMP      2   /* linear ramp                               */
#define FUE_ITV_SEASONAL  3   /* periodic seasonal dummy                   */
#define FUE_ITV_COS       4   /* cosine component: cos(2π·k/freq·j)        */
#define FUE_ITV_SIN       5   /* sine component:   sin(2π·k/freq·j)        */
#define FUE_ITV_ALTER     6   /* alternator: (-1)^j                        */

/* ── Single intervention with linear transfer function ω(B)/δ(B) ───────── */

typedef struct {
    int    type;                       /* FUE_ITV_* constant               */
    int    obs_index;                  /* 0-based index of the event        */
    double harmonic;                   /* harmonic k for cos/sin types      */

    int    nomega;                     /* degree of ω(B) numerator (≥0)    */
    double omega[FUE_MAX_POLYORD];     /* ω₀, ω₁, ..., ω_{nomega}          */
    int    omega_free[FUE_MAX_POLYORD];/* 1 = estimate, 0 = fix            */

    int    ndelta;                     /* degree of δ(B) denominator (≥0)  */
    double delta[FUE_MAX_POLYORD];     /* δ₁, δ₂, ..., δ_{ndelta}          */
    int    delta_free[FUE_MAX_POLYORD];/* 1 = estimate, 0 = fix            */
} FueIntervention;

/* ── Single AR or MA polynomial factor ─────────────────────────────────── */

typedef struct {
    int    order;                      /* polynomial degree (≥1)           */
    double coefs[FUE_MAX_POLYORD];     /* coefficients φ₁ … φ_p            */
    int    coef_free[FUE_MAX_POLYORD]; /* 1 = estimate, 0 = fix            */
} FueFactor;

/* ── Full model specification ───────────────────────────────────────────── */

typedef struct {
    /* ── Series ── */
    int     nobs;              /* number of observations                   */
    double *data;              /* pointer to data array [nobs]             */
    int     sper;              /* seasonal period: 1(A), 4(Q), 12(M)      */
    int     numbering;         /* 1 = plain numbering (no date labels)     */
    int     begyear;           /* first observation year                   */
    int     begtime;           /* first observation period within year     */

    /* ── Transformation ── */
    double  boxlam;            /* Box-Cox λ: 0.0 = log, 1.0 = levels      */
    double  refactor;          /* rescaling factor (usually 1.0)           */

    /* ── Differencing ── */
    int     nrdiff;            /* regular differences d                    */
    int     nadiff;            /* annual differences D                     */

    /* ── Mean ── */
    double  mu0;               /* initial value for μ                      */
    int     estimate_mu;       /* 1 = include μ in estimation              */

    /* ── Interventions ── */
    int            ninterventions;
    FueIntervention interventions[FUE_MAX_DETVARS];

    /* ── AR factors: regular φ(B) ── */
    int       nar1;
    FueFactor ar1[FUE_MAX_FACTORS];

    /* ── AR factors: seasonal Φ(Bˢ) ── */
    int       nar2;
    FueFactor ar2[FUE_MAX_FACTORS];

    /* ── MA factors: regular θ(B) ── */
    int       nma1;
    FueFactor ma1[FUE_MAX_FACTORS];

    /* ── MA factors: seasonal Θ(Bˢ) ── */
    int       nma2;
    FueFactor ma2[FUE_MAX_FACTORS];

    /* ── Fixed-frequency AR(2) factors: 1 − 2r·cos(ω)B + r²B² ── */
    /* phi2 = −r² < 0 is estimated; phi1 = 2·cos(2π·freq/sper)·√(−phi2) */
    int    nar1f;
    double ar1f_freq[FUE_MAX_FACTORS]; /* fixed frequency (pfre1, cycles/sper) */
    double ar1f_coef[FUE_MAX_FACTORS]; /* phi2 initial value (must be < 0)     */
    int    ar1f_free[FUE_MAX_FACTORS]; /* 1 = estimate, 0 = fix               */

    /* ── Fixed-frequency MA(2) factors ── */
    int    nma1f;
    double ma1f_freq[FUE_MAX_FACTORS];
    double ma1f_coef[FUE_MAX_FACTORS];
    int    ma1f_free[FUE_MAX_FACTORS];

    /* ── Optimizer settings ── */
    int    maxits;             /* max optimizer iterations (default 200)   */
    double grtol;              /* gradient tolerance (default 1e-5)        */
    double sptol;              /* step tolerance (default 1e-7)            */
    double xitol;              /* quick-recursion switch (default 1e-3)    */
    int    chkma;              /* 1 = enforce MA invertibility             */
    int    eml;                /* 1 = exact ML, 0 = approximate ML        */
} FueModelSpec;

/* ── Estimation results ────────────────────────────────────────────────── */

typedef struct {
    int     ifault;            /* 0 = ok; 1-6 = error (see fue_strerror)  */
    int     npar;              /* number of free parameters                */
    double *params;            /* estimated parameters    [npar]           */
    double *std_errors;        /* standard errors         [npar]           */
    double *cov_matrix;        /* covariance matrix       [npar*npar], row-major */
    double *residuals;         /* residuals               [nresiduals]     */
    int     nresiduals;        /* = nobs - (nrdiff + nadiff*sper)          */
    double  sigma2;            /* concentrated residual variance           */
    double  loglik;            /* exact log-likelihood                     */
    double  aic;               /* AIC = -2*loglik + 2*npar                */
    double  bic;               /* BIC = -2*loglik + npar*log(nresiduals)  */
} FueResult;

/* ── API entry points ───────────────────────────────────────────────────── */

/*
 * fue_estimate()  —  Estimate model parameters by exact ML.
 *
 * Returns a heap-allocated FueResult; caller must free with fue_result_free().
 * Returns NULL only on memory allocation failure.
 * Check result->ifault for estimation errors.
 */
FueResult *fue_estimate(const FueModelSpec *spec);

/*
 * fue_defaults()  —  Fill optimizer fields with safe defaults.
 * Call this before populating a FueModelSpec to avoid uninitialized fields.
 */
void fue_defaults(FueModelSpec *spec);

/*
 * fue_result_free()  —  Release all memory owned by a FueResult.
 */
void fue_result_free(FueResult *r);

/*
 * fue_strerror()  —  Human-readable description of an ifault code.
 */
const char *fue_strerror(int ifault);

#ifdef __cplusplus
}
#endif

#endif /* FUE_API_H */
