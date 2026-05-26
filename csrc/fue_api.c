/*
 * fue_api.c  —  Implementation of the public FUE estimation API.
 *
 * This file bridges fue_api.h (cffi-visible) and the internal estimation
 * engine (elfvarma.c, drvmlest.c, qnewtopt.c, usmelard.c, nlatools.c).
 *
 * The key step is converting FueModelSpec → Tusmodel → cast_us() → est().
 * cast_us() is extracted from fue.c and lives here.
 *
 * Copyright (C) 2009-2026 A.B Treadway, D.E. Guerrero & J.A. Mauricio
 * License: GPL v2 or later.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fue_api.h"
#include "internal/fue.h"
#include "internal/nlatools.h"

/* ── Globals required by the internal engine ───────────────────────────── */

real    macheps = 2.22e-16;   /* machine epsilon                           */
FILE   *outputv = NULL;       /* internal engine writes diagnostics here   */

/* ── fue_strerror ───────────────────────────────────────────────────────── */

const char *fue_strerror(int ifault)
{
    switch (ifault) {
    case 0: return "OK";
    case 1: return "bad initial estimate: Q not positive definite";
    case 2: return "bad initial estimate: AR has a unit root";
    case 3: return "bad initial estimate: AR is non-stationary";
    case 4: return "bad initial estimate: MA is non-invertible";
    case 5: return "bad initial estimate: unknown numerical problem";
    case 6: return "bad initial estimate: rejected by cast routine";
    default: return "unknown error";
    }
}

/* ── fue_defaults ───────────────────────────────────────────────────────── */

void fue_defaults(FueModelSpec *spec)
{
    memset(spec, 0, sizeof(*spec));
    spec->boxlam      = 1.0;   /* levels (no log)                          */
    spec->refactor    = 1.0;
    spec->estimate_mu = 0;
    spec->maxits      = 200;
    spec->grtol       = 1e-5;
    spec->sptol       = 1e-7;
    spec->xitol       = 1e-3;
    spec->chkma       = 1;
    spec->eml         = 1;
}

/* ── Internal helpers ───────────────────────────────────────────────────── */

static FueResult *alloc_result(int npar, int nresiduals)
{
    FueResult *r = calloc(1, sizeof(*r));
    if (!r) return NULL;
    r->npar       = npar;
    r->nresiduals = nresiduals;
    r->params     = calloc(npar,              sizeof(double));
    r->std_errors = calloc(npar,              sizeof(double));
    r->cov_matrix = calloc((size_t)npar*npar, sizeof(double));
    r->residuals  = calloc(nresiduals,        sizeof(double));
    if (!r->params || !r->std_errors || !r->cov_matrix || !r->residuals) {
        fue_result_free(r);
        return NULL;
    }
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

/* ── TODO: cast_us() extracted from fue.c ──────────────────────────────── */
/*
 * cast_us() converts the user model specification (Tusmodel) into the
 * internal VARMA structure (Tvarma) expected by est() / elf() / flikam().
 *
 * This function is currently in fue.c (~600 lines) and needs to be
 * extracted here as Phase 1 work.  The stub below marks the boundary.
 */

/* forward declaration — will be defined here after extraction from fue.c  */
static void cast_us(real *par, struct Tvarma *vm, int *ifault, int init, int free_flag);

/* placeholder Tusmodel used by cast_us (module-level, as in fue.c)        */
static struct Tusmodel usm;

/* ── fue_estimate ───────────────────────────────────────────────────────── */

FueResult *fue_estimate(const FueModelSpec *spec)
{
    /*
     * Phase 1 stub — full implementation pending cast_us() extraction.
     *
     * The flow will be:
     *   1. spec → usm  (populate Tusmodel from FueModelSpec)
     *   2. Count free parameters → npar
     *   3. Build par[] vector from initial estimates
     *   4. Allocate dev[], cov[][], a[]
     *   5. Call est(cast_us, npar, par, dev, cov, ...)
     *   6. Pack results into FueResult and return
     */

    FueResult *r = calloc(1, sizeof(*r));
    if (!r) return NULL;
    r->ifault = -1;  /* not yet implemented */
    return r;
}

/* ── cast_us stub ───────────────────────────────────────────────────────── */
/*
 * Full body to be filled in during Phase 1 (extraction from fue.c).
 * The signature matches what est() expects as its first argument.
 */
static void cast_us(real *par, struct Tvarma *vm, int *ifault, int init, int free_flag)
{
    (void)par; (void)vm; (void)ifault; (void)init; (void)free_flag;
    *ifault = 6;  /* signal: not yet implemented */
}
