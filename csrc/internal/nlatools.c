/***************************************************************************
 *   NLATOOLS -- Numerical linear algebra and dynamic memory allocation    *
 *   routines.                                                             *
 *                                                                         *
 *   Contributors:                                                         *
 *     Copyright (C) 1995  Jose Alberto Mauricio                           *
 *       LU and Cholesky decomposition, scalar and string utilities,       *
 *       Easter date calculation.                                          *
 *     Copyright (C) 2009  Arthur B. Treadway & David E. Guerrero          *
 *       davidesg@ucm.es                                                   *
 *       Eigenvalue solver (gsl_eigenqr) via GNU Scientific Library.       *
 *       Dynamic memory allocation routines and error handler.             *
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
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_matrix.h>

#include "nlatools.h"            /* Header file (prototype declarations)        */


/*****************************************************************************/
/* Jose Alberto Mauricio (1995) -- LU decomposition                          */
/*****************************************************************************/

void ludcp( double **a, int n, int *ip )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */

{
   int  i, j, k, kp1, m;
   double tmp;

   ip[n] = 1;
   for ( k = 1; k <= n; k++ )
       {
       if ( k != n )
          {
          kp1 = k + 1;
          m = k;
          for ( i = kp1; i <= n; i++ )
              if ( fabs( a[i][k] ) > fabs( a[m][k] ) ) m = i;
          ip[k] = m;
          if ( m != k ) ip[n] = -ip[n];
          tmp = a[m][k];
          a[m][k] = a[k][k];
          a[k][k] = tmp;
          if ( tmp != 0.0 )
             {
             for ( i = kp1; i <= n; i++ ) a[i][k] /= -tmp;
             for ( j = kp1; j <= n; j++ )
                 {
                 tmp = a[m][j];
                 a[m][j] = a[k][j];
                 a[k][j] = tmp;
                 if ( tmp != 0.0 )
                    for ( i = kp1; i <= n; i++ ) a[i][j] += a[i][k] * tmp;
                 }
             }
          }
       if ( a[k][k] == 0.0 ) ip[n] = 0;
       }
}

/****************************************************************************/

void lusol( double **a, double *b, int n, int *ip )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   int  i, k, km1, kp1, k1, m, nm1;
   double tmp;

   if ( n != 1 )
      {
      nm1 = n - 1;
      for ( k = 1; k <= nm1; k++ )
          {
          kp1 = k + 1;
          m = ip[k];
          tmp = b[m];
          b[m] = b[k];
          b[k] = tmp;
          for ( i = kp1; i <= n; i++ ) b[i] += a[i][k] * tmp;
          }
      for ( k1 = 1; k1 <= nm1; k1++ )
          {
          km1 = n - k1;
          k = km1 + 1;
          b[k] /= a[k][k];
          tmp = -b[k];
          for ( i = 1; i <= km1; i++ ) b[i] += a[i][k] * tmp;
          }
      }
   b[1] /= a[1][1];
}

/****************************************************************************/
/*****************************************************************************/
/* Jose Alberto Mauricio (1995) -- Cholesky decomposition                    */
/*****************************************************************************/

void choldcp( double **mat, int n, double *d1, double *d2, int *ifault )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   int   i, j, k;
   double  sum1, minl, maxoffl, minl2, maxadd, minljj, sqrteps, macheps;

   *ifault = 0;
   *d1     = 1.0;
   *d2     = 0.0;
   macheps = cmacheps();
   sqrteps = sqrt( macheps );

/* [1]: check wether mat is numerically different from zero (redundancy):   */

   for ( maxoffl = 0.0, j = 1; j <= n; j++ )
       if ( sqrt( fabs( mat[j][j] ) ) > maxoffl )
          maxoffl = sqrt( fabs( mat[j][j] ) );
   if ( maxoffl * maxoffl <= sqrteps )
      {
      for ( i = 1; i <= n; i++ ) for ( j = 1; j <= n; j++ ) mat[i][j] = 0.0;
      return;
      }

/* [2]: initialize finite-arithmetic constants:                             */

   minl   = 0.0;
   minl2  = sqrteps * maxoffl;
   maxadd = 0.0;

/* [3]: form the j-th column of the Cholesky factor of mat:                 */

   for ( j = 1; j <= n; j++ )
       {
       sum1 = mat[j][j];
       for ( i = 1; i <= j - 1; i++ ) sum1 -= mat[j][i] * mat[j][i];

       if ( (sum1 != fabs( sum1 )) && (fabs( sum1 ) > minl2) )
          {
          *ifault = 1;
          return;
          }
       else
          mat[j][j] = sum1;

       minljj = 0.0;
       for ( i = j + 1; i <= n; i++ )
           {
           sum1 = mat[j][i];
           for ( k = 1; k <= j - 1; k++ ) sum1 -= mat[i][k] * mat[j][k];
           mat[i][j] = sum1;
           if ( fabs( mat[i][j] ) > minljj ) minljj = fabs( mat[i][j] );
           }

       if ( (minljj / maxoffl) > minl )
          minljj /= maxoffl;
       else
          minljj = minl;

       if ( mat[j][j] > (minljj * minljj) )
          mat[j][j] = sqrt( mat[j][j] );
       else
          {
          if ( minljj < minl2 )
             minljj = minl2;
          if ( maxadd < (minljj * minljj - mat[j][j]) )
             maxadd = minljj * minljj - mat[j][j];
          mat[j][j] = minljj;
          }

       *d1 *= mat[j][j] * mat[j][j];
       while ( fabs( *d1 ) >= 1.0 )
             {
             *d1 *= 0.0625;
             *d2 += 4.0;
             }
       while ( fabs( *d1 ) < 0.0625 )
             {
             *d1 *= 16.0;
             *d2 -=  4.0;
             }

       for ( i = j + 1; i <= n; i++ ) mat[i][j] /= mat[j][j];
       }

   for ( j = 2; j <= n; j++ )
       for ( i = 1; i <= j - 1; i++ ) mat[i][j] = 0.0;

}

/****************************************************************************/

void cholfor( double **matl, int n, double *rhsol )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   int   i, j;
   double  tmp;

   rhsol[1] /= matl[1][1];

   for ( i = 2; i <= n; i++ )
       {
       tmp = 0.0;
       for ( j = 1; j <= i - 1; j++ )
           tmp += matl[i][j] * rhsol[j];
       rhsol[i] = (rhsol[i] - tmp) / matl[i][i];
       }
}

/****************************************************************************/

void cholbak( double **matl, int n, double *rhsol )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   int   i, j;
   double  tmp;

   rhsol[n] /= matl[n][n];

   for ( i = n - 1; i >= 1; i-- )
       {
       tmp = 0.0;
       for ( j = i + 1; j <= n; j++ )
           tmp += matl[j][i] * rhsol[j];
       rhsol[i] = (rhsol[i] - tmp) / matl[i][i];
       }
}

/****************************************************************************/

void cholsol( double **matl, int n, double *rhsol )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   cholfor( matl, n, rhsol );
   cholbak( matl, n, rhsol );
}

/*****************************************************************************/
/* Jose Alberto Mauricio (1995) -- scalar utilities                          */
/*****************************************************************************/

double rmax( double a, double b )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   return( ( a > b ) ? a : b );
}

/****************************************************************************/

double rmin( double a, double b )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   return( ( a < b ) ? a : b );
}

/****************************************************************************/

double cmacheps( void )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   double e;

   e = 1.0;
   do {
      e /= 2.0;
   } while ( (1.0 + e ) > 1.0 );
   return( 2.0 * e );
}

/****************************************************************************/

int iround( double num )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   double tmp1, tmp2;
   int  itmp;

   tmp1 = fabs( num );
   tmp2 = floor( tmp1 );
   if ( tmp1 - tmp2 >= 0.5 ) tmp2 = ceil( tmp1 );
   itmp = (int)tmp2;

   if ( num >= 0.0 )
      return( itmp );
   else
      return( -itmp );
}

/****************************************************************************/

double pythag( double a, double b )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   double at, bt;

   at = fabs( a );
   bt = fabs( b );
   at *= at;
   bt *= bt;
   return( sqrt( at + bt ) );
}

/*****************************************************************************/
/* Arthur B. Treadway & David E. Guerrero (2009) -- eigenvalue solver        */
/*                                                                           */
/* Computes the eigenvalues of a real unsymmetric n x n matrix a (1-based). */
/* Eigenvalues are returned in wr[1..n] (real parts) and wi[1..n] (imaginary */
/* parts).  Uses gsl_eigen_nonsymm from the GNU Scientific Library.          */
/*****************************************************************************/

void gsl_eigenqr( double **a, int n, double *wr, double *wi )
{
   int i, j;
   gsl_matrix *A = gsl_matrix_alloc( (size_t)n, (size_t)n );
   gsl_vector_complex *eval = gsl_vector_complex_alloc( (size_t)n );
   gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc( (size_t)n );
   gsl_complex ev;

   for ( i = 0; i < n; i++ )
      for ( j = 0; j < n; j++ )
         gsl_matrix_set( A, i, j, a[i+1][j+1] );

   /* balance=1: apply balancing similarity transform before eigendecomposition */
   gsl_eigen_nonsymm_params( 0, 1, w );
   gsl_eigen_nonsymm( A, eval, w );

   for ( i = 0; i < n; i++ )
      {
      ev = gsl_vector_complex_get( eval, i );
      wr[i+1] = GSL_REAL( ev );
      wi[i+1] = GSL_IMAG( ev );
      }

   /* copy modified matrix back; set diagonal to modulus of each eigenvalue */
   for ( i = 0; i < n; i++ )
      {
      for ( j = 0; j < n; j++ )
         a[i+1][j+1] = gsl_matrix_get( A, i, j );
      a[i+1][i+1] = pythag( wr[i+1], wi[i+1] );
      }

   gsl_eigen_nonsymm_free( w );
   gsl_vector_complex_free( eval );
   gsl_matrix_free( A );
}

/*****************************************************************************/
/* Arthur B. Treadway & David E. Guerrero (2009) -- dynamic memory           */
/*                                                                           */
/* All routines use 1-based (or nl-based) indexing compatible with the rest  */
/* of the codebase.  Allocations are sized to [0..nh] (or [0..nrh]/[0..ndh])*/
/* so that any nl >= 0 is valid without pointer arithmetic undefined         */
/* behavior.                                                                 */
/*****************************************************************************/

void nrerror( char error_text[] )
{
   fprintf( stderr, "Unrecoverable run-time error:\n%s\n... Exiting to system ...\n",
            error_text );
   exit( 1 );
}

/****************************************************************************/

double *vector( long nl, long nh )
{
   double *v = (double *)calloc( (size_t)(nh + 1), sizeof(double) );
   if ( !v ) nrerror( "ALLOCATION FAILURE in vector()" );
   return( v );
}

/****************************************************************************/

int *ivector( long nl, long nh )
{
   int *v = (int *)calloc( (size_t)(nh + 1), sizeof(int) );
   if ( !v ) nrerror( "ALLOCATION FAILURE in ivector()" );
   return( v );
}

/****************************************************************************/

double **matrix( long nrl, long nrh, long ncl, long nch )
{
   long i, nrow = nrh - nrl + 1;
   double **m;
   double *data;

   m = (double **)calloc( (size_t)(nrh + 1), sizeof(double *) );
   if ( !m ) nrerror( "ALLOCATION FAILURE 1 in matrix()" );

   data = (double *)calloc( (size_t)(nrow * (nch + 1)), sizeof(double) );
   if ( !data ) nrerror( "ALLOCATION FAILURE 2 in matrix()" );

   for ( i = nrl; i <= nrh; i++ )
      m[i] = data + (i - nrl) * (nch + 1);

   return( m );
}

/****************************************************************************/

int **imatrix( long nrl, long nrh, long ncl, long nch )
{
   long i, nrow = nrh - nrl + 1;
   int **m;
   int *data;

   m = (int **)calloc( (size_t)(nrh + 1), sizeof(int *) );
   if ( !m ) nrerror( "ALLOCATION FAILURE 1 in imatrix()" );

   data = (int *)calloc( (size_t)(nrow * (nch + 1)), sizeof(int) );
   if ( !data ) nrerror( "ALLOCATION FAILURE 2 in imatrix()" );

   for ( i = nrl; i <= nrh; i++ )
      m[i] = data + (i - nrl) * (nch + 1);

   return( m );
}

/****************************************************************************/

double ***tensor( long nrl, long nrh, long ncl, long nch, long ndl, long ndh )
{
   long i, j;
   long nrow = nrh - nrl + 1;
   long ncol = nch - ncl + 1;
   double ***t;
   double **planes;
   double *data;

   t = (double ***)calloc( (size_t)(nrh + 1), sizeof(double **) );
   if ( !t ) nrerror( "ALLOCATION FAILURE 1 in tensor()" );

   planes = (double **)calloc( (size_t)(nrow * (nch + 1)), sizeof(double *) );
   if ( !planes ) nrerror( "ALLOCATION FAILURE 2 in tensor()" );

   data = (double *)calloc( (size_t)(nrow * ncol * (ndh + 1)), sizeof(double) );
   if ( !data ) nrerror( "ALLOCATION FAILURE 3 in tensor()" );

   for ( i = nrl; i <= nrh; i++ )
      {
      t[i] = planes + (i - nrl) * (nch + 1);
      for ( j = ncl; j <= nch; j++ )
         t[i][j] = data + ( (i - nrl) * ncol + (j - ncl) ) * (ndh + 1);
      }

   return( t );
}

/****************************************************************************/

void free_vector( double *v, long nl, long nh )
{
   free( v );
}

/****************************************************************************/

void free_ivector( int *v, long nl, long nh )
{
   free( v );
}

/****************************************************************************/

void free_matrix( double **m, long nrl, long nrh, long ncl, long nch )
{
   if ( m ) { free( m[nrl] ); free( m ); }
}

/****************************************************************************/

void free_imatrix( int **m, long nrl, long nrh, long ncl, long nch )
{
   if ( m ) { free( m[nrl] ); free( m ); }
}

/****************************************************************************/

void free_tensor( double ***t, long nrl, long nrh, long ncl, long nch,
                  long ndl, long ndh )
{
   if ( t ) { free( t[nrl][ncl] ); free( t[nrl] ); free( t ); }
}

/*****************************************************************************/
/* Jose Alberto Mauricio (1995) -- string utilities                          */
/*****************************************************************************/

STRING NEW_STR( int size )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
// Returns the starting address of a (size + 1)-char string, allocates space
// and initializes the string to an empty string.

{
    STRING s;

    if ( (s = (STRING)malloc( (size_t)(size + 1) * sizeof( char ) )) == NULL )
       return( NULL );
    else
       {
       *s = '\0';
       return( s );
       }
}

/****************************************************************************/

void FREE_STR( STRING s )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
// Deallocates previously allocated space for string s.

{
    free( (char *)s );
}

/****************************************************************************/

int DELETE_STR( STRING s, int i, int n )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
// Removes n chars from string s starting from the i-th position of s.

{
    if ( (i >= 0) && (n > 0) && ((i + n) <= (int)strlen( s )) )
       {
       s[i] = '\0';
       strcat( s, s + i + n );
       return( OK );
       }
    else
       return( WRONG );
}

/****************************************************************************/

int COPY_STR( STRING source, int i, int n, STRING dest )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
// Fills string dest with a copy of a substring of source, starting from
// the i-th position of source and consisting of n chars.

{
    register int j;

    if ( (i >= 0) && (n > 0) && ((i + n) <= (int)strlen( source )) )
       {
       for ( j = 0; j <= n - 1; j++ )
           dest[j] = source[i + j ];
       dest[n] = '\0';
       return( OK );
       }
    else
       return( WRONG );
}

/****************************************************************************/

int INSERT_STR( STRING s1, STRING s, int i )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
// Inserts the string s1 into s at the i-th position of s;
// if i == strlen( s ), then s1 is added to s.

{
    register int j, k;

    j = strlen( s1 );
    if ( (i >= 0) && (i <= (int)strlen( s )) && (j > 0)  )
       {
       for ( k = strlen( s ) - i; k >= 0; k-- )
           s[i + j + k] = s[i + k];
       s += i;
       while ( *s1 ) *s++ = *s1++;
       return( OK );
       }
    else
       return( WRONG );
}

/****************************************************************************/

int POS_STR( STRING s1, STRING s2 )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
// Returns the position (in s2) of the first occurrence of s1 in s2,
// or -1 if s1 does not exist in s2.

{
    STRING aux1, aux2, prevs2 = s2;

    for ( aux1 = s1, aux2 = s2; *aux2; aux2++ )
        if ( *aux1 == *aux2 )
           {
           aux1++;
           if ( !*aux1 )
              return( aux2 - s2 - (aux1 - s1) + 1 );
           }
        else
          {
          aux1 = s1;
          aux2 = prevs2++;
          }
    return( -1 );
}

/****************************************************************************/

int CHANGE_STR( STRING s1, int i, STRING s2 )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
// Replaces strlen( s2 ) chars in s1, starting from the ith position of s1,
// with the chars from the string s2.

{
    int j, s1l, s2l;

    s1l = strlen( s1 );
    s2l = strlen( s2 );
    if ( (s1l == 0) || (i < 0) || (i > s1l) || ((i + s2l) > s1l) )
       return( WRONG );
    if ( s2l == 0 )
       return( OK );
    if ( s2l > (s1l - i + 1) )
       {
       s1[i] = '\0';
       strcat( s1, s2 );
       }
    else
       {
       s1 += i;
       for ( j = 0; j < s2l; j++ )
           *s1++ = *s2++;
       }
    return( OK );
}

/****************************************************************************/

void BLANKS_STR( STRING s )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   int k, len;

   len = strlen( s );
   while ( (len > 0) && (s[--len] == ' ') )     // Remove trailing blanks
      s[len] = '\0';
   k = 0;                                       // Remove leading blanks
   while ( (k <= len) && (s[k] == ' ') )
      k++;
   COPY_STR( s, k, len+1-k, s );
}

/****************************************************************************/

void UPCASE_STR( STRING s )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   int i, len;

   len = strlen( s );
   for ( i = 0; i < len; i++ )
       s[i] = toupper( s[i] );
}

/*****************************************************************************/
/* Jose Alberto Mauricio (1995) -- date utilities                            */
/*****************************************************************************/

void Easter( int *day, int *month, int year )
/*  Copyright (C) Jose Alberto Mauricio, 1995 */
{
   div_t a, b, c, d, e;

   a = div( year, 19 );
   b = div( year,  4 );
   c = div( year,  7 );
   d = div( 19 * a.rem + 24, 30 );
   e = div( 2 * b.rem + 4 * c.rem + 6 * d.rem + 5, 7 );

   *day = 22 + d.rem + e.rem;

   if ( *day <= 31 )
      *month = 3;
   else
      {
      *day  -= 31;
      *month =  4;
      }
}

/****************************************************************************/

#undef FREE_ARG
#undef SIGN
#undef SWAP
#undef RADIX
#undef IMIN

/****************************************************************************/
