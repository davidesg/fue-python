/***************************************************************************
 *   NLATOOLS -- Numerical linear algebra and dynamic memory allocation    *
 *   routines.                                                             *
 *                                                                         *
 *   Contributors:                                                         *
 *     Copyright (C) 1995  Jose Alberto Mauricio                           *
 *     Copyright (C) 2009  Arthur B. Treadway & David E. Guerrero          *
 *       abtreadway@telefonica.net & warriord@rocketmail.com               *
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


#define OK      1
#define WRONG   0
#define MAXSTR 90
#define FREE_ARG char*

/*****************************************************************************/

typedef char * STRING;

/*****************************************************************************/
/* Numerical linear algebra                                                  */
/*****************************************************************************/

void ludcp( double **a, int n, int *ip );
void lusol( double **a, double *b, int n, int *ip );
void choldcp( double **mat, int n, double *d1, double *d2, int *ifault );
void cholfor( double **matl, int n, double *rhsol );
void cholbak( double **matl, int n, double *rhsol );
void cholsol( double **matl, int n, double *rhsol );
void gsl_eigenqr( double **a, int n, double *wr, double *wi );

/*****************************************************************************/
/* Error handler and dynamic memory allocation                               */
/*****************************************************************************/

void nrerror( char error_text[] );
double *vector( long nl, long nh );
int  *ivector( long nl, long nh );
double **matrix( long nrl, long nrh, long ncl, long nch );
int  **imatrix( long nrl, long nrh, long ncl, long nch );
double ***tensor( long nrl, long nrh, long ncl, long nch, long ndl, long ndh );
void free_vector( double *v, long nl, long nh );
void free_ivector( int *v, long nl, long nh );
void free_matrix( double **m, long nrl, long nrh, long ncl, long nch );
void free_imatrix( int **m, long nrl, long nrh, long ncl, long nch );
void free_tensor( double ***t, long nrl, long nrh, long ncl, long nch,
                  long ndl, long ndh );

/*****************************************************************************/
/* Scalar utilities                                                          */
/*****************************************************************************/

double rmax( double a, double b );
double rmin( double a, double b );
double cmacheps( void );
int iround( double num );
double pythag( double a, double b );

/*****************************************************************************/
/* String utilities                                                          */
/*****************************************************************************/

STRING NEW_STR( int size );
void FREE_STR( STRING s );
int DELETE_STR( STRING s, int i, int n );
int COPY_STR( STRING source, int i, int n, STRING dest );
int INSERT_STR( STRING s1, STRING s, int i );
int POS_STR( STRING s1, STRING s2 );
int CHANGE_STR( STRING s1, int i, STRING s2 );
void BLANKS_STR( STRING s );
void UPCASE_STR( STRING s );

/*****************************************************************************/
/* Date utilities                                                            */
/*****************************************************************************/

void Easter( int *day, int *month, int year );
